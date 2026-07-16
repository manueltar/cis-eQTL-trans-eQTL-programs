"""
POSTAR3_preprocessing_chr.py

Intersects POSTAR3 CLIP-seq RBP binding peaks with Ensembl GRCh38.110 transcript
coordinates, for a single chromosome (intended for use in an LSF job array, one
job per chromosome, given the scale of the full genome-wide join OOM'd at 32GB).

Key decisions baked into this script, and why:

1. TRANSCRIPT-level intervals, not gene-body intervals. CLIP-seq measures RBP
   binding to transcribed RNA, not genomic DNA outside the transcript (promoters,
   untranscribed gene-body regions), so transcript coordinates are the correct
   unit. Peak-transcript overlaps are collapsed to gene level afterward.

2. QC filters applied to the GTF before joining:
   - gene_biotype == 'artifact' excluded (Ensembl's own tag for annotation
     artifacts, confirmed by inspecting a real example, ENSG00000280071).
   - transcript_support_level == 5 excluded (lowest support tier). TSL is NOT
     required to be present - ~86,706 transcripts have no TSL at all, dominated
     by lncRNA/pseudogene/miRNA biotypes where TSL assessment doesn't apply by
     Ensembl convention, not because they're unreliable. Only the transcripts
     that WERE assessed and scored worst (TSL5) are dropped; NaN is kept.

3. Strand handling: peaks and transcripts are SPLIT BY STRAND before joining,
   then each strand is joined separately with strandedness=False. This was
   necessary because this pyranges version's native strandedness='same' path
   throws a KeyError - PyRanges internally forces the Strand column into a
   3-level categorical ('.', '-', '+') regardless of input dtype, and its
   internal strand-comparison dictionary only has '+'/'-' as keys, causing a
   crash when pandas' groupby iterates the unused '.' category. Splitting by
   strand first avoids this entirely (single-strand data never triggers the
   buggy comparison code) AND avoids computing/discarding cross-strand overlaps
   (confirmed via direct comparison: identical row count to a
   strandedness=False + pandas post-filter approach, at ~half the raw output
   size and no wasted computation).

4. nb_cpu: kept as a genuine parameter (see --nb_cpu) but expected to always be
   1 in practice. Ray-based parallelism (nb_cpu>1) was tested and found to add
   pure startup overhead (~600x slower: 12.24s vs 0.02s on an identical small
   test) with zero speed benefit at the scale where it was tested, and was a
   likely contributor to the original genome-wide OOM crash (Ray duplicates
   data per worker). Single-threaded is both faster and safer here.

5. Chromosome naming mismatch: POSTAR3 uses 'M' for mitochondrial, Ensembl GTF
   uses 'MT'. Handled via an explicit mapping, not a shared variable name, to
   avoid silently matching zero rows on either side.

6. gene_name is NaN for many real, legitimate genes (novel/uncharacterized loci,
   confirmed e.g. for ENSG00000279669, a real lncRNA with gene_biotype='lncRNA',
   simply lacking a curated symbol). NaN gene_name is filled with gene_id as a
   fallback for readability - the row is NOT dropped. This also avoids a
   separate, non-obvious bug: pandas groupby() defaults to dropna=True, which
   silently discards entire groups where ANY grouping key is NaN. Since
   gene_name was one of our grouping keys, unfilled NaNs would have silently
   deleted real minus-strand results from a test window during development.

7. Peak-level detail retained per RBP-gene-transcript pair, for use in
   applications beyond this specific enrichment analysis: n_peaks, methods
   (deduplicated), and three POSITION-ALIGNED pipe-joined strings
   (peak_coords, confidence_scores, peak_methods) - built from a single sorted
   groupby-apply so that the Nth entry in each string always refers to the same
   individual peak. NOT deduplicated against each other, since deduplicating
   peak_coords independently of confidence_scores would desynchronize them.

8. confidence_score is NOT used for filtering. Investigated across all 18
   experiment_method values in this data: only 3 tools (CLIPper: Poisson
   p-value, lower=better; PureCLIP: log posterior-probability ratio,
   higher=better; PARalyzer: bounded ~0-1 score, higher=better) had score
   semantics confirmed against their own documentation AND consistent with the
   observed value ranges. For Piranha, CIMS, and bare "eCLIP", the documented
   native output did not match the value ranges actually present in this file,
   and the substituted statistic could not be identified - using raw scores
   across all methods risked applying filters in the wrong direction for those
   3. Scores are kept in the output for reference only.

9. Confidence filter actually applied: passes_method_filter flags (not drops)
   rows supported by >=2 distinct experiment_method values, with at least one
   from CONFIDENT_METHODS (CLIPper, PureCLIP, PARalyzer - the 3 with confirmed
   score semantics). This sidesteps the cross-method score incomparability
   problem entirely by counting independent methods rather than thresholding
   an incomparable score.
"""

import argparse
import pandas as pd
import pyranges as pr
import time

parser = argparse.ArgumentParser()
parser.add_argument('--chrom', type=str, required=True,
                     help='Chromosome to process, e.g. "1", "X", "MT" (GTF naming; mapped internally to "M" for POSTAR3 lookup)')
parser.add_argument('--nb_cpu', type=int, default=1,
                     help='Cores for pyranges join. Expected to always be 1 - see script docstring point 4.')
parser.add_argument('--mem_gb', type=int, default=8,
                     help='Memory allocated to this job in GB (for logging/record-keeping only; actual enforcement is by LSF, not this script)')
args = parser.parse_args()

POSTAR_PATH = '/nfs/team151/mt19/POSTAR3/human.txt.gz'
GTF_PATH = '/nfs/users/nfs_m/mt19/cardinal_analysis/ht/datasets/Flanders/trans_eQTLS_ERs/Homo_sapiens.GRCh38.110.gtf.gz'
OUT_DIR = '/nfs/team151/mt19/POSTAR3/per_chrom'

CONFIDENT_METHODS = ['CLIPper', 'PureCLIP', 'PARalyzer']

# GTF uses 'MT', POSTAR3 uses 'M' - see docstring point 5
POSTAR_CHROM_NAME = 'M' if args.chrom == 'MT' else args.chrom

print(f'Processing chromosome {args.chrom} (POSTAR3 chrom label: {POSTAR_CHROM_NAME})', flush=True)
print(f'Allocated: {args.mem_gb}GB RAM, nb_cpu={args.nb_cpu}', flush=True)

# ---------------------------------------------------------------------------
# 1. Load POSTAR3 peaks for this chromosome only
# ---------------------------------------------------------------------------
print('Loading POSTAR3 peaks (this chromosome only, chunked read to limit memory)...', flush=True)
postar_cols = ['chrom', 'start', 'end', 'peak_id', 'strand', 'RBP_name',
               'experiment_method', 'sample', 'accession', 'confidence_score']
chunks = []
for chunk in pd.read_csv(POSTAR_PATH, sep='\t', names=postar_cols, compression='gzip',
                          low_memory=False, chunksize=2_000_000):
    chunk = chunk[chunk.chrom == f'chr{POSTAR_CHROM_NAME}']
    if len(chunk) > 0:
        chunks.append(chunk)
postar_chr = pd.concat(chunks, ignore_index=True) if chunks else pd.DataFrame(columns=postar_cols)
postar_chr = postar_chr[postar_chr.RBP_name != 'RBP_occupancy'].copy()  # generic PIP-seq placeholder, not a named RBP
postar_chr['chrom_clean'] = postar_chr.chrom.str.replace('^chr', '', regex=True)
print(f'  {len(postar_chr)} peaks on chr{POSTAR_CHROM_NAME}', flush=True)

# ---------------------------------------------------------------------------
# 2. Load GTF transcripts for this chromosome, apply QC filters
# ---------------------------------------------------------------------------
print('Loading GTF transcripts (this chromosome only)...', flush=True)
gtf_cols = ['seqname','source','feature','start','end','score','strand','frame','attribute']
gtf = pd.read_csv(GTF_PATH, sep='\t', comment='#', names=gtf_cols, dtype={'seqname': str}, low_memory=False)
transcripts = gtf[(gtf.feature == 'transcript') & (gtf.seqname == args.chrom)].copy()
transcripts['gene_id'] = transcripts.attribute.str.extract(r'gene_id "([^"]+)"')
transcripts['gene_name'] = transcripts.attribute.str.extract(r'gene_name "([^"]+)"')
transcripts['transcript_id'] = transcripts.attribute.str.extract(r'transcript_id "([^"]+)"')
transcripts['gene_biotype'] = transcripts.attribute.str.extract(r'gene_biotype "([^"]+)"')
transcripts['tsl'] = transcripts.attribute.str.extract(r'transcript_support_level "([^"]+)"')
transcripts['tsl_numeric'] = transcripts.tsl.str.extract(r'^(\d)')  # strips version-reassignment suffixes, e.g. "3 (assigned to previous version 11)"
transcripts['start_0based'] = transcripts.start - 1
transcripts['end_0based'] = transcripts.end
print(f'  {len(transcripts)} raw transcript records', flush=True)

n_before_artifact = len(transcripts)
transcripts = transcripts[transcripts.gene_biotype != 'artifact'].copy()
print(f'  excluded {n_before_artifact - len(transcripts)} artifact-biotype transcripts', flush=True)

n_before_tsl = len(transcripts)
transcripts = transcripts[transcripts.tsl_numeric != '5'].copy()  # NaN kept deliberately - see docstring point 2
print(f'  excluded {n_before_tsl - len(transcripts)} TSL5 transcripts, {len(transcripts)} remain', flush=True)

if len(postar_chr) == 0 or len(transcripts) == 0:
    print('No data on this chromosome after filtering, writing empty output.', flush=True)
    pd.DataFrame(columns=['RBP_name','gene_id','gene_name','transcript_id','Strand','Strand_b',
                           'n_peaks','methods','peak_coords','confidence_scores','peak_methods',
                           'passes_method_filter']).to_csv(
        f'{OUT_DIR}/chr{args.chrom}_rbp_gene_transcript_pairs.csv', index=False)
    exit(0)

# ---------------------------------------------------------------------------
# 3. Build PyRanges objects, strand-split join - see docstring point 3
# ---------------------------------------------------------------------------
print('Building PyRanges objects...', flush=True)
postar_pr = pr.PyRanges(pd.DataFrame({
    'Chromosome': postar_chr.chrom_clean, 'Start': postar_chr.start, 'End': postar_chr.end,
    'Strand': postar_chr.strand.astype(str), 'RBP_name': postar_chr.RBP_name, 'peak_id': postar_chr.peak_id,
    'experiment_method': postar_chr.experiment_method, 'confidence_score': postar_chr.confidence_score,
}))
transcripts_pr = pr.PyRanges(pd.DataFrame({
    'Chromosome': transcripts.seqname, 'Start': transcripts.start_0based, 'End': transcripts.end_0based,
    'Strand': transcripts.strand.astype(str), 'gene_id': transcripts.gene_id, 'gene_name': transcripts.gene_name,
    'transcript_id': transcripts.transcript_id,
}))

print('Running strand-split join...', flush=True)
t0 = time.time()
results = []
for strand in ['+', '-']:
    p_sub = postar_pr.subset(lambda df: df.Strand == strand)
    t_sub = transcripts_pr.subset(lambda df: df.Strand == strand)
    if len(p_sub) == 0 or len(t_sub) == 0:
        continue
    ov = p_sub.join(t_sub, strandedness=False, nb_cpu=args.nb_cpu)
    results.append(ov.df)
overlaps_df = pd.concat(results, ignore_index=True) if results else pd.DataFrame()
print(f'  done in {time.time()-t0:.1f}s, {len(overlaps_df)} overlap rows', flush=True)

if len(overlaps_df) == 0:
    print('No overlaps found, writing empty output.', flush=True)
    pd.DataFrame(columns=['RBP_name','gene_id','gene_name','transcript_id','Strand','Strand_b',
                           'n_peaks','methods','peak_coords','confidence_scores','peak_methods',
                           'passes_method_filter']).to_csv(
        f'{OUT_DIR}/chr{args.chrom}_rbp_gene_transcript_pairs.csv', index=False)
    exit(0)

# ---------------------------------------------------------------------------
# 4. Collapse to RBP -> gene -> transcript pairs
# ---------------------------------------------------------------------------
print('Collapsing to RBP -> gene -> transcript pairs...', flush=True)
overlaps_df['gene_name'] = overlaps_df['gene_name'].fillna(overlaps_df['gene_id'])  # see docstring point 6
overlaps_df['peak_coord'] = overlaps_df.Chromosome.astype(str) + '_' + overlaps_df.Start.astype(str) + '_' + overlaps_df.End.astype(str)

def join_peaks_scores_methods(df):
    # single sorted pass so peak_coords/confidence_scores/peak_methods stay
    # positionally aligned - see docstring point 7
    df_sorted = df.sort_values('peak_coord')
    return pd.Series({
        'peak_coords': '|'.join(df_sorted.peak_coord.tolist()),
        'confidence_scores': '|'.join(df_sorted.confidence_score.astype(str).tolist()),
        'peak_methods': '|'.join(df_sorted.experiment_method.tolist()),
    })

group_cols = ['RBP_name', 'gene_id', 'gene_name', 'transcript_id', 'Strand', 'Strand_b']

peak_detail = (overlaps_df.groupby(group_cols).apply(join_peaks_scores_methods, include_groups=False).reset_index())

pair_summary = (overlaps_df.groupby(group_cols)
                 .agg(n_peaks=('peak_id', 'nunique'),
                      methods=('experiment_method', lambda s: '|'.join(sorted(set(s)))))
                 .reset_index())

pair_summary = pair_summary.merge(peak_detail, on=group_cols)

def passes_method_filter(methods_str):
    # >=2 distinct methods, >=1 from CONFIDENT_METHODS - see docstring point 9
    methods_list = methods_str.split('|')
    n_methods = len(set(methods_list))
    has_confident = any(any(cm in m for cm in CONFIDENT_METHODS) for m in methods_list)
    return n_methods >= 2 and has_confident

pair_summary['passes_method_filter'] = pair_summary.methods.apply(passes_method_filter)

print(f'  {len(pair_summary)} unique RBP-gene-transcript rows', flush=True)
print(f'  {pair_summary.passes_method_filter.sum()} pass the >=2-method / >=1-confident filter', flush=True)

pair_summary.to_csv(f'{OUT_DIR}/chr{args.chrom}_rbp_gene_transcript_pairs.csv', index=False)
print(f'Saved chr{args.chrom} results.', flush=True)
