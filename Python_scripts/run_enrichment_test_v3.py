"""
run_enrichment_test_v2.py

Runs one enrichment test: a single annotation resource, in a single direction,
against a single eQTL condition (edge.type x GWAS-colocalization subset), using
a degree-preserving curveball permutation null. Designed to be called once per
combination from an LSF array (see 404_enrichment_array_v2.sh) - each job handles
one (resource_key, direction, eqtl_condition) combination independently.

CHANGES FROM run_enrichment_test.py (v1), and why:

1. Single input file. v1 read a raw eQTL graph plus a separate GWAS-colocalized
   table and joined them via a 'graph' column built from matching Source_cell/
   Source_gene/Target_cell/Target_gene. That merge has already been done once,
   validated (0 colocalized edges missing from the base graph), and baked into
   whole_eqtl_annotated.tsv as a native boolean GWAS_colocalized column. This
   script reads that file directly - no re-merging, no 'graph' string column.

2. Direction is now an explicit, required argument (forward / reverse /
   undirected), rather than being implicit in how partner1/partner2 were laid
   out per resource. This is necessary because the edge-level classification
   this script is meant to reproduce (edge_roles$edge_role /
   edge_roles$edge_resource, built in R) tests directed resources in BOTH
   directions as two distinct, separately-labeled relationships
   (Source_TF->Target vs Target_TF->Source, etc.) - each direction is a
   separate biological claim and gets its own independent permutation test,
   not a single test where a hit in either direction counts. Undirected
   resources (STRING PPI variants) only ever get one direction: 'undirected'.

   - direction='forward':    pairs_set = {(partner1, partner2), ...}
                              tested against (Source_gene, Target_gene)
   - direction='reverse':    pairs_set = {(partner2, partner1), ...}
                              tested against (Source_gene, Target_gene)
                              (equivalent to testing (Target_gene, Source_gene)
                              against the original (partner1, partner2) set -
                              same logic as edge_roles$pair_key_rev in R)
   - direction='undirected': pairs_set = {frozenset((partner1, partner2)), ...}

3. --subtract_file allows building an incremental resource band (e.g. STRING
   experimental lenient MINUS strict, reproducing the 0.1<=score<0.4 band used
   in the R classification) by set-differencing two annotation files' pairs_set
   before testing. Only used for STRING_experimental_lenient_band in this
   round; kept general in case another incremental band is needed later.
   Gene universe (for assessability) is still drawn from the MAIN
   (--annotation_file) file only, not the subtracted one, consistent with
   "assessable = could this pair exist in the resource this test is actually
   about", not the resource being subtracted out.

4. Resource identity is read from the file's own resource_name column, exactly
   as in v1 - this script remains resource-agnostic and does not hardcode
   resource-specific logic beyond direction and the optional subtraction.

5. Zero-observed-hit handling, curveball trades, n_trades heuristic (5x source
   nodes), and the Z=NaN guard for null_sd==0 are UNCHANGED from v1 - these
   were already validated across the prior PPI/TF-target/L-R enrichment work
   and there is no reason to alter them here.

6. GWAS_colocalized is read as a string column from the tsv (R's fwrite writes
   TRUE/FALSE as literal text) and explicitly coerced to boolean - NOT assumed
   to auto-parse correctly, since pandas' automatic bool inference is version-
   and locale-dependent and silently getting this wrong would silently corrupt
   every GWAS_colocalized condition.
"""

import argparse
import os
import random
import time

import numpy as np
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('--annotation_file', type=str, required=True,
                     help='Path to a preprocessed annotation TSV (common schema)')
parser.add_argument('--subtract_file', type=str, default=None,
                     help='Optional: path to a second annotation TSV whose pairs_set is '
                          'subtracted from --annotation_file\'s pairs_set before testing '
                          '(e.g. lenient MINUS strict, to isolate an incremental band). '
                          'Only meaningful for direction=undirected in the current use case.')
parser.add_argument('--direction', type=str, required=True,
                     choices=['forward', 'reverse', 'undirected'],
                     help='forward: Source_gene tested as partner1 (regulator). '
                          'reverse: Target_gene tested as partner1 (regulator), i.e. '
                          'Source_gene is the regulated partner. '
                          'undirected: unordered pair match (STRING only).')
parser.add_argument('--eqtl_condition', type=str, required=True,
                     choices=['Across_full', 'Across_GWAS_colocalized', 'Within_full', 'Within_GWAS_colocalized'])
parser.add_argument('--n_perm', type=int, default=10000)
parser.add_argument('--seed', type=int, default=42)
parser.add_argument('--mem_gb', type=int, default=8,
                     help='Memory allocated to this job in GB (for logging/record-keeping only; actual enforcement is by LSF, not this script)')
args = parser.parse_args()

GTF_PATH = '/nfs/users/nfs_m/mt19/cardinal_analysis/ht/datasets/Flanders/trans_eQTLS_ERs/Homo_sapiens.GRCh38.110.gtf.gz'
EQTL_PATH = '/nfs/team151/mt19/overhaul_classification_factors_with_programs/whole_eqtl_annotated.tsv'
OVERLAPS_DIR = '/nfs/team151/mt19/overhaul_classification_factors_with_programs/overlaps_v2'
RESULTS_DIR = '/nfs/team151/mt19/overhaul_classification_factors_with_programs/results_v2'

os.makedirs(OVERLAPS_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)

print(f'Annotation file: {args.annotation_file}', flush=True)
print(f'Subtract file: {args.subtract_file}', flush=True)
print(f'Direction: {args.direction}', flush=True)
print(f'eQTL condition: {args.eqtl_condition}', flush=True)
print(f'Allocated: {args.mem_gb}GB RAM', flush=True)

# ---------------------------------------------------------------------------
# 0. GTF reverse mapping (ENSG -> symbol), used only for display in the saved
#    overlap detail table
# ---------------------------------------------------------------------------
print('Loading GTF for ENSG->symbol mapping...', flush=True)
gtf_cols = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
gtf = pd.read_csv(GTF_PATH, sep='\t', comment='#', names=gtf_cols, dtype={'seqname': str}, low_memory=False)
genes = gtf[gtf.feature == 'gene'].copy()
genes['gene_id'] = genes.attribute.str.extract(r'gene_id "([^"]+)"')
genes['gene_name'] = genes.attribute.str.extract(r'gene_name "([^"]+)"')
ens2sym = genes[['gene_name', 'gene_id']].dropna().drop_duplicates().set_index('gene_id')['gene_name'].to_dict()
print(f'  {len(ens2sym)} ENSG->symbol mappings loaded', flush=True)

# ---------------------------------------------------------------------------
# 1. Load whole_eqtl_annotated.tsv and filter to the requested condition
# ---------------------------------------------------------------------------
print('Loading whole_eqtl_annotated.tsv...', flush=True)
eqtl_full = pd.read_csv(EQTL_PATH, sep='\t')
print(f'  {len(eqtl_full)} rows total', flush=True)

# Explicit, non-assumed boolean coercion - see docstring point 6
eqtl_full['GWAS_colocalized'] = eqtl_full['GWAS_colocalized'].astype(str).str.upper() == 'TRUE'
n_coloc = eqtl_full['GWAS_colocalized'].sum()
print(f'  {n_coloc} rows flagged GWAS_colocalized (expect 11400 if this matches the validated R build)', flush=True)

edge_type, coloc_filter = args.eqtl_condition.split('_', 1)
if coloc_filter == 'full':
    eqtl_subset = eqtl_full[eqtl_full['edge.type'] == edge_type].copy()
else:
    eqtl_subset = eqtl_full[(eqtl_full['edge.type'] == edge_type) & (eqtl_full['GWAS_colocalized'])].copy()
print(f'  {args.eqtl_condition}: {len(eqtl_subset)} rows', flush=True)

pairs_full = eqtl_subset[['Source_gene', 'Target_gene']].drop_duplicates().reset_index(drop=True)
print(f'  {len(pairs_full)} unique gene pairs', flush=True)

# ---------------------------------------------------------------------------
# 2. Load the annotation file (and optional subtract_file)
# ---------------------------------------------------------------------------
print('Loading annotation file...', flush=True)
annot = pd.read_csv(args.annotation_file, sep='\t')
resource_name = annot.resource_name.iloc[0]
directed_flag = bool(annot.directed.iloc[0])
stringency_value = annot.stringency.iloc[0]
stringency_detail = annot.stringency_detail.iloc[0]
print(f'  resource: {resource_name}, directed={directed_flag}, {len(annot)} rows', flush=True)

if args.direction != 'undirected' and not directed_flag:
    raise ValueError(f'--direction {args.direction} requested but annotation file has directed=False. '
                      f'Undirected resources must use --direction undirected.')
if args.direction == 'undirected' and directed_flag:
    raise ValueError('--direction undirected requested but annotation file has directed=True. '
                      'Directed resources must use --direction forward or reverse.')

subtract_annot = None
if args.subtract_file is not None:
    print(f'Loading subtract_file...', flush=True)
    subtract_annot = pd.read_csv(args.subtract_file, sep='\t')
    print(f'  subtract resource: {subtract_annot.resource_name.iloc[0]}, {len(subtract_annot)} rows', flush=True)

# ---------------------------------------------------------------------------
# 3. Build pairs_set / gene universe (direction-aware)
# ---------------------------------------------------------------------------
def build_pairs_set(annot_df, direction):
    if direction == 'undirected':
        pairs = list(zip(annot_df.partner1_ensembl_gene_id, annot_df.partner2_ensembl_gene_id))
        return set(frozenset(p) for p in pairs)
    elif direction == 'forward':
        return set(zip(annot_df.partner1_ensembl_gene_id, annot_df.partner2_ensembl_gene_id))
    elif direction == 'reverse':
        return set(zip(annot_df.partner2_ensembl_gene_id, annot_df.partner1_ensembl_gene_id))
    else:
        raise ValueError(direction)

def build_gene_universe(annot_df):
    return set(annot_df.partner1_ensembl_gene_id) | set(annot_df.partner2_ensembl_gene_id)

pairs_set = build_pairs_set(annot, args.direction)
gene_universe = build_gene_universe(annot)

if subtract_annot is not None:
    subtract_pairs_set = build_pairs_set(subtract_annot, args.direction)
    n_before = len(pairs_set)
    pairs_set = pairs_set - subtract_pairs_set
    print(f'  pairs_set after subtraction: {n_before} -> {len(pairs_set)}', flush=True)
    stringency_detail = f'{stringency_detail} MINUS {subtract_annot.stringency_detail.iloc[0]}'

print(f'  pairs_set size: {len(pairs_set)}, gene universe size: {len(gene_universe)}', flush=True)

# ---------------------------------------------------------------------------
# 4. Assessable universe
# ---------------------------------------------------------------------------
assessable = pairs_full[pairs_full.Source_gene.isin(gene_universe) & pairs_full.Target_gene.isin(gene_universe)]
print(f'Assessable eQTL pairs: {len(assessable)}', flush=True)
print(f'  n unique source-side genes: {assessable.Source_gene.nunique()}', flush=True)
print(f'  n unique target-side genes: {assessable.Target_gene.nunique()}', flush=True)

# ---------------------------------------------------------------------------
# 5. Detailed overlap - always produced, even when empty
# ---------------------------------------------------------------------------
def observed_overlap_detail(assessable_df, pairs_set, direction):
    records = []
    for src, tgt in zip(assessable_df.Source_gene, assessable_df.Target_gene):
        key = frozenset((src, tgt)) if direction == 'undirected' else (src, tgt)
        if key in pairs_set:
            records.append({'Source_gene': src, 'Target_gene': tgt})
    if records:
        return pd.DataFrame(records).drop_duplicates()
    return pd.DataFrame(columns=['Source_gene', 'Target_gene'])

overlap = observed_overlap_detail(assessable, pairs_set, args.direction)
print(f'Observed hits: {len(overlap)}', flush=True)

overlap_columns = ['Source_gene', 'Target_gene', 'Source_symbol', 'Target_symbol'] + \
                   [c for c in eqtl_subset.columns if c not in ('Source_gene', 'Target_gene')] + \
                   ['resource_name', 'stringency', 'stringency_detail', 'direction', 'eqtl_condition']

if len(overlap) == 0:
    overlap_with_context = pd.DataFrame(columns=overlap_columns)
else:
    overlap_annotated = overlap.copy()
    overlap_annotated['Source_symbol'] = overlap_annotated.Source_gene.map(ens2sym)
    overlap_annotated['Target_symbol'] = overlap_annotated.Target_gene.map(ens2sym)
    overlap_with_context = overlap_annotated.merge(eqtl_subset, on=['Source_gene', 'Target_gene'], how='left')
    overlap_with_context['resource_name'] = resource_name
    overlap_with_context['stringency'] = stringency_value
    overlap_with_context['stringency_detail'] = stringency_detail
    overlap_with_context['direction'] = args.direction
    overlap_with_context['eqtl_condition'] = args.eqtl_condition

resource_label = f'{resource_name}_{stringency_value}_{args.direction}_{args.eqtl_condition}'
overlap_out_path = f'{OVERLAPS_DIR}/{resource_label}_overlap.tsv'
overlap_with_context.to_csv(overlap_out_path, sep='\t', index=False)
print(f'Saved overlap detail to {overlap_out_path}', flush=True)

# ---------------------------------------------------------------------------
# 6. Adjacency for the curveball permutation
# ---------------------------------------------------------------------------
def build_curveball_adjacency(assessable_df):
    all_genes = sorted(set(assessable_df.Source_gene) | set(assessable_df.Target_gene))
    gid = {g: i for i, g in enumerate(all_genes)}
    source_nodes = sorted(assessable_df.Source_gene.unique())
    adj0 = {gid[s]: set(gid[t] for t in assessable_df.loc[assessable_df.Source_gene == s, 'Target_gene']) for s in source_nodes}
    return adj0, gid, source_nodes

adj0, gid, source_nodes = build_curveball_adjacency(assessable)
n_edges = sum(len(v) for v in adj0.values())
print(f'Adjacency: {n_edges} edges, {len(source_nodes)} source nodes', flush=True)

# ---------------------------------------------------------------------------
# 7. Curveball permutation - UNCHANGED from v1
# ---------------------------------------------------------------------------
def curveball_trades(adj, n_trades, rng, src_list):
    for _ in range(n_trades):
        a, b = rng.sample(src_list, 2)
        A, B = adj[a], adj[b]
        only_A = A - B
        only_B = B - A
        if not only_A and not only_B:
            continue
        pool = list(only_A | only_B)
        rng.shuffle(pool)
        kA = len(only_A)
        new_A = set(pool[:kA])
        new_B = set(pool[kA:])
        if a in new_A or b in new_B:  # forbid self-loops
            continue
        common = A & B
        adj[a] = common | new_A
        adj[b] = common | new_B

def run_permutation(adj0, gid, pairs_set, direction, n_perm, seed):
    n_trades = 5 * max(len(adj0), 1)

    if direction == 'undirected':
        pairs_enc = set()
        for fs in pairs_set:
            a, b = tuple(fs)
            if a in gid and b in gid:
                pairs_enc.add(frozenset((gid[a], gid[b])))
    else:
        pairs_enc = set((gid[a], gid[b]) for a, b in pairs_set if a in gid and b in gid)

    rng = random.Random(seed)
    adj = {k: set(v) for k, v in adj0.items()}
    src_list = list(adj.keys())

    null = np.zeros(n_perm, dtype=int)
    t0 = time.time()
    for i in range(n_perm):
        curveball_trades(adj, n_trades, rng, src_list)
        if direction == 'undirected':
            n_hits = sum(1 for a, targets in adj.items() for t in targets if frozenset((a, t)) in pairs_enc)
        else:
            n_hits = sum(1 for a, targets in adj.items() for t in targets if (a, t) in pairs_enc)
        null[i] = n_hits
        if (i + 1) % 2000 == 0:
            print(f'  {i+1}/{n_perm} done ({time.time()-t0:.0f}s)', flush=True)
    print(f'  total permutation runtime: {time.time()-t0:.0f}s', flush=True)
    return null

if len(source_nodes) < 2:
    print('Fewer than 2 source nodes - curveball trades are impossible, skipping permutation.', flush=True)
    null = np.zeros(args.n_perm, dtype=int)
else:
    null = run_permutation(adj0, gid, pairs_set, args.direction, args.n_perm, args.seed)

# ---------------------------------------------------------------------------
# 8. Results row
# ---------------------------------------------------------------------------
def build_results_row(obs, null, resource_name, stringency_detail, direction, eqtl_condition):
    mean, sd = null.mean(), null.std(ddof=1)
    z = (obs - mean) / sd if sd > 0 else np.nan
    p = (np.sum(null >= obs) + 1) / (len(null) + 1)
    is_untestable = bool(sd == 0)
    return pd.DataFrame({
        'analysis': [resource_name],
        'direction': [direction],
        'condition': [f'{stringency_detail} | eQTL: {eqtl_condition}'],
        'observed_hits': [obs],
        'null_mean': [mean],
        'null_sd': [sd],
        'z_score': [z],
        'p_enrichment': [p],
        'is_untestable': [is_untestable],
    })

result_row = build_results_row(len(overlap), null, resource_name, stringency_detail, args.direction, args.eqtl_condition)
print(result_row.to_string(index=False), flush=True)

result_out_path = f'{RESULTS_DIR}/{resource_label}_result.tsv'
result_row.to_csv(result_out_path, sep='\t', index=False)
print(f'Saved result row to {result_out_path}', flush=True)
