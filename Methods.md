# Methods

## Data

 Cis-eQTL/trans-eQTL relationships were obtained from a cis-trans graph spanning multiple immune cell types, stratified by co-expression program membership (n=305,496 rows in the current build). A **factor**  is defined as the combination `(Source_module, Source_gene, program)`, which may span multiple `Source_cell`/`Target_cell` combinations; the same genomic locus/gene pair can therefore be split across several programs, each an independent factor for classification and testing purposes. Rows were classified by `edge.type` as *Within* (same source and target cell type) or *Across* (different source cell type from target cell type), and each factor was further classified by its edge-type composition as `Within-only`, `Across-only`, or `Mixed` (containing both). Under program stratification, no factors in the current build have `Mixed` composition — the finer granularity fully resolves what appeared as mixed composition under the coarser `(Source_module, Source_gene)` definition used previously, verified directly by tabulation rather than assumed. GWAS colocalization is tracked at the individual-edge level (`GWAS_colocalized`, matched on the full `Source_module`/`program`/`Source_cell`/`Target_gene`/`Target_cell` key against a colocalization table, with every colocalized edge confirmed present in the base graph before use), rather than collapsed across cell-type contexts as in the earlier design. Four eQTL subsets were tested independently: Within/full, Within/GWAS-colocalized, Across/full, Across/GWAS-colocalized. Gene symbols were mapped to Ensembl gene IDs using Ensembl GRCh38 release 110 gene annotations throughout, with an ALIAS-based fallback for symbols not resolved by the primary mapping; where the fallback produced two distinct Ensembl IDs sharing one display symbol within the same rendered network, the symbol was disambiguated for that graph only (Ensembl ID suffix appended) to prevent two genuinely distinct genes being merged into a single graph node.

## Reference annotations

Five resources were tested, each pre-processed into a standardized directed or undirected gene-pair schema at two stringency levels (except STRING, tested at three):

1. **Protein-protein interactions (PPI)**, from STRING (v12; undirected), using the *experimentally_determined_interaction* evidence channel alone (rescaled via STRING's documented noisy-OR formula, prior *p*=0.041, to recover a 0–1 score from that single channel), independent of text-mining, curated-database, and coexpression evidence — chosen to avoid circularity with trans-eQTL coexpression signal and to restrict to experimentally-supported physical interaction evidence specifically. Tested at two non-overlapping bands: strict (score ≥0.4) and an incremental lenient band (0.1≤score<0.4, obtained by set-subtracting the strict pair set from the lenient pair set, so no pair is double-counted across the two tests).


2. **Transcription factor (TF)–target relationships**, from CollecTRI, tested both filtered (≥2 independent supporting sub-resources; strict) and unfiltered (lenient).

3. **Ligand-receptor interactions**, from the Liana omni_resource (aggregating multiple curated ligand-receptor databases), with mouse-derived and multi-subunit complex entries excluded, tested using the CellPhoneDB sub-resource (strict) and the broader consensus sub-resource (lenient).

4. **RNA-binding protein (RBP)–target transcript relationships**, from POSTAR3 CLIP-seq peak data intersected with Ensembl transcript coordinates (strand-matched), tested requiring support from ≥2 independent experimental methods with ≥1 method having confirmed, documented score semantics (CLIPper, PureCLIP, or PARalyzer; strict) or without this filter (lenient).

5. **Peptidase-substrate relationships**, from MEROPS, tested restricted to physiological cleavage events (strict) or all annotated cleavage types including non-physiological and synthetic (lenient).

6. **Same dataset SCENIC+ regulons**,

A sixth and seventh reference annotation were added: TF→target relationships inferred by SCENIC+ from the same scRNA-seq dataset used to define the eQTL cis/trans modules. Regulons were provided as TF(sign)→target-gene-list pairs, where sign (`+`/`-`) denotes an activator or repressor regulon as called by SCENIC+. Activator and repressor regulons were split into two fully independent resources (`SCENIC_regulon_activator_same_dataset`,`SCENIC_regulon_repressor_same_dataset`) rather than one resource with a sign attribute. This allows either to be included, excluded, or weighted differently in any downstream analysis without needing to re-derive the split, and keeps the lower-confidence repressor calls from silently diluting or being conflated with the activator calls in enrichment testing, factor-level testing, or visualization. Both resources use a single stringency tier (presence in the regulon; no confidence/AUC threshold was available from the source data) and were mapped to Ensembl gene IDs via the same GRCh38.110 reference used throughout; a small number of gene symbols resolved to more than one Ensembl ID (segmental duplication/pseudogene pairs and antisense-transcript pairs) were retained via row explosion rather than arbitrarily resolved to a single ID. TF-autoregulation (a regulon target identical to its own TF) was retained and flagged (`is_autoregulation`).

## Addition — cis- and trans-eQTL effect sizes

Effect-size information was subsequently incorporated for both the cis-eQTL association and the trans-eQTL (Mendelian Randomization, MR) relationship: cis-eQTL effect (`most_likely_snp`, `most_likely_snp_beta`): unique per `(Source_module, program, Source_cell, Source_gene)`. trans-eQTL effect (`MR_effect`): merged onto `whole_eqtl_annotated` at the full edge key (`Source_module`, `program`, `Source_cell`,`Source_gene`, `Target_cell`, `Target_gene`). `MR_effect` takes a single magnitude shared across all of that hub's trans-eQTL targets, with only its sign varying per target (`MR_effect = program_beta × Target_gene_direction`. `MR_effect` should therefore be interpreted as encoding, per trans-eQTL target, whether that target's regulatory direction is concordant or discordant with the hub's single program-level effect magnitude — not as an independently estimated, continuously-varying effect size per target gene.

## Enrichment testing

For the six directed resources (CollecTRI, POSTAR3, Liana, MEROPS, SCENIC activator, SCENIC repressor), each gene pair was tested in **both directions independently** — forward (cis-eQTL gene as the regulatory/interacting partner, trans-eQTL gene as the target) and reverse (trans-eQTL gene as the regulatory/interacting partner, cis-eQTL gene as the target) — reflecting that these represent two distinct, non-equivalent biological hypotheses (e.g. "the cis-eQTL gene is a transcription factor for this trans-eQTL target" vs. "this trans-eQTL target's gene product regulates the cis-eQTL gene"). STRING (undirected) was tested with a single unordered-pair match. This yields 14 resource×direction combinations (12 directed + 2 undirected/STRING bands) × 4 eQTL conditions = 56 tests total in the final design (superseding the earlier "44 tests" and "5 resources × stringency" framing).  Statistical significance was assessed against a degree-preserving permutation null generated by the curveball algorithm, which randomizes bipartite source–target connectivity while exactly preserving each gene's observed degree within the assessable universe, with self-pairings forbidden to match the empirical data (10,000 permutations per test; trade count per permutation set to 5× the number of source-role nodes). Enrichment was quantified as Z = (observed − mean(null)) / SD(null) and as fold enrichment (observed / mean(null)), with one-sided empirical p-values (p = [1 + #(null ≥ observed)] / [1 + n_permutations]). Combinations in which the null distribution had zero variance (no permutation ever produced a hit, reflecting insufficient overlap between the eQTL and annotation gene sets rather than a true null result) were reported as untestable rather than assigned a spurious Z-score or p-value.

## Multi-level enrichment framework

Beyond the edge-level permutation test, resource-annotation enrichment was tested at three additional, progressively coarser grains, each with its own composition/Fisher's-exact-test logic:

1. **cis-eQTL (unique genes)** — `Source_gene` alone, collapsed across `Source_module` and `program`. The only grain at which `Mixed` (Within in some contexts, Across in others) survives rather than being collapsed away; tested with two independent Fisher families (`Across-only` vs. `Within-only`, and `Mixed` vs. `Within-only`), each BH-corrected separately.
2. **Unique cis-eQTL–trans-eQTL edges** — one row per unique `(Source_gene, Target_gene)` pair rather than per recurring eQTL-edge instance, addressing hub-gene concentration in the raw edge counts.
3. **cis/trans network** — the `(Source_module, Source_gene, program)` factor grain used throughout the rest of this project (previously referred to as "the factor level"; renamed for the combined multi-level figure to avoid a naming collision that had earlier caused this grain and the gene-level grain above to be silently tested as if identical).

All three grains, plus the edge-level permutation results, are combined into a single multi-panel figure; the edge-level panel uses fold-enrichment against a permutation null rather than a Fisher odds ratio and is not directly comparable in magnitude to the other three, only in direction/significance.

## cis/trans network-level enrichment (Fisher's exact test / logistic regression)

Complementary to the edge-level permutation enrichment, a factor-level ("cis/trans network") test was performed. For every factor `(Source_module, Source_gene, program)` with `Within-only`/`Across-only` composition (`Mixed` excluded — structurally absent under program stratification, see Data), binary presence/absence was determined for 14 resource×direction categories: CollecTRI, SCENIC activator, SCENIC repressor, POSTAR3 (RBP), Liana (Ligand-Receptor), and MEROPS, each split forward (cis-eQTL gene as regulator/ligand/peptidase) vs. reverse (cis-eQTL gene as the regulated/receptor/substrate partner), plus PPI (STRING-experimental) at its strict and lenient stringency bands (undirected, so no direction split). This is tested on the full graph only — no separate GWAS-colocalized subset at this level. TF-cascade is not tested as a factor-level category here (see the Pattern-finding section for cascade-specific analysis).

Two statistical approaches were compared for each category:

1. **Raw comparison**: Fisher's exact test, odds ratio and exact confidence interval, with explicit factor-level ordering (`Within-only` before `Across-only`) fixing OR = odds(Across)/odds(Within).
2. **Size-corrected comparison**: logistic regression    (`has_<category> ~ edge_type_composition + n_edges`, binomial family),   including factor size (total edge count) as a covariate, since factor   size is confounded with `Across-only`/`Within-only` status (Across-only   factors were, on average, ~5× larger by edge count than Within-only   factors in this dataset) and with raw annotation-presence rate more   generally (established separately in the edge-level enrichment work).   Odds ratios and 95% confidence intervals were computed via the Wald    approximation (estimate ± 1.96×SE, exponentiated).

Benjamini-Hochberg correction was applied once across all 14 categories within each statistical approach (no separate pre-specified-hypothesis/exploratory split at this level).

## Generic TF-cascade detection

A gene pair-level, resource-agnostic cascade detector was implemented and applied uniformly across every pattern and enrichment analysis in this project (not limited to a single pre-identified transcription factor). For a given set of trans-eQTL target genes belonging to one factor/cell-type sub-unit, the detector tests every pairwise combination of those target genes (excluding the cis-eQTL hub gene) against the CollecTRI strict TF→target table, flagging any pair where one target gene regulates another target gene within the same sub-unit. This surfaces cascading regulatory structure (e.g. a cis-eQTL's trans-eQTL target that is itself a transcription factor independently regulating a second, sibling trans-eQTL target of the same cis-eQTL) that would be missed by testing the cis-eQTL hub gene's direct relationships alone.


## Pattern-finding and visualization methodology

1. **TF_direct_cis_eQTL** cis-eQTL gene is a TF (CollecTRI, SCENIC activator, or SCENIC repressor) with ≥1 forward-direction (cis-eQTL as regulator) direct target among the factor's own trans-eQTL genes. (Note: with SCENIC activator/repressor now separate resources with their own edge-role labels — `Source_TF_activator->Target`, `Source_TF_repressor->Target`.

2. **PPI_TF_cascade** cis-eQTL gene is not itself a TF, is PPI-linked (STRING-experimental, either stringency band) to a trans-eQTL target that is generically a TF, and that specific TF partner independently regulates ≥1 sibling trans-eQTL target within the same `(Source_cell, Target_cell)` unit (i.e. the cascade criterion is now evaluated for the specific PPI-partner TF, not "any cascade present anywhere in the unit").

3. **Liana_LR_interaction** cis-eQTL gene participates in a Liana ligand-receptor pair with a trans-eQTL target (either direction).

4. **RBP_direct_cis_eQTL** cis-eQTL gene is an RBP (POSTAR3) with ≥1 forward-direction direct target among the factor's trans-eQTL genes.

5. **MEROPS_peptidase_substrate** cis-eQTL gene participates in a MEROPS peptidase-substrate pair with a trans-eQTL target (either direction).

6. **PPI_strict** searched only within factors unclassified by patterns 1–5; cis-eQTL gene has a strict-band (≥0.4) STRING-experimental PPI link to ≥1 trans-eQTL target, no TF/cascade requirement.

7. **PPI_lenient** searched only within factors unclassified by patterns 1–6; identical to pattern 6 but restricted to the lenient PPI band (0.1≤score<0.4).

For each pattern, an odds ratio for `Across-only` vs. `Within-only` composition was computed (Fisher's exact test), with the comparison background restricted to other classified factors only (excluding `Mixed`-composition and unclassified factors) — a narrower, pattern-focused complement to the full-factor-set Fisher's/logistic-regression comparison described in the "cis/trans network-level enrichment" section above, answering "is this pattern's composition distinctive relative to other structurally-annotated factors" rather than "relative to all factors."

## visualization: SCENIC sign, effect-size encoding, node ordering.

The description graph (all trans-eQTL targets, previously uniform grey) now colors each hub-spoke edge by the sign of `MR_effect` where available (with an explicit "no MR estimate" category for edges lacking a merged value). Node (spoke) ordering around the hub is now determined by `MR_effect` sign (negative-to-positive), with cascading TFs clustered within each sign block, applied consistently across the description graph and both topology panels so that spatial position remains comparable across all three views of a given factor. Each rendered factor's title additionally reports the cis-eQTL's SNP identifier and effect size (`most_likely_snp`, `most_likely_snp_beta`).

## Blood-exosome (healthy donor) membership

Gene-level presence in blood extracellular vesicles/particles (EVPs) of healthy donors was added as an additional enrichment category, from exoRBase 3.0 (`longRNAs_anno.csv`, exorbase.org), using the `Healthy frequency of EVPs(Sample number)` field. exoRBase stratifies EVP detection frequency by disease state (Healthy/Benign/Tumor) for blood specifically — other biofluids (urine, CSF, bile) carry a single undifferentiated frequency column — so this field is blood-specific by construction, without an additional sample-type filter. The parenthetical sample count in this field is a fixed cohort size (244 healthy donors), not a per-gene detected count; a per-gene detected-sample count was derived as `round(frequency × cohort size)`. A gene was scored present in healthy blood exosomes if it had a detected-sample count of ≥200 (out of the 244-donor healthy cohort) and a detection frequency of ≥0.9 among assessed healthy-donor samples; genes without any exoRBase healthy-EVP value were coded `NA` ("not assessed") and excluded from testing rather than treated as absent. Gene identifiers were matched to `whole_eqtl_annotated` by Ensembl ID with version suffixes stripped.

This category was tested identically to the 14 resource×direction categories above, as the 15th row in the cis-eQTL (unique genes) panel of the combined multi-level figure (gene-level presence, aggregated from the underlying `(Source_module, Source_gene, program)` annotation via `ANY()` per gene).

## Limitation

The permutation null controls for connectivity (degree) within the eQTL network itself but does not correct for literature/study ascertainment bias in the annotation resources — genes that are more extensively studied tend to accumulate more annotations of every type (PPI, TF-target, etc.) independent of true biological interaction, and this confound is not addressed by the present design.




---

*Notes: citations for STRING, CollecTRI, Liana, POSTAR3, MEROPS, SCENIC+, exoRBase 3.0, and the curveball algorithm (Strona et al. 2014) need to be added via Paperpile. Per-combination assessable-universe sizes (56 combinations) are best reported in a results table/supplement rather than in prose.*
