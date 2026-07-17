# cis-eQTL / trans-eQTL programs: annotation, enrichment, and pattern analysis

This repository documents the full workflow for annotating cis-eQTL/trans-eQTL
gene pairs against curated biological interaction resources (PPI, TF-target,
Ligand-Receptor, RBP-target, Peptidase-substrate), testing for enrichment at
both the edge level (permutation-based) and the factor level (Fisher's
exact/logistic regression), and identifying and visualizing specific
structural patterns (TF hubs, RBP hubs, Ligand-Receptor pairs, PPI-TF
cascades) within GWAS-colocalized cis/trans-eQTL "programs."

See [`Methods.md`](Methods.md) for the full methodological description
corresponding to each stage below.

## Repository structure

```
Bash_scripts/          LSF submission scripts (data download/preprocessing, permutation array)
Python_scripts/        Standalone preprocessing and enrichment-test scripts, called by the bash scripts
Jupyter_notebooks/     R-kernel notebooks for annotation, classification, enrichment, and visualization
Dependencies/          Conda environment specifications (exact package versions)
Methods.md             Full methods writeup
README.md              This file
```

## Environment setup

Four conda environments are used across the pipeline (see `Dependencies/`):

| Environment | Used for |
|---|---|
| `general_purpose` | Python preprocessing (STRING/CollecTRI/Liana/MEROPS parsing) |
| `pyranges_env` | POSTAR3 genomic-interval preprocessing |
| `mysql_db` | Standing up a local MEROPS MySQL instance |
| `R_general_purpose` | All Jupyter notebooks (R kernel) |

Recreate an environment with, e.g.:
```bash
conda env create -f Dependencies/R_general_purpose.yaml
```

## Pipeline, in execution order

### 1. Download and preprocess raw resources

- **Reference annotation**: `Homo_sapiens.GRCh38.110.gtf` (Ensembl GRCh38.110).
- **STRING**: `string_interactions.tsv` (STRING network export; node1/node2 as
  Ensembl gene IDs, with per-channel evidence scores including
  `experimentally_determined_interaction`).
- **CollecTRI**:
  ```bash
  wget https://zenodo.org/records/8192729/files/CollecTRI_regulons.csv
  ```
- **Liana** (Ligand-Receptor `omni_resource`):
  ```bash
  pip download liana==1.8.0 --no-deps -d .
  unzip liana-1.8.0-py3-none-any.whl liana/resource/omni_resource.csv -d .
  cp liana/resource/omni_resource.csv .
  ```
- **POSTAR3** (RBP CLIP-seq peaks): downloaded from http://postar.ncrnalab.org,
  then intersected with GRCh38.110 transcript coordinates per chromosome:
  ```bash
  bash Bash_scripts/403_POSTAR3_chr_array.sh normal
  ```
  (calls `Python_scripts/POSTAR3_preprocessing_v2.py` once per chromosome as
  an LSF array job; chromosomes that OOM at the default allocation are
  resubmitted at 24GB — see script comments for the retry list.)
- **MEROPS** (peptidase-substrate): a local MySQL instance is built from the
  MEROPS flat-file dump, then queried and exported.
  ```bash
  wget https://ftp.ebi.ac.uk/pub/databases/merops/current_release/meropsweb121.tar.gz
  conda create -y -p <env_path> -c conda-forge --override-channels mariadb
  bash Bash_scripts/400_Restarted_HT_outage_MEROPS_db_creation.sh normal
  bash Bash_scripts/401_Restarted_HT_outage_MEROPS_db_creation_extraction.sh normal
  ```
  → `Jupyter_notebooks/Merops_preprocessing.ipynb`

### 2. Standardize all resources into a common gene-pair schema

`Jupyter_notebooks/Preprocessing_datasets_of_pairs.ipynb`

Converts every resource above into a shared schema (`partner1_symbol`,
`partner1_ensembl_gene_id`, `partner2_symbol`, `partner2_ensembl_gene_id`,
`partner1_role`, `partner2_role`, `directed`, `stringency`,
`stringency_detail`, `resource_name`), at strict/lenient stringency levels
(two levels for STRING and MEROPS/Liana/CollecTRI/POSTAR3, since only these
subsets are carried forward into the final analysis).

**Output**: `annotation_preprocessed/{resource}_{stringency}.tsv` — six files
carried forward: `collectri_strict`, `postar3_strict`, `liana_lenient`,
`merops_strict`, `string_experimental_strict`, `string_experimental_lenient`.

### 3. Build the annotated cis/trans-eQTL edge table

`Jupyter_notebooks/Overhaul_classification_factors_v2.ipynb`

- Loads the base cis/trans-eQTL graph, stratified by *program*
  (`Whole_cis_tras_graph_with_programs.tsv`) — a **factor** is defined as the
  combination `(Source_module, Source_gene, program)`, which may span
  multiple `Source_cell`/`Target_cell` combinations.
- Classifies every factor's edge-type composition (`Within-only`,
  `Across-only`, `Mixed`) based on whether `Source_cell == Target_cell`
  per edge.
- Classifies every edge against all six annotation resources (both
  directions for directed resources), producing `edge_role`/`edge_resource`
  list-columns, plus a generic `add_cascade_edges_generic()`-style TF-cascade
  detector reused later in the pattern-finding notebook.
- Merges in the GWAS-colocalization flag at the individual-edge level.

**Output**: `whole_eqtl_annotated.rds` / `.tsv` — the single base table used
by every downstream notebook.

### 4. Edge-level enrichment (permutation-based)

```bash
bash Bash_scripts/405_enrichment_array_v2.sh normal
```

Submits one LSF job per (resource × direction × eQTL condition) combination,
each running `Python_scripts/run_enrichment_test_v3.py` — a degree-preserving
curveball permutation test (10,000 permutations) comparing observed
resource-annotation overlap against a randomized null.

`Jupyter_notebooks/Representation_ERs_v4.ipynb` aggregates the per-job result
files, applies Benjamini-Hochberg correction, computes Poisson confidence
intervals on fold enrichment, and produces the final forest plot
(`enrichment_forest_plot_v9`).

### 5. Factor-level enrichment (Fisher's exact / logistic regression)

`Jupyter_notebooks/Factor_ER_v2.ipynb`

A complementary, non-permutation-based test: for every factor
`(Source_module, Source_gene, program)`, tests whether presence of each
annotation resource (and TF-cascade) differs between `Within-only` and
`Across-only` factors, both on the full factor set and restricted to
GWAS-colocalized factors — size-corrected via logistic regression
(`n_edges` as covariate) alongside a raw Fisher's exact test for
transparency.

### 6. Pattern finding and visualization

`Jupyter_notebooks/Overhaul_finding_patterns_v5.ipynb`

Identifies and renders specific, mechanistically motivated structural
patterns among GWAS-colocalized factors (TF cis-eQTL hubs, RBP cis-eQTL
hubs, Ligand-Receptor pairs, PPI-reciprocal-TF, PPI-to-cascading-TF), each
rendered as a paired "description" graph (every trans-eQTL target, uniform
grey edges) and "topology" graph (only annotated/cascade edges, using an
identical node layout for direct visual cross-reference between the two).

### 7. Future steps

- Cascade pattern right now is only for TF but it can be extensible to the other annotations.
- Add other topologies
  * Add other resources of gene to gene relation: miRNAs-target genes, lncRNAs to target genes
  * Other same dataset topologies: cis-eQTL|trans-caQTL, cis-eQTL|trans-sQTL
  * Other same dataset exploration of cis-eQTL in cascading nodes
  * Exterior network topologies based on GO_BP or Kegg pathways to connect the nodes in the network

## Data provenance note

Every intermediate and final table used throughout this pipeline is derived
from primary sources (STRING, CollecTRI, Liana, POSTAR3, MEROPS, and the
project's own cis/trans-eQTL graph) with explicit sanity checks at each
merge/join step (see individual notebook cells for `stopifnot()` validation).
No values are fabricated, imputed, or assumed without an explicit note in the
corresponding notebook cell.
