#!/bin/bash
# 404_enrichment_array_v2.sh
#
# See before:
# ~/Jupyter_notebooks/HT_OUTAGE/CARDINAL/Preprocessing_cis-eQTL_trans-eQTL_pairs.ipynb
# ~/Jupyter_notebooks/HT_OUTAGE/CARDINAL/Preprocessing_datasets_of_pairs.ipynb
#
# CHANGES FROM 404_Permutation_test_v2.sh (v1), and why:
#
# 1. Single input file (whole_eqtl_annotated.tsv) - no eqtl_conditions logic
#    change needed there, GWAS_colocalized filtering now happens natively
#    inside run_enrichment_test_v2.py via the boolean column, not a separate
#    table join.
#
# 2. Resource list is now the 7 resources used to build edge_role /
#    edge_resource in the R classification (originally 6; SCENIC_regulon_
#    same_dataset added as a 7th - see point 6). This deliberately excludes
#    collectri_lenient, postar3_lenient, liana_strict, and the raw
#    combined-score string_ppi.tsv, none of which were part of the
#    classification this permutation test is meant to validate/extend.
#    If those are wanted later, they should be a separate, explicitly-scoped
#    run - not silently folded in here.
#
# 3. Each directed resource now gets TWO jobs per eqtl_condition (forward,
#    reverse), matching edge_role's Source_X->Target / Target_X->Source
#    distinction. Each undirected resource (STRING variants) gets ONE job
#    per eqtl_condition (direction=undirected).
#
# 4. STRING_experimental_lenient is run as an INCREMENTAL BAND (lenient MINUS
#    strict, via --subtract_file), reproducing the "PPI (STRING_experimental
#    >=0.1 <0.4)" label from the R classification - NOT the full lenient set,
#    which would double-count the strict pairs already tested separately.
#
# 5. Memory allocations carried over from v1 where the resource is unchanged
#    (POSTAR3 still by far the largest). Directed resources now run 2x as
#    many jobs (forward+reverse) at the SAME per-job memory as before, since
#    each direction's adjacency/permutation is no larger than the original
#    single-direction test.
#
# scenic_regulon_activator_same_dataset.tsv
# scenic_regulon_repressor_same_dataset.tsv


queue=$1

condaenv=/nfs/users/nfs_m/mt19/cardinal_analysis/ht/conda_envs/general_purpose
script=/nfs/users/nfs_m/mt19/Scripts/Python_scripts/run_enrichment_test_v3.py
annot_dir=/nfs/team151/mt19/annotation_preprocessed
outdir=/nfs/team151/mt19/overhaul_classification_factors_with_programs/results_v2
logdir=/nfs/team151/mt19/overhaul_classification_factors_with_programs/logs_v2

mkdir -p $outdir
mkdir -p $logdir

eqtl_conditions="Across_full Across_GWAS_colocalized Within_full Within_GWAS_colocalized"

# --- Resource definitions: file, subtract_file (or "none"), directions (space-separated), mem_gb ---
# Format: resource_key|annotation_file|subtract_file|directions|mem_gb



resources=(
  "collectri_strict|collectri_strict.tsv|none|forward reverse|8"
  "postar3_strict|postar3_strict.tsv|none|forward reverse|16"
  "liana_lenient|liana_lenient.tsv|none|forward reverse|10"
  "merops_strict|merops_strict.tsv|none|forward reverse|8"
  "string_experimental_strict|string_experimental_strict.tsv|none|undirected|8"
  "string_experimental_lenient_band|string_experimental_lenient.tsv|string_experimental_strict.tsv|undirected|8"
  "scenic_regulon_activator_same_dataset|scenic_regulon_activator_same_dataset.tsv|none|forward reverse|16"
  "scenic_regulon_repressor_same_dataset|scenic_regulon_repressor_same_dataset.tsv|none|forward reverse|16"
)








n_submitted=0

for resource_def in "${resources[@]}"; do
  IFS='|' read -r resource_key annot_file subtract_file directions mem <<< "$resource_def"
  mem_mb=$((mem * 1000))

  subtract_arg=""
  if [ "$subtract_file" != "none" ]; then
    subtract_arg="--subtract_file $annot_dir/${subtract_file}"
  fi

  for direction in $directions; do
    for eqtl_condition in $eqtl_conditions; do
      bsub -G team151 \
        -o $logdir/${resource_key}_${direction}_${eqtl_condition}.%J.out \
        -e $logdir/${resource_key}_${direction}_${eqtl_condition}.%J.err \
        -M $mem_mb \
        -R"select[mem>$mem_mb] rusage[mem=$mem_mb] span[hosts=1]" \
        -n1 \
        -q $queue -- \
        "module load ISG/conda && \
        conda activate $condaenv && \
        python $script --annotation_file $annot_dir/${annot_file} $subtract_arg \
          --direction $direction --eqtl_condition $eqtl_condition --mem_gb $mem"
      echo "Submitted ${resource_key} / ${direction} / ${eqtl_condition} (${mem}GB)"
      n_submitted=$((n_submitted + 1))
    done
  done
done

echo ""
echo "Total jobs submitted: $n_submitted"
echo "Expected: (5 directed resources x 2 directions + 2 undirected resources x 1 direction) x 4 eqtl_conditions = (10+2)*4 = 48"
