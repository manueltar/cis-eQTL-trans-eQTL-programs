# Workflow

## Preprocessing resources

Download Homo_sapiens.GRCh38.110.gtf

`string_interactions.tsv` — STRING network export (node1/node2 as Ensembl gene IDs, with `combined_score`)

wget https://zenodo.org/records/8192729/files/CollecTRI_regulons.csv

pip download liana==1.8.0 --no-deps -d .
unzip liana-1.8.0-py3-none-any.whl liana/resource/omni_resource.csv -d .
cp liana/resource/omni_resource.csv .

Download from URL: http://postar.ncrnalab.org

$ bash ~/Scripts/Wraper_scripts/403_POSTAR3_chr_array.sh normal

plus 

for chrom in 3 4 8 10 14 15 16 18 21 X; do
  bsub -G team151 -o /nfs/team151/mt19/POSTAR3/logs/chr${chrom}_retry.%J.out -e /nfs/team151/mt19/POSTAR3/logs/chr${chrom}_retry.%J.err \
    -M 24000 -R"select[mem>24000] rusage[mem=24000] span[hosts=1]" -n1 -q normal -- \
    "module load ISG/conda && conda activate /nfs/users/nfs_m/mt19/cardinal_analysis/ht/conda_envs/pyranges_env && python /nfs/users/nfs_m/mt19/Scripts/Python_scripts/POSTAR3_preprocessing_v2.py --chrom ${chrom} --nb_cpu 1 --mem_gb 24"
  echo "Resubmitted chr${chrom} at 24GB"
done

wget https://ftp.ebi.ac.uk/pub/databases/merops/current_release/meropsweb121.tar.gz

module load ISG/conda
conda create -y -p /nfs/users/nfs_m/mt19/cardinal_analysis/ht/conda_envs/merops_db -c conda-forge --override-channels mariadb
conda activate /nfs/users/nfs_m/mt19/cardinal_analysis/ht/conda_envs/merops_db
$ bash ~/Scripts/Wraper_scripts/400_Restarted_HT_outage_MEROPS_db_creation.sh normal
$ bash ~/Scripts/Wraper_scripts/401_Restarted_HT_outage_MEROPS_db_creation_extraction.sh normal

----------> Jupyter_notebooks/Merops_preprocessing.ipynb


----------> Jupyter_notebooks/Preprocessing_datasets_of_pairs.ipynb

## Classifying the cis-eQTL/trans-eQTL programs

----------> Jupyter_notebooks/Overhaul_classification_factors_v2.ipynb

## Enrichments

### edge-based

$ bash ~/Scripts/Wraper_scripts/405_enrichment_array_v2.sh normal

----------> Jupyter_notebooks/Representation_ERs_v4.ipynb
----------> Jupyter_notebooks/Factor_ER_v2.ipynb

## Pattern finding

----------> Jupyter_notebooks/Overhaul_finding_patterns_v5.ipynb

