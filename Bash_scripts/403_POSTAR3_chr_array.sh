#!/bin/bash
# 403_POSTAR3_chr_array.sh
queue=$1

condaenv=/nfs/users/nfs_m/mt19/cardinal_analysis/ht/conda_envs/pyranges_env
script=/nfs/users/nfs_m/mt19/Scripts/Python_scripts/POSTAR3_preprocessing_v2.py
outdir=/nfs/team151/mt19/POSTAR3/per_chrom
logdir=/nfs/team151/mt19/POSTAR3/logs

mkdir -p $outdir
mkdir -p $logdir

declare -A mem_gb=(
  [1]=16
  [2]=12  [17]=12
  [3]=10  [5]=10  [6]=10  [7]=10  [11]=10  [12]=10  [19]=10
  [4]=6   [8]=6   [9]=6   [10]=6  [13]=6   [14]=6  [15]=6  [16]=6  [20]=6  [22]=6  [X]=6
  [18]=4  [21]=4  [MT]=4  [Y]=4
)

for chrom in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT; do
  mem=${mem_gb[$chrom]}
  mem_mb=$((mem * 1000))

  bsub -G team151 \
    -o $logdir/chr${chrom}.%J.out \
    -e $logdir/chr${chrom}.%J.err \
    -M $mem_mb \
    -R"select[mem>$mem_mb] rusage[mem=$mem_mb] span[hosts=1]" \
    -n1 \
    -q $queue -- \
    "module load ISG/conda && \
    conda activate $condaenv && \
    python $script --chrom $chrom --nb_cpu 1 --mem_gb $mem"

  echo "Submitted chr${chrom} (${mem}GB)"
done
