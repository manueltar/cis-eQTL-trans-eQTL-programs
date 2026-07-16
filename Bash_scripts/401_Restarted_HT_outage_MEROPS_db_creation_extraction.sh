#!/bin/bash
queue=$1

condaenv=/nfs/users/nfs_m/mt19/cardinal_analysis/ht/conda_envs/mysql_db
dumpfile=/nfs/team151/mt19/merops_db/meropsweb121.tar.gz
dbhost=node-13-09
dbport=3307

bsub -G team151 -o /nfs/team151/mt19/merops_db/load_tables.%J.out -e /nfs/team151/mt19/merops_db/load_tables.%J.err -M 8000 -R"select[mem>8000] rusage[mem=8000] span[hosts=1]" -n1 -q $queue -- \
"module load ISG/conda && \
conda activate $condaenv && \
mysql -u root -h $dbhost -P $dbport -e 'CREATE DATABASE IF NOT EXISTS merops;' && \
for t in code domain sequence organism gene_name cleavage substrate; do \
  echo Loading \$t...; \
  tar -xzOf $dumpfile \${t}.sql | mysql -u root -h $dbhost -P $dbport merops; \
done"
