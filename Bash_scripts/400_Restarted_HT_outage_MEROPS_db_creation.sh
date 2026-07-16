#!/bin/bash
queue=$1

datadir=/nfs/team151/mt19/merops_db/data
socket=/nfs/team151/mt19/merops_db/mysql.sock
condaenv=/nfs/users/nfs_m/mt19/cardinal_analysis/ht/conda_envs/mysql_db
port=3307

bsub -G team151 -o /nfs/team151/mt19/merops_db/merops_db.%J.out -e /nfs/team151/mt19/merops_db/merops_db.%J.err -M 16000 -R"select[mem>16000] rusage[mem=16000] span[hosts=1]" -n2 -q $queue -- \
"module load ISG/conda && \
conda activate $condaenv && \
mysqld --datadir=$datadir --socket=$socket --port=$port --bind-address=0.0.0.0 --mysqlx=0 --innodb-buffer-pool-size=8G"
