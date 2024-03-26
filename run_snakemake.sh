#!/bin/bash
# yhbatch -N 1 -n 24 -p rhenv 

smk=/BIGDATA2/gzfezx_shhli_2/USER/luozhenyu/script/blood_viromes3/run.all.smk
echo ""
echo "start = $(date)"
echo "$(hostname)"
snakemake --keep-going -j 24 -pk -s ${smk} 2> snakemake.err.txt
echo "end = $(date)"
# echo "last status $? $(hostname)" |mail -s "job done: $(pwd)"  
