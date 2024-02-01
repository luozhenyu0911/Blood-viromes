#!/bin/bash
# yhbatch -N 1 -n 1 -p rhenv 

smk=./run.all.smk
echo ""
echo "start = $(date)"
echo "$(hostname)"
snakemake -j 24 -pk -s ${smk} 2> snakemake.err.txt
echo "end = $(date)"
# echo "last status $? $(hostname)" |mail -s "job done: $(pwd)"  
