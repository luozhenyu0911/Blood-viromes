#!/usr/bin/env bash

conda activate /home-02/zhenyuluo/anaconda3/envs/RNA-seq

smk=./run.all.smk
echo ""
echo "start = $(date)"
echo "$(hostname)"
snakemake -j 80 -pk -s ${smk} 2> snakemake.err.txt
echo "end = $(date)"
# echo "last status $? $(hostname)" |mail -s "job done: $(pwd)" luozhenyu@genomics.cn 
