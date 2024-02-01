def get_all_fq_fa(wildcards):
    all_fq_fa = []
    all_fq_fa.append("data/{}.R1.fq.gz".format(config['samples']['id']))
    all_fq_fa.append("01_align2hg38/{}.rm_repeats.R2.fq.gz".format(config['samples']['id']))
    all_fq_fa.append("02_align2T2T/{}.rm_T2T.R1.fq.gz".format(config['samples']['id']))
    all_fq_fa.append("03_align2pg/{}.rm_pg.R1.fq.gz".format(config['samples']['id']))
    all_fq_fa.append('04_annotation/{}.tid10239.R1.fa'.format(config['samples']['id']))
    all_fq_fa.append('04_annotation/{}.unclass.R1.fa'.format(config['samples']['id']))
    all_fq_fa.append('04_annotation/{}_viral_unclass.R1.fa'.format(config['samples']['id']))
    return all_fq_fa
    
rule stats_fq_fa:
    input:
        get_all_fq_fa
    output:
        "metrics/{id}_seq_count.txt"
    run:
        shell("seqkit stats {input} > {output}")
    

