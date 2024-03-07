def get_all_fq_fa(wildcards):
    all_fq_fa = []
    # all_fq_fa.append("data/{}.R1.fq.gz".format(config['samples']['id']))
    all_fq_fa.append("01_align2hg38/{}.rm_repeats.R2.fq.gz".format(config['samples']['id']))
    all_fq_fa.append("02_align2T2T/{}.rm_T2T.R1.fq.gz".format(config['samples']['id']))
    all_fq_fa.append("03_align2pg/{}.rm_pg.R1.fq.gz".format(config['samples']['id']))
    all_fq_fa.append('04_annotation/{}.tid10239.R1.fq'.format(config['samples']['id']))
    all_fq_fa.append('04_annotation/{}.unclass.R1.fq'.format(config['samples']['id']))
    all_fq_fa.append('04_annotation/{}_viral_unclass.R1.fq'.format(config['samples']['id']))
    return all_fq_fa
    
rule stats_fq_fa:
    input:
        get_all_fq_fa
    output:
        "metrics/{id}_seq_count.txt"
    run:
        shell("seqkit stats {input} > {output}")

rule get_singletons_id_hg38:
    input:
        "01_align2hg38/{id}.unmapped_hg38.bam"
    output:
        "metrics/{id}_singletons_hg38.id"
    threads:
        config['threads']['bwa']
    shell:
        "samtools view -@ {threads} {input} -f 4 -F 8 -o {output}"
        
rule get_singletons_id_T2T:
    input:
        "02_align2T2T/{id}.T2T.bam"
    output:
        "metrics/{id}_singletons_T2T.id"
    threads:
        config['threads']['bwa']
    shell:
        "samtools view -@ {threads} {input} -f 4 -F 8 -o {output}"

rule get_singletons_id_pg:
    input:
        "03_align2pg/{id}.pg.bam"
    output:
        "metrics/{id}_singletons_pg.id"
    threads:
        config['threads']['bwa']
    shell:
        "samtools view -@ {threads} {input} -f 4 -F 8 -o {output}"

def get_all_singletons_id(wildcards):
    all_singletons_id= []
    all_singletons_id.append("metrics/{}_singletons_hg38.id".format(config['samples']['id']))
    all_singletons_id.append("metrics/{}_singletons_T2T.id".format(config['samples']['id']))
    all_singletons_id.append("metrics/{}_singletons_pg.id".format(config['samples']['id']))
    return all_singletons_id

rule get_all_singletons:
    input:
        get_all_singletons_id
    output:
        "metrics/{id}_all_singleton.id"
    run:
        shell("cat {input} |cut -f1 | sort -u > {output}")

rule get_viral_seq_id:
    input:
        "05_viral_anno/{id}.viral.k2"
    output:
        "metrics/{id}.viral.seq.id"
    shell:
        "cat {input} |awk '$1~/C/' |cut -f2 |sort -u > {output}"


rule get_overlap_viral_singleton:
    input:
        singleton = "metrics/{id}_all_singleton.id",
        viral_id = "metrics/{id}.viral.seq.id"
    output:
        "metrics/{id}.overlap_viral_singleton.id"
    shell:
        "comm -1 -2 {input.singleton} {input.viral_id} > {output}"





