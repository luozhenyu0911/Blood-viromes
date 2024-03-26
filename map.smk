
from pathlib import Path

def get_fastqs_one(wildcards):
    fq_path = Path(config['samples']['fq_path'])
    name = config['samples']['id']
    fq_files = []
    for fp in fq_path.iterdir():
        # add the fastq file to the files list if it matches
        if fp.is_file and str(name) in str(fp) and str(fp).endswith("1.fq.gz"):
            print(fp)
            fq_files.append(str(fp))
    return fq_files    


def get_fastqs_two(wildcards):
    fq_path = Path(config['samples']['fq_path'])
    name = config['samples']['id']
    fq_files = []
    for fp in fq_path.iterdir():
        # add the fastq file to the files list if it matches
        if fp.is_file and str(name) in str(fp) and str(fp).endswith("2.fq.gz"):
            print(fp)
            fq_files.append(str(fp))
    return fq_files   

def get_cram(wildcards):
    fq_path = Path(config['samples']['fq_path'])
    name = config['samples']['id']
    cram_files = []
    for fp in fq_path.iterdir():
        # add the fastq file to the files list if it matches
        if fp.is_file and str(name).split('_')[0] in str(fp) and str(fp).endswith(".cram"):
            print(fp)
            cram_files.append(str(fp))
    return cram_files 

def get_cram_index(wildcards):
    fq_path = Path(config['samples']['fq_path'])
    name = config['samples']['id']
    cram_files = []
    for fp in fq_path.iterdir():
        # add the fastq file to the files list if it matches 
        if fp.is_file and str(name).split('_')[0] in str(fp) and str(fp).endswith(".crai"):
            print(fp)
            cram_files.append(str(fp))
    return cram_files

rule link_read_one:
    input:
        get_fastqs_one
    output:
        "data/{id}.R1.fq.gz"
    run:
        shell("ln -s {input} {output}")

rule link_read_two:
    input:
        get_fastqs_two
    output:
        "data/{id}.R2.fq.gz"
    run:
        shell("ln -s {input} {output}")

rule link_cram:
    input:
        get_cram
    output:
        "data/{id}.cram"
    run:
        shell("ln -s {input} {output}")

rule link_cram_index:
    input:
        get_cram_index
    output:
        "data/{id}.cram.crai"
    run:
        shell("ln -s {input} {output}")

# remove human seq. against hg38
rule rm_hg38seqs:
    input:
        "data/{id}.cram"
    output:
        # "01_align2hg38/{}.unmapped_hg38.bam".format(config['samples']['id'])
        "01_align2hg38/{id}.unmapped_hg38.bam"
    threads:
        config['threads']['bwa']
    shell:
        "samtools view -@ {threads} -bh -f 4 {input} > {output}"

# get fq1 after rm hg38 and repeats seq.
rule unmapped_seqid:
    input:
        out_bam = "01_align2hg38/{id}.unmapped_hg38.bam"
    output:
        seq_id = "01_align2hg38/{id}.unmapped.id"
    threads:
        config['threads']['bwa']
    shell:
        "samtools view -@ {threads} {input.out_bam} | cut -f1|sort -u > {output.seq_id}"

# get EBV id
rule mapped_EBV_seqid:
    input:
        cram = "data/{id}.cram",
        index = "data/{id}.cram.crai"
    output:
        "01_align2hg38/{id}.EBV.seqid"
    threads:
        config['threads']['bwa']
    shell:
        "samtools view {input.cram} chrEBV|cut -f1 > {output}"

# merge EBV and unmapped.id
rule merge_EBV_unmapped_id:
    input:
        repeatid="01_align2hg38/{id}.unmapped.id",
        EBVid = "01_align2hg38/{id}.EBV.seqid"
    output:
        mergeid = "01_align2hg38/{id}.unmapped_EBV.id"
    shell:
        "cat {input.repeatid} {input.EBVid} |sort -u > {output.mergeid}"

# get fq1 after rm hg38 seq.
rule unmapped_fq1:
    input:
        seq_id = "01_align2hg38/{id}.unmapped_EBV.id",
        rawfq1 = "data/{id}.R1.fq.gz"
    output:
        R1_id = "01_align2hg38/{id}.unmapped.R1.id",
        unmapped_fq1 = "01_align2hg38/{id}.unmapped.R1.fq.gz"
    threads:
        config['threads']['bwa']
    shell:
        "awk '{{print $1\"/1\"}}' {input.seq_id} > {output.R1_id} && "
        "seqtk subseq {input.rawfq1} {output.R1_id} |gzip > {output.unmapped_fq1}"

# get fq2 after rm hg38 and repeats seq.
rule unmapped_fq2:
    input:
        seq_id = "01_align2hg38/{id}.unmapped_EBV.id",
        rawfq2 = "data/{id}.R2.fq.gz"
    output:
        R2_id = "01_align2hg38/{id}.unmapped.R2.id",
        unmapped_fq2 = "01_align2hg38/{id}.unmapped.R2.fq.gz"
    threads:
        config['threads']['bwa']
    shell:
        "awk '{{print $1\"/2\"}}' {input.seq_id} > {output.R2_id} && "
        "seqtk subseq {input.rawfq2} {output.R2_id} |gzip > {output.unmapped_fq2} "

# rm pangenome + decoy seq.
rule Align2pangenome:
    input:
        ref_pg = config['params']['ref_pangenome'],
        rm_T2T_fq1 = "01_align2hg38/{id}.unmapped.R1.fq.gz",
        rm_T2T_fq2 = "01_align2hg38/{id}.unmapped.R2.fq.gz"
    output:
        pg_bam = "03_align2pg/{id}.pg.bam"
    threads:
        config['threads']['bwa']
    shell:
        "bwa mem -C -t {threads} {input.ref_pg} {input.rm_T2T_fq1} {input.rm_T2T_fq2} | samtools sort -o {output.pg_bam} -@ {threads} -O bam -"

# remove human seq. against T2T
rule rm_pg:
    input:
        pg_bam = "03_align2pg/{id}.pg.bam"
    output:
        unmapped_bam = "03_align2pg/{id}.unmapped_pg.bam",
        seq_id = "03_align2pg/{id}.rm_pg.id"
    threads:
        config['threads']['bwa']
    shell:
        "samtools view -@ {threads} -bh -f 4 {input.pg_bam} > {output.unmapped_bam} && "
        "samtools view -@ {threads} {output.unmapped_bam} | cut -f1|sort -u > {output.seq_id}"

# get fq1 after rm T2T.
rule pg_unmapped_fq1:
    input:
        seq_id = "03_align2pg/{id}.rm_pg.id",
        fq1 = "01_align2hg38/{id}.unmapped.R1.fq.gz"
    output:
        R1_id = "03_align2pg/{id}.rm_pg.R1.id",
        rm_pg_fq1 = "03_align2pg/{id}.rm_pg.R1.fq.gz"
    threads:
        config['threads']['bwa']
    shell:
        "awk '{{print $1\"/1\"}}' {input.seq_id} > {output.R1_id} && "
        "seqtk subseq {input.fq1} {output.R1_id} |gzip > {output.rm_pg_fq1}"

# get fq2 after rm T2T.
rule pg_unmapped_fq2:
    input:
        seq_id = "03_align2pg/{id}.rm_pg.id",
        fq2 = "01_align2hg38/{id}.unmapped.R2.fq.gz"
    output:
        R2_id = "03_align2pg/{id}.rm_pg.R2.id",
        rm_pg_fq2 = "03_align2pg/{id}.rm_pg.R2.fq.gz"
    threads:
        config['threads']['bwa']
    shell:
        "awk '{{print $1\"/2\"}}' {input.seq_id} > {output.R2_id} && "
        "seqtk subseq {input.fq2} {output.R2_id} |gzip > {output.rm_pg_fq2}"
		
