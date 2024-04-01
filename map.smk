
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

# get unmapped seq
rule rm_hg38seqs:
    input:
        "data/{id}.cram"
    output:
        # "Align2hg38/{}.unmapped_hg38.bam".format(config['samples']['id'])
        "Align2hg38/{id}.unmapped_hg38.bam"
    threads:
        config['threads']['bwa']
    shell:
        "samtools view -@ {threads} -bh -f 4 {input} > {output}"

# get unmapped_seqid
rule unmapped_seqid:
    input:
        out_bam = "Align2hg38/{id}.unmapped_hg38.bam"
    output:
        seq_id = "Align2hg38/{id}.unmapped_hg38.id"
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
        "Align2hg38/{id}.EBV.seqid"
    threads:
        config['threads']['bwa']
    shell:
        "samtools view {input.cram} chrEBV|cut -f1 > {output}"

# merge EBV and unmapped_seqid
rule merge_EBV_unmapped_seqid:
    input:
        unmappedid="Align2hg38/{id}.unmapped_hg38.id",
        EBVid = "Align2hg38/{id}.EBV.seqid"
    output:
        mergeid = "Align2hg38/{id}.EBV_unmapped.id"
    shell:
        "cat {input.unmappedid} {input.EBVid} |sort -u > {output.mergeid}"

# get fq1 of unmapped seqid.
rule unmapped_fq1:
    input:
        seq_id = "Align2hg38/{id}.EBV_unmapped.id",
        rawfq1 = "data/{id}.R1.fq.gz"
    output:
        R1_id = "Align2hg38/{id}.EBV_unmapped.R1.id",
        EBV_unmapped_fq1 = "Align2hg38/{id}.EBV_unmapped.R1.fq.gz"
    threads:
        config['threads']['bwa']
    shell:
        "awk '{{print $1\"/1\"}}' {input.seq_id} > {output.R1_id} && "
        "seqtk subseq {input.rawfq1} {output.R1_id} |gzip > {output.EBV_unmapped_fq1}"

# get fq2  of unmapped seqid.
rule unmapped_fq2:
    input:
        seq_id = "Align2hg38/{id}.EBV_unmapped.id",
        rawfq2 = "data/{id}.R2.fq.gz"
    output:
        R2_id = "Align2hg38/{id}.EBV_unmapped.R2.id",
        EBV_unmapped_fq2 = "Align2hg38/{id}.EBV_unmapped.R2.fq.gz"
    threads:
        config['threads']['bwa']
    shell:
        "awk '{{print $1\"/2\"}}' {input.seq_id} > {output.R2_id} && "
        "seqtk subseq {input.rawfq2} {output.R2_id} |gzip > {output.EBV_unmapped_fq2} "

# rm pangenome + decoy seq.
rule Align2pangenome:
    input:
        ref_pg = config['params']['ref_pangenome'],
        EBV_unmapped_fq1 = "Align2hg38/{id}.EBV_unmapped.R1.fq.gz",
        EBV_unmapped_fq2 = "Align2hg38/{id}.EBV_unmapped.R2.fq.gz"
    output:
        pg_bam = "Align2pg/{id}.pg.bam"
    threads:
        config['threads']['bwa']
    shell:
        "bwa mem -C -t {threads} {input.ref_pg} {input.EBV_unmapped_fq1} {input.EBV_unmapped_fq2} | samtools sort -o {output.pg_bam} -@ {threads} -O bam -"

# remove human seq. against T2T
rule rm_pg:
    input:
        pg_bam = "Align2pg/{id}.pg.bam"
    output:
        unmapped_bam = "Align2pg/{id}.unmapped_pg.bam",
        seq_id = "Align2pg/{id}.rm_pg.id"
    threads:
        config['threads']['bwa']
    shell:
        "samtools view -@ {threads} -bh -f 4 {input.pg_bam} > {output.unmapped_bam} && "
        "samtools view -@ {threads} {output.unmapped_bam} | cut -f1|sort -u > {output.seq_id}"

# get fq1 after rm pangenome.
rule pg_unmapped_fq1:
    input:
        seq_id = "Align2pg/{id}.rm_pg.id",
        fq1 = "Align2hg38/{id}.EBV_unmapped.R1.fq.gz"
    output:
        R1_id = "Align2pg/{id}.rm_pg.R1.id",
        rm_pg_fq1 = "Align2pg/{id}.rm_pg.R1.fq.gz"
    threads:
        config['threads']['bwa']
    shell:
        "awk '{{print $1\"/1\"}}' {input.seq_id} > {output.R1_id} && "
        "seqtk subseq {input.fq1} {output.R1_id} |gzip > {output.rm_pg_fq1}"

# get fq2 after rm pangenome.
rule pg_unmapped_fq2:
    input:
        seq_id = "Align2pg/{id}.rm_pg.id",
        fq2 = "Align2hg38/{id}.EBV_unmapped.R2.fq.gz"
    output:
        R2_id = "Align2pg/{id}.rm_pg.R2.id",
        rm_pg_fq2 = "Align2pg/{id}.rm_pg.R2.fq.gz"
    threads:
        config['threads']['bwa']
    shell:
        "awk '{{print $1\"/2\"}}' {input.seq_id} > {output.R2_id} && "
        "seqtk subseq {input.fq2} {output.R2_id} |gzip > {output.rm_pg_fq2}"
		
