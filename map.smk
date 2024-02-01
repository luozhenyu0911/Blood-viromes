
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

# remove LCR-hs38_rmsk seq
rule rm_repeats:
    input:
        bam = "01_align2hg38/{id}.unmapped_hg38.bam",
        region = config['params']['repeats']
    output:
        in_bam = "01_align2hg38/{id}.is_repeats.bam",
        out_bam = "01_align2hg38/{id}.rm_repeats.bam"
    threads:
        config['threads']['bwa']
    shell:
        "samtools view -@ {threads} {input.bam} -bh -o {output.in_bam} -U {output.out_bam} -L {input.region}"

# get fq1 after rm hg38 and repeats seq.
rule unmapped_seqid:
    input:
        out_bam = "01_align2hg38/{id}.rm_repeats.bam"
    output:
        seq_id = "01_align2hg38/{id}.rm_repeats.id"
    threads:
        config['threads']['bwa']
    shell:
        "samtools view -@ {threads} {input.out_bam} | cut -f1|sort -u > {output.seq_id}"

# get fq1 after rm hg38 and repeats seq.
rule unmapped_fq1:
    input:
        seq_id = "01_align2hg38/{id}.rm_repeats.id",
        rawfq1 = "data/{id}.R1.fq.gz"
    output:
        R1_id = "01_align2hg38/{id}.rm_repeats.R1.id",
        rmrepeats_fq1 = "01_align2hg38/{id}.rm_repeats.R1.fq.gz"
    threads:
        config['threads']['bwa']
    shell:
        "awk '{{print $1\"/1\"}}' {input.seq_id} > {output.R1_id} && "
        "seqtk subseq {input.rawfq1} {output.R1_id} |gzip > {output.rmrepeats_fq1}"

# get fq2 after rm hg38 and repeats seq.
rule unmapped_fq2:
    input:
        seq_id = "01_align2hg38/{id}.rm_repeats.id",
        rawfq2 = "data/{id}.R2.fq.gz"
    output:
        R2_id = "01_align2hg38/{id}.rm_repeats.R2.id",
        rmrepeats_fq2 = "01_align2hg38/{id}.rm_repeats.R2.fq.gz"
    threads:
        config['threads']['bwa']
    shell:
        "awk '{{print $1\"/2\"}}' {input.seq_id} > {output.R2_id} && "
        "seqtk subseq {input.rawfq2} {output.R2_id} |gzip > {output.rmrepeats_fq2} "

# remove Ref T2T seq
rule Align2T2T:
    input:
        ref_T2T = config['params']['ref_T2T'],
        rmrepeats_fq1 = "01_align2hg38/{}.rm_repeats.R1.fq.gz".format(config['samples']['id']),
        rmrepeats_fq2 = "01_align2hg38/{}.rm_repeats.R2.fq.gz".format(config['samples']['id'])
    output:
        T2T_bam = "02_align2T2T/{id}.T2T.bam"
    threads:
        config['threads']['bwa']
    shell:
        "bwa mem -C -t {threads} {input.ref_T2T} {input.rmrepeats_fq1} {input.rmrepeats_fq2} | samtools sort -o {output.T2T_bam} -@ 24 -O bam -"

# remove human seq. against T2T
rule rm_T2T:
    input:
        T2T_bam = "02_align2T2T/{id}.T2T.bam"
    output:
        unmapped_bam = "02_align2T2T/{id}.unmapped_T2T.bam",
        seq_id = "02_align2T2T/{id}.rm_T2T.id"
    threads:
        config['threads']['bwa']
    shell:
        "samtools view -@ {threads} -bh -f 4 {input.T2T_bam} > {output.unmapped_bam} && "
        "samtools view -@ {threads} {output.unmapped_bam} | cut -f1|sort -u > {output.seq_id}"

# get fq1 after rm T2T.
rule T2T_unmapped_fq1:
    input:
        seq_id = "02_align2T2T/{id}.rm_T2T.id",
        fq1 = "01_align2hg38/{}.rm_repeats.R1.fq.gz".format(config['samples']['id'])
    output:
        R1_id = "02_align2T2T/{id}.rm_T2T.R1.id",
        rm_T2T_fq1 = "02_align2T2T/{id}.rm_T2T.R1.fq.gz"
    threads:
        config['threads']['bwa']
    shell:
        "awk '{{print $1\"/1\"}}' {input.seq_id} > {output.R1_id} && "
        "seqtk subseq {input.fq1} {output.R1_id} |gzip > {output.rm_T2T_fq1}"

# get fq2 after rm T2T.
rule T2T_unmapped_fq2:
    input:
        seq_id = "02_align2T2T/{id}.rm_T2T.id",
        fq2 = "01_align2hg38/{}.rm_repeats.R2.fq.gz".format(config['samples']['id'])
    output:
        R2_id = "02_align2T2T/{id}.rm_T2T.R2.id",
        rm_T2T_fq2 = "02_align2T2T/{id}.rm_T2T.R2.fq.gz"
    threads:
        config['threads']['bwa']
    shell:
        "awk '{{print $1\"/2\"}}' {input.seq_id} > {output.R2_id} && "
        "seqtk subseq {input.fq2} {output.R2_id} |gzip > {output.rm_T2T_fq2}"

# rm pangenome + decoy seq.

rule Align2pangenome:
    input:
        ref_pg = config['params']['ref_pangenome'],
        rm_T2T_fq1 = "02_align2T2T/{id}.rm_T2T.R1.fq.gz",
        rm_T2T_fq2 = "02_align2T2T/{id}.rm_T2T.R2.fq.gz"
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
        fq1 = "02_align2T2T/{id}.rm_T2T.R1.fq.gz"
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
        fq2 = "02_align2T2T/{id}.rm_T2T.R2.fq.gz"
    output:
        R2_id = "03_align2pg/{id}.rm_pg.R2.id",
        rm_pg_fq2 = "03_align2pg/{id}.rm_pg.R2.fq.gz"
    threads:
        config['threads']['bwa']
    shell:
        "awk '{{print $1\"/2\"}}' {input.seq_id} > {output.R2_id} && "
        "seqtk subseq {input.fq2} {output.R2_id} |gzip > {output.rm_pg_fq2}"
		
