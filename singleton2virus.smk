# remove LCR-hs38_rmsk seq
rule get_singleton_seq:
    input:
        bam = "01_align2hg38/{}.unmapped_hg38.bam".format(config['samples']['id'])
    output:
        "07_single_virus/singleton.fq.gz"
    threads:
        config['threads']['bwa']
    shell:
        "samtools fastq -@ {threads} -f 4 -F 8 {input.bam} | gzip > {output}"
        
# remove Ref T2T seq
rule singleton_map2virus:
    input:
        ref_virus = config['params']['ref_viral'],
        fq = "07_single_virus/singleton.fq.gz"
    output:
        virus_bam = "07_single_virus/{id}.singleton_map2virus.bam"
    threads:
        config['threads']['bwa']
    shell:
        "bwa mem -C -t {threads} {input.ref_virus} {input.fq} | samtools sort -o {output.virus_bam} -@ 24 -O bam -"

# remove human seq. against T2T
rule get_sin_map2virus_id:
    input:
        "07_single_virus/{id}.singleton_map2virus.bam"
    output:
        mapped_bam = "07_single_virus/{id}.singleton_map2virus.mapped.bam",
        seq_id = "07_single_virus/{id}.singleton_map2virus.mapped.id"
    threads:
        config['threads']['bwa']
    shell:
        "samtools view -@ {threads} -bh -F 4 {input} > {output.mapped_bam} && "
        "samtools view -@ {threads} {output.mapped_bam} | cut -f1|sort -u > {output.seq_id}"