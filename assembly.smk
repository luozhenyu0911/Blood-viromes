# assembly for viral and unclassed seqs.
import os
import shutil

if os.path.exists('06_assembly/01_metaspades'):
    shutil.rmtree('06_assembly/01_metaspades', ignore_errors=True)
if os.path.exists('06_assembly/02_spades'):
    shutil.rmtree('06_assembly/02_spades', ignore_errors=True)
if os.path.exists('06_assembly/03_megahit'):
    shutil.rmtree('06_assembly/03_megahit', ignore_errors=True)

if os.path.exists('06_assembly/01_metaspades_quast'):
    shutil.rmtree('06_assembly/01_metaspades_quast', ignore_errors=True)
if os.path.exists('06_assembly/02_spades_quast'):
    shutil.rmtree('06_assembly/02_spades_quast', ignore_errors=True)
if os.path.exists('06_assembly/03_megahit_quast'):
    shutil.rmtree('06_assembly/03_megahit_quast', ignore_errors=True)
    
rule assembly_metaspades:
    input:
        fq1 = '04_annotation/{}_viral_unclass.R1.fq'.format(config['samples']['id']),
        fq2 = '04_annotation/{}_viral_unclass.R2.fq'.format(config['samples']['id'])
    output:
        directory('06_assembly/01_metaspades')
    threads:
        int(config['threads']['bwa'])
    shell:
        "metaspades.py -t {threads} -1 {input.fq1} -2 {input.fq2} -o {output}"

rule quast_metaspades:
    input:
        fa = '06_assembly/01_metaspades',
        viral_ref = "/BIGDATA2/gzfezx_shhli_2/USER/luozhenyu/database/viral_refseq/viral.1.1.genomic.fna"
    output:
        directory('06_assembly/01_metaspades_quast')
    threads:
        int(config['threads']['bwa'])
    shell:
        "quast.py -t {threads} -r {input.viral_ref} -o {output} {input.fa}/contigs.fasta"


rule assembly_spades:
    input:
        fq1 = '04_annotation/{}_viral_unclass.R1.fq'.format(config['samples']['id']),
        fq2 = '04_annotation/{}_viral_unclass.R2.fq'.format(config['samples']['id'])
    output:
        directory('06_assembly/02_spades')
    threads:
        int(config['threads']['bwa'])
    shell:
        "spades.py -t {threads} -1 {input.fq1} -2 {input.fq2} -o {output} --metaviral"

rule quast_spades:
    input:
        fa = '06_assembly/02_spades',
        viral_ref = "/BIGDATA2/gzfezx_shhli_2/USER/luozhenyu/database/viral_refseq/viral.1.1.genomic.fna"
    output:
        directory('06_assembly/02_spades_quast')
    threads:
        int(config['threads']['bwa'])
    shell:
        "quast.py -t {threads} -r {input.viral_ref} -o {output} {input.fa}/contigs.fasta"


rule assembly_megahit:
    input:
        fq1 = '04_annotation/{}_viral_unclass.R1.fq'.format(config['samples']['id']),
        fq2 = '04_annotation/{}_viral_unclass.R2.fq'.format(config['samples']['id'])
    output:
        directory('06_assembly/03_megahit')
    threads:
        int(config['threads']['bwa'])
    shell:
        "megahit -t {threads} -1 {input.fq1} -2 {input.fq2} -o {output}"

rule quast_megahit:
    input:
        fa = '06_assembly/03_megahit',
        viral_ref = "/BIGDATA2/gzfezx_shhli_2/USER/luozhenyu/database/viral_refseq/viral.1.1.genomic.fna"
    output:
        directory('06_assembly/03_megahit_quast')
    threads:
        int(config['threads']['bwa'])
    shell:
        "quast.py -t {threads} -r {input.viral_ref} -o {output} {input.fa}/final.contigs.fa"






















