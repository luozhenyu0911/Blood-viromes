rule cls_virus_kraken2:
    input:
        db = config['params']['k2_viral'],
        fq1 = '04_annotation/{id}_viral_unclass.R1.fq',
        fq2 = '04_annotation/{id}_viral_unclass.R2.fq'
    output:
        # unclassified = "04_annotation/{id}.Unclass#.fa",
        # classified = "04_annotation/{id}.Class#.fa",
        # unclassified = expand("04_annotation/{id}.Unclass_{i}.fa", i=range(1,3)),
        # classified = expand("04_annotation/{id}.Class_{i}.fa", i=range(1,3)),
        out_k2 = "05_viral_anno/{id}.viral.k2",
        report = "05_viral_anno/{id}.viral.k2_report.txt"
    threads:
        config['threads']['bwa']
    shell:
        "kraken2 --db {input.db} --paired {input.fq1} {input.fq2} "
        "--threads {threads} "
        "--output {output.out_k2} --report {output.report}"
        # " --unclassified-out 04_annotation/{id}.Unclass#.fa --classified-out 04_annotation/{id}.Class#.fa"

rule virus_bracken:
    input:
        db = config['params']['k2_viral'],
        report = "05_viral_anno/{id}.viral.k2_report.txt"
    output:
        bracken = "05_viral_anno/{id}.viral.bracken",
        breport = "05_viral_anno/{id}.viral.breport"
    threads:
        config['threads']['bwa']
    shell:
        "bracken -d {input.db} -i {input.report} -o {output.bracken} -w {output.breport} -t 0"
        
rule virus_diversity:
    input:
        bracken = "05_viral_anno/{id}.viral.bracken"
    output:
        "05_viral_anno/{id}.viral.beta_diversity.stats"
    shell:
        "alpha_diversity.py -f {input.bracken} -a BP >> {output} && echo '-----------' >> {output} && "
        "beta_diversity.py -i {input.bracken} --type bracken >> {output}"
        
rule virus_krona:
    input:
        breport = "05_viral_anno/{id}.viral.breport"
    output:
        krona_txt = '05_viral_anno/{id}.viral.krona.txt',
        krona_html = '05_viral_anno/{id}.viral.krona.html'
    shell:
        'kreport2krona.py -r {input.breport} -o {output.krona_txt} --no-intermediate-ranks &&  '
        'ktImportText {output.krona_txt} -o {output.krona_html}'


rule virus_extract_viral_reads:
    input:
        out_k2 = "05_viral_anno/{id}.viral.k2",
        report = "05_viral_anno/{id}.viral.k2_report.txt",
        fq1 = '04_annotation/{id}_viral_unclass.R1.fq',
        fq2 = '04_annotation/{id}_viral_unclass.R2.fq'
    output:
        fa1 = '05_viral_anno/{id}.tid10239.R1.fq',
        fa2 = '05_viral_anno/{id}.tid10239.R2.fq'
    shell:
        "extract_kraken_reads.py -k {input.out_k2} --fastq-output "
            "-s {input.fq1} -s2 {input.fq2} -t 10239 -r {input.report} "
            "-o {output.fa1} -o2 {output.fa2}"


