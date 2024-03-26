rule classify_kraken2:
    input:
        db = config['params']['k2_db'],
        fq1 = "03_align2pg/{id}.rm_pg.R1.fq.gz",
        fq2 = "03_align2pg/{id}.rm_pg.R2.fq.gz"
    output:
        # unclassified = "04_annotation/{id}.Unclass#.fa",
        # classified = "04_annotation/{id}.Class#.fa",
        # unclassified = expand("04_annotation/{id}.Unclass_{i}.fa", i=range(1,3)),
        # classified = expand("04_annotation/{id}.Class_{i}.fa", i=range(1,3)),
        out_k2 = "04_annotation/{id}.k2",
        report = "04_annotation/{id}.k2_report.txt"
    threads:
        config['threads']['bwa']
    shell:
        "kraken2 --db {input.db} --gzip-compressed --paired {input.fq1} {input.fq2} "
        "--threads {threads} "
        "--output {output.out_k2} --report {output.report}"
        # " --unclassified-out 04_annotation/{id}.Unclass#.fa --classified-out 04_annotation/{id}.Class#.fa"

rule bracken:
    input:
        db = config['params']['k2_db'],
        report = "04_annotation/{id}.k2_report.txt"
    output:
        bracken = "04_annotation/{id}.bracken",
        breport = "04_annotation/{id}.breport"
    threads:
        config['threads']['bwa']
    shell:
        "bracken -d {input.db} -i {input.report} -o {output.bracken} -w {output.breport} -t 0"
        
# rule diversity:
#     input:
#         bracken = "04_annotation/{id}.bracken"
#     output:
#         alpha_diversity = "04_annotation/{id}.alpha_beta_diversity.stats",
#         beta_diversity = "04_annotation/{id}.beta_diversity.stats"
#     shell:
#         "alpha_diversity.py -f {input.bracken} -a BP >> alpha_beta_diversity.stats && "
#         "beta_diversity.py -i {input.bracken} --type bracken >> alpha_beta_diversity.stats"
        
rule krona:
    input:
        breport = "04_annotation/{id}.breport"
    output:
        krona_txt = '04_annotation/{id}.krona.txt',
        krona_html = '04_annotation/{id}.krona.html'
    shell:
        'kreport2krona.py -r {input.breport} -o {output.krona_txt} --no-intermediate-ranks &&  '
        'ktImportText {output.krona_txt} -o {output.krona_html}'


rule extract_viral_reads:
    input:
        out_k2 = "04_annotation/{id}.k2",
        report = "04_annotation/{id}.k2_report.txt",
        fq1 = "03_align2pg/{id}.rm_pg.R1.fq.gz",
        fq2 = "03_align2pg/{id}.rm_pg.R2.fq.gz"
    output:
        fa1 = '04_annotation/{id}.tid10239.R1.fa',
        fa2 = '04_annotation/{id}.tid10239.R2.fa'
    shell:
        "extract_kraken_reads.py -k {input.out_k2} --include-children "
            "-s {input.fq1} -s2 {input.fq2} -t 10239 -r {input.report} "
            "-o {output.fa1} -o2 {output.fa2}"

rule extract_unclass_reads:
    input:
        out_k2 = "04_annotation/{id}.k2",
        report = "04_annotation/{id}.k2_report.txt",
        fq1 = "03_align2pg/{id}.rm_pg.R1.fq.gz",
        fq2 = "03_align2pg/{id}.rm_pg.R2.fq.gz"
    output:
        fa1 = '04_annotation/{id}.unclass.R1.fa',
        fa2 = '04_annotation/{id}.unclass.R2.fa'
    shell:
        #"extract_kraken_reads.py -k {input.out_k2} --include-children --fastq-output "
        "extract_kraken_reads.py -k {input.out_k2} --include-children "
            "-s {input.fq1} -s2 {input.fq2} -t 0 -r {input.report} "
            "-o {output.fa1} -o2 {output.fa2}"
            
rule merge_viral_unclass:
    input:
        viral_fa1 = '04_annotation/{id}.tid10239.R1.fa',
        viral_fa2 = '04_annotation/{id}.tid10239.R2.fa',
        unclass_fa1 = "04_annotation/{id}.unclass.R1.fa",
        unclass_fa2 = "04_annotation/{id}.unclass.R2.fa"
    output:
        fa1 = '04_annotation/{id}_viral_unclass.R1.fa',
        fa2 = '04_annotation/{id}_viral_unclass.R2.fa'
    shell:
        "cat {input.viral_fa1} {input.unclass_fa1} > {output.fa1} && "
        "cat {input.viral_fa2} {input.unclass_fa2} > {output.fa2} "

rule merge_fa:
    input:
        fa1 = '04_annotation/{id}_viral_unclass.R1.fa',
        fa2 = '04_annotation/{id}_viral_unclass.R2.fa'
    output:
        '04_annotation/{id}_viral_unclass.all.fa'
    shell:
        "cat {input.fa1} {input.fa2} > {output}"

rule blastn:
    input:
        all_fa = '04_annotation/{id}_viral_unclass.all.fa',
        db = config['params']['ref_viral']
    output:
        "04_annotation/{id}_virus_unclss.blastn.res.txt"
    threads:
        int(config['threads']['bwa'])-1
    shell:
        "blastn -query {input.all_fa} -db {input.db} -outfmt 7 -evalue 0.00001 "
            "-max_target_seqs 1 -num_threads {threads} -out {output}"
            
rule blastn_anno:
    input:
        res_blast = "04_annotation/{id}_virus_unclss.blastn.res.txt",
        seq_anno = config['params']['ref_viral'] # there is the seq. title (e.g.>NC_004066.1 Lactococcus phage ul36, complete genome) in the fasta file
    output:
        "04_annotation/{id}_virus_unclss.blastn.res.anno.txt"
    threads:
        int(config['threads']['bwa'])-1
    params:
        tools_dir = config['params']['tools']
    shell:
        "python {params.tools_dir}/get_annotation.3.py -b {input.res_blast} -a {input.seq_anno}_seqname.txt -o {output}"


