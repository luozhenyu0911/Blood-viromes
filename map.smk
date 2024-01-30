

# remove human seq. against hg38
rule rm_hg38seqs:
	input:
		"input/{}.bam".format(config['sample']['bam_id'])
	output:
		"01_align2hg38/{}.unmapped_hg38.bam".format(config['sample']['id']
	threads:
		config['threads']['bwa']
	shell:
		"samtools view -@ {threads} -bh -f 4 {input} > {output}"

# remove LCR-hs38_rmsk seq
rule rm_repeats:
	input:
		bam = "01_align2hg38/{id}.unmapped_hg38.bam"
		region = config['params']['repeats']
	output:
		in_bam = "01_align2hg38/{id}.is_repeats.bam"
		out_bam = "01_align2hg38/{id}.rm_repeats.bam"
	threads:
		config['threads']['bwa']
	shell:
		"samtools view -@ {threads} -bh -o {output.in_bam} -U {output.out_bam} -L {region}"

# get fq1 after rm hg38 and repeats seq.
rule fq1:
	input:
		out_bam = "01_align2hg38/{id}.rm_repeats.bam"
		rawfq1 = config['sample']['fq1']
		rawfq1 = config['sample']['fq2']
	output:
		seq_id = "01_align2hg38/{id}.rm_repeats.id"
		R1_id = "01_align2hg38/{id}.rm_repeats.R1.id"
		R2_id = "01_align2hg38/{id}.rm_repeats.R2.id"
		out_bam = "01_againsthg38/{id}.rm_repeats.bam"
		rmrepeats_fq1 = "01_againsthg38/{id}.rm_repeats.R1.fq.gz"
		rmrepeats_fq2 = "01_againsthg38/{id}.rm_repeats.R2.fq.gz"
	threads:
		config['threads']['bwa']
	shell:
		"samtools view -@ {threads} {input.out_bam} | cut -f1|sort -u > {output.seq_id} &&"
		"awk '{{print $1"/1"}}' {output.seq_id} > {output.R1_id} && "
		"awk '{{print $1"/2"}}' {output.seq_id} > {output.R2_id} && "
		"seqtk subseq {input.rawfq1} {output.R1_id} |gzip > {output.rmrepeats_fq1} && "
		"seqtk subseq {input.rawfq2} {output.R2_id} |gzip > {output.rmrepeats_fq1} "

# remove Ref T2T seq
rule Align2T2T:
	input:
		rmrepeats_fq1 = "01_againsthg38/{id}.rm_repeats.R1.fq.gz"
		rmrepeats_fq2 = "01_againsthg38/{id}.rm_repeats.R2.fq.gz"
	output:
		T2T_bam = "02_againsthg38/{id}.T2T.bam"
	threads:
		config['threads']['bwa']
	shell:
		"samtools view -@ {threads} -bh -o {output.in_bam} -U {output.out_bam} -L {region}"


















samtools view $id.in.bam |cut -f1|sort -u > $id.in.seqid
samtools view $i |cut -f1|sort -u > $i.total.seqid


samtools view -@ 24 $i -b -h -o $id.in.bam -U $id.out.bam -L LCR-hs38_rmsk.bed





# Rule for mapping
rule map_reads:
    input:
        ref = REF,
        fq_1 = expand("{sample}", sample=config['samples']['fq_1']),
        fq_2 = expand("{sample}", sample=config['samples']['fq_2'])
    output:
        sam = "Align/{}.aln_mem.sam".format(config['samples']['id'])
    threads:
        config['threads']['hisat2']
    shell:
        # map
        "hisat2 -p {threads} -x {input.ref} -1 {input.fq_1} -2 {input.fq_2} -S {output.sam} 2>Align/aln.err"

# generate raw bam 
rule sam2bam:
    input:
        "Align/{id}.aln_mem.sam"
    output:
        "Align/{id}.raw.bam"
    threads:
        config['threads']['hisat2']
    shell:
        "samtools view -@ {threads} -bhS {input} -o {output}"

# filter reads with mapping quality < INT and unmapped
rule get_uniq_bam:
    input:
        "Align/{id}.raw.bam"
    output:
        "Align/{id}.uniq.bam"
    threads:
        config['threads']['hisat2']
    params:
        reads_map_q = config['params']['reads_map_q']
    shell:
        "samtools view -@ {threads} -bhS -q {params.reads_map_q} -F 0x400 {input} -o {output}"

# sort uniq bam
rule sort_uniq_bam:
    input:
        "Align/{id}.uniq.bam"
    output:
        "Align/{id}.uniq.sort.bam"
    threads:
        config['threads']['hisat2']
    shell:
        "samtools sort -@ {threads} {input} -O bam -o {output} &&"
        "samtools index -@ {threads} {output}"

# samtools flagstat for uniq.bam
rule flagstat_uniq_bam:
    input:
        "Align/{id}.uniq.sort.bam"
    output:
        "Align/{id}.uniq.flagstat.txt"
    threads:
        config['threads']['hisat2']
    shell:
        "samtools flagstat -@ {threads} {input} >> {output}"
        
rule flagstat_raw_bam:
    input:
        "Align/{id}.raw.bam"
    output:
        "Align/{id}.raw.flagstat.txt"
    threads:
        config['threads']['hisat2']
    shell:
        "samtools flagstat -@ {threads} {input} >> {output}"


# get reads count of gene_id or transcript_id 
rule get_reads_count:
    input:
        "Align/{}.uniq.sort.bam".format(config['samples']['id'])
    output:
        gene_count = "Align/{}_gene_count.txt".format(config['samples']['id']),
        transcript_count = "Align/{}_transcript_count.txt".format(config['samples']['id'])
    threads:
        config['threads']['featureCounts']
    params:
        gtf = config['params']['gtf'],
        featureCounts = config['params']['featureCounts']
        
    shell:
        "{params.featureCounts} -T {threads} -p -a {params.gtf} -g gene_name -o {output.gene_count} {input} |"
        "{params.featureCounts} -T {threads} -p -a {params.gtf} -g transcript_id -o {output.transcript_count} {input}"

# This looks at the coverage across the genome, as well as percent coverage at particular depths (4X, 10X, 30X)
rule coverage_depth:
    input:
        "Align/{id}.uniq.sort.bam"
    output:
        "Align/{id}.uniq.bam_coverage_depth.txt"
    params:
        toolsdir = config['params']['toolsdir'],
        ref = config['params']['ref_fa']
    shell:
        "perl {params.toolsdir}/tools/depthV2.0.pl -l $({params.toolsdir}/tools/fasta_non_gapped_bases.py {params.ref}) {input} Align > {output}"
        



















