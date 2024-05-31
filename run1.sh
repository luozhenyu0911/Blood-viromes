#!/bin/sh
#############################################
#add your conda env path
# export PATH="/your/conda/env/path/bin/:$PATH"
refpath=/BIGDATA2/gzfezx_shhli_2/database/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa
Vrefpath=/BIGDATA2/gzfezx_shhli_2/database/viral_refseq/viral.1.1.genomic.fna
# fastqPATH=/your/sample/fastq/folder/path
n_threads=24
raw_bam=$1
sample=$(basename $raw_bam .cram)

#step 1 unmapped
echo "step 1 unmapped $samp"
echo "----------------------------------------------------"
echo ""
samtools view -@ ${n_threads} -f 4 -bh $raw_bam > $sample.f4.bam 
samtools view -@ ${n_threads} -f 4 -F 264 -bh $raw_bam > $sample.f4_264.bam

echo "step2 clip extraction $samp"
echo "----------------------------------------------------"
samtools view -@ ${n_threads} -H $raw_bam > header.txt
samtools view -@ ${n_threads} $raw_bam | awk '$6 ~ /S/{print}' > $sample.clip_tmp.bam
cat header.txt $sample.clip_tmp.bam | samtools view -@ ${n_threads} -bS > $sample.clip.bam
rm -f $sample.clip_tmp.bam

echo "step2 bam to fastq"
echo "----------------------------------------------------"
echo ""
samtools fastq -@ ${n_threads} $sample.f4.bam > $sample.f4.fq 
samtools fastq -@ ${n_threads} $sample.f4_264.bam > $sample.f4_264.fq
samtools fastq -@ ${n_threads} $sample.clip.bam >  $sample.clip.fq

echo "step3 trim SLIDINGWINDOW:4:20 MINLEN:50"
echo "-----------------did nothing-----------------------------------"
echo "step4  viral ref by bwa"
echo "----------------------------------------------------"
for i in $sample.*.fq
do
id=$(basename $i .fq)
bwa mem -C -t ${n_threads} ${Vrefpath} $i | samtools sort -o $id.V.bam -@ ${n_threads} -O bam -
done

echo "Step 4.1 count number of virus"
echo "----------------------------------------------------"
for vbam in *.V.bam
do
id=$(basename $vbam .V.bam)
samtools index -@ ${n_threads} ${n_threads} $vbam
samtools idxstats $vbam > $id.V.sort.tsv

echo "step 5 extract readname of Virus mapped"
echo "----------------------------------------------------"
samtools view $vbam |cut -f1 > $id.readnV.tsv

echo "step 6 extract alignment file only mapped V read in bam hg38 .."
echo "----------------------------------------------------"
# samtools view -bh -N $id.readnV.tsv $raw_bam -o $id.readnV.bam
samtools view -@ ${n_threads} $raw_bam -bh -o $sample.bam
samtools index -@ ${n_threads} $sample.bam
java -jar /BIGDATA2/gzfezx_shhli_2/software/picard.jar FilterSamReads \
       I=$sample.bam \
       O=$id.readnV.bam \
       READ_LIST_FILE=$id.readnV.tsv FILTER=includeReadList

samtools index $id.readnV.bam


#extract only mapped V read in sam hg38
echo "step 7 bedtools and cluster."
echo "----------------------------------------------------"
bedtools genomecov -bg -ibam $id.readnV.bam | bedtools cluster -d 200 -i |awk '$4>5'| awk '!arr[$5] {arr[$5]=$0; if(prevline) print prevline; print} {prevline=$0}' > $id.V.bed
bedtools genomecov -bg -ibam $id.readnV.bam | bedtools cluster -d 200 -i |awk '$4>5'| tail -1 >> $id.V.bed

###Group 1 clauster = 1 line
echo "step 8 group of cluster."
echo "----------------------------------------------------"
awk 'NR % 2 == 1'  $id.V.bed| cut -f1,2 >  $id.V.tmp1.bed
awk 'NR % 2 == 0'  $id.V.bed| cut -f3,4,5 >  $id.V.tmp2.bed
paste $id.V.tmp1.bed $id.V.tmp2.bed > $id.V.cluster.bed
rm -r *.tmp1.bed

echo "step 9 intersect position with Repeatmarker or knowngene from UCSC."
echo "----------------------------------------------------"

gene=/BIGDATA2/gzfezx_shhli_2/database/GRCh38/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.bed
repeat_bed=/WORK/gzfezx_shhli_3/BioDatahub/human_reference/GRCh38/rmsk.txt6col.bed

bedtools intersect -a $repeat_bed -b $id.V.cluster.bed > 9.1_repeat_${id}.V.cluster.bed
bedtools intersect -a $gene -b $id.V.cluster.bed > 9.2_gene_${id}.V.cluster.bed

done

echo "${samp} finished run"
done

echo "DONE DONE DONE"


