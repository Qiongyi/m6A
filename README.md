### Core commmands and custom Perl scripts used in m6A study (bioinformatic analysis)

## 1. Quality control of raw sequencing reads

# 1) using cutadapt to trim off low-quality nucleotides (Phred quality lower than 20) and Illumina adaptor sequences at the 3’ end

for i in "s_6_In_C13" "s_6_In_C15" "s_6_In_C18" "s_6_In_C20" "s_6_In_C3" "s_6_In_C8" "s_6_IP_C13" "s_6_IP_C15" "s_6_IP_C18" "s_6_IP_C20" "s_6_IP_C3" "s_6_IP_C8"
do
qsub -j y -o qsubout/cutadapt.$i.qsubout -b y -cwd -l h_vmem=50G -N $i /clusterdata/apps/cutadapt-1.3/bin/cutadapt -q 20 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC /illumina/Data/131212_7001408_0063_BC39MNACXX/TBJW_RNAMETHYL/fastq/Run/${i}_read1.fastq.gz -o /illumina/Data/131212_7001408_0063_BC39MNACXX/TBJW_RNAMETHYL/fastq/cutadapt_new/${i}_read1.fastq.gz
qsub -j y -o qsubout/cutadapt.$i.qsubout -b y -cwd -l h_vmem=50G -N $i /clusterdata/apps/cutadapt-1.3/bin/cutadapt -q 20 -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT /illumina/Data/131212_7001408_0063_BC39MNACXX/TBJW_RNAMETHYL/fastq/Run/${i}_read2.fastq.gz -o /illumina/Data/131212_7001408_0063_BC39MNACXX/TBJW_RNAMETHYL/fastq/cutadapt_new/${i}_read2.fastq.gz
done

# 2) using Bowtie2 to align against a custom contaminant list to filter contaminant reads

for i in "s_6_In_C13" "s_6_In_C15" "s_6_In_C18" "s_6_In_C20" "s_6_In_C3" "s_6_In_C8" "s_6_IP_C13" "s_6_IP_C15" "s_6_IP_C18" "s_6_IP_C20" "s_6_IP_C3" "s_6_IP_C8"
do
qsub -j y -o qsubout/$i.rRNA.PhiX.out -b y -cwd -pe onehost 10 -l h_vmem=3G -N $i.rRNA bowtie2 --local -q --phred33 -p 10 --no-hd --no-unal -x /clusterdata/hiseq_apps/resources/freeze001/mm10/mm10rRNA_PhiX -1 /illumina/Data/131212_7001408_0063_BC39MNACXX/TBJW_RNAMETHYL/fastq/cutadapt_new/${i}_read1.fastq.gz -2 /illumina/Data/131212_7001408_0063_BC39MNACXX/TBJW_RNAMETHYL/fastq/cutadapt_new/${i}_read2.fastq.gz -S ${i}.rRNAPhiX.sam
done

# 3) Custom Perl scripts to remove contaminant reads based on bowtie alignment (custom Perl scripts can be downloaded from https://github.com/Qiongyi/m6A/)
for i in "s_6_In_C3" "s_6_In_C8" "s_6_In_C13" "s_6_In_C15" "s_6_In_C18" "s_6_In_C20" "s_6_IP_C3" "s_6_IP_C8" "s_6_IP_C13" "s_6_IP_C15" "s_6_IP_C18" "s_6_IP_C20"
do
Discard_fastq_PE.pl /illumina/Data/131212_7001408_0063_BC39MNACXX/TBJW_RNAMETHYL/results/bowtie2rRNAPhiX/$i.rRNAPhiX.sam /illumina/Data/131212_7001408_0063_BC39MNACXX/TBJW_RNAMETHYL/fastq/cutadapt_new/${i}_read1.fastq.gz /illumina/Data/131212_7001408_0063_BC39MNACXX/TBJW_RNAMETHYL/fastq/cutadapt_new/${i}_read2.fastq.gz /illumina/Data/131212_7001408_0063_BC39MNACXX/TBJW_RNAMETHYL/fastq/cutadapt_new/filter/${i}_read1.fastq /illumina/Data/131212_7001408_0063_BC39MNACXX/TBJW_RNAMETHYL/fastq/cutadapt_new/filter/${i}_read2.fastq filter_stats.xls
done

## 2. Alignment of quality filtered reads

# 1) Using Tophat2 to align the quality filtered paired-end reads to the mouse genome (build mm10)

for i in "s_6_In_C3" "s_6_In_C8" "s_6_In_C13" "s_6_In_C15" "s_6_In_C18" "s_6_In_C20" "s_6_IP_C3" "s_6_IP_C8" "s_6_IP_C13" "s_6_IP_C15" "s_6_IP_C18" "s_6_IP_C20"
do
qsub -j y -o qsubout/$i.tophat.cutadapt.out -b y -cwd -pe onehost 8 -l h_vmem=5G -N cutadapt.$i.tophat tophat --mate-inner-dist 50 --mate-std-dev 50 --microexon-search --min-anchor 10 -p 8 -o ${i}_tophat_filter -G /clusterdata/hiseq_apps/resources/freeze001/mm10/mm10.gtf --transcriptome-index /illumina/Data/131212_7001408_0063_BC39MNACXX/TBJW_RNAMETHYL/results/tophat2/transcriptome_data/known /clusterdata/hiseq_apps/resources/freeze001/mm10/bowtie2_index/mm10 /illumina/Data/131212_7001408_0063_BC39MNACXX/TBJW_RNAMETHYL/fastq/cutadapt_new/filter/${i}_read1.fastq /illumina/Data/131212_7001408_0063_BC39MNACXX/TBJW_RNAMETHYL/fastq/cutadapt_new/filter/${i}_read2.fastq
done

# 2) Using Samtools to keep concordantly aligned reads

for i in "s_6_In_C3" "s_6_In_C8" "s_6_In_C13" "s_6_In_C15" "s_6_In_C18" "s_6_In_C20" "s_6_IP_C3" "s_6_IP_C8" "s_6_IP_C13" "s_6_IP_C15" "s_6_IP_C18" "s_6_IP_C20"
do
samtools view -h -bu -f2 -q 4 /illumina/Data/131212_7001408_0063_BC39MNACXX/TBJW_RNAMETHYL/results/tophat2/${i}_tophat_filter/accepted_hits.bam |samtools sort - /illumina/Data/131212_7001408_0063_BC39MNACXX/TBJW_RNAMETHYL/results/tophat2/${i}_tophat_filter/Properly_accepted_hits.sorted
samtools index /illumina/Data/131212_7001408_0063_BC39MNACXX/TBJW_RNAMETHYL/results/tophat2/${i}_tophat_filter/Properly_accepted_hits.sorted.bam
done

# 3) Using Picard to mark duplicates

for i in "s_6_In_C3" "s_6_In_C8" "s_6_In_C13" "s_6_In_C15" "s_6_In_C18" "s_6_In_C20" "s_6_IP_C3" "s_6_IP_C8" "s_6_IP_C13" "s_6_IP_C15" "s_6_IP_C18" "s_6_IP_C20"
do
qsub -j y -o qsubout/MarkDuplicates.$i -b y -cwd -pe onehost 8 -l h_vmem=5G -N $i java -Xmx5g -XX:ParallelGCThreads=8 -jar /clusterdata/hiseq_apps/bin/freeze001/picard/picard-tools-1.72/MarkDuplicates.jar INPUT=/illumina/Data/131212_7001408_0063_BC39MNACXX/TBJW_RNAMETHYL/results/tophat2/${i}_tophat_filter/Properly_accepted_hits.sorted.bam OUTPUT=/illumina/Data/131212_7001408_0063_BC39MNACXX/TBJW_RNAMETHYL/results/tophat2/${i}_tophat_filter/Properly_accepted_hits.asd.bam METRICS_FILE=/illumina/Data/131212_7001408_0063_BC39MNACXX/TBJW_RNAMETHYL/results/tophat2/${i}_tophat_filter/Properly_accepted_hits.asd.bam.dupl AS=true VALIDATION_STRINGENCY=LENIENT
done

# 4) Using Samtools to keep uniquely aligned reads
for i in "s_6_In_C3" "s_6_In_C8" "s_6_In_C13" "s_6_In_C15" "s_6_In_C18" "s_6_In_C20" "s_6_IP_C3" "s_6_IP_C8" "s_6_IP_C13" "s_6_IP_C15" "s_6_IP_C18" "s_6_IP_C20"
do
samtools view -b -o /illumina/Data/131212_7001408_0063_BC39MNACXX/TBJW_RNAMETHYL/results/tophat2/${i}_tophat_filter/Properly_accepted_hits.asd.rmdup.bam -F 1024 /illumina/Data/131212_7001408_0063_BC39MNACXX/TBJW_RNAMETHYL/results/tophat2/${i}_tophat_filter/Properly_accepted_hits.asd.bam 
samtools index /illumina/Data/131212_7001408_0063_BC39MNACXX/TBJW_RNAMETHYL/results/tophat2/${i}_tophat_filter/Properly_accepted_hits.asd.rmdup.bam


## 3. Detection of m6A sites (winscore method)

for i in "C3" "C8" "C13" "C15" "C18" "C20"
do
qsub -j y -o qsubout/$i.winscore.out -b y -cwd -l h_vmem=20G -N $i.winscore RNAMethy_winscore.pl /clusterdata/hiseq_apps/resources/freeze001/mm10/mm10.gtf /illumina/Data/131212_7001408_0063_BC39MNACXX/TBJW_RNAMETHYL/results/tophat2/ReadsNum_for_EachGene_normalized_to_1kb.xls /illumina/Data/131212_7001408_0063_BC39MNACXX/TBJW_RNAMETHYL/results/tophat2/s_6_IP_${i}_tophat_filter/Properly_accepted_hits.asd.rmdup.bam /illumina/Data/131212_7001408_0063_BC39MNACXX/TBJW_RNAMETHYL/results/tophat2/s_6_In_${i}_tophat_filter/Properly_accepted_hits.asd.rmdup.bam ${i}.winscore.xls ${i}.winscore.all.xls
done

## 4. Detection of differentially methylated m6A sites

# 1) get the unique peak summit for each condition

RNAMethy_GetUniqSummit.pl Naive_SummitPosition.xls Naive_UniqSummitPosition.xls
RNAMethy_GetUniqSummit.pl Context_SummitPosition.xls Context_UniqSummitPosition.xls
RNAMethy_GetUniqSummit.pl FC_SummitPosition.xls FC_UniqSummitPosition.xls

# 2) combine the summit position for three conditions
RNAMethy_CombineUniqSummit.pl Naive_UniqSummitPosition.xls Context_UniqSummitPosition.xls FC_UniqSummitPosition.xls Naive_Context_FC_SummitPosition.xls

# 3) remove redundant positions (select highest peak)
RNAMethy_GetUniqSummit2.pl 150 Naive_Context_FC_SummitPosition.xls Naive_Context_FC_SummitPosition.group.xls Naive_Context_FC_UniqSummitPosition.xls

# 4) student’s t-test 
RNAMethy_Summit_TTest.pl Naive_Context_FC_UniqSummitPosition.xls Naive_Context_FC.ttest.xls


## 5. Analysis of m6A location (shown is an example for Naive group)

# 1) for unique summit positions
RNAMethy_AddGenomicLocation_Case.pl /clusterdata/hiseq_apps/resources/freeze001/mm10/mm10.gtf Naive_UniqSummitPosition.xls Naive.final.GenomicLocation.xls

# 2) for all summit positions
RNAMethy_GenomicLocation.pl /clusterdata/hiseq_apps/resources/freeze001/mm10/mm10.gtf Naive_SummitPosition.xls Naive_SummitPosition_TranscriptSegments.xls Naive_SummitPosition_TranscriptSegments_stat.xls

## 6. Identification and clustering of enriched motifs (shown is an example for Naive group)

# 1) to count occurrences of each 4–6-nucleotide k-mer
for i in {4..6}
do
RNAMethy_MotifFinder1_Kmer.pl $i 25 /clusterdata/hiseq_apps/resources/freeze001/mm10/mm10_gene.fasta /illumina/Data/131212_7001408_0063_BC39MNACXX/TBJW_RNAMETHYL/results/winscore/Naive.final.GenomicLocation.xls ../Naive_control_GenomicLocation.xls m6A_Kmers_${i}bp.xls
done

# 2) to select motifs enriched more than 1.5 fold with the corrected P value < 0.05
for i in {4..6}
do
awk '$1=="Kmer" || ($4>=$5*1.5 && $3<0.05)' m6A_Kmers_${i}bp.xls > m6A_Kmers_${i}bp_1.5fold.xls
done

# 3) to cluster kmers
RNAMethy_MotifFinder2_Cluster.pl 1 m6A_Kmers_4bp_1.5fold.xls m6A_Kmers_5bp_1.5fold.xls m6A_Kmers_6bp_1.5fold.xls Kmer_cluster.xls

# 4) to calculate the position-specific scoring matrix for detected motif (shown is an example for the first seed)
RNAMethy_MotifFinder3_PSSM.pl Kmer_sequences_seed1.fa seed1.PSSM.xls
