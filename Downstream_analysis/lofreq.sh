#!/bin/bash

# load modules
module load htslib lofreq/2.1.5 samtools

bam=$1
# bam=Simulations_3/var/jobs/Mapping/Sample0/Sample0.sorted.chrR.f2F2308q20.wo_XA.bam
sample_id=$(echo $bam | awk -F'/' '{print $5}')

# paths
outpath=$(echo $bam | awk -F'/' '{print $1}') #E.g.: Simulations_1
outpath=${outpath}/var/jobs/Lofreq_all/
mkdir -p $outpath


echo "Starting with $bam"
# infiles
reference=/gpfs/projects/bsc83/Data/assemblies/T2T_CHM13/chrR/Human_hs1-rDNA_genome_v1.0/hs1-rDNA_v1.0.fa

#First step: add base- and indel-alignment qualities (BAQ, IDAQ) to BAM file
lofreq alnqual $bam $reference > ${outpath}/${sample_id}.qual.sam
samtools view -S -b ${outpath}/${sample_id}.qual.sam > ${outpath}/${sample_id}.qual.bam

#Run variant calling with lofreq if only SNPs
#lofreq call --call-indels -f $reference -l /gpfs/projects/bsc83/Projects/ribosomal_RNAs/Jose/07_RepeatMasker/data/pre-rRNA_47S.regions.bed --verbose ${TMPDIR}/${sample_id}.qual.bam -o ${outpath}/${sample_id}.output.vcf --force-overwrite

lofreq indelqual --dindel -f $reference -o ${outpath}/${sample_id}.realigned.bam ${outpath}/${sample_id}.qual.bam

lofreq call --call-indels -f $reference -l /gpfs/projects/bsc83/Projects/ribosomal_RNAs/Jose/07_RepeatMasker/data/pre-rRNA_47S.regions.bed --verbose ${outpath}/${sample_id}.realigned.bam -o ${outpath}/${sample_id}.output_indel.vcf --force-overwrite

rm ${outpath}/${sample_id}.qual.sam ${outpath}/${sample_id}.qual.bam ${outpath}/${sample_id}.realigned.bam
echo "Finished with ${outpath}/${sample_id}.output_indel.vcf"


