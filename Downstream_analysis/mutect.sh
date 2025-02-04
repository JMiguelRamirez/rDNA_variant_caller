#!/bin/bash

## GATk version
#gatk=4.5.0.0

# load modules
module load  java-openjdk gatk bcftools

bam=$1
sample_id=$(echo $bam | awk -F'/' '{print $5}')

# paths
outpath=$(echo $bam | awk -F'/' '{print $1}') #E.g.: Simulations_1
outpath=${outpath}/var/jobs/Mutect_all/
mkdir -p $outpath


# infiles
reference=/gpfs/projects/bsc83/Data/assemblies/T2T_CHM13/chrR/Human_hs1-rDNA_genome_v1.0/hs1-rDNA_v1.0.fa

#I will not be saving in TMPDIR because it is giving problems

# Single-sample gVCF calling
gatk --java-options "-Xmx115g -Xms100g" Mutect2 \
	-I ${bam} \
	-R ${reference} \
	-O ${TMPDIR}/${sample_id}.mutect_1.g.vcf.gz \
	--mitochondria-mode \
	--native-pair-hmm-threads $(nproc) \
	-L /gpfs/projects/bsc83/Projects/ribosomal_RNAs/Jose/07_RepeatMasker/data/pre-rRNA_47S.regions.bed

bcftools index -t -f ${TMPDIR}/${sample_id}.mutect_1.g.vcf.gz	

#Split multiallelic positions into biallelic positions. If we added the reference, reference alleles could to be adjusted but we don't need that
bcftools annotate \
    	-x ^INFO/DP,^FORMAT/GT,^FORMAT/AD,^FORMAT/DP,^FORMAT/GQ \
    	${TMPDIR}/${sample_id}.mutect_1.g.vcf.gz > ${TMPDIR}/${sample_id}_filtering.g_wo_PL.vcf.gz

#Remove alternative alleles that have not been called from g.vcf. Sometime GATK adds a variant in the ALT column when it has not been called
#-O z is to get the output as .gz
bcftools view --trim-alt-alleles ${TMPDIR}/${sample_id}_filtering.g_wo_PL.vcf.gz -O z -o ${TMPDIR}/${sample_id}_trim.vcf.gz 

bcftools norm \
        -m -any \
        --old-rec-tag INFO \
	--fasta-ref ${reference} \
        ${TMPDIR}/${sample_id}_trim.vcf.gz > ${outpath}/${sample_id}.mutect.g.vcf.gz
        	
bcftools index -t -f ${outpath}/${sample_id}.mutect.g.vcf.gz

echo "${outpath}/${sample_id}.mutect.g.vcf.gz created"

#Creating index in .tbi format (-t), we will need index if we plan to use bcftools merge
