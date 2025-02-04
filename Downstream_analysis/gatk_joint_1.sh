#!/bin/bash

# load modules
module load java-openjdk gatk

# input data
bam=$1
ploidy=$2

#Data modified from the input
genotypes=$(python3 -c "import math; num_alleles=6; ploidy=$ploidy; maxgen=math.factorial(num_alleles + ploidy - 1) // (math.factorial(num_alleles - 1) * math.factorial(ploidy)); print(maxgen)")
sample=$(echo $bam | awk -F'/' '{print $5}')

#Reference genome
reference=/gpfs/projects/bsc83/Data/assemblies/T2T_CHM13/chrR/Human_hs1-rDNA_genome_v1.0/hs1-rDNA_v1.0.fa


# paths
outpath=$(echo $bam | awk -F'/' '{print $1}') #E.g.: Simulations_1
outpath=${outpath}/var/jobs/Joint_all/
mkdir -p $outpath

#echo "$sample in $outpath for ploidy ${ploidy} and sample ${bam} has started at $(date)"

gatk --java-options "-Xmx800g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" HaplotypeCaller \
	-I ${bam} \
	-R ${reference} \
	-O ${outpath}/${sample}.${ploidy}.g.vcf.gz \
	-ERC BP_RESOLUTION \
	--sample-ploidy ${ploidy} \
	--intervals /gpfs/projects/bsc83/Projects/ribosomal_RNAs/Jose/07_RepeatMasker/data/pre-rRNA_47S.regions.bed \
	--max-reads-per-alignment-start 0 \
	--dont-use-soft-clipped-bases true \
	--native-pair-hmm-threads $(nproc) \
	--max-genotype-count ${genotypes} 


echo "$outpath for ploidy ${ploidy} and sample ${bam} finished after gatk at $(date)"

#I don't specify number of threads in the function, as it detects it automatically (except snakemake):
#	--native-pair-hmm-threads ${threads} \
