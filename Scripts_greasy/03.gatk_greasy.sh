#!/bin/bash

# load modules
module load java-openjdk gatk bcftools R
#module load java-openjdk/17.0.11+9
#module load gatk/4.5.0.0
       
sample=$1
input_folder=$2
output_folder=$3
ploidy=$4
type=$5

genotypes=$(python3 -c "import math; num_alleles=6; ploidy=$ploidy; maxgen=math.factorial(num_alleles + ploidy - 1) // (math.factorial(num_alleles - 1) * math.factorial(ploidy)); print(maxgen)")

echo $ploidy
echo $genotypes
reference=/gpfs/projects/bsc83/Data/assemblies/T2T_CHM13/chrR/Human_hs1-rDNA_genome_v1.0/hs1-rDNA_v1.0.fa
bam=${input_folder}/${sample}.sorted.chrR.f2F2308q20.wo_XA.bam
intervals=/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Jose/07_RepeatMasker/data/pre-rRNA_47S.included_intervals.bed
	
if [[ $type == "RNA" ]]; then
	gatk --java-options "-Xmx115g -Xms100g" HaplotypeCaller \
		-I ${bam} \
		-R ${reference} \
		-O ${TMPDIR}/${sample}.${ploidy}.g.vcf.gz \
		--sample-ploidy ${ploidy} \
		--intervals $intervals \
		--max-reads-per-alignment-start 0 \
		--dont-use-soft-clipped-bases true \
		--native-pair-hmm-threads $(nproc) \
		--max-genotype-count ${genotypes} \
		--min-base-quality-score 0 \
		-ERC BP_RESOLUTION
else #DNA
#The differences are that in RNA we want to include -ERC BP_RESOLUTION to have information on all positions and not only variant positions. In this way, we can know if a DNA variant is not called in RNA because there are not enough reads or because it is not expressed
	gatk --java-options "-Xmx115g -Xms100g" HaplotypeCaller \
		-I ${bam} \
		-R ${reference} \
		-O ${TMPDIR}/${sample}.${ploidy}.g.vcf.gz \
		--sample-ploidy ${ploidy} \
		--intervals $intervals \
		--max-reads-per-alignment-start 0 \
		--dont-use-soft-clipped-bases true \
		--native-pair-hmm-threads $(nproc) \
		--max-genotype-count ${genotypes} \
		--min-base-quality-score 0
fi

gatk --java-options "-Xmx3g -Xms3g" VariantFiltration \
	-R ${reference} \
	-V ${TMPDIR}/${sample}.${ploidy}.g.vcf.gz \
	-O ${TMPDIR}/${sample}.${ploidy}_filtering_1.g.vcf.gz \
	--filter-expression "QUAL < 1000.0"  \
	--filter-name "QUAL1000"  

#Filter by the ones that PASS the filter of QUAL > 1000 on DNA, but not in RNA
if [[ $type == "RNA" ]]; then
#If RNA we don't need the filter on quality because we will consider only the variants that were also found in DNA
	bcftools view ${TMPDIR}/${sample}.${ploidy}_filtering_1.g.vcf.gz > ${TMPDIR}/${sample}.${ploidy}_filtering.g.vcf.gz
else #DNA
#If DNA we use a filter on quality.
	bcftools view -f PASS ${TMPDIR}/${sample}.${ploidy}_filtering_1.g.vcf.gz > ${TMPDIR}/${sample}.${ploidy}_filtering.g.vcf.gz
fi


#Removing PL from the vcf, this is necessary for bcftools norm
bcftools annotate \
    	-x ^INFO/DP,^FORMAT/GT,^FORMAT/AD,^FORMAT/DP,^FORMAT/GQ \
    	${output_folder}/${sample}.${ploidy}_filtering.g.vcf.gz > ${TMPDIR}/${sample}.${ploidy}_filtering.g_wo_PL.vcf

#Remove alternative alleles that have not been called from g.vcf. Sometime GATK adds a variant in the ALT column when it has not been called
#-O z is to get the output as .gz
bcftools view --trim-alt-alleles ${TMPDIR}/${sample}.${ploidy}_filtering.g_wo_PL.vcf -O z -o ${TMPDIR}/${sample}.${ploidy}_trim.vcf.gz 

#Split multiallelic positions into different rows of "biallelic" positions. If we add the reference fasta, reference alleles are adjusted so all donors call them in the same way
#-O z is to get the output as .gz.
bcftools norm \
        -m -any \
        --old-rec-tag INFO \
	--fasta-ref $reference \
        ${TMPDIR}/${sample}.${ploidy}_trim.vcf.gz -O z --threads $(nproc) > ${output_folder}/${sample}.vcf.gz

bcftools index -t ${output_folder}/${sample}.vcf.gz #Creating index in .tbi format (-t), we will need index if we plan to use bcftools merge
        	

echo "$bam has finished at $(date)"
