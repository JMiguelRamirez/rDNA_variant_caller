#!/bin/bash

#SBATCH --job-name=gatk
#SBATCH --output=./out/gatk.%A_%a.out
#SBATCH --error=./out/gatk.%A_%a.err
#SBATCH --cpus-per-task=112
#SBATCH --account=bsc83

# load modules
module load java-openjdk gatk bcftools R
#module load java-openjdk/17.0.11+9
#module load gatk/4.5.0.0
       
file_tab=$1
input_folder=$2
output_folder=$3
ploidy=$4
export sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${file_tab} | cut -f1)

genotypes=$(python3 -c "import math; num_alleles=6; ploidy=$ploidy; maxgen=math.factorial(num_alleles + ploidy - 1) // (math.factorial(num_alleles - 1) * math.factorial(ploidy)); print(maxgen)")


echo $ploidy
echo $genotypes
reference=/gpfs/projects/bsc83/Data/assemblies/Mm39/Mouse_mm39-rDNA_genome_v1.0/mm39-rDNA_v1.0.fa
bam=${input_folder}/${sample}.sorted.chrR.f2F2308q20.wo_XA.bam
echo "$bam has started at $(date)"
threads=112

# The first time we need to create the .dict file for the reference
# gatk CreateSequenceDictionary R=mm39-rDNA_v1.0.fa O=mm39-rDNA_v1.0.dict

gatk --java-options "-Xmx115g -Xms100g" HaplotypeCaller \
	-I ${bam} \
	-R ${reference} \
	-O ${TMPDIR}/${sample}.${ploidy}.g.vcf.gz \
	--sample-ploidy ${ploidy} \
	--intervals ./data/pre-rRNA_47S.included_intervals.bed \
	--max-reads-per-alignment-start 0 \
	--dont-use-soft-clipped-bases true \
	--native-pair-hmm-threads ${threads} \
	--max-genotype-count ${genotypes} \
	--min-base-quality-score 0 #No needed for DNA, but needed for RNA
	 #by default gatk gets the maximum number of threads available, but we are adding it just in case
	
gatk --java-options "-Xmx3g -Xms3g" VariantFiltration \
	-R ${reference} \
	-V ${TMPDIR}/${sample}.${ploidy}.g.vcf.gz \
	-O ${TMPDIR}/${sample}.${ploidy}_filtering.g.vcf.gz \
	--filter-expression "QUAL < 1000.0"  \
	--filter-name "QUAL1000"   
	

#Only get the ones that PASS the filter of QUAL > 1000. Don't use it for RNA
bcftools view -f PASS ${TMPDIR}/${sample}.${ploidy}_filtering.g.vcf.gz > ${TMPDIR}/${sample}.${ploidy}_filtering_2.g.vcf.gz

#Removing PL from the vcf, this is necessary for bcftools norm
bcftools annotate \
    	-x ^INFO/DP,^FORMAT/GT,^FORMAT/AD,^FORMAT/DP,^FORMAT/GQ \
    	${TMPDIR}/${sample}.${ploidy}_filtering_2.g.vcf.gz > ${TMPDIR}/${sample}.${ploidy}_filtering.g_wo_PL.vcf

#Remove alternative alleles that have not been called from g.vcf. Sometime GATK adds a variant in the ALT column when it has not been called
#-O z is to get the output as .gz
bcftools view --trim-alt-alleles ${TMPDIR}/${sample}.${ploidy}_filtering.g_wo_PL.vcf -O z -o ${TMPDIR}/${sample}.${ploidy}_trim.vcf.gz 

#Split multiallelic positions into different rows of "biallelic" positions. If we add the reference fasta, reference alleles are adjusted so all donors call them in the same way
#-O z is to get the output as .gz.
bcftools norm \
        -m -any \
        --old-rec-tag INFO \
	--fasta-ref $reference \
        ${TMPDIR}/${sample}.${ploidy}_trim.vcf.gz -O z --threads ${threads} > ${output_folder}/${sample}.vcf.gz

bcftools index -t ${output_folder}/${sample}.vcf.gz #Creating index in .tbi format (-t), we will need index if we plan to use bcftools merge
        	

echo "$bam has finished at $(date)"
