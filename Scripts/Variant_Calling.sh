#!/bin/bash

# load modules
###Everything with ### was used in my hpc environment, but everything would already be loaded with the singularity image
###module load java-openjdk gatk bcftools R
       
sample=$1
input_folder=$2
output_folder=$3
ploidy=$4
type=$5
reference=$6

genotypes=$(python3 -c "import math; num_alleles=6; ploidy=$ploidy; maxgen=math.factorial(num_alleles + ploidy - 1) // (math.factorial(num_alleles - 1) * math.factorial(ploidy)); print(maxgen)")

threads=$(cat /proc/cpuinfo | grep "processor"| wc -l)

#Create a temporal folder:
TMPDIR=$output_folder/temporal
mkdir -p $TMPDIR

echo $ploidy
echo $genotypes
bam=${input_folder}/${sample}.sorted.chrR.f2F2308q20.wo_XA.bam

#GATK will use as many threads as available by default
	
if [[ $type == "RNA" ]]; then
#If RNA we output all positions, as we can have the counts for variants found in DNA but not RNA. -ERC BP_RESOLUTION
gatk --java-options "-Xmx20g -Xms20g" HaplotypeCaller \
	-I ${bam} \
	-R ${reference} \
	-O ${TMPDIR}/${sample}.${ploidy}.g.vcf.gz \
	--sample-ploidy ${ploidy} \
	--intervals Reference/Intervals/pre-rRNA_47S.included_intervals.bed \
	--max-reads-per-alignment-start 0 \
	--dont-use-soft-clipped-bases true \
	--max-genotype-count ${genotypes} \
	--min-base-quality-score 0 \
	-ERC BP_RESOLUTION
else
#If DNA we only output the variants
gatk --java-options "-Xmx20g -Xms20g" HaplotypeCaller \
	-I ${bam} \
	-R ${reference} \
	-O ${TMPDIR}/${sample}.${ploidy}.g.vcf.gz \
	--sample-ploidy ${ploidy} \
	--intervals Reference/Intervals/pre-rRNA_47S.included_intervals.bed \
	--max-reads-per-alignment-start 0 \
	--dont-use-soft-clipped-bases true \
	--max-genotype-count ${genotypes} \
	--min-base-quality-score 0 \
	--native-pair-hmm-threads ${threads}
fi

gatk --java-options "-Xmx3g -Xms3g" VariantFiltration \
	-R ${reference} \
	-V ${TMPDIR}/${sample}.${ploidy}.g.vcf.gz \
	-O ${TMPDIR}/${sample}.${ploidy}_filtering_p.g.vcf.gz  \
	--filter-expression "QUAL < 1000.0"  \
	--filter-name "QUAL1000"   

        
if [[ $type == "RNA" ]]; then
#If RNA we don't need the filter on quality because we will consider only the variants that were also found in DNA
	bcftools view ${TMPDIR}/${sample}.${ploidy}_filtering_p.g.vcf.gz > ${TMPDIR}/${sample}.${ploidy}_filtering.g.vcf
else
#If DNA we use a filter on quality.
	bcftools view -f PASS ${TMPDIR}/${sample}.${ploidy}_filtering_p.g.vcf.gz > ${TMPDIR}/${sample}.${ploidy}_filtering.g.vcf
fi

#Split multiallelic positions into biallelic positions. By adding the reference fasta, reference alleles are adjusted to be consistent across callings in different donors
bcftools annotate \
    	-x ^INFO/DP,^FORMAT/GT,^FORMAT/AD,^FORMAT/DP,^FORMAT/GQ \
    	${TMPDIR}/${sample}.${ploidy}_filtering.g.vcf > ${TMPDIR}/${sample}.${ploidy}_filtering.g_wo_PL.vcf

#Remove alternative alleles that have not been called from g.vcf. Sometime GATK adds a variant in the ALT column when it has not been called
#-O z is to get the output as .gz
bcftools view --trim-alt-alleles ${TMPDIR}/${sample}.${ploidy}_filtering.g_wo_PL.vcf -O z -o ${TMPDIR}/${sample}.${ploidy}_trim.vcf.gz 

#Split multiallelic positions into different rows of "biallelic" positions. If we add the reference fasta, reference alleles are adjusted so all donors call them in the same way
#-O z is to get the output as .gz.
bcftools norm \
        -m -any \
        --old-rec-tag INFO \
	--fasta-ref ${reference} \
        ${TMPDIR}/${sample}.${ploidy}_trim.vcf.gz -O z --threads ${threads} > ${output_folder}/${sample%.bam}.vcf.gz

bcftools index -t ${output_folder}/${sample%.bam}.vcf.gz #Creating index in .tbi format (-t), we will need index if we plan to use bcftools merge to merge different samples
        	
#rm -r $TMPDIR
