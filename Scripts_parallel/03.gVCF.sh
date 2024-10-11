#!/bin/bash

#SBATCH --job-name=gatk
#SBATCH --output=./out/gatk.%A_%a.out
#SBATCH --error=./out/gatk.%A_%a.err
#SBATCH --cpus-per-task=56
#SBATCH --account=bsc83
#SBATCH --qos=gp_bscls
#SBATCH --time=00:30:00

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
reference=/gpfs/projects/bsc83/Data/assemblies/T2T_CHM13/chrR/Human_hs1-rDNA_genome_v1.0/hs1-rDNA_v1.0.fa
bam=${input_folder}/${sample}.sorted.chrR.f2F2308q20.wo_XA.bam

#We iterate over regions in case we need to parallelize. Usually one one whole bed including all regions should be fine. But in some cases we need to send different jobs per region. This is an intermediate way easily adaptable to both cases
for region in ITS1 ITS2 18S 5.8S 28S 5_ETS 3_ETS;do
	echo $region
	gatk --java-options "-Xmx115g -Xms100g" HaplotypeCaller \
	-I ${bam} \
	-R ${reference} \
	-O ${TMPDIR}/${sample}.${ploidy}.${region}.g.vcf.gz \
	--sample-ploidy ${ploidy} \
	--intervals /gpfs/projects/bsc83/Projects/ribosomal_RNAs/Alba/data/pre-rRNA_47S.included_intervals.${region}.bed \
	--max-reads-per-alignment-start 0 \
	--dont-use-soft-clipped-bases true \
	--native-pair-hmm-threads 112 \
	--max-genotype-count ${genotypes} 
	
	gatk --java-options "-Xmx3g -Xms3g" VariantFiltration \
	-R ${reference} \
	-V ${TMPDIR}/${sample}.${ploidy}.${region}.g.vcf.gz \
	-O ${TMPDIR}/${sample}.${ploidy}.${region}_filtering.g.vcf.gz \
	--filter-expression "QUAL < 1000.0"  \
	--filter-name "QUAL1000"   
done

#Merging the .vcf of all regions
bcftools concat \
        ${TMPDIR}/${sample}.${ploidy}.5_ETS_filtering.g.vcf.gz \
        ${TMPDIR}/${sample}.${ploidy}.18S_filtering.g.vcf.gz \
        ${TMPDIR}/${sample}.${ploidy}.ITS1_filtering.g.vcf.gz \
        ${TMPDIR}/${sample}.${ploidy}.5.8S_filtering.g.vcf.gz \
        ${TMPDIR}/${sample}.${ploidy}.ITS2_filtering.g.vcf.gz \
        ${TMPDIR}/${sample}.${ploidy}.28S_filtering.g.vcf.gz \
        ${TMPDIR}/${sample}.${ploidy}.3_ETS_filtering.g.vcf.gz \
        -o ${TMPDIR}/${sample}.${ploidy}_filtering_p.g.vcf.gz 
bcftools view -f PASS ${TMPDIR}/${sample}.${ploidy}_filtering_p.g.vcf.gz > ${TMPDIR}/${sample}.${ploidy}_filtering.g.vcf.gz


#Split multiallelic positions into biallelic positions. By adding the reference fasta, reference alleles are adjusted to be consistent across callings in different donors
bcftools annotate \
    	-x ^INFO/DP,^FORMAT/GT,^FORMAT/AD,^FORMAT/DP,^FORMAT/GQ \
    	${TMPDIR}/${sample}.${ploidy}_filtering.g.vcf.gz > ${TMPDIR}/${sample}.${ploidy}_filtering.g_wo_PL.vcf.gz

#This normalization splits multiallelic into different rows of biallelic. It also "normalizes" indes to be comparable across samples, thanks to adding the reference genome, the alternative alleles will be less interpretable, as they now need to be paired with the reference to be interpretable. For example this: 

#22309	.	CCGCGCGCG	C,CCGCG,CCGCGCG,CCGCGCGCGCGCG will now be

# chrR  22309	.	CCGCGCGCG	C	16376.2	PASS	DP=2639;INFO=chrR|22309|CCGCGCGCG|C,CCGCG,CCGCGCG,CCGCGCGCGCGCG|1	GT:AD:DP:GQ	0/0/0/0/0/0/0/0/0/0/0/0/0/1/1/1/0/0/0/0:1129,288:1777:58
# chrR	22309	.	CCGCG	C	16376.2	PASS	DP=2639;INFO=chrR|22309|CCGCGCGCG|C,CCGCG,CCGCGCG,CCGCGCGCGCGCG|2	GT:AD:DP:GQ	0/0/0/0/0/0/0/0/0/0/0/0/0/0/0/0/1/0/0/0:1129,129:1777:58
# chrR	22309	.	CCG	C	16376.2	PASS	DP=2639;INFO=chrR|22309|CCGCGCGCG|C,CCGCG,CCGCGCG,CCGCGCGCGCGCG|3	GT:AD:DP:GQ	0/0/0/0/0/0/0/0/0/0/0/0/0/0/0/0/0/1/0/0:1129,51:1777:58
# chrR	22309	.	C	CCGCG

bcftools norm \
        -m -any \
        --old-rec-tag INFO \
	--fasta-ref /gpfs/projects/bsc83/Data/assemblies/T2T_CHM13/chrR/Human_hs1-rDNA_genome_v1.0/hs1-rDNA_v1.0.fa \
        ${TMPDIR}/${sample}.${ploidy}_filtering.g_wo_PL.vcf.gz > ${TMPDIR}/${sample}.${ploidy}_filtering_norm.vcf




#If -ERC base_resolution, which maybe I would revisit to merge and know the covarage of non-called variants, but for now this is only relevant for RNA
#awk -F'\t' '$5!="<NON_REF>"' ${output_folder}/${sample}.${ploidy}_filtering_norm.vcf.gz > ${output_folder}/${sample}.${ploidy}_filtering_norm.vcf

#Remove alternative alleles that have not been called
Rscript /gpfs/projects/bsc83/Projects/ribosomal_RNAs/Jose/04_Pipeline/05.Remove_not_called_alleles_from_norm.R ${TMPDIR}/${sample}.${ploidy}_filtering_norm.vcf ${output_folder}/${sample}.${ploidy}_filtering_norm_called.vcf.gz


#Option 2. Better but I didn't run. This could then be used with bcftools merge
#bcftools view --trim-alt-alleles ${TMPDIR}/${sample}.${ploidy}_filtering_norm.vcf -O z -o ${output_folder}/${sample}.${ploidy}_filtering_norm_called.vcf.gz 
#bcftools index -t ${output_folder}/${sample}.${ploidy}_filtering_norm_called.vcf.gz
#In RNA I don't use this either because I want to keep the information on quality per donor
	
#Other relevant lines:
#These lines are to merge the .vcfs
#bcftools annotate -x ^FORMAT/GT ${TMPDIR}/${sample}.${ploidy}_filtering_norm.vcf > ${TMPDIR}/${sample}.${ploidy}_filtering_norm_GT.vcf
#bcftools view -I ${output_folder}/${sample}.${ploidy}_filtering_norm.vcf -O z -o ${output_folder}/${sample}.${ploidy}_filtering_norm.vcf.gz

