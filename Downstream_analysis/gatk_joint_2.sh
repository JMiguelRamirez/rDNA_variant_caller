#!/bin/bash

## GATk version
#gatk=4.5.0.0

# load modules
module load java-openjdk gatk bcftools

input_folder_og=$1
ploidy=$2

genotypes=$(python3 -c "import math; num_alleles=6; ploidy=$ploidy; maxgen=math.factorial(num_alleles + ploidy - 1) // (math.factorial(num_alleles - 1) * math.factorial(ploidy)); print(maxgen)")

reference=/gpfs/projects/bsc83/Data/assemblies/T2T_CHM13/chrR/Human_hs1-rDNA_genome_v1.0/hs1-rDNA_v1.0.fa



# GenomicsDBImport https://gatk.broadinstitute.org/hc/en-us/articles/5358869876891-GenomicsDBImport
# Import single-sample GVCFs into GenomicsDB before joint genotyping
# creating sample map file

input_folder="${input_folder_og}/var/jobs/Joint_all/"
output_folder="${input_folder_og}/var/jobs/Joint_called_all/"
mkdir -p $output_folder

threads=112


echo "Starting $input_folder_og at $(date) for ${ploidy}"
	
	#I save with input_folder_og because in parallel maybe different jobs access the same memory and hence, the same object. This is to make sure the files are not edited
for i in `find ${input_folder}/ -size +0c -name *.${ploidy}.g.vcf.gz`;do file=`basename $i`; sample="$(cut -d'.' -f1 <<<"$file")" ;echo -e ${sample}"\t"${i} >> ${output_folder}/simulated_WGS.${input_folder_og}.${ploidy}.g.vcf.map;done
	
# Provide sample map with sample IDs and paths to each g.vcf to consolidate all samples before performing joint genotyping

	# NOTE: for big ploidy values (such as 100), you might have to increase the buffer size. 
gatk --java-options "-Xmx200g -Xms100g"  GenomicsDBImport \
	--sample-name-map ${output_folder}/simulated_WGS.${input_folder_og}.${ploidy}.g.vcf.map \
	--genomicsdb-workspace-path ${output_folder}/simulated_WGS.${input_folder_og}.${ploidy}.database \
	--max-num-intervals-to-import-in-parallel ${threads} \
	--batch-size 50 \
	--intervals /gpfs/projects/bsc83/Projects/ribosomal_RNAs/Jose/07_RepeatMasker/data/pre-rRNA_47S.regions.bed \
	--genomicsdb-vcf-buffer-size 20000000 \
	--merge-input-intervals
	
#For ploidy 40 and 50 we need to increase buffer size of GenomicsDBImport: from --genomicsdb-vcf-buffer-size  16384000 to --genomicsdb-vcf-buffer-size  163840000 and from -Xmx32g -Xms30g to -Xmx200g -Xms100g and add --merge-input-intervals. But the next step does not work with ploidy 50 because of buffer overflow but I could not fix that
gatk --java-options "-Xmx32g -Xms30g" GenotypeGVCFs \
	-R ${reference} \
	-V gendb://${output_folder}/simulated_WGS.${input_folder_og}.${ploidy}.database \
	-O ${output_folder}/simulated_WGS.${input_folder_og}.${ploidy}.vcf.gz \
	--sample-ploidy ${ploidy} \
	--max-genotype-count ${genotypes} \
	--max-alternate-alleles 6

rm ${output_folder}/simulated_WGS.${input_folder_og}.${ploidy}.g.vcf.map 
rm -r ${output_folder}/simulated_WGS.${input_folder_og}.${ploidy}.database

#Separating SNPs
gatk SelectVariants \
   	-V ${output_folder}/simulated_WGS.${input_folder_og}.${ploidy}.vcf.gz \
    	-select-type SNP \
    	-O ${output_folder}/snps.${input_folder_og}.${ploidy}.vcf.gz
    	
#Separating indels and mixed:   
gatk SelectVariants \
   	-V ${output_folder}/simulated_WGS.${input_folder_og}.${ploidy}.vcf.gz \
   	--select-type-to-exclude SNP \
   	-O ${output_folder}/indels.${input_folder_og}.${ploidy}.vcf.gz
   
   	
#Following GATK recommendations
gatk --java-options "-Xmx3g -Xms3g" VariantFiltration \
	-R ${reference} \
	-V ${output_folder}/snps.${input_folder_og}.${ploidy}.vcf.gz \
	-O ${output_folder}/snps_filtered.${input_folder_og}.${ploidy}.vcf.gz \
	-filter "QD < 2.0" --filter-name "QD2" \
   	-filter "QUAL < 30.0" --filter-name "QUAL30" \
    	-filter "SOR > 3.0" --filter-name "SOR3" \
    	-filter "FS > 60.0" --filter-name "FS60" \
    	-filter "MQ < 40.0" --filter-name "MQ40" \
    	-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    	-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8"
	
gatk --java-options "-Xmx3g -Xms3g" VariantFiltration \
	-R ${reference} \
	-V ${output_folder}/indels.${input_folder_og}.${ploidy}.vcf.gz \
	-O ${output_folder}/indels_filtered.${input_folder_og}.${ploidy}.vcf.gz \
    	-filter "QD < 2.0" --filter-name "QD2" \
    	-filter "QUAL < 30.0" --filter-name "QUAL30" \
    	-filter "FS > 200.0" --filter-name "FS200" \
    	-filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" 
    	
bcftools concat \
        ${output_folder}/snps_filtered.${input_folder_og}.${ploidy}.vcf.gz \
        ${output_folder}/indels_filtered.${input_folder_og}.${ploidy}.vcf.gz \
        -o ${output_folder}/Simulated.${ploidy}.vcf.gz


rm ${output_folder}/simulated_WGS.${input_folder_og}.${ploidy}.vcf.gz ${output_folder}/snps.${input_folder_og}.${ploidy}.vcf.gz ${output_folder}/indels.${input_folder_og}.${ploidy}.vcf.gz ${output_folder}/snps_filtered.${input_folder_og}.${ploidy}.vcf.gz ${output_folder}/indels_filtered.${input_folder_og}.${ploidy}.vcf.gz

bcftools annotate \
    	-x ^INFO/DP,^FORMAT/GT,^FORMAT/AD,^FORMAT/DP,^FORMAT/GQ \
    	${output_folder}/Simulated.${ploidy}.vcf.gz > ${output_folder}/Simulated.${ploidy}_filtering.wo_PL.vcf

#Remove alternative alleles that have not been called from g.vcf. Sometime GATK adds a variant in the ALT column when it has not been called
#-O z is to get the output as .gz. Maybe not necessary for joint genotyping
bcftools view --trim-alt-alleles ${output_folder}/Simulated.${ploidy}_filtering.wo_PL.vcf -O z -o ${output_folder}/Simulated.${ploidy}_trim.vcf.gz 

bcftools norm \
        -m -any \
        --old-rec-tag INFO \
	--fasta-ref /gpfs/projects/bsc83/Data/assemblies/T2T_CHM13/chrR/Human_hs1-rDNA_genome_v1.0/hs1-rDNA_v1.0.fa \
        ${output_folder}/Simulated.${ploidy}_trim.vcf.gz  > ${output_folder}/Simulated.${ploidy}_filtering_norm.vcf

bcftools view -I ${output_folder}/Simulated.${ploidy}_filtering_norm.vcf -O z -o ${output_folder}/Simulated.${ploidy}_filtering_norm.vcf.gz

rm ${output_folder}/Simulated.${ploidy}.vcf.gz ${output_folder}/Simulated.${ploidy}_filtering.wo_PL.vcf ${output_folder}/Simulated.${ploidy}_trim.vcf.gz ${output_folder}/Simulated.${ploidy}_filtering_norm.vcf

echo "${output_folder}/Simulated.${ploidy}_filtering_norm.vcf.gz"
echo "Finished $input_folder_og at $(date) for ploidy ${ploidy}"
