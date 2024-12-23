#!/bin/bash


# load modules
###Everything with ### was used in my hpc environment, else, everything would already be loaded with the singularity image
###module load samtools anaconda
###source activate seqkit

# sample (GTEx tissue sample )
#file=$2
#donor=$1
#tissue=$3
sample=$1
file=$sample
#sample=${donor}_${tissue}
input=$2
#input=/gpfs/scratch/bsc83/MN4/bsc83/bsc83535/GTEx/v8/bam_files/${tissue}/cram_files/
output=$3
#output=GTEx/


echo "Retrieving paired reads mapped to rDNA-like regions for sample $sample"
#threads=2
threads=4

# infiles
###reference=/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Jose/04_Pipeline/data/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta 
reference=Reference/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta 
###bed_file=/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Jose/04_Pipeline/data/hg38_rDNA-like_regions.bed
bed_file=Reference/hg38_rDNA-like_regions.bed

# paths
# path to RNA-seq cram files
cram_file=${input}/${file}.Aligned.sortedByCoord.out.patched.md.cram

#inpath=$output/${donor}/${tissue}
inpath=$input
#outpath=$output/${donor}/${tissue}
outpath=$output
# Check directory if it does not exits
if [[ ! -d "${outpath}" ]]; then
  echo "Create directory: ${outpath}"
  mkdir -p "${outpath}"
fi

# 1. cram2fastq rDNA-like genomic regions
samtools view -@ ${threads} -h -T ${reference} -L ${bed_file} ${cram_file} | samtools fastq -@ ${threads} -N --reference ${reference} -1 ${TMPDIR}/${sample}.rDNA-like.R1.fastq.gz -2 ${TMPDIR}/${sample}.rDNA-like.R2.fastq.gz 

#Sorting fastqs:
seqkit sort ${TMPDIR}/${sample}.rDNA-like.R1.fastq.gz | gzip > ${TMPDIR}/${sample}.R1.sorted.fastq.gz
seqkit sort ${TMPDIR}/${sample}.rDNA-like.R2.fastq.gz | gzip > ${TMPDIR}/${sample}.R2.sorted.fastq.gz


echo "Mapping to chrR $sample"
# 6. map rDNA reads
###module load bwa/0.7.17

# reference genome info
###reference=/gpfs/projects/bsc83/Data/assemblies/T2T_CHM13/chrR/Human_hs1-rDNA_genome_v1.0/hs1-rDNA_v1.0.fa
reference=Reference/Human_hs1-rDNA_genome_v1.0/hs1-rDNA_v1.0.fa
###bwa_index=/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Jose/04_Pipeline/new_index
bwa_index=Reference/new_index

#I only want paired reads. So I will create .out files with the list of IDs
seqkit fx2tab ${TMPDIR}/${sample}.R1.sorted.fastq.gz | cut -f 1 | sort > ${TMPDIR}/${sample}.R1.out
seqkit fx2tab ${TMPDIR}/${sample}.R2.sorted.fastq.gz | cut -f 1 | sort > ${TMPDIR}/${sample}.R2.out

#Get list of paired reads
# GTEx read IDs include /1 and /2 -> remove to keep read-pair ids
awk 'NR==FNR{a[$1]=$1; next } {if (a[$1]) {print} }' <(awk '{OFS="\t"; split($1,a,"/");print a[1]}' ${TMPDIR}/${sample}.R1.out) <(awk '{OFS="\t"; split($1,a,"/");print a[1]}' ${TMPDIR}/${sample}.R2.out) > ${TMPDIR}/${sample}.uniq_read_pair_IDs.out
sed 's/$/\/1/'  ${TMPDIR}/${sample}.uniq_read_pair_IDs.out > ${TMPDIR}/${sample}.R1_read_ids.out
sed 's/$/\/2/'  ${TMPDIR}/${sample}.uniq_read_pair_IDs.out > ${TMPDIR}/${sample}.R2_read_ids.out

#Getting fastqs with only paired reads
seqkit grep -f ${TMPDIR}/${sample}.R1_read_ids.out ${TMPDIR}/${sample}.R1.sorted.fastq.gz -j 16 | gzip > ${outpath}/${sample}_1.rDNA_reads.fastq.gz
seqkit grep -f ${TMPDIR}/${sample}.R2_read_ids.out ${TMPDIR}/${sample}.R2.sorted.fastq.gz -j 16 | gzip > ${outpath}/${sample}_2.rDNA_reads.fastq.gz

fastq1=${outpath}/${sample}_1.rDNA_reads.fastq.gz
fastq2=${outpath}/${sample}_2.rDNA_reads.fastq.gz

bwa mem -t ${threads} -h 1000 -R $(echo "@RG\tID:${sample}\tSM:${sample}\tLB:lib1\tPL:ILLUMINA\tPU:unit1") ${bwa_index} ${fastq1} ${fastq2} | samtools sort -T ${outpath}/${sample} -o ${outpath}/${sample}.sorted.bam -
samtools index ${outpath}/${sample}.sorted.bam 

# filter bam
samtools view -hb -f 2 -F 2308 -q 20  ${outpath}/${sample}.sorted.bam chrR > ${outpath}/${sample}.sorted.chrR.f2F2308q20.bam
samtools index ${outpath}/${sample}.sorted.chrR.f2F2308q20.bam

# exclude multimapped reads with primary alignment to chrR
samtools view -h -f 2 ${outpath}/${sample}.sorted.chrR.f2F2308q20.bam | grep -v "XA:" | samtools view -hb > ${outpath}/${sample}.sorted.chrR.f2F2308q20.wo_XA.bam
samtools index ${outpath}/${sample}.sorted.chrR.f2F2308q20.wo_XA.bam

# samtools flagstat
samtools flagstat ${outpath}/${sample}.sorted.bam -O tsv > ${outpath}/${sample}.sorted.flagstat
samtools flagstat ${outpath}/${sample}.sorted.chrR.f2F2308q20.bam -O tsv > ${outpath}/${sample}.sorted.chrR.f2F2308q20.flagstat
samtools flagstat ${outpath}/${sample}.sorted.chrR.f2F2308q20.wo_XA.bam -O tsv > ${outpath}/${sample}.sorted.chrR.f2F2308q20.wo_XA.flagstat

echo "The mapping has finished for $sample"
#We could add here the gatk step. We still need to decide
