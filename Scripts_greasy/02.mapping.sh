#!/bin/bash


# load modules
module load bwa/0.7.17 samtools

sample_id=$1
input_folder=$2
output_folder=$3
species=$4
threads=112


# reference genome info
if [[ $species == "mouse" ]]; then
	reference=/gpfs/projects/bsc83/Data/assemblies/Mm39/Mouse_mm39-rDNA_genome_v1.0/mm39-rDNA_v1.0.fa
	bwa_index=/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Jose/11_Mice/new_index
elif [[ $species == "human" ]]; then
	reference=/gpfs/projects/bsc83/Data/assemblies/T2T_CHM13/chrR/Human_hs1-rDNA_genome_v1.0/hs1-rDNA_v1.0.fa
	bwa_index=/gpfs/projects/bsc83/Projects/ribosomal_RNAs/Jose/04_Pipeline/new_index
	#The reference comes from this paper: https://www.jbc.org/article/S0021-9258(23)01794-5/fulltext 
	#the bwa index comes from running bwa mem index on the reference fasta
else
	echo "Species not available. Try mouse or human"
fi


# fq files
fastq1=${input_folder}/${sample_id}_1.rDNA_reads.fastq.gz
fastq2=${input_folder}/${sample_id}_2.rDNA_reads.fastq.gz

# for MAGE
#fastq1=${input_folder}/${sample_id}_r1.fq.gz
#fastq2=${input_folder}/${sample_id}_r2.fq.gz

# bwa mapping 
#-t is number of threads
#-h INT[,INT]  if there are <INT hits with score >80% of the max score, output all in XA [0,0]

#the final - from samtools sort means that it takes as input the standard input instead of a specified input, so the result from the pipe, this is only necessary for samtools
bwa mem -t ${threads} -h 1000 -R $(echo "@RG\tID:${sample_id}\tSM:${sample_id}\tLB:lib1\tPL:ILLUMINA\tPU:unit1") ${bwa_index} ${fastq1} ${fastq2} | samtools sort -T ${TMPDIR}/${sample_id} -@ ${threads} -o ${TMPDIR}/${sample_id}.sorted.bam -
#-@ for samtools means the number of threads
samtools index ${TMPDIR}/${sample_id}.sorted.bam -@ ${threads}
                                                   
                                                   
# filter bam
samtools view -hb -f 2 -F 2308 -q 20 ${TMPDIR}/${sample_id}.sorted.bam chrR > ${TMPDIR}/${sample_id}.sorted.chrR.f2F2308q20.bam
samtools index ${TMPDIR}/${sample_id}.sorted.chrR.f2F2308q20.bam

# exclude multimapped reads with primary alignment to chrR
samtools view -h -f 2 ${TMPDIR}/${sample_id}.sorted.chrR.f2F2308q20.bam | grep -v "XA:" | samtools view -hb > ${output_folder}/${sample_id}.sorted.chrR.f2F2308q20.wo_XA.bam
samtools index ${output_folder}/${sample_id}.sorted.chrR.f2F2308q20.wo_XA.bam


#Flagstat to get number of reads
samtools flagstat ${output_folder}/${sample_id}.sorted.chrR.f2F2308q20.wo_XA.bam -O tsv > ${output_folder}/${sample_id}.sorted.chrR.f2F2308q20.wo_XA.flagstat                                           
                                                                                         
