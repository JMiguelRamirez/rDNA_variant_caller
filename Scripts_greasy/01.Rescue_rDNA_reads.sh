#!/bin/bash


#INPUT : fastq
#Step 1   (using seqkit)        : Convert fastq file to fasta (nucmer takes fasta as input) 
#Step 2   (using nucmer)        : Extract candidate rDNA reads (reads with a 30nt exact match with a provided reference sequence)
#Step 2.2   (using show-coords)   : Convert delta file generated by step 2 to tab-seperated readable coords file
#Step 3   (using bash)          : coord file column 13 contains read name, extract all read names for all reads with a 30nt match and keep PE for which both pairs have a 3nt match
#Step 4   (using seqkit)        : use readIDS from step 4 to extract those reads from original fastq files

# load modules
module load mummer
module load anaconda
source activate seqkit

sample_id=$1
input_folder=$2
output_folder=$3 
species=$4
length=30

#This takes less than 1 hour in general. With this mice data it takes 1:30
echo $sample_id

# input fq files
fastq1=${input_folder}/${sample_id}_1.fastq.gz
fastq2=${input_folder}/${sample_id}_2.fastq.gz

# output fq file
out_fastq1=${output_folder}/${sample_id}_1.rDNA_reads.fastq.gz
out_fastq2=${output_folder}/${sample_id}_2.rDNA_reads.fastq.gz

# Step 1. We store the fasta in a temporal directory
#-j is for numbers of threads, increasing this value might not be faster
#This takes around 15 minutes
echo "start: converting fastq to fasta - $(date)"
seqkit fq2fa ${fastq1} -j 16 > ${TMPDIR}/${sample_id}_1.fasta
seqkit fq2fa ${fastq2} -j 16 > ${TMPDIR}/${sample_id}_2.fasta
echo "end: converting fastq to fasta - $(date)"
echo ""

# Step 2
if [[ $species == "mouse" ]]; then
	echo "Running nucmer to keep candidate rDNA reads"
	rDNA_reference=data/chrR.fa
elif [[ $species == "human" ]]; then
	rDNA_reference=/gpfs/projects/bsc83/Data/assemblies/T2T_CHM13/chrR/Human_hs1-rDNA_genome_v1.0/chroms/hs1-rDNA_v1.0.chrR.47S_pre-rRNA.500padded.fa
	#The reference is adapted from this paper: https://www.jbc.org/article/S0021-9258(23)01794-5/fulltext
else
	echo "Species not available. Try mouse or human"
fi



echo "start: finding 30 nt matches: $(date)"
#This takes 15 minutes
# 2.1
#FLAGS used :
# --maxmatch      Use all anchor matches regardless of their uniqueness
# -l  minimum length of a single exact match
# --delta=PATH    Output delta file to PATH (instead of PREFIX.delta)
nucmer --maxmatch -l ${length} --delta=${TMPDIR}/${sample_id}_1.nucmer.rDNA.delta --threads=112 ${rDNA_reference} ${TMPDIR}/${sample_id}_1.fasta
nucmer --maxmatch -l ${length} --delta=${TMPDIR}/${sample_id}_2.nucmer.rDNA.delta --threads=112 ${rDNA_reference} ${TMPDIR}/${sample_id}_2.fasta

# 2.2
# delta file can be used to extract coordinates and alignment data table, and also to see alignments, but for analysis main file is .coord file
#FLAGS used :
#-r      Sort output lines by reference
#-c      Include percent coverage columns in the output
#-l      Include sequence length columns in the output
#-T      Switch output to tab-delimited for
show-coords -r -c -l -T ${TMPDIR}/${sample_id}_1.nucmer.rDNA.delta > ${TMPDIR}/${sample_id}_1.nucmer.rDNA.coords
show-coords -r -c -l -T ${TMPDIR}/${sample_id}_2.nucmer.rDNA.delta > ${TMPDIR}/${sample_id}_2.nucmer.rDNA.coords
echo "end: finding  reads with 30nt exact match - $(date)"
echo ""

# Step 3
echo "getting read ids - $(date)"
#It takes 1 min
# first five lines are headers
tail -n +5 ${TMPDIR}/${sample_id}_1.nucmer.rDNA.coords | awk '{print $13}' > ${TMPDIR}/${sample_id}.R1.out
tail -n +5 ${TMPDIR}/${sample_id}_2.nucmer.rDNA.coords | awk '{print $13}' > ${TMPDIR}/${sample_id}.R2.out

# only keep reads for which both end of the read have a 30nt match
# hs1 (T2T-CHM13) WGS data does not include /1 and /2 in read names
# GTEx read IDs include /1 and /2 -> remove to keep read-pair ids
awk 'NR==FNR{a[$1]=$1; next } {if (a[$1]) {print} }' ${TMPDIR}/${sample_id}.R1.out ${TMPDIR}/${sample_id}.R2.out > ${TMPDIR}/${sample_id}.uniq_read_pair_IDs.out

# Step 4
echo "start: extracting reads with match from fasta - $(date)"
#It takes 10 min
seqkit grep -f ${TMPDIR}/${sample_id}.uniq_read_pair_IDs.out ${fastq1} -j 30 | gzip > ${out_fastq1}
seqkit grep -f ${TMPDIR}/${sample_id}.uniq_read_pair_IDs.out ${fastq2} -j 30 | gzip > ${out_fastq2}
echo "end: extracting reads with match from fasta - $(date)"
echo ""

