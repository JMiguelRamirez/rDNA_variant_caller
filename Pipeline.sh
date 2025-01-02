#!/bin/bash

#Help message
Help()
{
   # Display Help
   echo "This code runs a variant caller for rDNA reads"
   echo
   echo "Syntax: Pipeline.sh [-h|-t|-n|-f|-i|-o]"
   echo "parameters:"
   echo "h     Print the help message."
   echo "t     Type of input data: Either DNA or RNA."
   echo "n     Basename of the input data: (e.g., for a bam file named sample.bam: -n sample, for fastqs named sample_1.fastq.gz and sample_2.fasq.gz: -n sample. This means that paired fastqs should end by _1.fastq.gz and _2.fastq.gz)"
   echo "f     Type of input data format: Either bam (for bam/cram/sam) or fastq."
   echo "i     Input folder."
   echo "o     Output folder."
   echo
}



# process the parameters used as input for the pipeline:
while getopts 'ht:n:f:i:o:' OPTION; do
   case "$OPTION" in
	h) # Echo Help
	Help
	exit
	;;
	t) # Input data type: Either DNA or RNA
	if [[ "$OPTARG" != "DNA" && "$OPTARG" != "RNA" ]]; then
		echo "Error: -t option must be either 'DNA' or 'RNA'."
		exit 1
	fi
        type="$OPTARG"
        ;;
        n) # base name of the input files
        basename="$OPTARG"
        ;;
	f) # Input data format: Either BAM or FASTQ
	if [[ "$OPTARG" != "bam" && "$OPTARG" != "fastq" ]]; then
		echo "Error: -t option must be either 'bam' or 'fastq'. Option bam works for bam/cram/sam files."
		exit 1
	fi
        format="$OPTARG"
        ;;
	i) # base name of the input files
        input="$OPTARG"
        ;;
	o) # base name of the input files
        output="$OPTARG"
        ;;
	?) # Invalid option
	echo "Error: Invalid parameters. Check -h for help"
	exit
	;;
   esac
done

echo "Preprocessing the data and mapping to our custom genome reference"
#I have another folder named Reference with the following information:
reference="Reference/Human_hs1-rDNA_genome_v1.0/hs1-rDNA_v1.0.fa" #reference from https://github.com/vikramparalkar/rDNA-Mapping-Genomes: Human_hs1. + .fa.fai needed for gatk (samtools faidx) + .fa.dict (gatk-launch CreateSequenceDictionary -R reference.fa)
index="Reference/new_index" #bwa index of the reference sequence (takes ~2 hours to obtain from the reference)
	
if [[ $type == "DNA" ]]; then
	echo "Running nucmer to keep candidate rDNA reads"
	if [[ "$format" == "fastq" ]]; then
		echo "Note: make sure the fastq file does not contain adaptors. If so, use trimmomtatric to remove them and fastqc to check it"
		#If we are using DNA we should run mummer to retrieve candidate reads, this speeds up the process and removes reads coming from pseudogenes
		./Scripts/Rescue_rDNA_reads_from_fastq_sequential.sh $basename $input $output
	elif [[ "$format" == "bam" ]]; then
		echo "Retrieving candidate rDNA reads from the bam file. The bam must be mapped to a reference genome, and not transcriptome."
		./Scripts/Rescue_rDNA_reads_from_bam_sequential.sh $basename $input $output
	fi
	echo "Mapping the fastq file to the reference genome (human + chrR)"
	./Scripts/Mapping_sequential.sh $basename $output $output $reference $index
	#We now have the output as: ${output}/${basename}.sorted.chrR.f2F2308q20.wo_XA.bam
	
elif [[ $type == "RNA" ]]; then
	#The steps are the same as with DNA but removing the nucmer step, as we can assume that the pseudogenes are not expressed
	if [[ "$format" == "fastq" ]]; then
		echo "Note: make sure the fastq file does not contain adaptors. If so, use trimmomtatric to remove them and fastqc to check it"
		echo "Mapping the fastq file to the reference genome (human + chrR)"
		./Scripts/Mapping_sequential.sh $basename $input $output $reference $index
		#We now have the output as: ${output}/${basename}.sorted.chrR.f2F2308q20.wo_XA.bam

	elif [[ "$format" == "bam" ]]; then
		echo "Retrieving candidate rDNA reads from the bam file and mapping to our reference. The bam must be mapped to a reference genome, and not transcriptome."

		./Scripts/Rescue_rDNA_reads_from_bam_RNA_sequential.sh $basename $input $output
		#We now have the output as: ${output}/${basename}.sorted.chrR.f2F2308q20.wo_XA.bam
	fi
fi


echo "Running the variant caller"
#At this point we have a .bam file, we can do the variant calling
./Scripts/Variant_Calling.sh $basename $output $output 20 $type $reference #20 is the ploidy




