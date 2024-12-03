# rDNA_variant_caller
Repository to call rDNA variants from short read sequencing in parallel using the software: greasy
Greasy saves X nodes and keeps submitting "jobs" when the X nodes have space. This does better use of resources when there are many samples and tends to be faster than using arrays.

# Software dependencies

| Package | Version | 
| -------- | ------- | 
| MUMMER | 4.0.0rc1 | 
| ANACONDA | 2023.07 | 
| seqkit | 2.8.0 |
| Samtools | 1.19.2 |
| bwa | 0.7.17 |
| gatk | 4.5.0.0 |
| java-openjdk | 17.0.11+9|
| bcftools | 1.19 |
| R | 4.3.2 |

# Running the code:

To run the code in parallel using greasy, we just need a .txt file with all the sample IDs (file names): IDs.tab. If a sample is called SRA123.cram, the ID should be SRA123. All .out and .err files will be stored in a folder named out/ in the source file location

Step 1A: Retrieving candidate rDNA reads from fastq files. If the input is RNA, we can skip this script that calls MUMmer. Available for human and mouse
```
rm 01.greasy_file_rescue_rDNA_reads.txt
for item in $(cat IDs.tab); do
	echo -e "Scripts_greasy/01.rescue_rDNA_reads_greasy.sh\t$item\t$input_folder\t$output_folder\thuman" >> 01.greasy_file_rescue_rDNA_reads.txt
done
sbatch --ntasks=160 --cpus-per-task=56 01.to_submit_rescue_rDNA_reads_greasy.sh
```
Step 1B: Retrieving candidate rDNA reads from cram/bam files. If the input is RNA, we can skip this script that calls MUMmer:
```
rm 01.greasy_file_rescue_rDNA_reads_bam.txt
for item in $(cat IDs.tab); do
	echo -e "Scripts_greasy/01.Rescue_rDNA_reads_bam.sh\t$item\t$input_folder\t$output_folder\thuman" >> 01.greasy_file_rescue_rDNA_reads_bam.txt
done
sbatch --ntasks=160 --cpus-per-task=56 01.to_submit_rescue_rDNA_reads_bam_greasy.sh
```
Step 2: Mapping candidate rDNA reads to our custom reference (it has a bwa index):
```
rm 02.greasy_file_mapping.txt
for item in $(cat IDs.tab); do
	echo -e "Scripts_greasy/02.mapping.sh\t$item\t$input_folder\t$output_folder\thuman" >> 02.greasy_file_mapping.txt
done
sbatch --ntasks=160 --cpus-per-task=56 02.to_submit_mapping_greasy.sh
```

Step 3: Variant calling. 20 refers to the ploidy we want to use:
```
rm 03.greasy_file_gatk.txt
for item in $(cat IDs.tab); do
	echo -e "Scripts_greasy/03.gatk_greasy.sh\t$item\t$input_folder\t$output_folder\t20" >> 03.greasy_file_gatk.txt
done
sbatch --ntasks=160 --cpus-per-task=56 03.to_submit_gatk_greasy.sh

#For RNA:
rm 03.greasy_file_gatk.txt
for item in $(cat IDs.tab); do
	echo -e "Scripts_greasy/03.gatk_greasy.sh\t$item\t$input_folder\t$output_folder\t20\tRNA" >> 03.greasy_file_gatk.txt
done
sbatch --ntasks=160 --cpus-per-task=56 03.to_submit_gatk_greasy.sh
```

Step 4: merging the vcfs from the different samples into one:
```
Rscript Scripts_greasy/04.merge_gVCF.R IDs.tab output/
```

# Running the code for RNA: