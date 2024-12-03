# rDNA_variant_caller
Repository to call rDNA variants from short read sequencing in parallel

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

To run the code in parallel using SLURM arrays, we just need a .txt file with all the sample IDs (file names): IDs.tab. If a sample is called SRA123.cram, the ID should be SRA123. All .out and .err files will be stored in a folder named out/ in the source file location

Step 1A: Retrieving candidate rDNA reads from fastq files. If the input is RNA, we can skip this script that calls MUMmer. Available for human and mouse
```
sbatch -a 1-10 Scripts_parallel/01.Rescue_rDNA_reads.sh IDs.tab input/ output/ human
```
Step 1B: Retrieving candidate rDNA reads from cram/bam files. If the input is RNA, we can skip this script that calls MUMmer:
```
sbatch -a 1-10 Scripts_parallel/01.Rescue_rDNA_reads_bam.sh IDs.tab input/ output/ human
```
Step 2: Mapping candidate rDNA reads to our custom reference (it has a bwa index):
```
sbatch -a 1-10 --time=00:30:00 Scripts_parallel/02.mapping.sh IDs.tab input/ output/ human
#DNA is the default, if using RNA:
sbatch -a 1-10 --time=00:30:00 Scripts_parallel/02.mapping.sh IDs.tab input/ output/ human RNA
```

Step 3: Variant calling
```
sbatch -a 1-10 Scripts_parallel/03.gVCF.sh IDs.tab input/ output/ human
#DNA is the default, if using RNA:
sbatch -a 1-10 Scripts_parallel/03.gVCF.sh IDs.tab input/ output/ human RNA

```

Step 4: merging the vcfs from the different samples into one:
```
#Original version
Rscript Scripts_parallel/04.merge_gVCF.R IDs.tab output/
#New version that outputs read depth (DP) and is more efficient
Rscript Scripts_parallel/04.merge_gVCF_RNA.R IDs.tab output/
```
