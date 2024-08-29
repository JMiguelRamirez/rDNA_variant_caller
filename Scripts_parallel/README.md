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

To run the code in parallel using SLURM arrays, we just need a .txt file with all the sample IDs (file names): IDs.tab

Step 1: Retrieving candidate rDNA reads from cram files. If the input is RNA, we can skip the part of this script that calls MUMmer:
```
sbatch -a 1-10 Parallel_scripts/01.Rescue_rDNA_reads.sh IDs.tab input/ output/
```
Step 2: Mapping candidate rDNA reads to our custom reference (it has a bwa index):
```
sbatch -a 1-10 Parallel_scripts/02.mapping.sh IDs.tab input/ output/
```

Step 3: Variant calling. 20 refers to the ploidy we want to use:
```
sbatch -a 1-10 Parallel_scripts/03.gVCF.sh IDs.tab input/ output/ 20
```

Step 4: merging the vcfs from the different samples into one:
```
Rscript Parallel_scripts/04.merge_gVCF.R IDs.tab output/
```

# Running the code for RNA: