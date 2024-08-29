# rDNA_variant_caller
Repository to call rDNA variants from short read sequencing

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

# Running the code

To run the code sequentially
```
./Pipeline.sh -h
This code runs a variant caller for rDNA reads

Syntax: Pipeline.sh [-h|-t|-n|-f|-i|-o]
parameters:
h     Print the help message.
t     Type of input data: Either DNA or RNA.
n     Basename of the input data: (e.g., for a bam file named sample.bam: -n sample, for fastqs named sample_1.fastq.gz and sample_2.fasq.gz: -n sample. This means that paired fastqs should end by _1.fastq.gz and _2.fastq.gz)
f     Type of input data format: Either bam (for bam/cram/sam) or fastq.
i     Input folder.
o     Output folder.
```

To run the code in a parallel manner, check inside Scripts_parallel/
