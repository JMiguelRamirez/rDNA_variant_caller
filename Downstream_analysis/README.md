# Code to replicate results of the paper "rDNAcaller: a fast and robust workflow to call ribosomal DNA variants"


Scripts to run different variant callers
```
./lofreq.sh $bam 
./mutect.sh $bam
./gatk_joint_1 $bam $sample_ploidy
./gatk_joint_2 $folder_with_multiple_bam_files $sample_ploidy
```

Scripts to compute precision, recall and F1 score
```
Rscript X
```

Script to compute nucleotide diversity score:
```
Rscript NucleotideDiversity.R
```
