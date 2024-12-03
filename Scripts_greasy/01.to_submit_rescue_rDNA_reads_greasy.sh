#!/bin/bash

#SBATCH --job-name=greasy_gatk
#SBATCH --output=./out/greasy_gatk.%A.out
#SBATCH --error=./out/greasy_gatk.%A.err
#SBATCH --account=bsc83
#SBATCH --qos=gp_bscls
#SBATCH --time=04:00:00

# load modules
module load ucx greasy


GREASY_LOGFILE="./out/greasy_gatk.log"
export GREASY_LOGFILE=$GREASY_LOGFILE
GREASY_RESTARTFILE="./out/greasy_gatk.rst"
export GREASY_RESTARTFILE=$GREASY_RESTARTFILE


greasy 01.greasy_file_rescue_rDNA_reads.txt 
#5000 rows with 200 tasks and 56 cpus per tasks, so 100 nodes, we finish in 2 hours, maybe ploidy 50 takes a bit more, but 4 hours or so


