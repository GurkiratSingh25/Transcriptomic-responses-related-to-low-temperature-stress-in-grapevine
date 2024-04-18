#!/bin/bash

#SBATCH --job-name=RNASeqcopy_V2_20230905
#SBATCH --time=6:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=16gb
#SBATCH --mail-type=FAIL
#SBATCH --account=PAS1444


/users/PAS1444/li10917/miniconda2/envs/rnaseq/bin/fastqc -t 4 -f fastq *.fq

mkdir fastQC_output_before

mv *fastqc.html fastQC_output_before
mv *fastqc.zip fastQC_output_before
