#!/bin/bash

#SBATCH --job-name=RNAseq_pipeline_pt1_20230912
#SBATCH --time=24:00:00
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=4
#SBATCH --mem=16gb
#SBATCH --mail-type=FAIL
#SBATCH --account=PAS1444


cd $SLURM_SUBMIT_DIR
perl RNAseq_pipeline_pt1V15.pl