#!/bin/bash

#SBATCH --job-name=RNAseq_multiQC_20230906
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --mail-type=FAIL
#SBATCH --account=PAS1444


cd /fs/ess/PAS1444/Freezing_tolerance_RNA_Seq/work_dir/fastQC_output_before/multiqcoutput

/users/PAS1444/li10917/miniconda2/envs/rnaseq/bin/multiqc ../

multiqc ../
