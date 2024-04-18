#!/bin/bash

#SBATCH --job-name=CoCo_annotation_try7.sh_20230922
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100gb
#SBATCH --mail-type=FAIL
#SBATCH --account=PAS1444

#### The reason we have to even run this is so that CoCo converts the .gft file into a version of .gtf
#		that is compatible with the 'correct_count' step, which is the one we really care about.
            

## Load the module you intend run (need for samtools) (with this specific version of conda which is very necessary) 
# samtools is used within the CoCo program
module load python/3.6-conda5.2

### module load python/3.9-2022.05 (This is the new version GS was previously using)


# Activate the RNAseq environment which contains samtools
source activate /users/PAS2250/singh1815/.conda/envs/RNAseq


## Set the placeholder variables
gtf_file=/fs/ess/PAS1444/Freezing_tolerance_RNA_Seq/work_dir/CoCo/Vitis_labrusca_annoV2.plusgenesplusbiotypeminusexquot.gtf
output_location=/fs/ess/PAS1444/Freezing_tolerance_RNA_Seq/work_dir/CoCo/Vitis_labrusca_annoV2.cococorrectannoout.gtf


## Run CoCo (correct_annotation)
/fs/ess/PAS1444/Freezing_tolerance_RNA_Seq/work_dir/CoCo/Software/coco/bin/coco correct_annotation \
	-o "$output_location" \
	-b FILLER \
	-V \
	"$gtf_file"


# Options:
# -h = show this help message and exit
# -o = Name of the correct_annotation output gtf. Default: add 'correct_annotation' suffix to original gtf name.
# -b = List of gene biotypes to correct for. Must be given in a comma seperated list (no spaces). Default:
#                       snoRNA,scaRNA,tRNA,miRNA,snRNA
#						I am giving it a biotype to search for of FILLER because I don't want it to do any correcting
# -V = Print the progression of the gtf file reading
# -f = Minimal reciprocal overlap fraction to merge embedded  genes to avoid conflict. 
#						To disable, set to -1. Default: 0.85.

