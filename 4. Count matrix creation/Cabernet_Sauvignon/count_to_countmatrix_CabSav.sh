#!/bin/bash

#SBATCH --job-name=count_to_countmatrix_CabSav.sh
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100gb
#SBATCH --mail-type=FAIL
#SBATCH --account=PAS1444

###   Creating a Count Matrix from All of the Individual .count Files output by 'coco correct_count'   ###

# All of the following was run at location '/fs/ess/PAS1444/Freezing_tolerance_RNA_Seq/work_dir/CoCo' on 11/21/22 - CWD 11/30/22
# Sort the file by the first column (gene name) which is important because for some odd reason all of the count files don't have their
#       genes sorted in the first column in the same order.


## GS- This file was modified on 10/11/2023

for file in *.count
do
  (tail -n +3 ${file} | sort -k 1 ) > ${file}.sorted
done
#          ^ skip the first two lines (its a reverse pyschology head (tail?) command)

# Add a header to the column over the count column for use later when we cut them out
for file in *.count.sorted
do
  awk 'NR==1 {print "\t" "\t" FILENAME "\t" "\t"}; 1' ${file} > temp && mv temp ${file}.headed # the "; 1" portion ensures that it just doesn't rewrite the entire file as just the filename.  Without this, this it deletes all text within and just prints the tab + file name.
done

# Combine the various sorted count files into one large file
paste *.sorted.headed | # this says take all of the .headed files just made above and combine them into one file, moving from left to right, with each file being placed in their own successive columns

# Cut out the columns we want in the count matrix file ('gene name/id' & 'count')
#       IMPORTANT: THIS \/ IS THE ORDER OF THE INFORMATION IN THE .count FILES:
#             gene_id gene_name count cpm tpm
cut -f1,3,8,13,18,23,28,33,38,43,48,53,58,63,68,73,78,83,88,93,98,103,108,113,118 \
    > countmatrix_CabSav.txt