#!/bin/bash

#SBATCH --job-name=diamond_blastp_blast_v5.sh_20240124
#SBATCH --time=0:30:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=85gb
#SBATCH --mail-type=FAIL
#SBATCH --account=PAS1444


### Change the directory to where you need to be
cd /fs/ess/PAS1444/Freezing_tolerance_RNA_Seq/diamond


### Run the blastp of the personally created GREM4 and CabSav concatentated protein database to produce a .blast output file
~/.conda/envs/diamond/bin/diamond blastp \
-d /fs/ess/PAS1444/Freezing_tolerance_RNA_Seq/diamond/database/concat_Vl_CabSav_proteins \
-q concat_Vl_CabSav_proteins.fasta \
--threads 20 \
--ultra-sensitive \
--no-self-hits \
-e 1e-2 \
--max-target-seqs 10  \
-f 6 \
-o /fs/ess/PAS1444/Freezing_tolerance_RNA_Seq/diamond/diamondout4.0.blast

# -d = database to use (just use the name of the database file you wish to use, not the file extension too)
# -q = your query input (this will just be the concatenated protein fasta file again containing both the GREM4 and PN proteins)
# --threads = dictates how many threads it can spread the job over
# --iterate = Run multiple rounds of searches with increasing sensitivity
# --sensitive = Enable the sensitive mode designed for full sensitivity for hits of >40% identity.
# --no-self-hits = prevents diamond from reporting alignments that are just the same identical gene as the query back at you.
# -e = dicates the e-value cutoff (Bo recommends 1e-10)
# --max-target-seq = dicate the maximum number of matches you wish to report
# -f = dictates the format of the output file
# -o = dicates the name of the output file