#!/bin/bash

#SBATCH --job-name=PullgenesforKEGG_GREM_Chill_Unique_20240123
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --mail-type=FAIL
#SBATCH --account=PAS1444

for file in DEGs_Unique_to_Chill_Stress_in_GREM.txt

do
   /fs/ess/PAS1444/Freezing_tolerance_RNA_Seq/KEGG/seqtk/seqtk subseq Vitis_labrusca.cds-transcripts_noT.fa ${file} > ${file}.cdsseqs.fa
done