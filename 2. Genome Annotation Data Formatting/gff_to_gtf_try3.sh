#Originally prepared by Cullen Dixon

#### Initial Notes - 
# '~' = /users/PAS1444/singh1815
# All output files will be placed at the location from which the job was submitted of the submission file

#### .gff3 File to .gtf File Conversion - 
# When quantifying gene expression, stringtie needs a .gtf file - if you do not have .gft but rather a more conventional .gff3 file, you need to convert the .gff file to a .gft file first

#Turn on the python module
module load python/3.7-2019.10

#Turn on the RNAseq environment within which 'gffread' resides (activate the environment)
source activate RNAseq

#Run the command which creates a .gtf file from the .gff file
gffread -T /fs/ess/PAS1444/Freezing_tolerance_RNA_Seq/work_dir/Cab08.chr.gff3 -o /fs/ess/PAS1444/Freezing_tolerance_RNA_Seq/work_dir/CabSav_annoV3.gtf


#Remember to now back on out of the RNAseq environment and miniconda and python