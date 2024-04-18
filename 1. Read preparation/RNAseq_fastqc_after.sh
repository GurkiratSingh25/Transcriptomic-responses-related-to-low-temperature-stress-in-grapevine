/users/PAS1444/li10917/miniconda2/envs/rnaseq/bin/fastqc -t 4 -f fastq *clean.fq

mkdir fastQC_output_after

mv *fastqc.html fastQC_output_after
mv *fastqc.zip fastQC_output_after