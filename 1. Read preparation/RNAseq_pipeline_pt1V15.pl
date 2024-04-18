#!usr/bin/perl

use strict;
use warnings;

### Step 1: Quality Check Before Cleanup with Trimmomatic ###

## After you download RNA-seq and convert it into fastq file format, you can start with a quality control step to have a preliminary (and control dataset to reference) look at the input files
#run the 'RNAseq_fastQC.sh' script here


### Step 2: Run Trimmomatic to Cleanup Reads (remove adapters, poor sequencing run files) ###

## We use trimomatic to remove adapters and low quality reads
## single-end and paired-end reads need to be run in different ways, thus the hashing out of lines below

# run Trimmomatic
my @allfastq=glob("*.fq"); # read in all fastq files; must have extension be exactly as input is
foreach my $each (@allfastq){
	my $name;
	## pair end reads
	next if($each=~ /(.*)_1.fq/); # specifies demarcating this as one read file; must match input file name exactly; "for pair end reads, only deal this pair once" - Bo
	if ($each =~ /(.*)_2.fq/) {   # specifies demarcating this as the other read file; must match input file name exactly; "get fastq file name, eg, SRA0000001_1.fastq" - Bo
		$name=$1;                 # "now $name is SRA0000001" - Bo
		my $read1=$name."_1";      
		my $read2=$name."_2";
		system "trimmomatic PE -phred33 $read1.fq $read2.fq  $read1.clean.fq $read1.forward_unpaired.fastq  $read2.clean.fq $read2.reverse_unpaired.fastq  ILLUMINACLIP:Novogene-PE.fa:2:30:10 SLIDINGWINDOW:4:5 LEADING:5 TRAILING:5 MINLEN:36";
	    #            ^must use YOUR OWN file path to the program, or else it won't have permission to go into others' folders and run it from that, thus, it won't run
	}
	## single end reads

	#elsif ($each=~/(.*).fastq/){
		#$name=$1;
		#system "~/miniconda3/envs/RNAseq/bin/trimmomatic SE -phred33 $name.fastq  $name.clean.fq   ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:5 LEADING:5 TRAILING:5 MINLEN:36";
	#}


}


## Remove files containing unpaired reads - ie a paired end read without a pair

`rm *_unpaired.fastq`;


## We now have 'clean' RNA-seq data which is ready for downstream analysis


### Step 3: Quality Check Post-Cleanup ###
# This is the quality control step which matters most as this is the data that will be used in the downstream analysis.  Be sure to scrutinize the fastQC output html files closely.

#run the below to call the RNAseq_fastQC_after.sh script to run another round of quality control on the newly cleaned up read files
system "/fs/ess/PAS1444/Freezing_tolerance_RNA_Seq/work_dir/RNAseq_fastQC_after.sh";