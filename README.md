# RNAseq-container
The programs run in this pipeline are:
* fastp - does quality control filtering on paired end reads  
* ParDRE - removes artifact duplicates from the paired end reads  
* kallisto - does pseudoalignment on the reads to a transcriptome  
* DESeq - does statistical analysis on the differential expressions between case and controls  

The files are split into the programs run in bash using a python script and a R script to run the R program.

### Python script  
This script takes in the input files and output folder. 

#### Input files  
The input files are provided to the program in a single text file.  
The fasta files are provided on a single line per paired ends in the form of `Input: path/file_R1 \t path/file_R2`.  
The full path of the transcriptome file is provided on its own line.   
The input files are the paired end fasta files and a transcriptome file. 
The input files are paired end fasta reads. Each paired end set have the same file name but ends in R1 and R2.  
The transcriptome file can either be a transcriptome file containing transcripts for a genome or a kallisto transcript file that was previously created from another run of the kallisto program.  

#### Output folder
The script will make directories within the provided output folder to sort and write the output files  

The script will call other python scripts to create histograms of read lengths for pre and post fastp filtering and the transcripts in the transcriptome.  

The script will also use the transcriptome file to generate a conversion table to convert the transcripts into genes for the final expression levels.  

### R script  
The R script gets called internally at the end of the python script to run the DESeq.  
The script will run DESeq and will calculate gene rankings based on the -log(pvalue).  
The output of this program is the differential expressions of the genes as well as a subset of the top 100 differentially expressed genes.  





