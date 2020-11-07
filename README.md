# NGS-paired-end-reads-fastq-generator
A python program used to generate a perfect NGS paired-end reads fastq file. 
It's created for simulation purpose. I used it to evaluate the performance of alignments, calling accuracy, etc.  
Since it aims to create a "perfect" data, the phred quality score for each base is setted to "I".
Reads are assumed to be uniformly distributed with a normal(0,1) shift.

# Requirement

pyhon3

pysam (pip install pysam)



# Arguments

Fasta file at first place 

-t1   Specified your target reference name. (must be identical with the name in fasta file)")

-t2   (potional) Specified another target reference name if you want to generate heterozigous data

-o    Output file name without extention.

-r    The average length of reads. (Default: 150)

-f    The average length of DNA fragments. (Default: 450)

-d    The average alignment depth. (Default: 30)


# usage
**Example for homozygous** 
```
python fastq_generator.py file.fasta -t1 HLA:HLA02169 -o output -r 150 -f 450 -d 30
```
It generates two fq files (output1.fq and output2.fq) for first read and second read. 


**Example for heterozygous** 
```
python fastq_generator.py file.fasta -t1 HLA:HLA02169 -t1 HLA:HLA02169 -t2 HLA:HLA15760 -o output_name -r 150 -f 450 -d 30
```
It generates two fq files (output1.fq and output2.fq)  for first read and second read.

# output
**Output format**
```
@simulate_targetname_1_read0/1
GATTGGGGAGTCCCAGCCTTGGGGATTCCCCAACTCC...
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII...
```
