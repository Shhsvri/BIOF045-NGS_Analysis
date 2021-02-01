---
title: Alignment with BWA
author: Shahin Shahsavari
date: October 2020
duration: 90
---

## Learning Objectives:

* Exploring the genomics files
* Exploring the variant calling workflow
* Hands-on alignment 

## Variant Calling Workflow

The variant calling workflow begins with quality control and alignment, similar to the other NGS applications. Alignment is followed by alignment clean-up to prepare data for variant calling. Then, variant calling is performed, followed by filtering and annotation of the variant calls.

![var\_calling\_workflow](../img/variant_calling_workflow.png)


## Directory Set-up

To start with variant calling, we need to set-up our directory structure, and ensure we have all the softwares that we need.

Create the following directory structure for the variant calling project in your home directory:

```bash
~/var_calling/
    ├── logs/
    ├── meta/
    ├── raw_data/
    ├── scripts/
    ├── quality_control
    ├── results/
        ├── bwa/
```

```bash
$ mkdir ~/var_calling
$ cd ~/var_calling

$ mkdir -p quality_control raw_data scripts logs meta results/bwa
$ ls -l
```

Now that we have the directory structure created, let's copy over the raw data (fastq) to perform quality control and alignment:

```bash
$ cp /data/DNAseq/human/*.fq ~/var_calling/raw_data
$ ls -l ~/var_calling/raw_data
```

## Part0: executables and tools Setup

We will be using the following softwares. These have been already installed on our server.

1. [bwa](https://github.com/lh3/bwa)
2. [samtools](https://github.com/samtools/samtools)
3. [htslib](https://github.com/samtools/htslib)
5. [freebayes](https://github.com/ekg/freebayes)
6. [vcflib](https://github.com/ekg/vcflib/)
7. [sambamba](https://github.com/lomereiter/sambamba)
8. [seqtk](https://github.com/lh3/seqtk)
9. [mutatrix](https://github.com/ekg/mutatrix)
10. [sra-tools](https://github.com/ncbi/sra-tools/wiki/HowTo:-Binary-Installation)

Depending on the software, you can download and compile the source code using this kind of pattern:

> **DO NOT RUN THIS**

```bash
$ git clone https://github.com/lh3/bwa
$ cd bwa
$ make
```

## Part1: Aligning raw sequencing data with `bwa mem`

### 1.1 QC

> **NOTE:** It is a good practice to perform a short Quality Control step to ensure there are no
serious problems with our samples such as the presence of vector contamination or low quality reads.
In general, bwa mem doesn't need trimming.

```bash
$ cd ~/var_calling/raw_data
$ ls -l
$ fastqc 08008_r1.fq
$ mv *fastqc* ~/var_calling/quality_control
```

#### Viewing the result of fastqc

The fastqc command creates an html file that could be opened using your web browser.
On the server, we have _firefox_ installed.

```bash
cd ~/var_calling/quality_control
firefox 08008_r1_fastqc.html
```

### 1.2 Alignment

Choice of alignment tool is often determined by the type of NGS application being conducted.
For variant calling we will use [BWA (Burrows-Wheeler Aligner)](http://bio-bwa.sourceforge.net).
It may be slower than Bowtie2 or some other alternatives, though it is generally considered to be more accurate.

#### BWA Modes

Depending on read length, BWA has different modes optimized for different sequence lengths:

- BWA-backtrack: designed for Illumina sequence reads up to 100bp (3-step)
- BWA-SW: designed for longer sequences ranging from 70bp to 1Mbp, long-read support and split alignment
- BWA-MEM: shares similar features to BWA-SW, but BWA-MEM is the latest, and is generally recommended for
high-quality queries as it is faster and more accurate. BWA-MEM also has better performance than
BWA-backtrack for ≥70 bp Illumina reads.

### 1.2.1 Creating BWA-MEM index

The BWA aligner needs to create an index of our reference genome to process our data faster.
We are going to copy our hg38 human reference genome and its bwa index in our reference directory

> **NOTE** The human reference genome contains 3 billion base pairs. Indexing is a computationally intensive process.
I already generated the BWA index which took about an hour.
Your institution's High Performance Computing server admin may have most likely generating these for you.
On NIH Biowulf, these references are in `/fdb/igenomes/Homo_Sapiens/UCSC/hg38/Sequence/BWAIndex/genome.fa*`

```bash
$ cd /data/DNAseq/BWAIndex/
$ ls -l genome.fa*
```

Let's explore this fasta reference file before alignment.

---
**Exercise**
1. How large (Bytes) is this fasta reference?
2. View the top 20 lines of your fasta reference
3. Count the number of lines in your fasta reference
4. Count the number of characters in your fasta reference
4. Find all the lines that contain ">"
5. Find all the lines that contain ">" and in addition to 3 more lines after each hit.
---

### 1.2.2 Aligning reads with BWA-MEM

Now that we have our indexes, let's perform alignment in our paired-end reads for sample 08008.

> We will find out what disease this individual has later today.

Before we start, let's check the [user manual](http://bio-bwa.sourceforge.net/bwa.shtml) for bwa.
I encourage you to peruse the document for bwa. Any other bioinformatics tool you plan to use
in the future should have a similar manual telling you how to choose the correct options for your data.

One of the most used options for aligning reads to a reference genome using BWA-MEM is:

- `-t`: number of threads/cores

Additionally, we will specify:

- the path to genome indexes including prefix
- FASTQ files for paired-end reads (r1 and r2)
- `2>`: save standard error to file
- `>`: save alignment output to a SAM(Sequence Alignment Map) file. aka standard output

```bash
$ cd ~/var_calling
$ bwa mem -t 2 \
	/data/DNAseq/BWAIndex/genome.fa \
	raw_data/08008_r1.fq raw_data/08008_r2.fq \
	2> logs/bwa.err \
	> results/bwa/08008.sam
```

It may take some time for the process to complete. When alignment is over, you can view the content of your sam file.

```bash
$ cd results/bwa
$ head -n 200 08008.sam
$ grep -v "@" 08008.sam | head
```

> Every sam file starts with a header. It contains the reference chromosomes, command that was used to generate it and more. Sam header lines start with "@"

> -v option for `grep` indicates that we want to exclude all the lines that containt our character.

**Questions**

1. Is the output of our aligner sorted?
2. What are the columns in a sam file?

### 1.3 Convert your sam file to bam

SAM files are pure text files which take too much space. Common practice is to compress these files
using samtools into bam files.

```bash
$ samtools view -b 08008.sam > 08008.bam
```

### 1.4 Sorting your alignment bam file

The next step is to sort all the aligned records by their chromosomal location.
There are a variety of tools for this task and they all perform the same task. Some are faster than
others. We will use `sambamba` which is very fast and generates all the needed downstream files.

```bash
$ sambamba sort 08008.bam
```

This will generate the sorted bam file in the same directory.

**Exercise**

1. How much smaller is that BAM file compared with the SAM?
2. How does the size of the sorted BAM file compare with our unsorted BAM file?

### 1.5 Marking duplicates

The last stage of the alignment phase is marking duplicates, and it is usually only required for variant calling.
We need to find reads that are likely artifacts from the PCR amplification as they can bias variant calls.

![align\_cleanup](../img/workflow_cleanup.png)

If duplicates aren't marked, then the PCR-based errors will be picked up again and again as false
positive variant calls. Duplicates are easy to detect: since they have the same mapping information
and CIGAR string:

![dedup1](../img/dedup_begin.png)

Marking duplicates with tools such as sambamba will result in the variant caller ignoring these PCR-based errors, and instead seeing:

![dedup1](../img/dedup_end.png)

The variant caller will be more likely to discard the error, instead of calling it as a variant.

```bash
$ sambamba markdup -t 2 08008.sorted.bam 08008.sorted.markdup.bam
```


### 1.6 Creating an index for the final bam file

Now that we have a sorted BAM file that has duplicates marked, we need to ensure the index file
for it exist. sambamba creates the index by default.

This file is now ready for visualization in IGV or any other visualization tool.

### 1.7 Viewing the files

[IGV](http://software.broadinstitute.org/software/igv/) (Integrative Genomics Viewer) is a graphical
tool for the visualiztion of sorted bam files.

> **NOTE** BAM files must be indexed before viewing in IGV. Ensure a bam index file exists in the same directory alongside your BAM file.

```
$ la -l

-rw-rw-r-- 1 shahin shahin 56382104 Oct 20 17:21 08008.sorted.markdup.bam
-rw-rw-r-- 1 shahin shahin  1534200 Oct 20 17:54 08008.sorted.markdup.bam.bai
```

##Create a script

For the final step, we would like to create a shell script that can run all commands we have used today in one go. As always, it is good practice to leave comments in your shell scripts.

Use vim to write the following into a shell script named "DNA\_alignment.sh":

```
#!/bin/bash

# This is for DNA alignment, sorting, and indexing
  
i=08008
## Set up the variables for reference and raw data 

### Assign the full path of your reference genome
### the input fastq files
### and the output directory you would like to save the end result in

ref=/data/DNAseq/BWAIndex/genome.fa
fastq1=~/var\_calling/raw\_data/"$i"\_r1.fq
fastq2=~/var\_calling/raw\_data/"$i"\_r2.fq
outdir=~/var\_calling/results/


## 1. BWA MEM alignment with 2 threads

bwa mem -t 2 $ref $fastq1 $fastq2 2> $outdir/bwa.err > $outdir/bwa/"$i".sam


## 2. Convert sam to bam

samtools view -b $outdir/bwa/"$i".sam > $outdir/bwa/"$i".bam


## 3. Sort your bam file using sambamba

sambamba sort $outdir/bwa/"$i".bam


## 4. markduplicates with 2 threads

sambamba markdup -t 2 $outdir/bwa/"$i".sorted.bam "$i".sorted.marked.bam
```

For future alignments, you could change the variables (ref, fastq, and outdir), and
run this alignment script using `source DNA_alignment.sh`


---

This lesson has been developed by Shahin Shahsavari, adapted from the teaching team at the Harvard Chan Bioinformatics Core (HBC). These are open access materials distributed under the terms of the Creative Commons Attribution license (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.

---
