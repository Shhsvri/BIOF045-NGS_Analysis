---
title: Alignment with BWA
author: Shahin Shahsavari
date: March 2021
duration: 120
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

In order to remain organized, I always prepare my directory as follows

```bash
~/Day2/
    ├── raw_data/
    ├── quality_control/
    ├── results/
    ├── genome/
```

This has already been prepared for you using the `mkdir` function.


Then I always copy the raw fastq data into this folder using the `cp` command


## Part0: executables and tools Setup

We will be using the following softwares. These have been already installed on our server.

1. [bwa](https://github.com/lh3/bwa)
2. [samtools](https://github.com/samtools/samtools)
3. [PICARD](https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard-)
4. [freebayes](https://github.com/ekg/freebayes)
5. [IGV](https://software.broadinstitute.org/software/igv/download)

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
$ cd ~/Day2/raw_data
$ ls -l
$ fastqc ptA_R1.fastq
$ fastqc ptA_R2.fastq
$ ls -l
```

#### Viewing the result of fastqc

The fastqc command creates an html file that could be opened using your web browser.
On the server, we have _firefox_ installed.

```bash
$ firefox ptA_R1_fastqc.html
```

When you are done viewing the fastqc results, move those to the quality\_control
folder.

```bash
$ mv *fastqc* ~/Day2/quality_control
```

### 1.2 Alignment

Choice of alignment tool is often determined by the type of NGS application being conducted.
For variant calling we will use [BWA (Burrows-Wheeler Aligner)](http://bio-bwa.sourceforge.net).
It may be slower than `Bowtie2` or some other alternatives, though it is generally considered to be
more accurate.

#### BWA Modes

Depending on your read length, BWA has different modes optimized for different sequence lengths:

- BWA-backtrack: designed for Illumina sequence reads up to 100bp (3-step)
- BWA-SW: designed for longer sequences ranging from 70bp to 1Mbp, long-read support and split alignment
- BWA-MEM: shares similar features to BWA-SW, but BWA-MEM is the latest, and is generally recommended for
high-quality queries as it is faster and more accurate. BWA-MEM also has better performance than
BWA-backtrack for ≥70 bp Illumina reads.

### 1.2.1 Creating BWA-MEM index

The BWA aligner needs to create an index of our reference genome to process our data faster.
I have obtained the latest major assembly of the human reference genome (GRCh38).

In order to practice and see how this works, I created a small reference into the `reference_test` folder.
Let's view this file first

```bash
$ cd ~/Day2/reference_test
$ cat test_genome.fa
$ bwa index test_genome.fa
$ ls -l -h
```

> **NOTE** The human reference genome contains 3 billion base pairs. Indexing is a computationally intensive process.
I already generated the BWA index which took about an hour.
Your institution's High Performance Computing server admin may have most likely generating these for you.
On NIH Biowulf, these references are in `/fdb/igenomes/Homo_Sapiens/UCSC/hg38/Sequence/BWAIndex/genome.fa*`


Let's explore the hg38 fasta reference file at `~/Day2/genome/` before alignment.

---
**Exercise**
1. How large (Bytes) is this fasta reference? `ls`
2. View the top 20 lines of your fasta reference. `head`
3. Count the number of lines in your fasta reference `wc`
4. Find all the lines that contain ">" `grep '>'`
5. Find all the lines that contain ">" and in addition to 3 more lines after each hit. `grep -A 3 '>'`
---

### 1.2.2 Aligning reads with BWA-MEM

Now that we have our indexes, let's perform alignment on our paired-end reads for ptA (patient A).

> We will find out what disease this individual has later today.

Before we start, let's check the [user manual](http://bio-bwa.sourceforge.net/bwa.shtml) for bwa.
I encourage you to peruse the document for bwa. Any other bioinformatics tool you plan to use
in the future should have a similar manual telling you how to choose the correct options for your data.

One of the most used options for aligning reads to a reference genome using BWA-MEM is:

- `-t`: number of threads/cores

Additionally, we will specify:

- the path to genome indexes including prefix
- FASTQ files for paired-end reads (R1 and R2)
- `>`: save alignment output to a SAM(Sequence Alignment Map) file. aka standard output

```bash
$ cd ~/Day2
$ bwa mem -t 2 \
	genome/hg38.fa \
	raw_data/ptA_R1.fastq raw_data/ptA_R2.fastq \
	> results/ptA.sam
```

It may take some time for the process to complete. When alignment is over, you can view the content of your sam file.

```bash
$ cd ~/Day2/results
$ head -n 500 ptA.sam
```

> Every sam file starts with a header. It contains the reference chromosomes, command that was used to generate it and more. Sam header lines start with "@"

```bash
$ grep -v "@" ptA.sam | head
```

> -v option for `grep` indicates that we want to exclude all the lines that containt our character.

**Questions**

1. Is the output of our aligner sorted?
2. What are the columns in a sam file?

<img src="../img/sam_bam.png" size=300>

### 1.3 Convert your sam file to bam

SAM (Sequence Alignment Map) files are pure text files which take too much space. Common practice
is to compress these files using `samtools` into BAM (Binary Alignment Map) files.

```bash
$ samtools view -b ptA.sam > ptA.bam
```

### 1.4 Sorting your alignment bam file

The next step is to sort all the aligned records by their chromosomal location.
There are a variety of tools for this task and they all perform the same task. Some are faster than
others. We will use `samtools` which is very fast and generates all the needed downstream files.

```bash
$ samtools sort ptA.bam > ptA.sorted.bam
```

This will generate the sorted bam file in the same directory.

---
**Exercise**

We currently have 3 files. `ptA.sam`, `ptA.bam`, and `ptA.sorted.bam`

1. How much smaller is that BAM file compared with the SAM?
2. How does the size of the sorted BAM file compare with our unsorted BAM file?
---

### 1.5 Marking duplicates

The last stage of the alignment phase is marking duplicates, and it is usually only required for variant calling.
We need to find reads that are likely artifacts from the PCR amplification as they can bias variant calls.

![align\_cleanup](../img/workflow_cleanup.png)

If duplicates aren't marked, then the PCR-based errors will be picked up again and again as false
positive variant calls. Duplicates are easy to detect: since they have the same mapping information
and CIGAR string:

![dedup1](../img/dedup_begin.png)

Marking duplicates with tools such as samtools will result in the variant caller ignoring these PCR-based errors, and instead seeing:

![dedup1](../img/dedup_end.png)

The variant caller will be more likely to discard the error, instead of calling it as a variant.

```bash
$ PicardCommandLine MarkDuplicates I=ptA.sorted.bam O=ptA.markdup.sorted.bam M=ptA_markdup_metrics.txt
```


### 1.6 Creating an index for the final bam file

Now that we have a sorted BAM file that has duplicates marked, we need to ensure the index file
for it exist. Just like a book that needs a table of contents, a bam file needs an index.

```bash
samtools index ptA.markdup.sorted.bam
```

This file is now ready for visualization in IGV or any other visualization tool.

### 1.7 Viewing the files

[IGV](http://software.broadinstitute.org/software/igv/) (Integrative Genomics Viewer) is a graphical
tool for the visualiztion of sorted bam files.

> **NOTE** BAM files must be indexed before viewing in IGV. Ensure a bam index file exists in the same directory alongside your BAM file.

```
$ la -l

-rw-rw-r-- 1 shahin shahin 56382104 Mar 20 17:21 ptA.markdup.sorted.bam
-rw-rw-r-- 1 shahin shahin  1534200 Mar 20 17:54 ptA.markdup.sorted.bam.bai
```

##Create a script

For the final step, we would like to create a shell script that can run all commands we have used today in one go. As always, it is good practice to leave comments in your shell scripts.

Use vim to write the following into a shell script named "DNA\_alignment.sh":

```
#!/bin/bash

# This is for DNA alignment, sorting, and indexing
  
## 1. BWA MEM alignment with 2 threads

bwa mem -t 2 genome/hg38.fa raw_data/ptA_R1.fastq raw_data/ptA_R2.fastq > results/ptA.sam


## 2. Convert sam to bam

samtools view -b results/ptA.sam > results/ptA.bam


## 3. Sort your bam file using samtools

samtools sort results/ptA.bam > results/ptA.sorted.bam


## 4. markduplicates with PICARD

PicardCommandLine MarkDuplicates I=ptA.sorted.bam O=ptA.markdup.sorted.bam M=ptA_md_metrics.txt
```

---

This lesson has been developed by Shahin Shahsavari, adapted from the teaching team at the Harvard Chan Bioinformatics Core (HBC). These are open access materials distributed under the terms of the Creative Commons Attribution license (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.

---
