---
title: Day2 Review
author: Shahin Shahsavari
date: October 2020
duration: 60 minutes
---

## DNAseq Alignment Review

We introduced the various genomics file formats yesterday. Let's quickly review those:


### File Formats

1. [fasta](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=BlastHelp)
<img src="../img/fasta.png" width="600">

2. [fastq](https://support.illumina.com/bulletins/2016/04/fastq-files-explained.html)
<img src="../img/fastq.png" width="600">

3. [SAM/BAM](https://gatk.broadinstitute.org/hc/en-us/articles/360035890791-SAM-or-BAM-or-CRAM-Mapped-sequence-data-formats)
<img src="../img/sam.png" width="600">


4. [VCF/BCF](https://www.internationalgenome.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-40/)
<img src="../img/vcf.png" width="900">

### Alignment Steps

1. Index your reference genome

```bash
$ bwa index genome.fa
```

2. Align your raw data

```bash
bwa mem -t $threads $ref $fastq_r1 $fastq_r2 > aln.sam
```

3. Convert (compress) your sam file to bam

```bash
samtools view -b aln.sam > aln.bam
```

4. Sort your bam file

```bash
sambamba sort aln.bam
```

5. mark PCR duplicates (recommended for variant calling)

```bash
sambamba markdup -t 2 aln.sorted.bam aln.sorted.markdup.bam
```

> You could view this in IGV. Make sure the index is present in the same directory.

6. Variant Calling

##### freebayes

```bash
freebayes -f ref.fa aln.sorted.markdup.bam > var.vcf
```

OR

##### bcftools

```bash
bcftools mpileup -Ou -f ref.fa aln.bam | bcftools call -mv -Ob -o var.bcf
bcftools index var.bcf
```

> **NOTE** You can replace `-Ob` with `-Ov` if you want to produce a text file (var.vcf)

OR

##### gatk

```bash
gatk HaplotypeCaller -f ref.fa -I aln.bam -O var.vcf
```

## BASH

Let's look at a few text file:

```bash
cd /data/review
```
---

**Exercise**

- print the content of `file1.tsv` to your terminal
- print the values in the second column of `files.tsv`
- view the content of `file2.csv`
- can you print the values in the second column of `file2.csv` using the same command?
- does the file extention matter bash?
---

## Text File Delimiters

We are working with 4 text files in the `/data/review` directory. What are the data columns separated by in each file?

The most common delimiter in a text file would be `\t` or <kbd>tab</kbd>. These files are frequently given the extention `.tsv`
However, it is also common to see `, ; :` used as delimiters.

If you want to view the second column of `file4.txt`, you can specify a delimiter using cut -d:

```bash
cut -d ";" -f 2 file4.txt 
```

**Exercise1**

- View the content of `cars.csv`?
- What is the delimiter?
- Find all the lines that contain "US" and print only the first and fourth and ninth columns

****
