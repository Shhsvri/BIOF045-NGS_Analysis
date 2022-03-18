---
title: Day3 Review
author: Shahin Shahsavari
date: March 2022
duration: 60 minutes
---

## RNAseq Alignment Review

### Alignment Steps

1. Index your reference genome

```bash
$ STAR	--runThreadN $threads \
	--runMode genomeGenerate \
	--genomeDir $genome_dir \
	--genomeFastaFiles $genome_fa \
	--sjdbGTFfile $genome_gtf \
	--sjdbOverhang 149
```

2. Align your raw data

```bash
$ STAR	--runThreadN 2 \
	--genomeDir STAR_chr2_genome/ \
	--readFilesIn raw_data/sample1.fastq \
	--quantMode GeneCounts \
	--outFileNamePrefix results/sample1_
```

3. Convert (compress) your sam file to bam

```bash
$ samtools view -b sample1_Aligned.out.sam  > sample1.bam
```

4. Sort your bam file

```bash
$ samtools sort sample1.bam > sample1.sorted.bam
```

6. generate the index

```bash
$ samtools index sample1.sorted.bam
```

> This generates the index file that ends with `.bai`. You could then view the bam file in IGV.
Make sure the index is present in the same directory.


## Create a script

You could use `gedit` through *X2go* to create the RNAseq script. Alternatively, if you feel
comfortable using `vim`, you could use your personal computer's terminal to connect to the
server and use `vim` for scripting and text editing.

In order to connect to the server over `ssh`, open the terminal on your local computer and
run:

```bash
$ ssh USERNAME@SERVER_IP
```

```bash
$ gedit RNAseq.sh
```

OR

```bash
$ vim RNAseq.sh
```

---

```bash
#!/bin/bash

# BIOF045: 03/17/2022
# This script for RNA alignment, sorting, and indexing

# 0. set up the file structure and change into the review directory
cd ~/Day4/Day3_review
genome_dir=~/Day3/STAR_chr2_genome/
threads=2

for fastq_file in *.fastq
 do
	# 1. create a file prefix
	## "sample1.fastq" is stored in $i
	## "sample1"       is stored in $sample as follows
	sample=${fastq_file%.*}

	# 2. STAR alignment using 2 threads
	STAR	--runThreadN $threads \
		--genomeDir $genome_dir \
		--readFilesIn $fastq_file \
		--quantMode GeneCounts \
		--outFileNamePrefix $sample

	# 3. Convert your sam file into bam
	samtools view -b $sample"Aligned.out.sam" > $sample".bam"

	# 4. Sort your bam file
	samtools sort $sample".bam" > $sample".sorted.bam"

	# 5. index the bam files
	samtools index $sample".sorted.bam"
 done

```


After you create the script, you could run it using `source RNAseq.sh`.


For the last step, we will extract our counts from the STAR output with the
*ReadsPerGene.out.tab* suffix.
Per the [manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf),
column 2 contains counts for unstranded RNA-seq.

```bash
$ cut -f 1 sample1ReadsPerGene.out.tab > genes.txt
$ cut -f 2 sample1ReadsPerGene.out.tab > sample1.txt
$ cut -f 2 sample2ReadsPerGene.out.tab > sample2.txt
$ paste genes.txt sample1.txt sample2.txt > counts.tsv

```

or you could use redirection as follows:

```bash
$ paste <(cut -f 1,2 sample1ReadsPerGene.out.tab) <(cut -f 1 sample2ReadsPerGene.out.tab) > counts.tsv
```


---
