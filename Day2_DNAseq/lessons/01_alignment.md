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
    ├── reference_data/
    ├── scripts/
    ├── quality_control
    ├── results/
        ├── bwa/
```

```bash
$ mkdir ~/var_calling
$ cd ~/var_calling

$ mkdir -p quality_control raw_data reference_data scripts logs meta results/bwa
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

**DO NOT RUN THIS**

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
firefox 08008_r1.html
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
On NIH Biowulf, these references are in `/fdb/igenomes/Homo_Sapiens/UCSC/hg38/BWAIndex/genome.fa*`

```bash
cd ~/var_calling/reference_data
cp /data/DNAseq/human/genome.fa* .
```

Let's explore this fasta reference file before alignment.

---
**Exercise**
1. How large (Bytes) is this fasta reference?
2. View the top 20 lines of your fasta reference
3. Count the number of lines in your fasta reference
4. find all the lines that contain ">"
5. find all the lines that contain ">" and in addition to 3 more lines after each hit.
---

### 1.2.2 Aligning reads with BWA-MEM

Now that we have our indeces, let's perform alignment in our paired-end reads for sample 08008.

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
$ bwa mem -t 2 \
	reference_data/genome \
	raw_data/08008_r1.fq raw_data/08008_r2.fq \
	2> logs/bwa.err \
	> results/bwa/08008.sam
```

It may take some time for the process to complete. When alignment is over, you can view the content of your sam file.

```bash
$ cd results/bwa
$ head 08008.sam
```

**Question**

Is this file sorted?

### 1.3 Convert your sam file into bam

SAM files are pure text files which take too much space. Common practice is to compress these files
using samtools into bam files.

```bash
$ samtools view -b 08008.sam > 08008.bam
```

### 1.4 Sorting your alignment bam file

The next step is to sort all the aligned records by their chromosomal location.
There are a variety of tools for this task and they all do the same task. Some are faster than others. sambamba is very fast.

```bash
$ sambamba sort 08008.bam
```

This will generated the sorted bam file in the same directory.

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
$ sambamba markdup -t 2 08008.sorted.bam 08008.sorted.dedup.bam
```

### 1.6 Creating an index for final bam file

Now that we have a sorted BAM file that has duplicates marked, let's index it for visualization with
IGV.


This lesson has been developed by Shahin Shahsavari, adapted from the teaching team at the Harvard Chan Bioinformatics Core (HBC). These are open access materials distributed under the terms of the Creative Commons Attribution license (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.






---
# Ecoli Alignment

Go into the `` directory and download the FASTA reference there:

```bash
$ cd ~/DNAseq/Reference
$ wget http://hypervolu.me/~erik/genomes/E.coli_K12_MG1655.fa
```

> Any other organism you are working with, you need to download its assembled genome from online resources

### 1.2 Setting up our reference Indexes

#### BWA's FM-index

BWA uses the [FM-index](https://en.wikipedia.org/wiki/FM-index), which a compressed full-text substring index based around the [Burrows-Wheeler transform](https://en.wikipedia.org/wiki/Burrows%E2%80%93Wheeler_transform).
To use this index, we first need to build it:

```bash
bwa index E.coli_K12_MG1655.fa
```

You should see `bwa` generate some information about the build process:

```text
[bwa_index] Pack FASTA... 0.04 sec
[bwa_index] Construct BWT for the packed sequence...
[bwa_index] 2.26 seconds elapse.
[bwa_index] Update BWT... 0.04 sec
[bwa_index] Pack forward-only FASTA... 0.03 sec
[bwa_index] Construct SA from BWT and Occ... 0.72 sec
[main] Version: 0.7.8-r455
[main] CMD: bwa index E.coli_K12_MG1655.fa
[main] Real time: 3.204 sec; CPU: 3.121 sec
```

And, you should notice a new index file which has been made using the FASTA file name as prefix:

```bash
$ ls -rt1 E.coli_K12_MG1655.fa*
# -->
E.coli_K12_MG1655.fa
E.coli_K12_MG1655.fa.fai
E.coli_K12_MG1655.fa.bwt
E.coli_K12_MG1655.fa.pac
E.coli_K12_MG1655.fa.ann
E.coli_K12_MG1655.fa.amb
E.coli_K12_MG1655.fa.sa
```

When done, go back to the parent directory:

```bash
$ cd ../
```

### 1.3 Aligning our data against the E. Coli K12 reference

Here's an outline of the steps we'll follow to align our K12 strain against the K12 reference:

1. use bwa to generate SAM - Sequence Alignment Map - records for each read
2. convert the output to BAM - Binary Alignment Map
3. sort the output based on position
4. mark PCR duplicates that result from exact duplication of a template during amplification

We could the steps one-by-one, generating an intermediate file for each step.
However, this isn't really necessary unless we want to debug the process, and it will make a lot of excess files which will do nothing but confuse us when we come to work with the data later.
Thankfully, it's easy to use [unix pipes](https://en.wikiepdia.org/wiki/Pipeline_%28Unix%29) to stream most of these tools together (see this [nice thread about piping bwa and samtools together on biostar](https://www.biostars.org/p/43677/) for a discussion of the benefits and possible drawbacks of this).

> Details on BWA and its functionality can be found in the [user manual](http://bio-bwa.sourceforge.net/bwa.shtml); I encourage you to peruse through to get familiar with all available options.

You can now run the alignment using a piped approach. _Replace `$threads` with the number of CPUs you would like to use for alignment._ Not all steps in `bwa` run in parallel, but the alignment, which is the most time-consuming step, does. You'll need to set this given the available resources you have.

```bash
bwa mem -t 2 -R '@RG\tID:K12\tSM:K12' \
    ~/DNAseq/Reference/E.coli_K12_MG1655.fa ~/DNAseq/Raw_data/SRR1770413_1.fastq.gz ~/DNAseq/Raw_data/SRR1770413_2.fastq.gz \
    | samtools view -b - > SRR1770413.raw.bam
smatools sort SRR1770413.raw.bam
samtools markdup SRR1770413.raw.sorted.bam SRR1770413_fin.bam
```

Breaking it down by line:

- *alignment with bwa*: `bwa mem -t $threads -R '@RG\tID:K12\tSM:K12'` --- this says "align using so many threads" and also "give the reads the read group K12 and the sample name K12"
- *reference and FASTQs* `E.coli_K12_MG1655.fa SRR1770413_1.fastq.gz SRR1770413_2.fastq.gz` --- this just specifies the base reference file name (`bwa` finds the indexes using this) and the input alignment files. The first file should contain the first mate, the second file the second mate.
- *conversion to BAM*: `samtools view -b -` --- this reads SAM from stdin (the `-` specifier in place of the file name indicates this) and converts to BAM.
- *sorting the BAM file*: `sambamba sort SRR1770413.raw.bam` --- sort the BAM file, writing it to `.sorted.bam`.
- *marking PCR duplicates*: `sambamba markdup SRR1770413.raw.sorted.bam SRR1770413.bam` --- this marks reads which appear to be redundant PCR duplicates based on their read mapping position. It [uses the same criteria for marking duplicates as picard](http://lomereiter.github.io/sambamba/docs/sambamba-markdup.html).

Now, run the same alignment process for the O104:H4 strain's data. Make sure to specify a different sample name via the `-R '@RG...` flag incantation to specify the identity of the data in the BAM file header and in the alignment records themselves:

```bash
bwa mem -t 2 -R '@RG\tID:O104_H4\tSM:O104_H4' \
    ~/DNAseq/Reference/E.coli_K12_MG1655.fa ~/DNAseq/Raw_data/SRR341549_1.fastq.gz  ~/DNAseq/Raw_data/SRR341549_2.fastq.gz \
    | samtools view -b - >SRR341549.raw.bam
samtools sort SRR341549.raw.bam
samtools markdup SRR341549.raw.sorted.bam SRR341549.bam
```

As a standard post-processing step, it's helpful to add a BAM index to the files. This let's us jump around in them quickly using BAM compatible tools that can read the index. `sambamba` does this for us by default, but if it hadn't or we had used a different process to generate the BAM files, we could use samtools to achieve exactly the same index.

```bash
samtools index SRR1770413.bam
samtools index SRR341549.bam
```

## Part 2: Calling variants

Now that we have our alignments sorted, we can quickly determine variation against the reference by scanning through them using a variant caller.
There are many options, including [samtools mpileup](http://samtools.sourceforge.net/samtools.shtml), [platypus](http://www.well.ox.ac.uk/platypus), and the [GATK](https://www.broadinstitute.org/gatk/).

For this tutorial, we'll keep things simple and use [freebayes](https://github.com/ekg/freebayes). It has a number of advantages in this context (bacterial genomes), such as long-term support for haploid (and polyploid) genomes. However, the best reason to use it is that it's very easy to set up and run, and it produces a very well-annotated VCF output that is suitable for immediate downstream filtering.

### 2.1 Variant calls with `freebayes`

It's quite easy to use `freebayes` provided you have your BAM file completed. We use `--ploidy 1` to indicate that the sample should be genotyped as haploid.

```bash
freebayes -f ~/DNAseq/Reference/E.coli_K12_MG1655.fa --ploidy 1 SRR1770413.bam >SRR1770413.vcf
```

### Joint calling

We can put the samples together if we want to find differences between them. Calling them jointly can help if we have a population of samples to use to help remove calls from paralogous regions. The Bayesian model in freebayes combines the data likelihoods from sequencing data with an estimate of the probability of observing a given set of genotypes under assumptions of neutral evolution and a [panmictic](https://en.wikipedia.org/wiki/Panmixia) population. For instance, [it would be very unusual to find a locus at which all the samples are heterozygous](https://en.wikipedia.org/wiki/Hardy%E2%80%93Weinberg_principle). It also helps improve statistics about observational biases (like strand bias, read placement bias, and allele balance in heterozygotes) by bringing more data into the algorithm.

However, in this context, we only have two samples and the best reason to call them jointly is to make sure we have a genotype for each one at every locus where a non-reference allele passes the caller's thresholds in either sample.

We would run a joint call by dropping in both BAMs on the command line to freebayes:

```bash
freebayes -f ~/DNAseq/Reference/E.coli_K12_MG1655.fa --ploidy 1 SRR1770413.bam SRR341549.bam >e_colis.vcf
```

As long as we've added the read group (@RG) flags when we aligned (or did so after with [bamaddrg](https://github.com/ekg/bamaddrg), that's all we need to do to run the joint calling. (NB: due to the amount of data in SRR341549, this should take about 20 minutes.)

### `bgzip` and `tabix`

We can speed up random access to VCF files by compressing them with `bgzip`, in the [htslib](https://github.com/samtools/htslib) package.
`bgzip` is a "block-based GZIP", which compresses files in chunks of lines. This chunking let's us quickly seek to a particular part of the file, and support indexes to do so. The default one to use is tabix. It generates indexes of the file with the default name `.tbi`.

```bash
bgzip SRR1770413.vcf # makes SRR1770413.vcf.gz
tabix -p vcf SRR1770413.vcf.gz
```

Now you can pick up a single part of the file. For instance, we could count the variants in a particular region:
