---
title: "Variant calling with Freebayes"
author: "Shahin Shahsavari"
date: July 2021
duration: 45
---

## Learning Objectives:

* Call variants with Freebayes
* Get familiar with the Variant Call Format (VCF)
* Use vcftools to perform some simple filtering on the variants in the VCF file


## Variant Calling

We have the aligned and cleaned up the data, and have a BAM file ready for calling variants. 

<img src="../img/variant_calling_workflow_2.png" width="450">

Some of the more popular tools for calling variants include [bcftools mpileup](http://samtools.sourceforge.net/mpileup.shtml), [the GATK suite](https://www.broadinstitute.org/gatk/about/) and [FreeBayes](https://github.com/ekg/freebayes#freebayes-a-haplotype-based-variant-detector) ([Garrison and Marth, 2012](http://arxiv.org/abs/1207.3907)). While it can be useful to work through the [GATK Best Practices](https://software.broadinstitute.org/gatk/best-practices/) we will be using FreeBayes in this module as it is just as sensitive and precise (maybe even better), but has java requirements. Also, it is much more straightforward to use. After calling variants, we will filter out low quality variants using *[vcftools](https://vcftools.github.io/index.html)*, a toolkit designed to work with Variant Call Format or VCF files.

## Freebayes

*FreeBayes* is a **haplotype-based** variant detector and is a great tool for calling variants from a population. 

> "FreeBayes is a Bayesian genetic variant detector designed to find small polymorphisms, specifically SNPs (single-nucleotide polymorphisms), indels (insertions and deletions), MNPs (multi-nucleotide polymorphisms), and complex events (composite insertion and substitution events) smaller than the length of a short-read sequencing alignment."

> "FreeBayes is haplotype-based, in the sense that it calls variants based on the literal sequences of reads aligned to a particular target, not their precise alignment. This model is a straightforward generalization of previous ones (e.g. PolyBayes, samtools, GATK) which detect or report variants based on alignments. This method avoids one of the core problems with alignment-based variant detection--- that identical sequences may have multiple possible alignments:"

<img src="../img/freebayes_2.png" width="600">
---
<img src="../img/freebayes_1.png" width="200">

> "FreeBayes uses short-read alignments (BAM files with Phred+33 encoded quality scores, now standard) for any number of individuals from a population and a reference genome (in FASTA format) to determine the most-likely combination of genotypes for the population at each position in the reference. It reports positions which it finds putatively polymorphic in variant call file (VCF) format."

### Running the FreeBayes Variant Caller

```bash
$ cd ~/Day2/results/
```

Let's check the freebayes options:

```bash
$ freebayes -h
```

> **DO NOT RUN THIS**

```bash
$ freebayes -f ~/Day2/genome/hg38.fa ptA.markdup.sorted.bam > ptA.vcf
```
Freebayes is a very popular variant caller, however, it takes a long time for it to create vcf files.

#### Running the bcftools Variant Caller (Alternative)

```bash
$ cd ~/Day2/results
$ bcftools mpileup -f ~/Day2/genome/hg38.fa ptA.markdup.sorted.bam | bcftools call -mv -Ov -o ptA.vcf
```

## _Optional Steps_

### Variant Call Format (VCF)

VCF is a text format. It usually has several header lines before the actual data; the header lines start with `##`. There is usually only 1 VCF file generated for all the samples in an experiment. Variants are represented in the rows, and each sample has a column with the status of a given variant:

```
##format=VCFv4.0
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=1000GenomesPilot-NCBI36
##phasing=partial
#CHROM  POS     ID        REF   ALT    QUAL  FILTER  INFO                                 FORMAT       NA00001         NA00002         
20      14370   rs6054257 G     A      29    0       NS=55;DP=255;AF=0.768;DB;H2          GT:GQ:DP:HQ  0|0:48:1:51,51  1|0:48:8:51,51  
20      13330   .         T     A      3     q10     NS=55;DP=202;AF=0.024                GT:GQ:DP:HQ  0|0:49:3:58,50  0|1:3:5:65,3    
20      1110696 rs6040355 A     G,T    67    0       NS=55;DP=276;AF=0.421,0.579;AA=T;DB  GT:GQ:DP:HQ  1|2:21:6:23,27  2|1:2:0:18,2    
20      10237   .         T     .      47    0       NS=57;DP=257;AA=T                    GT:GQ:DP:HQ  0|0:54:7:56,60  0|0:48:4:51,51  
20      123456  microsat1 G     D4,IGA 50    0       NS=55;DP=250;AA=G                    GT:GQ:DP     0/1:35:4        0/2:17:2        
```

Often the header lines will have some explanation about the various columns in the VCF, including the confusing looking INFO column. Here's an explanation of the INFO column for the first entry in the example above (the example below is representing the same variant as above, "rs6054257", but the VCF was excerpted from a much larger experiment):

<img src="../img/vcf_3.png" width="600">

Below is another example with slightly different fields in the INFO column:

<img src="../img/vcf_2.png" width="600">

Now, let's take a look at the one we just generated:

```bash
$ cat ptA.vcf
```

How does this compare to the 2 examples we have seen so far? How does the ID column compare?

---
# Optional Steps
## Generating a summary of VCF files

```bash
bcftools stats ptA.vcf
```
