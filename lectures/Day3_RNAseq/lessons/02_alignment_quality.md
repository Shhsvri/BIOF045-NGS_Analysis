---
title: "Assessing the output of STAR"
author: "Shahin Shahsavari"
date: "July 2021"
---

Approximate time: 30 minutes

## Learning objectives

* Evaluating the STAR aligner output files
* Understanding the standard alignment file (SAM/BAM) structure
* Using `samtools` to evaluate alignment quality 
* Visualizing alignment quality using IGV (genome browser)  


## Assessing alignment quality

After running our single FASTQ file through the STAR aligner, you should have noticed a number of output files in the `~/Day3/results/` directory. Let's turn this **sam** file into **bam**, then sort, and lastly index it.

```bash
$ cd ~/Day3/results
$ samtools view -b sample1_Aligned.out.sam  > sample1.bam
$ samtools sort sample1.bam > sample1.sorted.bam
$ samtools index sample1.sorted.bam
```

Then let's take a quick look at some of the files that were generated and explore the content of some of them. 

```bash
$ cd ~/Day3/results	
$ ls -l -h
```

What you should see, is that for each FASTQ file you have **5 output files** and a single tmp directory. Briefly, these files are described below:

* `Log.final.out` - a summary of mapping statistics for the sample
* `Aligned.sortedByCoord.out.bam` - the aligned reads, sorted by coordinate, in BAM format
* `Log.out` - a running log from STAR, with information about the run 
* `Log.progress.out` -  job progress with the number of processed reads, % of mapped reads etc., updated every ~1 minute
* `SJ.out.tab` - high confidence collapsed splice junctions in tab-delimited format. Only junctions supported by uniquely mapping reads are reported

## Mapping statistics

Having completed the alignment, the first thing we want to know is how well did our reads align to the reference. Rather than looking at each read alignment, it can be more useful to evaluate statistics that give a general overview for the sample. One of the output files from the STAR aligner contains mapping statistics, let's take a closer look at one of those files. We'll use the `head` command which allows us to scroll through it easily: 

```bash
$ head sample1_Log.final.out
```	

The log file provides information on reads that 1) mapped uniquely, 2) reads that mapped to mutliple locations and 3) reads that are unmapped. Additionally, we get details on splicing, insertion and deletion. From this file the most informative statistics include the **mapping rate and the number of multimappers**.

* As an example, a good quality sample will have **at least 75% of the reads uniquely mapped**. Once values start to drop lower than 60% it's advisable to start troubleshooting. The lower the number of uniquely mapping reads means the higher the number of reads that are mapping to multiple locations. It is best to keep this number low because multi-mappers are not included when we start counting reads

> NOTE: The thresholds suggested above will vary depending on the organism that you are working with. Much of what is discussed here is in the context of working with human or mouse data. For example, 75% of mapped reads holds true only if the genome is good or mature. For badly assembled genomes we may not observe a high mapping rate, even if the actual sequence sample is good.


## Other quality checks

In addition to the aligner-specific summary we can also obtain quality metrics using tools like [Qualimap](http://qualimap.bioinfo.cipf.es/doc_html/intro.html#what-is-qualimap) or [RNASeQC](https://www.broadinstitute.org/publications/broad4133). These tools examine sequencing alignment data according to the features of the mapped reads and their genomic properties and **provide an overall view of the data that helps to to the detect biases in the sequencing and/or mapping of the data**.The input can be one or more BAM files and the output consists of HTML or PDF reports with useful figures and tab delimited files of metrics data.

We will not be performing this step in the workshop, but we describe some of the features below to point out things to look for when assessing alignment quality of RNA-seq data:

* **Reads genomic origin**: Even if you have high genomic mapping rate for all samples, check to see where the reads are mapping. Ensure that there is not an unusually high number of **reads mapping to intronic regions** (~30% expected) and fewer than normally observed **mapping to exons** (~55%). A high intronic mapping suggests possible genomic DNA contamination and/or pre-mRNA. 
* **Ribosomal RNA (rRNA)** constitutes a large majority of the RNA species in any total RNA preparation. Despite depletion methods, you can never achieve complete rRNA removal. Even with Poly-A enrichment a small percentage of ribosomal RNA can stick to the enrichment beads non-specifically. **Excess ribosomal content (> 2%)** will normally have to be filtered out so that differences in rRNA mapped reads across samples do not affect alignment rates and skew subsequent normalization of the data.
* **Transcript coverage and 5'-3' bias**: assesing the affect on expression level and on levels of transcript GC content
* **Junction analysis**: analysis of junction positions in spliced alignments (i.e known, partly known, novel) 
* **Strand specificity:** assess the performance of strand-specific library construction methods. The percentage of sense-derived reads is given for each end of the read pair. A non-strand-specific protocol would give values of 50%/50%, whereas strand-specific protocols typically yield 99%/1% or 1%/99% for this metric.


### Viewing the SAM file

By default, `samtools view` excludes the header lines.`samtools view -h` writes the header alongside the read entries to the terminal.
``

```
$ samtools view -H
$ samtools view -h sample1.sorted.bam | head -n 20

``` 

> ### Filtering the SAM file
> 
> Now we know that we have all of this information for each of the reads -- wouldn't it be useful to summarize and filter based on selected criteria? Suppose we wanted to set a **threshold on mapping quality**. For example, we want to know how many reads aligned with a quality score higher than 30. To do this, we can combine the `view` command with additional flags `q 30` and `-c` (to count):
> 
> ```
> $ samtools view -q 30 sample1.sorted.bam
> 
> ```
> *How many of reads have a mapping quality of 30 or higher?*
> 
> We can also **apply filters to select reads based on where they fall within the `FLAG` categories**. Remember that the bitwise flags are like boolean values. If the flag exists, the statement is true. Similar to when filtering by quality we need to use the `samtools view` command, however this time use the `-F` or `-f` flags.
> 
> * `-f` - to find the reads that agree with the flag statement 
> * `-F`  - to find the reads that do not agree with the flag statement
> 
> ```
> ## This will tell us how many reads are unmapped
> $ samtools view -f 4 -c sample1.sorted.bam
> 
> ## This should give us the remaining reads that do not have this flag set (i.e reads that are mapped)
> $ samtools view -F 4 -c sample1.sorted.bam
> ```


****


*This lesson has been developed by Shahin Shahsavari adapted from the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
