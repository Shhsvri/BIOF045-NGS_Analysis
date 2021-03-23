---
title: "Variant annotation with snpEff"
author: "Shahin Shahsavari"
date: "March 2021"
---



## Learning Objectives:

* Adding information on known SNPs to our VCF
* Adding functional information to variants in the VCF
* Understanding where the annotation is added in the VCF format



## Annotating variants

Variant annotation is a crucial step in linking sequence variants with changes in phenotype. Annotation results can have a strong influence on the ultimate conclusions of disease studies. Incorrect or incomplete annotations can cause researchers both to overlook potentially disease-relevant DNA variants and to dilute interesting variants in a pool of false positives. 

<img src="../img/variant_calling_workflow_3.png" width="450">

At this stage, we have a large tab-delimited file containing loci at which a variation was found in the sample DNA sequence relative to the reference. We have filtered out these variations (also referred to as 'variant calls') to keep only those we are highly confident in, and now need to find out more. We can do this by **comparing our variants against known variants, and also use genome annotations to help predict information about our variants.** 

<img src="../img/prioritize.png" width="200">


### Setting up

Let's create a new directory for the results of our annotation steps:

```
$ cd ~/var_calling
$ mkdir results/annotation
```


## Annotation with known variants 

Variant annotation is the process of assigning information to DNA variants. There are many different types of information that can be associated with variants, and a first commonly used resource is using databases which contain variants that have previously been described. One popular example is [dbSNP](http://www.ncbi.nlm.nih.gov/SNP/),a free, public archive for genetic variation within and across different species. It is hosted by NCBI in collaboration with NHGRI and although the name implies SNPs; it actually includes range of molecular variation.

<img src="../img/dbsnp.png" width="500">

To add dbSNP information you need to download the organism specific data using their FTP download. **We have already done this for you** and was the zip file that you copied over into your `reference_data` folder. 

> For the full dbSNP dataset that corresponds with our data you can download it via the FTP site: ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/ . You may notice is that there are alot of files, and the README does not provide enough information. To find out more on how the datasets differ you can access the [NCBI human variation docs](http://www.ncbi.nlm.nih.gov/variation/docs/human_variation_vcf/).  

To annotate our data with dbSNP information we wil be using [`bcftools`](https://samtools.github.io/bcftools/), a command-line utility for variant calling and manipulating VCF files and its binary counterpart BCF files. It is a part of the `samtools` project, a tool that we are by now pretty familiar with. 

The `bcftools annotate` command allows the user to **add or remove annotations**. 

```bash
$ bcftools annotate --help
```

The annotation we wish to add and the file we are annotating must be a Bgzip-compressed and tabix-indexed file (usually VCF or BED format). Tabix indexes a tab-delimited genome position file and creates an index file (.tbi), which facilitates quick retrieval of data lines overlapping regions. *NOTE: this has already been done for our dbSNP file*

```bash
$ bgzip ~/var_calling/results/variants/08008_q20.recode.vcf
$ tabix ~/var_calling/results/variants/08008_q20.recode.vcf.gz
```

When running `bcftools annotate`, we also need to specify the column(s) to carry over from the annotation file, which in our case is ID.

```bash
$ bcftools annotate -c ID \
-a /data/DNAseq/dbSNP/dbSNP.151.38.gz \
~/var_calling/results/variants/08008_q20.recode.vcf.gz \
> ~/var_calling/results/annotation/08008_q20_annot.vcf
```

Take a quick peek at the new VCF file that was generated using `less`. You should now see in the ID column `rs` ids which correspond to identifiers from dbSNP. For the variants that are not known, you will find the `.` in place of an ID indicating novelty of that sequence change.

## Functional annotation with SnpEff

One fundamental level of variant annotation involves categorising each variant based on its relationship to coding sequences in the genome and how it may change the coding sequence and affect the gene product. To do this we will be using a tool called [SnpEff](http://snpeff.sourceforge.net/), a **variant effect predictor program**. 

Our understanding of the protein-coding sequences in the genome is summarised in the set of transcripts we believe to exist. Thus, **variant annotation depends on the set of transcripts used as the basis for annotation**. The widely used annotation databases and browsers – ENSEMBL, RefSeq, UCSC – contain sets of transcripts that can be used for variant annotation, as well as a wealth of information of many other kinds as well, such as ENCODE data about the function of non-coding regions of the genome. SnpEff will take information from the provided annotation database and populate our VCF file by adding it into the `INFO` field name `ANN`. Data fields are encoded separated by pipe sign "\|"; and the order of fields is written in the VCF header.

<img src="../img/snpeff.png" width="700">

Some **common annotations** are listed below, but [the manual](http://snpeff.sourceforge.net/SnpEff_manual.html#input)
provides a more comprehensive list.

* Putative_impact/impact: A simple estimation of putative impact / deleteriousness : (HIGH, MODERATE, LOW, MODIFIER)
* Gene Name: Common gene name (HGNC). Optional: use closest gene when the variant is “intergenic”.
* Feature type: Which type of feature (e.g. transcript, motif, miRNA, etc.). It is preferred to use Sequence Ontology (SO) terms, but ‘custom’ (user defined) are allowed. 
* Feature ID: Depending on the annotation sources, this may be: Transcript ID (preferably using version number), Motif ID, miRNA, ChipSeq peak, Histone mark, etc. Note: Some features may not have ID (e.g. histone marks from custom Chip-Seq experiments may not have a unique ID).
* Biotype: The bare minimum is at least a description on whether the transcript is (“Coding”, “Noncoding”). Whenever possible, use ENSEMBL biotypes.

Take a look at the options available. We will be using `snpEff` to annotate our variants.

```bash

$ cd results/annotation

$ java -jar $SNPEFF/snpEff.jar -h
```

> `snpEff` is a java based tool, similar to `picard` and we have to use a similar syntax to run it.

To run the snpEff command we will need to specify two things:

1. The appropriate genome
2. The VCF file we want to annotate

An additional parameter to add to our command is `Xmx8G`, a Java parameter to define available memory. Since we are in an interactive session with 8GB, if we had requested more before starting the session we could increase the number here.

The final command will look like this:

```bash
## DO NOT RUN THIS CODE

$ java -Xmx8G -jar $SNPEFF/snpEff.jar eff hg19 \
na12878_q20_annot.vcf \
> na12878_q20_annot_snpEff.vcf
```	
