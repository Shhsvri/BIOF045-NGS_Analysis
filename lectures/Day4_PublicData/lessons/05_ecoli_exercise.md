---
title: "Extras"
author: "Shahin Shahsavari"
date: "October, 2020"
---

## Mapping a public Ecoli fastq

#### 1. Download the reference genome

```bash
cd ~/Day5/ecoli_exercise
ls
cd genome
wget http://hypervolu.me/%7Eerik/genomes/E.coli_K12_MG1655.fa
```

> Use bwa to index your genome


#### 2. Download the SRR1770413 Ecoli data

```bash
cd ~/Day5/ecoli_exercise/raw_data
fastq-dump --split-files SRR1770413
```

#### 3. Align your reads using bwa mem

```bash
cd ~/Day5/ecoli_exercise
```
> Use `bwa` to align your raw data
> Use `samtools` to compress, sort, and index your sam file

#### 4. Perform variant calling

> Use `freebayes` to call variants
> **NOTE:** with freebayes, you need to set the ploidy option to 1. You can find the documentation [here](https://github.com/freebayes/freebayes)


---

* This exercise was developed by Shahin Shahsavari adapted from Erik Garrison's Tutorial on Variant Calling available on [GitHub](https://github.com/ekg/alignment-and-variant-calling-tutorial)

