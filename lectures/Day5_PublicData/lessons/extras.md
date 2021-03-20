---
title: "Extras"
author: "Shahin Shahsavari"
date: "October, 2020"
---

## Mapping a public Ecoli fastq

#### 1. Download the reference genome

```bash
mkdir -p ~/var_calling/ecoli/{genome,results,raw_data}
cd ~/var_calling/ecoli/genome
wget http://igenomes.illumina.com.s3-website-us-east-1.amazonaws.com/Escherichia_coli_K_12_MG1655/NCBI/2001-10-15/Escherichia_coli_K_12_MG1655_NCBI_2001-10-15.tar.gz
tar -xvf Escherichia_coli_K_12_MG1655_NCBI_2001-10-15.tar.gz
ref=~/var_calling/ecoli/genome/Escherichia_coli_K_12_MG1655/NCBI/2001-10-15/Sequence/BWAIndex/genome.fa
```

#### 2. Download the SRR1770413 Ecoli data

```bash
cd ~/var_calling/ecoli/raw_data
fastq-dump --split-files SRR1770413
fq1=~/var_calling/ecoli/raw_data/SRR1770413_1.fastq
fq2=~/var_calling/ecoli/raw_data/SRR1770413_2.fastq
```

#### 3. Align your reads using bwa mem

```bash
cd ~/var_calling/ecoli/results
threads=2
bwa mem -t $threads -R '@RG\tID:K12\tSM:K12' $ref $fq1 $fq2 | samtools view -b > SRR1770413.bam
sambamba sort SRR1770413.bam
Overwrite the initial bam file
sambamba markdup -t 2 SRR1770413.sorted.bam SRR1770413.bam
```

> **NOTE** In the alignment step, we added the RG flag (Read Group). If we decide to run joint calling, it is required to keep track of our samples within the fastq file.

#### 4. Perform variant calling

```bash
freebayes -f $genome --ploidy 1 SRR1770413.bam > SRR1770413.vcf
```

## 5. Download the data to your local computer

If you decide to transfer the data to your local machine, you could use the scp command from your local computer:

```bash
scp username@server_ip:~/var_calling/ecoli/results /path/to/local_computer
```

---

* This exercise was developed by Shahin Shahsavari adapted from Erik Garrison's Tutorial on Variant Calling available on [GitHub](https://github.com/ekg/alignment-and-variant-calling-tutorial)

