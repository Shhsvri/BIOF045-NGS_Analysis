—
title: Differential gene expression
author: Jonathan Vi Perrie
date: 03/23/2021
duration: ~2-3h 
—


### Analysis of RNA-seq in R

### Libraries
```R
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("DESeq2")) BiocManager::install("DESeq2")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("colorspace")) install.packages("pheatmap")
if (!require("colorspace")) install.packages("colorspace")
if (!require("org.Hs.eg.db")) BiocManager::install("org.Hs.eg.db")
```

### DESeq2
DESeq2 uses the negative binomial distribution to model counts after correcting for between sample variability. 

It requires a counts matrix and a design matrix. 

### Negative binomial distribution


### Mean ratio normalization 


### Count matrix
Here, we're not working directly with the full count matrix. Instead, we're taking a short-cut and using a publicly avaible count matrix (GSE153311).

As such, it has some extra columns. 
```R
cts <- read.table("GSE153310_Raw_gene_counts_matrix.txt",header=TRUE)
```
Output
```R
> colnames(cts)
 [1] "st_gene_id"            "gene_id"               "gene_symbol"           "exp_AS1"               "exp_AS2"               "exp_AS3"              
 [7] "exp_NC1"               "exp_NC2"               "exp_NC3"               "read_count_AS1"        "read_count_AS2"        "read_count_AS3"       
[13] "read_count_NC1"        "read_count_NC2"        "read_count_NC3"        "average_exp_AS"        "average_exp_NC"        "average_read_count_AS"
[19] "average_read_count_NC"
```
We want the read counts, as DESeq2 works with count data, so lets select columns
```
# headers 
exp_headers <- colnames(cts)[4:9]
# row names
row_data <- cts[,2]
# raw count matrix
cts <- round(cts[,10:15])
rownames(cts) <-row_data
```
### Design matrix 
We need to specify where our replicates and treatments/controls are.
```
# There are six experiments
experiments<-seq(6)
# Three replicates per condition
replicates<-c(1,2,3,1,2,3)
# Two conditions (
treatment<-c(rep(1,3),rep(2,3))

meta$experiment<-as.factor(meta$experiment)
meta$treatment<-as.factor(meta$treatment)
meta$replicate<-as.factor(meta$replicate)

meta
```
Output
```
> meta
  experiment treatment replicate
1          1         1         1
2          2         1         2
3          3         1         3
4          4         2         1
5          5         2         2
6          6         2         3
```
Not every gene will be high quality and by removing some genes with low counts, we can improve DESeq2's statistical power. Here, we filter out some genes. This is a much bigger problem in scRNA-seq where dropout is common. 

First, we look at the distribution of gene counts summed across all experiments. 
```
hist(log(rowSums(cts)),breaks=100)
quantile(rowSums(cts),probs=c(0.01,0.05,0.10,0.20))
cts <- cts[rowSums(cts)>=36,]
```
![alt text](tmp.png)

Okay, now we load the data into DESeq2. The design of this experiment is pretty simply as we are just comparing two populations. But more complex designs could be considered in temporal experiments. e.g. `design= ~ treatment + time`
```
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = meta,
                              design= ~ treatment)
```

We can then run DESeq2 and order the results by adjusted p-value. 
```
dds <- DESeq(dds)
res <- results(dds)
res <- res[order(res$padj),]
head(res)
```
Output, the DESeq2 data structure has 6 columns with each row corresponding to a gene (ENTREZID)
1. baseMean: average normalized count / size factors over all samples
2. log2FoldChange: effect size estimate of treatment
3. lfcSE: standard error estimate of fold change
4. stat: Wald statistic
5. pvalue: Wald test p-value
6. padj: Behnjamini-Hochberg p-value = p-value rank / # of tests * FDR and find largest p-value smaller than critical value
```
> head(res)
log2 fold change (MLE): treatment 2 vs 1 
Wald test p-value: treatment 2 vs 1 
DataFrame with 6 rows and 6 columns
            baseMean log2FoldChange     lfcSE      stat      pvalue        padj
           <numeric>      <numeric> <numeric> <numeric>   <numeric>   <numeric>
64129       365.5798      -1.212309  0.144442  -8.39306 4.73631e-17 6.31161e-13
4046        832.9958       1.553207  0.210053   7.39437 1.42085e-13 9.46711e-10
3371      11609.6847       1.642778  0.272421   6.03029 1.63662e-09 7.26986e-06
8324       1041.9838      -0.693306  0.119186  -5.81700 5.99134e-09 1.99601e-05
57419       158.1467      -2.411077  0.467300  -5.15959 2.47487e-07 6.59603e-04
101928589    11.1888       6.922034  1.367017   5.06360 4.11405e-07 9.13731e-04
```






