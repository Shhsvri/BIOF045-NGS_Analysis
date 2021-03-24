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
row_data <- cts[,3]
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

Okay, now we load the data into DESeq2
```
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = meta,
                              design= ~ treatment+experiment)
```







