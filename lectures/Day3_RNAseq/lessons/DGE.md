### Analysis of RNA-seq in R

### Libraries
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("DESeq2")) install.packages("DESeq2")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("colorspace")) install.packages("colorspace")

### DESeq2
DESeq2 uses the negative binomial distribution to model counts after correcting for between sample variability. 

It requires a counts matrix and a design matrix. 

### Negative binomial distribution


### Mean ratio normalization 


### Count matrix
Here, we're not working directly with the full count matrix. Instead, we're taking a short-cut and using a publicly avaible count matrix (GSE153311).

As such, it has some extra columns. 
```
cts <- read.table("GSE153310_Raw_gene_counts_matrix.txt",header=TRUE)
```
Output
```
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
