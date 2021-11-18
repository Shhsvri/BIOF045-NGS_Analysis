# load in the libraries
library(DESeq2)
library(ggplot2)
library(colorspace)
library(pheatmap)
library(tidyverse)

# Change Directory
setwd("/data/Day4/combined_counts/")

# Read in the Count Matrix and Metadata file
counts <- read.table("counts.tsv")
metadata <- read.table("metadata.tsv")

# Filter Genes with Low Counts (across all samples)
rowSums(counts) %>% log() %>% hist(breaks=100, main="Sum of Read Counts Across all Samples")
rowSums(counts) %>% quantile(probs=c(0.01,0.05,0.10,0.20))
## remove the bottom 10th percentile based on the sum of gene expression
counts_filtered <- counts[rowSums(counts)>=5,]

# How many genes were filtered?
nrow(counts)
nrow(counts_filtered)
nrow(counts) - nrow(counts_filtered)

View(counts_filtered)

# Create a DESeq object
## These are the steps the DESeq2 author recommends
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata,
                              design= ~ treatment)
dds <- DESeq(dds)
res <- results(dds)

# Order genes based on adjusted p-Values (the p-values are calculated by DESeq)
res_ordered <- res[order(res$padj),]
head(res_ordered)

# Generate a PCA
rld <- rlog(dds, blind=FALSE)
plotPCA(rld, intgroup="treatment")

# Generate a heatmap of the genes with the highest variance
topVarGenes <- rowVars(assay(rld)) %>% order() %>% tail(20)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
pheatmap(mat,annotation_col=metadata,cluster_rows=FALSE,cluster_cols=FALSE)
