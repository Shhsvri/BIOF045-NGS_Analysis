---
title: DESeq2 Normalization Steps
author: Shahin Shahsavari
date: 03/16/2022
duration: 60 minutes 
---

## DESeq2 Normalization Outline

We will use a large gene count matrix to create plots of out counts data before and after
DESeq2 normalization.

### Load the libraries that you need

```R
# load in the libraries
library(DESeq2)
library(ggplot2)
library(colorspace)
library(pheatmap)
library(tidyverse)
library(affy)
library(gplots)
```
### Read in the combined count matrix

```R
# Read the counts matrix and its associated metadata file
setwd("~/Day4/DESeq2_normalization")
countTable <- read.table("yeast_counts.tsv", header=TRUE, row.names=1)
phenoTable <- read.table("yeast_metadata.tsv", header=TRUE, row.names=1)
```

Then look at the content of each dataframe

```R
# Min, Max, median (...). 
# Here on the first 4 samples

countTable[1:6,1:4]
summary(countTable[1:6,1:4])
```

### Distributions

The summary only displays a few milestone values (mean, median, quartiles). In order to get
a better intuition of the data, we can draw an histogram of all values.

```R
# create a histogram of the raw counts data, and the log2 of raw counts
par(mfrow=c(2,1))
hist(as.matrix(countTable), col="blue", border="white",
     breaks=20000, xlim=c(0,2000), main="Counts per gene",
     xlab="Counts (truncated axis)", ylab="Number of genes", 
     las=1, cex.axis=0.7)

epsilon <- 1 # pseudo-count to avoid problems with log(0)
hist(as.matrix(log2(countTable + epsilon)), breaks=100, col="blue", border="white",
     main="Log2-transformed counts per gene", xlab="log2(counts+1)", ylab="Number of genes", 
     las=1, cex.axis=0.7)
```


### Interpretation

* The top histogram is not very informative so far, apparently due to the presence of a few
very high count values, that impose a very large scale on the X axis.
* The middle histogram shows the representative range. Note the height of the first bin,
which includes the zero counts.
* The logarithmic transformation (bottom histogram) improves the readability. Note that we
added a pseudo-count of 1 to avoid problems with the log transformation of zero counts (which gives −∞).

To get better insights into the distribution per sample, the boxplot offer a good perspective.

```R
# create a boxplot of the log2(conutsTable)
boxplot(log2(countTable + epsilon), col=phenoTable$color, pch=".", 
        horizontal=TRUE, cex.axis=0.5,
        las=1, ylab="Samples", xlab="log2(counts +1)")
```

Another way to get an intuition of the value distributions is to use the plotDensity function, which draws one density curve for each sample.

```R
# create a plotDensity of the log2(countsTable)
plotDensity(log2(countTable + epsilon), lty=1, col=phenoTable$color, lwd=2)
grid()
legend("topright", legend=names(strainColor), col=strainColor, lwd=2)
```


## Normalization

```R
# Load the data as a DESeq2 object and estimate SizeFactors (aka normalize)
dds0 <- DESeqDataSetFromMatrix(countData = countTable, colData = phenoTable, design = ~ strain)
dds.norm <-  estimateSizeFactors(dds0)

# Pick a sample of 10 WT samples and 10 Snf2 samples
nb.replicates = 10
samples.WT <- sample(1:48, size=nb.replicates, replace=FALSE)
# Random sampling of the Snf2 replicates (columns 49 to 96)
samples.Snf2 <- sample(49:96, size=nb.replicates, replace=FALSE)
selected.samples <- c(samples.WT, samples.Snf2)
# Don't forget to update colors
col.pheno.selected <- phenoTable$color[selected.samples]


# Plot the counts before and after normalization to compare
par(mfrow=c(2,2),cex.lab=0.7)
boxplot(log2(counts(dds.norm)+epsilon),  col=col.pheno.selected, cex.axis=0.7, 
        las=1, xlab="log2(counts+1)", horizontal=TRUE, main="Raw counts")
boxplot(log2(counts(dds.norm, normalized=TRUE)+epsilon),  col=col.pheno.selected, cex.axis=0.7, 
        las=1, xlab="log2(normalized counts)", horizontal=TRUE, main="Normalized counts") 
plotDensity(log2(counts(dds.norm)+epsilon),  col=col.pheno.selected, 
            xlab="log2(counts+1)", cex.lab=0.7, panel.first=grid()) 
plotDensity(log2(counts(dds.norm, normalized=TRUE)+epsilon), col=col.pheno.selected, 
            xlab="log2(normalized counts)", cex.lab=0.7, panel.first=grid()) 
```



### Performing differential expression call
```R
# run estimateDispersions on the normalized data and run the
# negative binomial test to obrain p values
dds.disp <- estimateDispersions(dds.norm)
alpha <- 0.0001
waldTestResult <- nbinomWaldTest(dds.disp)
resultDESeq2 <- results(waldTestResult, alpha=alpha, pAdjustMethod="BH")
head(resultDESeq2)

# Order the table by decreasing p-values
resultDESeq2 <- resultDESeq2[order(resultDESeq2$padj),]
head(resultDESeq2)
```

### Volcano Plot
```R
alpha <- 0.01 # Threshold on the p-value

# Compute significance, with a maximum of 320 for the p-values set to 0 due to limitation of computation precision
resultDESeq2$sig <- -log10(resultDESeq2$padj)
sum(is.infinite(resultDESeq2$sig))
resultDESeq2[is.infinite(resultDESeq2$sig),"sig"] <- 350
genes.to.plot <- !is.na(resultDESeq2$pvalue)
## Volcano plot of adjusted p-values
cols <- densCols(resultDESeq2$log2FoldChange, resultDESeq2$sig)
cols[resultDESeq2$pvalue ==0] <- "purple"
resultDESeq2$pch <- 19
resultDESeq2$pch[resultDESeq2$pvalue ==0] <- 6
plot(resultDESeq2$log2FoldChange, 
     resultDESeq2$sig, 
     col=cols, panel.first=grid(),
     main="Volcano plot", 
     xlab="Effect size: log2(fold-change)",
     ylab="-log10(adjusted p-value)",
     pch=resultDESeq2$pch, cex=0.4)
abline(v=0)
abline(v=c(-1,1), col="brown")
abline(h=-log10(alpha), col="brown")

## Plot the names of a reasonable number of genes, by selecting those begin not only significant but also having a strong effect size
gn.selected <- abs(resultDESeq2$log2FoldChange) > 2 & resultDESeq2$padj < alpha 
text(resultDESeq2$log2FoldChange[gn.selected],
     -log10(resultDESeq2$padj)[gn.selected],
     lab=rownames(resultDESeq2)[gn.selected ], cex=0.6)
```

### Check the expression levels of the most differentially expressed gene

```R
# Plot counts for most diff genes
selectedGenes <- c(
  "Most significant" =  rownames(resultDESeq2)[which.max(resultDESeq2$sig)])

## Select a gene with small fold change but high significance
sel1 <- subset(
  na.omit(resultDESeq2), 
  sig >= 50 & log2FoldChange > 0 & log2FoldChange < 0.5)
selectedGenes <- append(selectedGenes, 
                         c("Small FC yet significant"=rownames(sel1)[1]))

## Select the non-significant gene with the highest fold change
sel2 <- subset(x = na.omit(resultDESeq2), padj > alpha & log2FoldChange > 0 & baseMean > 1000 & baseMean < 10000)
sel2 <- sel2[order(sel2$log2FoldChange, decreasing = TRUE),][1,]
selectedGenes <- append(
  selectedGenes, 
  c("Non-significant"=rownames(sel2)[1]))

## Plot the counts for the selected genes
par(mfrow=c(length(selectedGenes),1))
for (g in selectedGenes) {
  barplot(counts(dds.norm, normalized=TRUE)[g,], 
          col=phenoTable$color, 
          main=g, las=2, cex.names=0.5)
}
```

### Hierarchical clustering
```R
## generate a heatmap
gene.kept <- rownames(resultDESeq2)[resultDESeq2$padj <= alpha & !is.na(resultDESeq2$padj)]
countTable.kept <- log2(countTable + epsilon)[gene.kept, ]
heatmap.2(as.matrix(countTable.kept), 
          scale="row", 
          hclust=function(x) hclust(x,method="average"), 
          distfun=function(x) as.dist((1-cor(t(x)))/2), 
          trace="none", 
          density="none", 
          labRow="",
          cexCol=0.7)
```

---
This tutorial was adapted from Nicolas Terrapon, Denis Puthier and Jacques van Helden (adapted from tutorial from Hugo Varet, Julie Aubert and J. van Helden) which is available on [github](http://dputhier.github.io/jgb71e-polytech-bioinfo-app/practical/rna-seq_R/rnaseq_diff_Snf2.html)
