# DGE Review

<img src="../img/jq1_bet.png" width="800">

Here, we will be reviewing the commands we used yesterday, but we will be using a cleaner data set that looks at how BET inhibition affects cancer cells. 

## Data

### Libraries 
```R
# load in the libraries
library(DESeq2)
library(ggplot2)
library(colorspace)
library(pheatmap)
library(tidyverse)
```
### Read in the data
```R
setwd("~/Day4/combined_counts2/")
counts <- read.table("counts2.tsv",header=TRUE,sep="\t",row.names=1)
metadata <- read.table("metadata2.tsv",sep="\t",header=TRUE,row.names=1)
View(metadata)
View(counts)
# Change column names of counts to match rownames of metadata
colnames(counts) = rownames(metadata)
```
### Look at distribution
```R
hist(log(rowSums(counts)),breaks=100)
quantile(rowSums(counts),probs=c(0.01,0.05,0.10,0.50))
counts_filtered <- counts[rowSums(counts)>=11,]
```

### Load data in DESeq2 with appropriate design 
```R
> dds <- DESeqDataSetFromMatrix(countData = counts_filtered,
				colData = metadata,
				design= ~ treatment)
```

### Run DESeq2
```R
dds <- DESeq(dds)
res <- results(dds)
res_ordered <- res[order(res$padj),]
```

### Transform data for quality control
```R
rld <- rlog(dds, blind=FALSE)
plotPCA(rld, intgroup="treatment")

colnames(rld) <- rownames(metadata)
pheatmap(cor(assay(rld),method="spearman"),display_numbers=TRUE,annotation_col=metadata,
         number_format='%.4f',cluster_rows=FALSE,cluster_cols=FALSE)
```

### Evaluate some genes related to BET and cancer
Note that we have made it simpler for all of you, so you don't have to query the data base to find gene names. In reality however you will need to make a conversion between gene name formats. Also, this is mouse data, so if one were to `mapIds`, one would have to load db as `db <- org.Mm.eg.db`. See https://bioconductor.org/packages/release/data/annotation/html/org.Mm.eg.db.html for installation. The suffix number `_#`, so don't worry.
```R
par(mfrow=c(2,2))
plotCounts(dds,gene="Apc_37324",intgroup="treatment",main="Apc")
plotCounts(dds,gene="Brd4_35665",intgroup="treatment",main="BRD4")
plotCounts(dds,gene="Myc_33094",intgroup="treatment",main="MYC")
plotCounts(dds,gene="Pms2_14724",intgroup="treatment",main="PMS2")
```

### Volcano plot 
```R
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj')
```

### Heatmap of top 20 most variable genes
```R
topVarGenes <- head(order(-rowVars(assay(rld))),20)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
pheatmap(mat,annotation_col=metadata,cluster_rows=FALSE,cluster_cols=FALSE, show_rownames=FALSE)
```

### Heatmap of top 20 most differentially expressed genes
```R
topDiffGenes <- head(order(assay(rld)),20)
mat <- assay(rld)[ topDiffGenes, ]
mat <- mat - rowMeans(mat)
pheatmap(mat,annotation_col=metadata,cluster_rows=FALSE,cluster_cols=FALSE, show_rownames=FALSE)
```

### correlation plot
```R
pheatmap(cor(assay(rld), method="spearman"),
         annotation_col=metadata,
         cluster_rows=FALSE,
         cluster_cols=FALSE,
         number_format='%.4f',
         display_numbers=TRUE)
```
