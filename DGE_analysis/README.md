Differential Gene Expression (DGE) Analysis 
================
Maria Tocora
June, 2022

The following pipeline is modified from "Analyzing RNA-seq data with DESeq2" by Love et al (2022) available here: http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#differential-expression-analysis

### Installing and loading packages
```{r Install and Load Packages}
### DESeq2 (https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("DESeq2")

### EnhacedVolcano (https://bioconductor.org/packages/release/bioc/html/EnhancedVolcano.html) and (https://bioconductor.org/packages/devel/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html)
devtools::install_github('kevinblighe/EnhancedVolcano')
#OR
if (!requireNamespace('BiocManager', quietly = TRUE))
    install.packages('BiocManager')
  BiocManager::install('EnhancedVolcano')
  
### tidyverse (https://www.tidyverse.org/)
install.packages("tidyverse")

### RColorBrewer (https://cran.r-project.org/web/packages/RColorBrewer/RColorBrewer.pdf)
install.packages("RColorBrewer")


library("DESeq2")
library(tidyverse)
library(RColorBrewer)
library(ggplot2)
library(EnhancedVolcano)
library("pheatmap")

  
````
### Imput data
We manually obtained a matrix of read counts from the files obtained in Salmon and imported this matrix under "countdata". I also imported the sample information table "coldata". The "countdata" matrix has the samples information in columns organized in the same order as the rows in the "coldata" table.  then, construct a DESeqDataSet with the DESeqDataSetFromMatrix() function.  

```{r Install and Load Packages}
countdata <- read.csv("activity_count_matrix.csv", header = TRUE, row.names = 1)
coldata <- read.delim("activity_samples.txt", header = TRUE, row.names = 1)
all(rownames(coldata) == colnames(countdata))
### Need TRUE as output 
ddsFullCountTable <- DESeqDataSetFromMatrix(countData=round(countdata), 
                                  colData=coldata, 
                                  design=~activity.level)
ddsFullCountTable
```

The DESeqDataSet got as output is described as follows: 

- class: DESeqDataSet 
- dim: 67613 14 
- metadata(1): version
- assays(1): counts
- rownames(67613): g0 g1 ... g9998 g9999
- rowData names(0):
- colnames(14): sample1 sample2 ... sample17 sample19
- colData names(5): plant.id code treatment attack.avg activity.level

### Differential expression analysis
Differential gene expression analysis based on the negative binomial distribution described in Love et al (2014). For mor information check the Methods section of the paper. 

```{r Install and Load Packages}
dds <- ddsFullCountTable
dds$activity.level <- relevel(dds$activity.level, "low")
as.data.frame(colData(dds))
dds <- DESeq(dds)
res <- results(dds)
res
````
The output obtained "res" containes: 
- log2 fold change (MLE): activity.level high vs low 
- Wald test p-value: activity.level high vs low 
- DataFrame with 67613 rows and 6 columns

## Extracting significant differentially expressed genes
The following pipeline modifies code from https://hbctraining.github.io/Intro-to-R-with-DGE/lessons/B1_DGE_visualizing_results.html
Please consider that thet The lfc.cutoff is set to 0.58, so since we are working with log2 fold changes so this translates to an actual fold change of 1.5; and the padj.cutoff is set to 0.01. The final database "Activity_sigOE.csv" is available in the folder "Data"  

```{r Install and Load Packages}
res_ordered = res[order(res$pvalue),]
head(res_ordered)
summary(res_ordered)
sum(res_ordered$padj < 0.01, na.rm=TRUE)
### Setting up parameters: 
padj.cutoff <- 0.01
lfc.cutoff <- 0.58
threshold <- res_ordered$padj < padj.cutoff & abs(res_ordered$log2FoldChange) > lfc.cutoff ### vector containing genes that meet our criteria (parameters stated above). 
length(which(threshold)) ### the vector has a length equal to the total number of significant genes in the dataset.
res_ordered$threshold <- threshold 
sigOE <- data.frame(subset(res_ordered, threshold==TRUE)) ## subsetting significant genes
write.csv(sigOE, file="Activity_sigOE.csv")
normalized_counts <- counts(dds, normalized=T) ###normalized counts
norm_OEsig <- normalized_counts[rownames(sigOE),] ###Extract normalized expression for significant genes
```

## Plotting Sig. genes in a heatmap 
In this section I extract the normalized values of all the significant genes and plot a heatmap of their expression using pheatmap()

```{r Install and Load Packages}
### Run pheatmap
heat.colors <- brewer.pal(6, "Greys") ### Set a color palette

annotation <- data.frame(activity.level=coldata[,'activity.level'], 
                     row.names=rownames(coldata))
my_colour = list(activity.level = c(high = "#BA55D3", low = "#FFD700"))

ActivityHeatmap <- pheatmap(norm_OEsig, color = heat.colors, cluster_rows = T, show_rownames=T,
annotation= annotation, border_color=NA, fontsize = 10, scale="row",
     fontsize_row = 10, height=20, annotation_colors = my_colour)

```

![alt text](https://github.com/mariatocora/Transcriptomic-analysis-ant-plant/blob/main/DGE_analysis/Figures/Activity_heatmap.png) 

Figure 1. Heatmap with significant genes. Please consider the figure was edited. 

## Visualization of top significant genes
In this section I plotted the top differentially expressed genes considering a log2FoldChange equal to 2 and padj value of 10e-6 as cutoff parameters. I used the R package "EnhacedVolcano" (Blighe, et al., 2018). Please consider that the log transformed adjusted p-values are plotted on the y-axis and log2 fold change values on the x-axis. 

```{r Install and Load Packages}
keyvals <- ifelse(
    res$log2FoldChange < -2 & res$padj < 10e-6, 'royalblue',
      ifelse(res$log2FoldChange > 2 & res$padj < 10e-6, 'gold',
        'black')) 
  keyvals[is.na(keyvals)] <- 'black'
  names(keyvals)[keyvals == 'gold'] <- 'Upregulated high active bodyguards'
  names(keyvals)[keyvals == 'black'] <- 'NS/Only Log2 FC/Only p-value'
  names(keyvals)[keyvals == 'royalblue'] <- 'Upregulated low active bodyguards'
  
 ActivityVolcanoPlot <- EnhancedVolcano(res,
    lab = rownames(res),
    x = 'log2FoldChange',
    y = 'pvalue',
    selectLab = rownames(res)[which(names(keyvals) %in% c('Upregulated high active bodyguards', 'Upregulated low active bodyguards'))],
    xlab = bquote(~Log[2]~ 'fold change'),
    title = 'Activity Analysis',
    pCutoff = 10e-6,
    FCcutoff = 2.0,
    pointSize = 4,
    labSize = 4.5,
    boxedLabels = TRUE,
    colCustom = keyvals,
    colAlpha = 1,
    legendPosition = 'top',
    legendLabSize = 15,
    legendIconSize = 5.0,
    drawConnectors = TRUE,
    widthConnectors = 1.0,
    colConnectors = 'black',
    arrowheads = TRUE,
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
    border = 'partial',
    borderWidth = 1.5,
    borderColour = 'black')
ActivityVolcanoPlot
```

![alt text] (https://github.com/mariatocora/Transcriptomic-analysis-ant-plant/blob/main/DGE_analysis/Figures/Activity_Volcano_plot.png)


## __References__

1. DESeq2: Love MI, Huber W, Anders S (2014). “Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.” Genome Biology, 15, 550. doi: 10.1186/s13059-014-0550-8.
2. EnhacedVolcano: Blighe, K, S Rana, and M Lewis. 2018. “EnhancedVolcano: Publication-ready volcano plots with enhanced colouring and labeling.” https://github.com/kevinblighe/EnhancedVolcano.
