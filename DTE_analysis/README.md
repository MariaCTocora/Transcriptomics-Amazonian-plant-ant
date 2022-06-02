Differential Transcript Expression 
================
Maria Tocora
June, 2022

The following analyses are modified from the Swish method for differential expression analysis of bulk or single-cell RNA-seq data using inferential replicate counts is described in the following reference: Zhu et al. (2019) doi: 10.1093/nar/gkz622, and available here https://bioconductor.org/packages/release/bioc/vignettes/fishpond/inst/doc/swish.html#Differential_transcript_expression

### Installing and loading packages
```{r Install and Load Packages}
###tximeta (https://bioconductor.org/packages/release/bioc/html/tximeta.html)
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("tximeta")

###Fishpond (https://bioconductor.org/packages/release/bioc/html/fishpond.html)
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("fishpond")

###SummarizedExperiment (https://bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html)
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("SummarizedExperiment")

###Load packages
library(tximeta)
library(fishpond)
library(SummarizedExperiment)
```

## __Activity Analysis__

### Importing data

In this step we are going to create a SummarizedExperiment object that is a matrix-like container where rows represent features of interest (e.g. genes, transcripts, exons, etc...) and columns represent samples (with sample data summarized as a DataFrame). This object contains one or more assays, each represented by a matrix-like object of numeric or other mode. The following code describes a DTE Analysis for the activity comparison (high vs low-activity bodyguards), modifications of the code for the caste analysis is stated at the end. 

```{r Install and Load Packages}
## We start by reading in a CSV with the column data, that is, information about the samples, which are represented as columns of the SummarizedExperiment object we will construct containing the counts of reads per gene or transcript.

coldata <- read.csv(file.path("C:/Users/Paula/Desktop/BODYGUARD PROJECT/Transcriptome analysis/DESeq2_analysis/AllTranscripts/QuantFiles/Activity_level", "samplesActivity.csv"))
head(coldata)
names(coldata) <- c("names","plant.id","code","treatment", "attack.avg", "condition")
head(coldata)

## coldata needs to have a column files which specifies the path to the quantification files.

coldata$files <- file.path("C:/Users/Paula/Desktop/BODYGUARD PROJECT/Transcriptome analysis/DESeq2_analysis/AllTranscripts/QuantFiles/Activity_level", coldata$names, "quant.sf")
all(file.exists(coldata$files))

### We use the tximeta (Love et al. 2020) package to read in the data:
se <- tximeta(coldata)
y <- se
y <- y[,y$condition %in% c("low","high")]
y$condition <- factor(y$condition, levels=c("low","high"))
```

### Running Swish at the transcript level
Running swish has three steps: scaling the inferential replicates, labeling the rows with sufficient counts for running differential expression, and then calculating the statistics.
```{r Install and Load Packages}
y <- scaleInfReps(y)
y <- labelKeep(y) ###labelKeep by default  keep features with minN=3 samples with a minimal count of 10
y <- y[mcols(y)$keep,]
set.seed(1)
y <- swish(y, x="condition") ###The default number of permutations for computing p-values is nperms=100
```
We can see how many transcripts are in a 5% FDR set:

```{r Install and Load Packages}
table(mcols(y)$qvalue < .05)
### With the following code chunk, we construct two vectors that give the significant genes with the lowest (most negative) and highest (most positive) log2 fold changes: 
with(mcols(y),
     table(sig=qvalue < .05, sign.lfc=sign(log2FC))
     )
```

| FALSE | TRUE |
| --- | --- |
| 39646 | 4 | 

Table 1. Transcripts in the 5% FDR set.

|     |    | sign.lfc |
| --- | --- | ---|
| sig | -1 | 1 |
| FALSE | 19700 | 19946 |
| TRUE | 0 | 4 |
  
 Table 2. Significant genes with the lowest (most negative) and highest (most positive) log2 fold changes

### Get significant genes list
```{r Install and Load Packages}
sig <- mcols(y)$qvalue < .05
sig2 <- mcols(y)[(sig),c("log2FC","qvalue")] ###to get just the significant genes
sig3 <- print(as.data.frame(sig2)) ##Print the significant genes as a list in a dataframe (write.csv(sig3, file="Activity_DTESigGenes.csv"))
```

| gene ID | log2FC | qvalue |
| --- | --- | --- |
| NODE_29252_length_1645_cov_6183.978342_g7216_i2 | 3.188244 | 0.0025 |
| NODE_39386_length_1080_cov_3627.980019_g13731_i1 | 4.734310 | 0.0025 |
| NODE_46557_length_810_cov_19810.451985_g6591_i2 | 5.379874 | 0.0025 |
| NODE_47477_length_782_cov_2782.565737_g9305_i1 | 3.402328 | 0.0025 |

### Plot of distributtion of p-values
There is not an enrichment of transcripts with p-values near 0 (Figure 1). 
```{r Install and Load Packages}
hist(mcols(y)$pvalue, col="grey")
```

![alt text](https://github.com/mariatocora/Transcriptomic-analysis-ant-plant/blob/main/DTE_analysis/README_Figures/Activity_Distribution_of_p-values.png)

### Plots
Plot the scaled counts for the inferential replicates, and also group the samples by a covariate. The analysis was paired, so the statistic assessed if the change within pairs was consistent (Figure 2). 
```{r Install and Load Packages}
hi <- order(-mcols(y)$log2FC * sig)
plotInfReps(y, idx=hi[1], x="condition")
plotInfReps(y, idx=hi[2], x="condition")
plotInfReps(y, idx=hi[3], x="condition")
plotInfReps(y, idx=hi[4], x="condition")
```
![alt text](https://github.com/mariatocora/Transcriptomic-analysis-ant-plant/blob/main/DTE_analysis/README_Figures/DTE_g6591.png)

We can make an MA plot, where the transcripts in our FDR set are colored (Figure 3.):
```{r Install and Load Packages}
We can make an MA plot, where the transcripts in our FDR set are colored (Figure 3.): 
```

![alt text](https://github.com/mariatocora/Transcriptomic-analysis-ant-plant/blob/main/DTE_analysis/README_Figures/MAplot_DTE_Activity.png)

## __Caste Analysis__


## __References__
1. Tximeta: Love MI, Soneson C, Hickey PF, Johnson LK, Pierce NT, Shepherd L, Morgan M, Patro R (2020). “Tximeta: Reference sequence checksums for provenance identification in RNA-seq.” PLOS Computational Biology, 16, e1007664. doi: 10.1371/journal.pcbi.1007664.
2. fishpond: Zhu A, Srivastava A, Ibrahim JG, Patro R, Love MI (2019). “Nonparametric expression analysis using inferential replicate counts.” Nucleic Acids Research, 47, e105. doi: 10.1093/nar/gkz622, https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6765120.
3. SummarizedExperiment: Morgan M, Obenchain V, Hester J, Pagès H (2022). SummarizedExperiment: SummarizedExperiment container. R package version 1.26.1, https://bioconductor.org/packages/SummarizedExperiment.
