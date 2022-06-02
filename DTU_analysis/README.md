Differential Transcript Usage 
================
Maria Tocora
June, 2022

The following analyses are based on Love MI, Soneson C and Patro R. Swimming downstream: statistical analysis of differential transcript usage following Salmon quantification [version 3; peer review: 3 approved]. F1000Research 2018, 7:952 (https://doi.org/10.12688/f1000research.15398.3)

### Installing and loading packages
```{r Install and Load Packages}
###DRIMSeq (https://bioconductor.org/packages/release/bioc/html/DRIMSeq.html)
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("DRIMSeq")
### rnaseqDTU (https://bioconductor.org/packages/release/workflows/html/rnaseqDTU.html)
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("rnaseqDTU")
###tximport (https://bioconductor.org/packages/release/bioc/html/tximport.html)
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("tximport")
###stageR (https://bioconductor.org/packages/release/bioc/html/stageR.html)
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("stageR")

###Load packages
library(DRIMSeq)
library(rnaseqDTU)
library(tximport)
library(stageR)
```

### Importing counts into R/Bioconductor
This example in going to provide code for the activity (high vs low-activity bodyguards) but files for the Caste analysis are under the Caste section at the end. 
```{r Install and Load Packages}
samps <- read.csv(file.path("C:/Users/Paula/Desktop/BODYGUARD PROJECT/Transcriptome analysis/DESeq2_analysis/AllTranscripts/QuantFiles/Activity_level", "samplesActivity.csv"))
head(samps)
names(samps) <- c("sample_id","plant.id","code","treatment", "attack.avg", "condition")
samps$condition <- factor(samps$condition)
table(samps$condition)
files <- file.path("C:/Users/Paula/Desktop/BODYGUARD PROJECT/Transcriptome analysis/DESeq2_analysis/AllTranscripts/QuantFiles/Activity_level", samps$sample_id, "quant.sf")
names(files) <- samps$sample_id
head(files)
### We used the quant.sf files from Salmon (Samples' names i.e. NODE_1_length_27414_cov_892.910316_g0_i0)
```

### Workflow DRIMSeq
```{r Install and Load Packages}
txi <- tximport(files, type="salmon", txOut=TRUE,
 countsFromAbundance="scaledTPM")
cts <- txi$counts
cts <- cts[rowSums(cts) > 0,]
### This cts file has the original samples names from Salmon, nevertheless the file must match the samples' names from the txdf file obtained in the Transcript-to-gene mapping step in the pipeline. For more information check Love et al (2018). I recommend downloading the file and changing the names (i.e. write.csv(cts, file="cts.csv").

cts <- read.csv(file.path("C:/Users/Paula/Desktop/BODYGUARD PROJECT/Transcriptome analysis/DESeq2_analysis/AllTranscripts/QuantFiles/Activity_level", "cts.csv"), row.names = 1) ###Newfile with new names. 
txdf <- read.csv(file.path("C:/Users/Paula/Desktop/BODYGUARD PROJECT/Transcriptome analysis/DESeq2_analysis/AllTranscripts/QuantFiles/Activity_level", "txdf.csv")) ###File with gene names. 
all(rownames(cts) %in% txdf$TXNAME)
txdf <- txdf[match(rownames(cts),txdf$TXNAME),]
all(rownames(cts) == txdf$TXNAME)

### In order to run DRIMSeq, we build a data.frame with the gene ID, the feature (transcript) ID, and then columns for each of the samples:
counts <- data.frame(gene_id=txdf$GENEID,
 feature_id=txdf$TXNAME,
 cts)

### Create a dmDSdata object, with our counts and samps data.frames. Typing in the object name and pressing return will give information about the number of genes:

d <- dmDSdata(counts=counts, samples=samps)
d

### The dmDSdata object has a number of specific methods. Note that the rows of the object are gene-oriented, so pulling out the first row corresponds to all of the transcripts of the first gene:

methods(class=class(d))
counts(d[1,])[,1:4]

### Filtering the object, before running procedures to estimate model parameters

n <- 14 ###total number of samples
n.small <- 7 ###sample size of the smallest group

### The following filters follow the default parameters in Love et al (2018)

d <- dmFilter(d,
 min_samps_feature_expr=n.small, min_feature_expr=10,
 min_samps_feature_prop=n.small, min_feature_prop=0.1,
 min_samps_gene_expr=n, min_gene_expr=10)
d

### We can find out how many of the remaining genes have N isoforms by tabulating the number of times we see a gene ID, then tabulating the output again:

table(table(counts(d)$gene_id))

### We create a design matrix, using a design formula and the sample information:

design_full <- model.matrix(~condition, data=DRIMSeq::samples(d))
colnames(design_full)

### We then use the following three functions to estimate the model parameters and test for DTU: 

set.seed(1)
system.time({
 d <- dmPrecision(d, design=design_full)
 d <- dmFit(d, design=design_full)
 d <- dmTest(d, coef="conditionlow")
})

### To build a results table, we run the results function. We can generate a single p-value per gene, which tests whether there is any differential transcript usage within the gene, or a single p-value per transcript, which tests whether the proportions for this transcript changed within the gene:
  
res <- DRIMSeq::results(d)
head(res)
res.txp <- DRIMSeq::results(d, level="feature")
head(res.txp)
write.csv(res, file="DTU.DRIMSeq.withDecoy.csv")


### Because the pvalue column may contain NA values, we use the following function to turn these into 1’s. The NA values would otherwise cause problems for the stage-wise analysis: 
  
no.na <- function(x) ifelse(is.na(x), 1, x)
res$pvalue <- no.na(res$pvalue)
res.txp$pvalue <- no.na(res.txp$pvalue)

### Plot the estimated proportions for the significant genes, where we can see evidence of switching
idx <- which(res$adj_pvalue < 0.05)[1]
res[idx,]
plotProportions(d, res$gene_id[idx], "condition") 
### You can change the idx value to check for different genes.
```
![alt text](https://github.com/mariatocora/Transcriptomic-analysis-ant-plant/blob/main/DTU_analysis/README_figures/g7229_DTU.png)

### stageR following DRIMSeq Analysis
The stage-wise analysis has been adopted from (Heller et al. 2009) and consists of a screening stage and a confirmation stage. In the screening stage, genes are screened by calculating p-values that aggregate evidence across the different hypotheses of interest for the gene. The screening p-values are then adjusted for FDR control after which significance of the screening hypothesis is assessed. In the confirmation stage, only genes passing the screening stage are considered for analysis. 

```{r Install and Load Packages}
### We show below how stageR is used to detect DTU and how to interpret its output.
### We first construct a vector of p-values for the screening stage: 

pScreen <- res$pvalue
strp <- function(x) substr(x,1,15)
names(pScreen) <- strp(res$gene_id)

### We construct a one column matrix of the confirmation p-values:

pConfirmation <- matrix(res.txp$pvalue, ncol=1)
rownames(pConfirmation) <- strp(res.txp$feature_id)

### We arrange a two column data.frame with the transcript and gene identifiers.

tx2gene <- res.txp[,c("feature_id", "gene_id")]
for (i in 1:2) tx2gene[,i] <- strp(tx2gene[,i])

### The following functions then perform the stageR analysis. We must specify an alpha, which will be the overall false discovery rate target for the analysis, defined below. Unlike typical adjusted p-values or q-values, we cannot choose an arbitrary threshold later: after specifying alpha=0.05, we need to use 5% as the target in downstream steps.

stageRObj <- stageRTx(pScreen=pScreen, pConfirmation=pConfirmation,
 pScreenAdjusted=FALSE, tx2gene=tx2gene)
stageRObj <- stageWiseAdjustment(stageRObj, method="dtu", alpha=0.05)
suppressWarnings({
 drim.padj <- getAdjustedPValues(stageRObj, order=FALSE,
 onlySignificantGenes=TRUE)
})
head(drim.padj)

### The final table with adjusted p-values summarizes the information from the two-stage analysis, you can download the database with  write.csv(drim.padj, file="Activity_drim.padj.csv")

###The significant genes can be returned with the getSignificantGenes function (https://bioconductor.org/packages/devel/bioc/vignettes/stageR/inst/doc/stageRVignette.html#references)

StageR_Act_DTUSigGenes <- (getSignificantGenes(stageRObj))
head(StageR_Act_DTUSigGenes)
### Download the final database with write.csv(StageR_Act_DTUSigGenes, file="StageR_Activity_DTUSigGenes.csv")
```
The final database summarizes the list of significant genes from the two-stage stageR analysis (screening and confirmation). Only genes that
passed those filters are included in the table. Please consider that the returned adjusted p-values are based on a stage-wise testing approach and are only valid for the provided target  overall false discovery rate  (OFDR) level of 5% (alpha = 0.05).

|gen ID|FDR adjusted p-value|
| --- | --- |
| g15 |0.004044799
| g3338 | 1.31E-07 |
| g4350	| 0.042032382 |
| g4830	| 9.13E-05 |
| g4886	| 0.025726378 |
| g5024	| 0.042032382 |
| g5260	| 0.001595937 |
| g5456	| 1.41E-15 |
| g5511	| 1.53E-05 |
| g5823	| 0.041709246 |
| g7229	| 0.037695361 |
| g8013	| 0.001589396 |
| g8428	| 0.025726378 |

### OPTIONAL: Post-hoc filtering on the standard deviation in proportions
Love et al (2018) describe the following Post-hoc, non-specific filtering of the DRIMSeq transcript p-values and adjusted p-values, that can improve the FDR and OFDR control considered the standard deviation (SD) of the per-sample proportions as a filtering statistic: 

```{r Install and Load Packages}
res.txp.filt <- DRIMSeq::results(d, level="feature")
smallProportionSD <- function(d, filter=0.1) {
 cts <- as.matrix(subset(counts(d), select=-c(gene_id, feature_id)))
 gene.cts <- rowsum(cts, counts(d)$gene_id)
 total.cts <- gene.cts[match(counts(d)$gene_id, rownames(gene.cts)),]
 props <- cts/total.cts
 propSD <- sqrt(rowVars(props))
 propSD < filter
}
filt <- smallProportionSD(d)
res.txp.filt$pvalue[filt] <- 1
res.txp.filt$adj_pvalue[filt] <- 1
```

### Caste
The following are modifications to the code we need to consider in order to perform the analysis between broodcare workers and bodyguards.

```{r Install and Load Packages}
samps <- read.csv(file.path("C:/Users/Paula/Desktop/BODYGUARD PROJECT/Transcriptome analysis/DESeq2_analysis/AllTranscripts/QuantFiles/Caste_level", "samplesCaste.csv"))
names(samps) <- c("sample_id","plant.id","treatment.code","condition", "attack.avg", "activity.level", "date", "time")
files <- file.path("C:/Users/Paula/Desktop/BODYGUARD PROJECT/Transcriptome analysis/DESeq2_analysis/AllTranscripts/QuantFiles/Caste_level", samps$sample_id, "quant.sf")
cts <- read.csv(file.path("C:/Users/Paula/Desktop/BODYGUARD PROJECT/Transcriptome analysis/DESeq2_analysis/AllTranscripts/QuantFiles/Caste_level", "cts_caste.csv"), row.names = 1)
txdf <- read.csv(file.path("C:/Users/Paula/Desktop/BODYGUARD PROJECT/Transcriptome analysis/DESeq2_analysis/AllTranscripts/QuantFiles/Caste_level", "txdf.csv"))
n <- 12 ### The number of samples is different from the Activity analysis
n.small <- 6 ### The n.small to be the sample size of the smallest group is different from the Activity analysis. 
set.seed(1)
system.time({
 d <- dmPrecision(d, design=design_full)
 d <- dmFit(d, design=design_full)
 d <- dmTest(d, coef="conditionguard") ###The condition label changes. 
})
```

## __References__
1. DRIMSeq: Nowicka M, Robinson MD (2016). “DRIMSeq: a Dirichlet-multinomial framework for multivariate count outcomes in genomics [version 2; referees: 2 approved].” F1000Research, 5(1356). doi: 10.12688/f1000research.8900.2, https://f1000research.com/articles/5-1356/v2. 
2. tximport: Charlotte Soneson, Michael I. Love, Mark D. Robinson. Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences, F1000Research, 4:1521, December 2015. doi: 10.12688/f1000research.7563.1
3. rnaseqDTU:Love MI, Soneson C, Patro R (2018). “Swimming downstream: statistical analysis of differential transcript usage following Salmon quantification.” F1000Research. doi: 10.12688/f1000research.15398.3.
4. stageR:Van den Berge K, Clement L (2022). stageR: stageR: stage-wise analysis of high throughput gene expression data in R. R package version 1.18.0
5. stageR: Heller, Ruth, Elisabetta Manduchi, Gregory R Grant, and Warren J Ewens. 2009. “A flexible two-stage procedure for identifying gene sets that are differentially expressed.” Bioinformatics (Oxford, England) 25 (8): 1019–25. https://doi.org/10.1093/bioinformatics/btp076.
