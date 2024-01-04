---
DGE Activity comparison results
---
```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```
Installing and loading packages

```{r echo=TRUE}
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
  
### Pheatmap (https://cran.r-project.org/web/packages/pheatmap/index.html)
install.packages("pheatmap")

### tidyverse (https://www.tidyverse.org/)
install.packages("tidyverse")

### RColorBrewer (https://cran.r-project.org/web/packages/RColorBrewer/RColorBrewer.pdf)
install.packages("RColorBrewer")

###BigPint
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("bigPint")

library("DESeq2")
library(tidyverse)
library(RColorBrewer)
library(EnhancedVolcano)
library("pheatmap") 
```

##Imput data
###ACTIVITY
```{r echo=TRUE}
setwd("~/Documents/Trancriptomics A. octoarticulatus/DEAnalyses/DGE")
countdata <- read.csv("ACTIVITYgene_count_matrix_collapsedTranscriptomeWithDecoys.csv", header = TRUE, row.names = 1)
coldata <- read.delim("activitylevel_modified.txt", header = TRUE, row.names = 1)
all(rownames(coldata) == colnames(countdata))
### Need TRUE as output 
ddsFullCountTable <- DESeqDataSetFromMatrix(countData=round(countdata), 
                                  colData=coldata, 
                                  design=~activity.level)
ddsFullCountTable
```

```{r echo=TRUE}
dds <- ddsFullCountTable
dds$activity.level <- relevel(dds$activity.level, "low")
as.data.frame(colData(dds))
dds <- DESeq(dds)
res <- results(dds)
res
```

##Extracting significant differentially expressed genes
###ACTIVITY
```{r echo=TRUE}
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
###write.csv(sigOE, file="Activity_sigOE.csv")
normalized_counts <- counts(dds, normalized=T) ###normalized counts
norm_OEsig <- normalized_counts[rownames(sigOE),] ###Extract normalized expression for significant genes
###write.csv(norm_OEsig, file="Activity_norm_OEsig.csv")
```

```{r echo=TRUE}
### Run pheatmap
heat.colors <- brewer.pal(6, "Greys") ### Set a color palette

annotation <- data.frame(activity.level=coldata[,'activity.level'], 
                     row.names=rownames(coldata))
my_colour = list(activity.level = c(high = "#0000ff", low = "#ff1493"))

ActivityHeatmap <- pheatmap(norm_OEsig, color = heat.colors, cluster_rows = T, show_rownames=F,annotation= annotation, border_color=NA, fontsize = 10, scale="row", fontsize_row = 10, height=20, annotation_colors = my_colour, legend_labels = F, annotation_names_row = F, annotation_names_col = F)

###Figure modified in Adobe
```

```{r echo=TRUE}
NES <- read.csv("~/Documents/METAANALYSIS/Experimental/Meta-analysis Graphs/Funct.Enrich/NES(Act).csv", header = TRUE, row.names = 1) ###First row as columns
NES

####GO NAME
PlotNES <- ggplot(data = NES, aes(x = GO.Category, y = GO.Name, 
                        color = -log10(`FDR`), size = Size)) + geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  theme_bw() + 
  ylab("") + 
  xlab("") + theme_cowplot(12) 
```

```{r echo=TRUE}
res2 <- results(dds, tidy=TRUE, contrast=c("activity.level", "high", "low")) %>%
  arrange(padj, pvalue) %>%
  tbl_df()
res2

#Metalloproteinase (g16223)
##Where is the gene we want to plot? 
which(grepl("g16223", res2$row))
##plot the gene normalized counts
library(ggpubr)
goi <- res2$row[45]
tcounts <- t(log2((counts(dds[goi, ], normalized=TRUE, replaced=FALSE)+.5))) %>%
    merge(colData(dds), ., by="row.names") %>%
    gather(gene, expression, (ncol(.)-length(goi)+1):ncol(.))
g16223 <- ggboxplot(tcounts, y = "expression", x = "activity.level", fill = "activity.level") +  labs(x="labour", y="Expression (log normalized counts)", title = "Metalloproteinase", fill = "Activity level") + theme(plot.title = element_text(color="black", size=15, hjust = 0.5)) +  scale_fill_manual(values = c("#ff1493", "#0000ff")) + theme(legend.position="none") + annotate(geom="text", x =1.5, y=10, label="p = 0.006", color="black")
#check ylim in the tcounts dataframe
###For no legend: theme(legend.position="none")
##Where is the significant value? 
sigOE["g16223",]

#Chymotrypsin (g7229)
##Where is the gene we want to plot? 
which(grepl("g7229", res2$row))
##plot the gene normalized counts
library(ggpubr)
goi <- res2$row[23]
tcounts <- t(log2((counts(dds[goi, ], normalized=TRUE, replaced=FALSE)+.5))) %>%
    merge(colData(dds), ., by="row.names") %>%
    gather(gene, expression, (ncol(.)-length(goi)+1):ncol(.))
g7229 <- ggboxplot(tcounts, y = "expression", x = "activity.level", fill = "activity.level") +  labs(x="labour", y="Expression (log normalized counts)", title = "Chymotrypsin 1-like", fill = "Activity level") + theme(plot.title = element_text(color="black", size=15, hjust = 0.5)) +  scale_y_continuous(limits=c(0, 13)) + scale_fill_manual(values = c("#ff1493", "#0000ff")) + theme(legend.position="none") + annotate(geom="text", x =1.5, y=13, label="p = 0.001", color="black")
##Where is the significant value? 
sigOE["g7229",]
```

```{r echo=TRUE}
#Read in png figure
image <- readPNG("venn_diagram.png")

#Set png figure up in ggplot
venn <- ggdraw() +
  draw_image(
    image, scale = 1, x = 1, y = 0.1,
    hjust = 1, halign = 1, valign = 0
  )  

#Read in png figure
image2 <- readPNG("ActivityHeatmap1.png")

#Set png figure up in ggplot
system_fig <- ggdraw() +
  draw_image(
    image2, scale = 1, x = 1, y = 0.1,
    hjust = 1, halign = 1, valign = 0
  )  

##List of plots
##system_fig
##PlotNes
#g16223
#g7229 
##venn

#Assemble multi-panel figure
up_row <- plot_grid(system_fig, PlotNES, labels = c('A', 'B'), label_size = 12, ncol = 2, rel_widths = c(1, 2))
bottom_row <- plot_grid(g16223, g7229, labels = c('C', 'D'), label_size = 12)
bottom_row
bottom_row <- plot_grid(bottom_row, venn, labels = c('', 'E'), label_size = 12, ncol = 2, rel_widths = c(2, 1))
Figure_3 <- plot_grid(up_row, bottom_row, label_size = 12, ncol = 1, rel_heights = c(1.3,1.1))
#Save plot
save_plot("Figure3.pdf", Figure_3, base_width=10, base_height=7)
```
