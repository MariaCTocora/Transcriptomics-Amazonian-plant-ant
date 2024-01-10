### DGE comparison across broodcare workers and bodyguards 

###Imput data 
###CASTE(Broodcare workers vs high-quality)
```{r echo=TRUE}
setwd("~/Documents/Trancriptomics A. octoarticulatus/DEAnalyses/DGE")
countdata <- read.csv("caste(BROODvsHIGH)withDecoy.csv", header = TRUE, row.names = 1)
coldata <- read.delim("caste(BROODvsHIGH)_modified.txt", header = TRUE, row.names = 1)
ddsFullCountTable <- DESeqDataSetFromMatrix(countData=round(countdata), 
                                  colData=coldata, 
                                  design=~treatment)
ddsFullCountTable
```

###CASTE(Broodcare workers vs low-quality)
```{r echo=TRUE}
countdata <- read.csv("caste(BROODvsLOW)withDecoy.csv", header = TRUE, row.names = 1)
coldata <- read.delim("caste(BROODvsLOW)_modified.txt", header = TRUE, row.names = 1)
ddsFullCountTable <- DESeqDataSetFromMatrix(countData=round(countdata), 
                                  colData=coldata, 
                                  design=~treatment)
ddsFullCountTable
```

##Differential expression analysis
###CASTE
```{r echo=TRUE}
dds <- ddsFullCountTable
dds$treatment <- relevel(dds$treatment, "broodcare")
as.data.frame(colData(dds))
dds <- DESeq(dds)
res <- results(dds)
res
```

##Extracting significant differentially expressed genes
###CASTE(Broodcare workers vs high-quality)
```{r echo=TRUE}
### Extracting significant differentially expressed genes
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
###write.csv(sigOE, file="Caste(BroodcareVSBodyguards(high))_sigOE.csv")
normalized_counts <- counts(dds, normalized=T) ###normalized counts
norm_OEsig <- normalized_counts[rownames(sigOE),] ###Extract normalized expression for significant genes
###write.csv(norm_OEsig, file="Caste(BroodcareVSBodyguards(high))_norm_OEsig.csv")
```

###CASTE(Broodcare workers vs low-quality)
```{r echo=TRUE}
### Extracting significant differentially expressed genes
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
###write.csv(sigOE, file="Caste(BroodcareVSBodyguards(low)_sigOE.csv")
normalized_counts <- counts(dds, normalized=T) ###normalized counts
norm_OEsig <- normalized_counts[rownames(sigOE),] ###Extract normalized expression for significant genes
###write.csv(norm_OEsig, file="Caste(BroodcareVSBodyguards(low)_norm_OEsig.csv")
```

###Plotting sig genes
###CASTE(Broodcare workers vs high-quality)
```{r echo=TRUE}
### Plot significant genes: 
heat.colors <- brewer.pal(6, "Greys") ### Set a color palette
annotation <- data.frame(treatment=coldata[,'treatment'], 
                     row.names=rownames(coldata))
my_colour = list(treatment = c("bodyguard(high)" = "#0000ff", broodcare = "#ffd500"))

HeatmapHvsB <- pheatmap(norm_OEsig, color = heat.colors, cluster_rows = T, show_rownames=F, annotation= annotation, border_color=NA, fontsize = 10, scale="row", fontsize_row = 10, height=20, legend_labels = F, annotation_names_row = F, annotation_colors = my_colour, annotation_names_col = F)
```

###CASTE(Broodcare workers vs low-quality)
```{r echo=TRUE}
### Plot significant genes: 
heat.colors <- brewer.pal(6, "Greys") ### Set a color palette
annotation <- data.frame(treatment=coldata[,'treatment'], 
                     row.names=rownames(coldata))
my_colour = list(treatment = c("bodyguard(low)" = "#ff1493", broodcare = "#ffd500"))

HeatmapLvsB <- pheatmap(norm_OEsig, color = heat.colors, cluster_rows = T, show_rownames=F, annotation= annotation, border_color=NA, fontsize = 10, scale="row", fontsize_row = 10, height=20, legend_labels = F, annotation_names_row = F, annotation_colors = my_colour, annotation_names_col = F)
```

### Enrichment Analysis 
```{r echo=TRUE}
###GSEA-Analysis
NES1 <- read.csv("~/Documents/METAANALYSIS/Experimental/Meta-analysis Graphs/Funct.Enrich/NES(Caste).csv")
NES1

NES1 <- ggplot(data = NES1, aes(x = GO.Category, y = GO.Name, 
                        color = -log10(`FDR`), size = Size)) + geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  ggtitle("BroodvsBodyguards - GSEA") + theme_cowplot(12) 
NES1
```

### Getting significant genes for all three treatments 

```{r echo=TRUE}
countdata <- read.csv("CASTEgene_count_matrix_collapsedTranscriptomeWithDecoys.csv", header = TRUE, row.names = 1)
coldata <- read.delim("caste_modified.txt", header = TRUE, row.names = 1)
ddsFullCountTable <- DESeqDataSetFromMatrix(countData=round(countdata), 
                                  colData=coldata, 
                                  design=~treatment)
ddsFullCountTable
dds <- ddsFullCountTable
dds$treatment <- relevel(dds$treatment, "brood")
as.data.frame(colData(dds))
dds <- DESeq(dds)
res <- results(dds)
res
### Extracting significant differentially expressed genes
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
####write.csv(sigOE, file="Caste_sigOE.csv")
normalized_counts <- counts(dds, normalized=T) ###normalized counts
norm_OEsig <- normalized_counts[rownames(sigOE),] ###Extract normalized expression for significant genes
```

### Extracting significant genes

```{r echo=TRUE}
###Candidate genes
res2 <- results(dds, tidy=TRUE, contrast=c("treatment", "brood", "guard")) %>%
  arrange(padj, pvalue) %>%
  tbl_df()
res2

##Vitellogenin 
####write.csv(res2, file="caste(Activity)_res2.csv")
which(grepl("g8459", res2$row))

goi <- res2$row[12]
tcounts <- t(log2((counts(dds[goi, ], normalized=TRUE, replaced=FALSE)+.5))) %>%
  merge(colData(dds), ., by="row.names") %>%
  gather(gene, expression, (ncol(.)-length(goi)+1):ncol(.))

g8459 <- ggboxplot(tcounts, y = "expression", x = "activity.level", fill = "activity.level") +  labs(x="labour", y="Expression (log normalized counts)", title = "Vitellogenin", fill = "Activity level") + theme(plot.title = element_text(color="black", size=15, hjust = 0.5)) +  scale_y_continuous(limits=c(0, 13)) + theme(legend.position="none") + annotate(geom="text", x =1.5, y=8, label="p = 2.14E-06", color="black") + scale_x_discrete(limits = c("brood", "high", "low")) + scale_fill_manual(values = c("#ff1493", "#ffd500", "#0000ff")) + annotate(geom="text", x =2, y=10, label="p = 0.0001", color="black") + annotate(geom="text", x =2.5, y=8, label="p = n.s.", color="black") ##BvsH(2.14E-06)##BvsL(0.0001)
g8459

##Where is the significant value? 
sigOE["g8459",]

##MRJP 1-like(g8960) 
####write.csv(res2, file="caste(Activity)_res2.csv")
which(grepl("g8960", res2$row))

goi <- res2$row[88]
tcounts <- t(log2((counts(dds[goi, ], normalized=TRUE, replaced=FALSE)+.5))) %>%
  merge(colData(dds), ., by="row.names") %>%
  gather(gene, expression, (ncol(.)-length(goi)+1):ncol(.))

g8960 <- ggboxplot(tcounts, y = "expression", x = "activity.level", fill = "activity.level") +  labs(x="labour", y="Expression (log normalized counts)", title = "MRJP 1-like", fill = "Activity level") + theme(plot.title = element_text(color="black", size=15, hjust = 0.5)) +  scale_y_continuous(limits=c(0, 11)) + theme(legend.position="none") + annotate(geom="text", x =1.5, y=10, label="p = 0.003", color="black") + scale_x_discrete(limits = c("brood", "high", "low")) + scale_fill_manual(values = c("#ff1493", "#ffd500", "#0000ff")) + annotate(geom="text", x =2, y=11, label="p = 1.76E-15", color="black") + annotate(geom="text", x =2.5, y=9, label="p = n.s.", color="black") ##BvsH(0.003)##BvsL(1.76E-15)
g8960

##JH biosynthesis (g5161) 
####write.csv(res2, file="caste(Activity)_res2.csv")
which(grepl("g5161", res2$row))

goi <- res2$row[144]
tcounts <- t(log2((counts(dds[goi, ], normalized=TRUE, replaced=FALSE)+.5))) %>%
  merge(colData(dds), ., by="row.names") %>%
  gather(gene, expression, (ncol(.)-length(goi)+1):ncol(.))

g5161 <- ggboxplot(tcounts, y = "expression", x = "activity.level", fill = "activity.level") +  labs(x="labour", y="Expression (log normalized counts)", title = "JH biosynthesis", fill = "Activity level") + theme(plot.title = element_text(color="black", size=15, hjust = 0.5)) +  scale_y_continuous(limits=c(0, 10)) + theme(legend.position="none") + annotate(geom="text", x =1.5, y=9, label="p = n.s.", color="black") + scale_x_discrete(limits = c("brood", "high", "low")) + scale_fill_manual(values = c("#ff1493", "#ffd500", "#0000ff")) + annotate(geom="text", x =2, y=10, label="p = 3.30E-05", color="black") + annotate(geom="text", x =2.5, y=9, label="p = n.s.", color="black") ##BvsH(0.003)##BvsL(3.3015613963483E-05)
g5161
```

### Serine-threonine protein kinase genes

```{r echo=TRUE}
###GATA
which(grepl("g508", res2$row))

goi <- res2$row[109]
tcounts <- t(log2((counts(dds[goi, ], normalized=TRUE, replaced=FALSE)+.5))) %>%
  merge(colData(dds), ., by="row.names") %>%
  gather(gene, expression, (ncol(.)-length(goi)+1):ncol(.))

g508 <- ggboxplot(tcounts, y = "expression", x = "activity.level", fill = "activity.level") +  labs(x="labour", y="Expression (log normalized counts)", title = "GATA", fill = "Activity level") + theme(plot.title = element_text(color="black", size=15, hjust = 0.5)) +  scale_y_continuous(limits=c(0, 13)) + theme(legend.position="none") + annotate(geom="text", x =1.5, y=11.5, label="p = 0.0004", color="black") + scale_x_discrete(limits = c("brood", "high", "low")) + scale_fill_manual(values = c("#ff1493", "#ffd500", "#0000ff")) + annotate(geom="text", x =2, y=12.5, label="p = n.s.", color="black") + annotate(geom="text", x =2.5, y=11.5, label="p = n.s.", color="black") ##BvsH(0.0004)##BvsL(n.s.)
g508

###DYRK1A (g5480)
which(grepl("g5480", res2$row))

goi <- res2$row[43]
tcounts <- t(log2((counts(dds[goi, ], normalized=TRUE, replaced=FALSE)+.5))) %>%
  merge(colData(dds), ., by="row.names") %>%
  gather(gene, expression, (ncol(.)-length(goi)+1):ncol(.))

g5480 <- ggboxplot(tcounts, y = "expression", x = "activity.level", fill = "activity.level") +  labs(x="labour", y="Expression (log normalized counts)", title = "DYRK1A", fill = "Activity level") + theme(plot.title = element_text(color="black", size=15, hjust = 0.5)) +  scale_y_continuous(limits=c(0, 11.5)) + theme(legend.position="none") + annotate(geom="text", x =1.5, y=10.5, label="p = 2.59E-12", color="black") + scale_x_discrete(limits = c("brood", "high", "low")) + scale_fill_manual(values = c("#ff1493", "#ffd500", "#0000ff")) + annotate(geom="text", x =2, y=11.5, label="p = 0.002", color="black") + annotate(geom="text", x =2.5, y=7, label="p = n.s.", color="black") ##BvsH(2.5864791371359E-12)##BvsL(0.002)
g5480
```




