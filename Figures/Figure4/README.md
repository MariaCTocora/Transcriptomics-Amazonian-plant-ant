### Inmune-related gene expression among workers (broodcare and bodyguard ants displaying different levels of aggressiveness)

##Getting significant genes for all three treatments 

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

```{r echo=TRUE}
###g68084_hymenoptaecin 
which(grepl("g68084", res2$row))

goi <- res2$row[91]
tcounts <- t(log2((counts(dds[goi, ], normalized=TRUE, replaced=FALSE)+.5))) %>%
  merge(colData(dds), ., by="row.names") %>%
  gather(gene, expression, (ncol(.)-length(goi)+1):ncol(.))

g68084 <- ggboxplot(tcounts, y = "expression", x = "activity.level", fill = "activity.level") +  labs(x="labour", y="Expression (log normalized counts)", title = "Hymenoptaecin", fill = "Activity level") + theme(plot.title = element_text(color="black", size=15, hjust = 0.5)) +  scale_y_continuous(limits=c(0, 10)) + theme(legend.position="none") + annotate(geom="text", x =1.5, y=8, label="p = n.s.", color="black") + scale_x_discrete(limits = c("brood", "high", "low")) + scale_fill_manual(values = c("#ff1493", "#ffd500", "#0000ff")) + annotate(geom="text", x =2, y=10, label="p = 2.67E-07", color="black") + annotate(geom="text", x =2.5, y=8, label="p = ns.s", color="black") ##BvsH(p = 9.89E-07)##BvsL(ns)#HvsL(ns)
g68084


###g9196_metallopeptidase(3)
which(grepl("g9196", res2$row))

goi <- res2$row[75]
tcounts <- t(log2((counts(dds[goi, ], normalized=TRUE, replaced=FALSE)+.5))) %>%
  merge(colData(dds), ., by="row.names") %>%
  gather(gene, expression, (ncol(.)-length(goi)+1):ncol(.))

g9196 <- ggboxplot(tcounts, y = "expression", x = "activity.level", fill = "activity.level") +  labs(x="labour", y="Expression (log normalized counts)", title = "Metallopeptidase", fill = "Activity level") + theme(plot.title = element_text(color="black", size=15, hjust = 0.5)) +  scale_y_continuous(limits=c(0, 14)) + theme(legend.position="none") + annotate(geom="text", x =1.5, y=8, label="p = 8.018E-06", color="black") + scale_x_discrete(limits = c("brood", "high", "low")) + scale_fill_manual(values = c("#ff1493", "#ffd500", "#0000ff")) + annotate(geom="text", x =2, y=10, label="p = 2.675E-07", color="black") + annotate(geom="text", x =2.5, y=8, label="p = n.s.", color="black") ##BvsH(p = 9.89E-07)##BvsL(ns)#HvsL(ns)
g9196

plot_grid(g68084, g9196, labels = c('C', ''), label_size = 12)

```
