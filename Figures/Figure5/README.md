
## Figure 4 - Compared viral expression between high and low activity bodyguards. 

#Import data
```{r echo=TRUE}
countdata <- read.csv("ACTIVITY_gene_count_matrix.csv", header = TRUE, row.names = 1)
coldata <- read.delim("activity.txt", header = TRUE, row.names = 1)
all(rownames(coldata) == colnames(countdata))
### Need TRUE as output 
ddsFullCountTable <- DESeqDataSetFromMatrix(countData=round(countdata), 
                                  colData=coldata, 
                                  design=~activity.level)
ddsFullCountTable
dds <- ddsFullCountTable
dds$activity.level <- relevel(dds$activity.level, "low")
as.data.frame(colData(dds))
dds <- DESeq(dds)
res <- results(dds)
res
```

#Find virues
```{r echo=TRUE}
res2 <- results(dds, tidy=TRUE, contrast=c("activity.level", "high", "low")) %>%
  arrange(padj, pvalue) %>%
  tbl_df()
res2

##Narnavirus (g9141) -(AOV-1)
which(grepl("g9141", res2$row))
##plot the gene normalized counts
library(ggpubr)
goi <- res2$row[7]
tcounts <- t(log2((counts(dds[goi, ], normalized=TRUE, replaced=FALSE)+.5))) %>%
    merge(colData(dds), ., by="row.names") %>%
    gather(gene, expression, (ncol(.)-length(goi)+1):ncol(.))
g9141<- ggboxplot(tcounts, y = "expression", x = "activity.level", fill = "activity.level") +  labs(x="labour", y="Expression (log normalized counts)", title = "Narnavirus", fill = "Activity level") + theme(plot.title = element_text(color="black", size=15, hjust = 0.5, face = "italic")) +  scale_fill_manual(values = c("#ff1493", "#0000ff")) + theme(legend.position="left") + annotate(geom="text", x =1.5, y=9, label="p = 6.45E-11", color="black")
g9141

##With no label 
g9141<- ggboxplot(tcounts, y = "expression", x = "activity.level", fill = "activity.level") +  labs(x="labour", y="Expression (log normalized counts)", title = "Narnavirus", fill = "Activity level") + theme(plot.title = element_text(color="black", size=15, hjust = 0.5, face = "italic")) +  scale_fill_manual(values = c("#ff1493", "#0000ff")) + theme(legend.position="none") + annotate(geom="text", x =1.5, y=9, label="p = 6.45E-11", color="black")
g9141

#Putative A. octoarticulatus virus 2 (g2275_i0) (AOV-2)
##Where is the gene we want to plot? 
which(grepl("g2275", res2$row))
##plot the gene normalized counts
library(ggpubr)
goi <- res2$row[45]
tcounts <- t(log2((counts(dds[goi, ], normalized=TRUE, replaced=FALSE)+.5))) %>%
    merge(colData(dds), ., by="row.names") %>%
    gather(gene, expression, (ncol(.)-length(goi)+1):ncol(.))
g2275 <- ggboxplot(tcounts, y = "expression", x = "activity.level", fill = "activity.level") +  labs(x="labour", y="Expression (log normalized counts)", title = "Putative A. octoarticulatus virus 1", fill = "Activity level") + theme(plot.title = element_text(color="black", size=15, hjust = 0.5, face = "italic")) +  scale_fill_manual(values = c("#ff1493", "#0000ff")) + theme(legend.position="none") + annotate(geom="text", x =1.5, y=10, label="p = 0.0011", color="black")
g2275
#check ylim in the tcounts dataframe
###For no legend: theme(legend.position="none")

#Putative A. octoarticulatus virus 3 (g5994_i0) (AOV-3)
##Where is the gene we want to plot? 
which(grepl("g5994", res2$row))
##plot the gene normalized counts
library(ggpubr)
goi <- res2$row[45]
tcounts <- t(log2((counts(dds[goi, ], normalized=TRUE, replaced=FALSE)+.5))) %>%
    merge(colData(dds), ., by="row.names") %>%
    gather(gene, expression, (ncol(.)-length(goi)+1):ncol(.))
g5994<- ggboxplot(tcounts, y = "expression", x = "activity.level", fill = "activity.level") +  labs(x="labour", y="Expression (log normalized counts)", title = "Putative A. octoarticulatus virus 2", fill = "Activity level") + theme(plot.title = element_text(color="black", size=15, hjust = 0.5, face = "italic")) +  scale_fill_manual(values = c("#ff1493", "#0000ff")) + theme(legend.position="none") + annotate(geom="text", x =1.5, y=10, label="p = 0.0001", color="black")
g5994

#Putative A. octoarticulatus virus 4 (g6982_i4) (AOV-4)
##Where is the gene we want to plot? 
which(grepl("g6982", res2$row))
##plot the gene normalized counts
library(ggpubr)
goi <- res2$row[4]
tcounts <- t(log2((counts(dds[goi, ], normalized=TRUE, replaced=FALSE)+.5))) %>%
    merge(colData(dds), ., by="row.names") %>%
    gather(gene, expression, (ncol(.)-length(goi)+1):ncol(.))
g6982<- ggboxplot(tcounts, y = "expression", x = "activity.level", fill = "activity.level") +  labs(x="labour", y="Expression (log normalized counts)", title = "Putative A. octoarticulatus virus 3", fill = "Activity level") + theme(plot.title = element_text(color="black", size=15, hjust = 0.5, face = "italic")) +  scale_fill_manual(values = c("#ff1493", "#0000ff")) + theme(legend.position="none") + annotate(geom="text", x =1.5, y=10, label="p = 7.85E-13", color="black")
g6982

##7.84660589394401E-13

#Putative A. octoarticulatus virus 5 (g1894_i0) (AOV-5)
##Where is the gene we want to plot? 
which(grepl("g1894", res2$row))
##plot the gene normalized counts
library(ggpubr)
goi <- res2$row[11]
tcounts <- t(log2((counts(dds[goi, ], normalized=TRUE, replaced=FALSE)+.5))) %>%
    merge(colData(dds), ., by="row.names") %>%
    gather(gene, expression, (ncol(.)-length(goi)+1):ncol(.))
g1894<- ggboxplot(tcounts, y = "expression", x = "activity.level", fill = "activity.level") +  labs(x="labour", y="Expression (log normalized counts)", title = "Putative A. octoarticulatus virus 4", fill = "Activity level") + theme(plot.title = element_text(color="black", size=15, hjust = 0.5, face = "italic")) +  scale_fill_manual(values = c("#ff1493", "#0000ff")) + theme(legend.position="none") + annotate(geom="text", x =1.5, y=11, label="p = 1.48E-06", color="black")
g1894
###1.47704778358419E-06

###Putative A. octoarticulatus virus 6 (g2855_i0) (AOV-6)
which(grepl("g2855", res2$row))
##plot the gene normalized counts
library(ggpubr)
goi <- res2$row[15]
tcounts <- t(log2((counts(dds[goi, ], normalized=TRUE, replaced=FALSE)+.5))) %>%
    merge(colData(dds), ., by="row.names") %>%
    gather(gene, expression, (ncol(.)-length(goi)+1):ncol(.))
g2855<- ggboxplot(tcounts, y = "expression", x = "activity.level", fill = "activity.level") +  labs(x="labour", y="Expression (log normalized counts)", title = "Putative A. octoarticulatus virus 5", fill = "Activity level") + theme(plot.title = element_text(color="black", size=15, hjust = 0.5, face = "italic")) +  scale_fill_manual(values = c("#ff1493", "#0000ff")) + theme(legend.position="none") + annotate(geom="text", x =1.5, y=9, label="p = 8.65E-05", color="black")
g2855
###8.65254052277891E-05

###Putative A. octoarticulatus virus 7 (g3684_i0) (AOV-7)
which(grepl("g3684", res2$row))
##plot the gene normalized counts
library(ggpubr)
goi <- res2$row[24]
tcounts <- t(log2((counts(dds[goi, ], normalized=TRUE, replaced=FALSE)+.5))) %>%
    merge(colData(dds), ., by="row.names") %>%
    gather(gene, expression, (ncol(.)-length(goi)+1):ncol(.))
g3684<- ggboxplot(tcounts, y = "expression", x = "activity.level", fill = "activity.level") +  labs(x="labour", y="Expression (log normalized counts)", title = "Putative A. octoarticulatus virus 6", fill = "Activity level") + theme(plot.title = element_text(color="black", size=15, hjust = 0.5, face = "italic")) +  scale_fill_manual(values = c("#ff1493", "#0000ff")) + theme(legend.position="none") + annotate(geom="text", x =1.5, y=10, label="p = 0.001", color="black")
g3684

##Putative A. octoarticulatus virus 8 (g3814_i0) (AOV-8)
which(grepl("g3814", res2$row))
##plot the gene normalized counts
library(ggpubr)
goi <- res2$row[38]
tcounts <- t(log2((counts(dds[goi, ], normalized=TRUE, replaced=FALSE)+.5))) %>%
    merge(colData(dds), ., by="row.names") %>%
    gather(gene, expression, (ncol(.)-length(goi)+1):ncol(.))
g3814<- ggboxplot(tcounts, y = "expression", x = "activity.level", fill = "activity.level") +  labs(x="labour", y="Expression (log normalized counts)", title = "Putative A. octoarticulatus virus 7", fill = "Activity level") + theme(plot.title = element_text(color="black", size=15, hjust = 0.5, face = "italic")) +  scale_fill_manual(values = c("#ff1493", "#0000ff")) + theme(legend.position="none") + annotate(geom="text", x =1.5, y=10, label="p = 0.004", color="black")
g3814

