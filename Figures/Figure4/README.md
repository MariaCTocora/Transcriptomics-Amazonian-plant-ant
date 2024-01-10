


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
