
Differential Expression Analysis (DGE & DTE & DTU)
==================================================
Maria Tocora 
Junio, 2022

The following code is modified from https://jokergoo.github.io/ComplexHeatmap-reference/book/upset-plot.html#upset-plots-as-heatmaps. Tables with the number of genes per every differential expression analysis category are provided below under every comparison. The UpSet plots group co-occurring variables into sets and shows you a bar chart of their frequency (Lex et al., 2014; Conway, et al., 2017). 

## Install and load packages 
```{r Install and Load Packages}
### UpSetR (https://github.com/hms-dbmi/UpSetR)
install.packages("UpSetR")

### ComplexHeatmap (https://www.bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html)
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")

library(UpSetR)
library(ComplexHeatmap)
```

## __Activity differential expression analysis__

| Diferential Expression Analysis |	No. Genes |
| --- | --- | 
| DGE	| 34 |
| DTU |	10 |
| DTE |	0 |
| DGE & DTE | 4 |
| DTE & DTU |	0 |
| DGE & DTU |	2 |
| DGE & DTE & DTU |	0 |
| Total No. genes | 51 |

### UpSet Plot
```
CasteDGE_DTE_DTU <- read_csv("C:/Users/Paula/Desktop/BODYGUARD PROJECT/Transcriptome analysis/DESeq2_analysis/AllTranscripts/Caste_DGE_DTE_DTU_analysis.csv")
head(CasteDGE_DTE_DTU)
m = make_comb_mat(CasteDGE_DTE_DTU, top_n_sets = 3)
UpSet(m)
comb_elements = lapply(comb_name(m), function(nm) extract_comb(m, nm))
Interactions = lapply(comb_elements, function(ind) CasteDGE_DTE_DTU$Gene[ind])
CasteUpSetPlot <- UpSet(t(m)) + rowAnnotation(Interactions = anno_boxplot(Interactions))
CasteUpSetPlot
```

![alt text](https://github.com/mariatocora/Transcriptomic-analysis-ant-plant/blob/main/DGE_DTE_DTU/Figures/Activity_UpSetPlot.png)

Figure 1. Significant genes in the activity comparison per differential expression analysis. 

## __Caste differential expression analysis__

| Diferential Expression Analysis |	No. Genes |
| --- | --- | 
| DGE | 358 |
| DTU |	74 |
| DTE |	400 |
| DGE & DTE |	251 |
| DTE & DTU |	30 |
| DGE & DTU |	11 |
| DGE & DTE & DTU |	9 |
| Total No. genes | 553 |

### UpSet Plot

```
CasteDGE_DTE_DTU <- read_csv("C:/Users/Paula/Desktop/BODYGUARD PROJECT/Transcriptome analysis/DESeq2_analysis/AllTranscripts/Caste_DGE_DTE_DTU_analysis.csv")
head(CasteDGE_DTE_DTU)
m = make_comb_mat(CasteDGE_DTE_DTU, top_n_sets = 3)
UpSet(m)
comb_elements = lapply(comb_name(m), function(nm) extract_comb(m, nm))
Interactions = lapply(comb_elements, function(ind) CasteDGE_DTE_DTU$Gene[ind])
CasteUpSetPlot <- UpSet(t(m)) + rowAnnotation(Interactions = anno_boxplot(Interactions))
CasteUpSetPlot
```

![alt text](https://github.com/mariatocora/Transcriptomic-analysis-ant-plant/blob/main/DGE_DTE_DTU/Figures/Caste_UpSetPlot.png)

Figure 2. Significant genes in the caste comparison per differential expression analysis. 

## __References__

1. UpSetR: Conway, J. R., Lex, A., & Gehlenborg, N. (2017). UpSetR: an R package for the visualization of intersecting sets and their properties. Bioinformatics (Oxford, England), 33(18), 2938–2940. https://doi.org/10.1093/bioinformatics/btx364
2. UpSetR: Lex, A., Gehlenborg, N., Strobelt, H., Vuillemot, R., & Pfister, H. (2014). UpSet: Visualization of Intersecting Sets. IEEE transactions on visualization and computer graphics, 20(12), 1983–1992. https://doi.org/10.1109/TVCG.2014.2346248
3. ComplexHeatmap: Gu Z, Eils R, Schlesner M (2016). “Complex heatmaps reveal patterns and correlations in multidimensional genomic data.” Bioinformatics.
