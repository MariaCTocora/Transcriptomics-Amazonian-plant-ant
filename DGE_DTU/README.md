Differential Expression Analysis (DGE & DTU)
==================================================
Maria Tocora 

The following code is modified from https://jokergoo.github.io/ComplexHeatmap-reference/book/upset-plot.html#upset-plots-as-heatmaps. Tables with the number of genes per every differential expression analysis category are provided below under every comparison. The UpSet plots group co-occurring variables into sets and shows you a bar chart of their frequency (Lex et al., 2014; Conway, et al., 2017). The databases are available in the "Data" folder. 

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

### UpSet Plot
```
mydata = read.csv("Final_geneSummary.csv")
Analysis = colnames(mydata)[2:4]

upset(
    mydata,
    Analysis,
    base_annotations=list(
        'Intersection size'=intersection_size(
            counts=FALSE,
            mapping=aes(fill=A)
        )
    ),
    width_ratio=0.1
)
```

## __References__

1. UpSetR: Conway, J. R., Lex, A., & Gehlenborg, N. (2017). UpSetR: an R package for the visualization of intersecting sets and their properties. Bioinformatics (Oxford, England), 33(18), 2938–2940. https://doi.org/10.1093/bioinformatics/btx364
2. UpSetR: Lex, A., Gehlenborg, N., Strobelt, H., Vuillemot, R., & Pfister, H. (2014). UpSet: Visualization of Intersecting Sets. IEEE transactions on visualization and computer graphics, 20(12), 1983–1992. https://doi.org/10.1109/TVCG.2014.2346248
3. ComplexHeatmap: Gu Z, Eils R, Schlesner M (2016). “Complex heatmaps reveal patterns and correlations in multidimensional genomic data.” Bioinformatics.

