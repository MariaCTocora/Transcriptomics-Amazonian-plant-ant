Differential Gene Expression (DGE) Analysis 
================
Maria Tocora
June, 2022

```{r Install and Load Packages}


library("DESeq2")
library(tidyverse)


countdata <- read.csv("activity_count_matrix.csv", header = TRUE, row.names = 1)
coldata <- read.delim("activity_samples.txt", header = TRUE, row.names = 1)
all(rownames(coldata) == colnames(countdata))
### TRUE 
