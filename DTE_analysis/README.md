Differential Transcript Expression 
================
Maria Tocora
June, 2022

The following analyses are based the Swish method for differential expression analysis of bulk or single-cell RNA-seq data using inferential replicate counts is described in the following reference: Zhu et al. (2019) doi: 10.1093/nar/gkz622, and available here https://bioconductor.org/packages/release/bioc/vignettes/fishpond/inst/doc/swish.html#Differential_transcript_expression

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

###Load packages
library(tximeta)
library(fishpond)
```


## __References__
1. Tximeta: Love MI, Soneson C, Hickey PF, Johnson LK, Pierce NT, Shepherd L, Morgan M, Patro R (2020). “Tximeta: Reference sequence checksums for provenance identification in RNA-seq.” PLOS Computational Biology, 16, e1007664. doi: 10.1371/journal.pcbi.1007664.
2. fishpond: Zhu A, Srivastava A, Ibrahim JG, Patro R, Love MI (2019). “Nonparametric expression analysis using inferential replicate counts.” Nucleic Acids Research, 47, e105. doi: 10.1093/nar/gkz622, https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6765120.
