Salmon pipeline for samples alignment and counts quatification
================
Haoran Xue
June, 2022

salmon_indexing.sh was used to build a salmon index with a decoy-aware transcriptome file. The decoy we used was the genome of Wasmannia auropuctata (available here https://www.ncbi.nlm.nih.gov/genome/?term=txid64793[orgn]), a species closely related to Allomerus octoarticulatus, as the genome of A. octoarticulatus was not available.

salmon_quant.sh was used to quantify the reads against the index for each sample.

## __References__
1. Patro, R., Duggal, G., Love, M. I., Irizarry, R. A., & Kingsford, C. (2017). Salmon provides fast and bias-aware quantification of transcript expression. Nature Methods.
