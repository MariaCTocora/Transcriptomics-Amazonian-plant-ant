Assembly validation
=================
Maria Tocora 
Junio, 2022

## __BUSCO Validation__

run_BUSCO.sh describes code to assess genomic data quality by using BUSCO (Benchmarking Universal Single-Copy Orthologs) scores, which look for the presence or absence of highly conserved genes in an assembly. More information available at https://busco.ezlab.org/

![alt text](https://github.com/mariatocora/Transcriptomic-analysis-ant-plant/blob/main/Assembly_validation/Figures/BUSCO_results.png)

Figure 1. BUSCO Results

## __rnaQUAST Validation__

run_rnaQUAST.sh describes code to access quast, a tool for evaluating RNA-Seq assemblies (Bushmanova, et al., 2016). Please considered I self-installed quast in order to run the program. For more information visit:  https://cab.spbu.ru/files/rnaquast/release2.2.1/manual.html

![alt text](https://github.com/mariatocora/Transcriptomic-analysis-ant-plant/blob/main/Assembly_validation/Figures/AllseqsTranscript_length.png)

Figure 2. Assembled transcripts length distribution obtained with rnaQUAST

## __References__

1. BUSCO: Mosè Manni, Matthew R Berkeley, Mathieu Seppey, Felipe A Simão, Evgeny M Zdobnov, BUSCO Update: Novel and Streamlined Workflows along with Broader and Deeper Phylogenetic Coverage for Scoring of Eukaryotic, Prokaryotic, and Viral Genomes. Molecular Biology and Evolution, Volume 38, Issue 10, October 2021, Pages 4647–4654
2. BUSCO: Manni, M., Berkeley, M. R., Seppey, M., & Zdobnov, E. M. (2021). BUSCO: Assessing genomic data quality and beyond. Current Protocols, 1, e323. doi: 10.1002/cpz1.323
3. rnaQUAST: Bushmanova, E., Antipov, D., Lapidus, A., Suvorov, V. and Prjibelski, A.D., 2016. rnaQUAST: a quality assessment tool for de novo transcriptome assemblies. Bioinformatics, 32(14), pp.2210-2212.
