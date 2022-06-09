Data analyses for Tocora et al (in preparation)
================
Maria Tocora 
June, 2022

## Data analyses for Tocora et al (in preparation): Transcriptomic analysis of cooperative behavior in a mutualistic ant-plant interaction

### Abstract: 
Some of the most well-studied mutualisms are ant-plant interactions, where partners benefit from each other with nutrition, protection, or dispersal. Even though we have knowledge of the ecology and evolution of these interactions, how mutualistic traits are molecularly regulated remains extensively overlooked.  Here, we focused on investigating the molecular basis of bodyguard behavior in the ant species Allomerus octoarticulatus, which aggressively defends its host myrmecophyte, Cordia nodosa, against herbivores. We performed a de novo transcriptome assembly high-throughput RNA sequences from 14 ant colonies collected from the Peruvian Amazon. A total of 93,122 transcripts and 67,613 unigenes were generated and functionally annotated using BLAST with the UniProtKB dataset as the query. We obtained 136 Gene Ontology (GO) functional sub-groups in the three main GO categories: ‘biological process (86)’, ‘cell component (17)’ and ‘molecular function (31)’. We then analyzed genes that are differentially expressed between high- and low-quality bodyguards and between bodyguards and brood-care workers. We identified 11 differentially expressed genes between the best and worst bodyguards. Most of those genes were upregulated (9/11) in active bodyguards, including a gene with putative functions in cuticle formation, and a protein kinase gene. In the behavioral task comparison, we found 59 differentially expressed genes, most of which were upregulated in brood care workers (58/59), including a putative immunity-related gene and a vitellogenin gene. Our results provide insights into the molecular basis of ant-plant mutualisms and their molecular evolution.

NOTE: Abstract for Evolution (2022) needs to be updated. 

### Description of pipeline, relevant folders in parentheses:
 1.  De novo transcriptome assembly (./De_novo_assembly/): 
 2.  De novo transcriptome assembly validation (./Assembly_validation): 
 3.  De novo transcriptome annotation (./Transcriptome_annotation):
 4.  Samples alignment to a Decoy-aware transcriptome and counts quantification (./Salmon_pipeline): Based on Patro et al (2017)
 5.  Differential gene expression analysis (./DGE_analysis): Based on Love et al (2014) 
 6.  Differential transcript expression (./DTE_analysis): Based on Love et al (2018)  
 7.  Differential transcript usage (./DTU_analysis): Based on Love et al (2018) 
 8.  Differential experssion analysis (./DGE_DTE_DTU)
 9.  Alternative splicing analysis (./DS_analysis): Based on Li et al (2018)
 10.  Phylogenetic analysis of virus (./Virus_phylo): 
 11.  References (./references): References for methods used in the analyses. 
