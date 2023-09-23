Data analyses for Tocora et al (in preparation)
================
Maria Tocora 
Sept, 2023

## Data analyses for Tocora et al (in preparation): Transcriptomic analysis of cooperative behavior in a mutualistic ant-plant interaction

### Abstract: 
Ant-plant symbioses are a classic example of mutualism in which ant ‘bodyguards’ defend myrmecophytic plants against enemies in exchange for nest sites and often food. Although the ecology of these interactions has been well studied, how this mutualism is molecularly regulated has been largely overlooked. Here, we investigated the molecular basis of bodyguard behavior in the ant species Allomerus octoarticulatus, which aggressively defends its host myrmecophyte, Cordia nodosa, against herbivores. Field observations in the Peruvian Amazon show variation among colonies in their symbiotic effectiveness. RNA-seq data unravel inter-colony variation in plant protection and mutualist quality is linked to gene and isoform expression. Transcriptomic analysis suggests signatures of task specialization between brood care (nurse) and bodyguard workers providing evidence of labor-dependent molecular expression. As well, our results suggest a condition-specific molecular landscape underlying this interaction as less cooperative partners also differ in immune-related genetic expression from ants actively engaging in cooperation. Finally, we find potential viral infections may explain some of the variations in ants’ willingness to reciprocate the mutualistic interaction. Altogether, our study provides insights into the molecular basis of defense mutualisms, enabling the identification of mechanisms driving cooperative investment. 

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
 10. Phylogenetic analysis of virus (./Virus_phylo): 
 11. References (./references): References for methods used in the analyses. 
