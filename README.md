Data analyses for Tocora et al (in preparation)
================
Maria Tocora 

## Data analyses for Tocora et al:  In sickness and in health: viral infections correlate with the quality of ant bodyguarding behavior in an Amazonian ant-plant symbiosis

### Abstract: 
Ant-plant symbioses are classic examples of mutualism in which ant ‘bodyguards’ defend myrmecophytic plants against enemies in exchange for nest sites and often food. We used RNA-Seq to profile the transcriptomes of Allomerus octoarticulatus ant workers, which aggressively defend the Amazonian plant Cordia nodosa against herbivores, but to varying degrees. Field behavioral assays with herbivores in the Peruvian Amazon showed striking variation among colonies in the relative zeal with which A. octoarticulatus workers defend their host plant. Highly effective and ineffective bodyguards differed in their gene expression profiles, which revealed viral infections significantly associated with ant bodyguarding behavior. Transcripts from eight new positive-sense single-stranded RNA viruses were differentially expressed between colonies with high- or low-quality bodyguards. Colonies of ‘good’ and ‘bad’ bodyguards were infected by distinct viruses, including viruses clustering phylogenetically with viruses known to cause aggression or reduced locomotion, respectively, in bees. Gene expression, including of immunity-related genes, also differed between broodcare workers and bodyguard ants, suggesting bodyguarding is a distinct worker task. Ant colony health and viral infections may influence ant cooperation with plants in ant-plant mutualisms.


### Description of pipeline, relevant folders in parentheses:
 1.  De novo transcriptome assembly (./De_novo_assembly/) 
 2.  De novo transcriptome assembly validation (./Assembly_validation) 
 3.  De novo transcriptome annotation (./Transcriptome_annotation)
 4.  Samples alignment to a Decoy-aware transcriptome and counts quantification (./Salmon_pipeline). Based on Patro et al (2017)
 5.  Differential gene expression analysis (./DGE_analysis). Based on Love et al (2014) 
 6.  Differential transcript usage (./DTU_analysis). Based on Love et al (2018) 
 7.  Differential experssion analysis comparison (./DGE_DTU)
 8. Phylogenetic analysis of virus (./Virus_phylo)
 9. Field data analysis (.FieldData)
 10. References (./references): References for methods used in the analyses. 
