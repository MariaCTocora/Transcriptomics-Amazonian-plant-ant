De novo transcriptome assembly 
===========================
Maria Tocora 
June, 2022

### Assess raw read quality with FastQC
run_fastqc.sh describes code to assess the intrinsic quality of the raw sequences (Andrews, S., 2010). More information available at https://dnacore.missouri.edu/PDF/FastQC_Manual.pdf

### Quality control with Trimmomatic
run_trimmomatic.sh describes code to access a trimming tool for Illumina NGS data. Trimmomatic allows users to remove adapters and low-quality sequences that might affect the quality of our reads (Bolger et al., 2014). For more information visit the Trimmomatic Manual (http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf).

### De novo transcriptome assembly with rnaSPAdes
run_spades.sh describes code to run rnaSPAdes – a de novo transcriptome assembler from RNA-Seq data (Bushmanova et al., 2019). Manual available at https://cab.spbu.ru/files/release3.12.0/rnaspades_manual.html. 

### Reduce transcript redundancy with CD-HIT
run_CD-HIT.sh describes code to run the clustering program CD-HIT. CD-HIT clusters highly similar sequences by sorting and grouping them into clusters defined by a given threshold (Huang, et al, 2010; Li & Godzik, 2006). I used a  sequence identity threshold equal to 0.98. More information at http://www.bioinformatics.org/cd-hit/cd-hit-user-guide.pdf.

## __References__
1. FastQC: Andrews, S. (2010). FastQC: A Quality Control Tool for High Throughput Sequence Data [Online]. Available
online at: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
2. Trimmomatic:  Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: A flexible trimmer for Illumina Sequence Data.
Bioinformatics, btu170
3. rnaSPAdes: Bushmanova, E., Antipov, D., Lapidus, A., & Prjibelski, A. D. (2019). rnaSPAdes: a de novo transcriptome assembler and its application to RNA-Seq data. GigaScience, 8(9), giz100. https://doi.org/10.1093/gigascience/giz100.
4. CD-HIT: Huang, Y., Niu, B., Gao, Y., Fu, L., & Li, W. (2010). CD-HIT Suite: a web server for clustering and comparing biological sequences. Bioinformatics (Oxford, England), 26(5), 680–682. https://doi.org/10.1093/bioinformatics/btq003.
5. CD-HIT: Li, W., & Godzik, A. (2006). Cd-hit: a fast program for clustering and comparing large sets of protein or nucleotide sequences. Bioinformatics (Oxford, England), 22(13), 1658–1659. https://doi.org/10.1093/bioinformatics/btl158

