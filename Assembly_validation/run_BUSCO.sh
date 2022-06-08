#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=80
#SBATCH --time=24:00:00

module load NiaEnv/2019b intelpython3
source activate myPythonEnv

busco -m transcriptome -i AllSeqtranscriptsCollapsed.fasta -o busco_Allseqs_insecta -l insecta_odb10
busco -m transcriptome -i AllSeqtranscriptsCollapsed.fasta -o busco_Allseqs_eukaryota -l eukaryota_odb10
i
