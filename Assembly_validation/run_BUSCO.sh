#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=80
#SBATCH --time=24:00:00

module load NiaEnv/2019b intelpython3
source activate myPythonEnv

busco -m transcriptome -i AllSeqtranscriptsCollapsed.fasta -o busco_Allseqs_insectaV2 -l insecta_odb10
i
