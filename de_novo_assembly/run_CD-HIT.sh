#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=80
#SBATCH --time=24:00:00

module load NiaEnv/2019b intelpython3
source activate myPythonEnv

cd-hit-est -i transcripts.fasta -o transcriptsCollapsed.fasta -c 0.98 -n 8
i
