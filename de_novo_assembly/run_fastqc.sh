#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --time=18:00:00

module load NiaEnv/2019b intelpython3
source activate myPythonEnv
fastqc *.gz -o FastQC_Results
i
