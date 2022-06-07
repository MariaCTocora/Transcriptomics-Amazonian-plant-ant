#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=80
#SBATCH --time=24:00:00

spades --threads 80 --rna -1 R1_seq_cat.fastq -2 R2_seq_cat.fastq -o Spades_output
i

