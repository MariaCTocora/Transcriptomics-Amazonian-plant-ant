#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=80
#SBATCH --time=24:00:00


module load NiaEnv/2019b intelpython3
source activate myPythonEnv

blastx -query AllSeqtranscriptsCollapsed.fasta -db uniprot_sprot.pep -num_threads 80 -max_target_seqs 1 -outfmt 5 -evalue 1e-5 -out blastxAllSeqs_uniprot5.xml

