#Generate the list of decoys (genomic sequences)
grep "^>" Wa_genome.fasta | cut -d " " -f 1 > decoys.txt
sed -i.bak -e 's/>//g' decoys.txt

#Concatenate the transcriptome and the genome
cat AllSeqtranscriptsCollapsed.fasta Wa_genome.fasta | gzip >  gentrome_collapsedTranscriptome.fa.gz

#Build Salmon index
salmon index -t gentrome_collapsedTranscriptome.fa.gz -d decoys.txt -p 12 -i salmon_index_collapsedTranscriptomeWithDecoys --gencode
