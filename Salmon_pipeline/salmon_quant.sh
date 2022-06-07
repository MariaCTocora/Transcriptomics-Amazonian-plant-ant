#Read quantification with Salmon
#We used Gibbs sampling over the fragment equivalence classes to estimate the variance in abundance estimates: --numGibbsSamples 20 --thinningFactor 100
#We passed the --gcBias flag so Salmon will learn and correct for fragment-level GC biases in the input data

while read i
do

echo "Processing sample $i ..."i

gunzip $i\_R1_output_forward_paired.fq.gz

gunzip $i\_R2_output_reverse_paired.fq.gz

/ohta1/haoran.xue/programs/salmon-1.8.0_linux_x86_64/bin/salmon quant --numGibbsSamples 20 --thinningFactor 100  --gcBias -i ../salmon_index_collapsedTranscriptomeWithoutDecoys -l A -1 $i\_R1_output_forward_paired.fq -2 $i\_R2_output_reverse_paired.fq --validateMappings -o $i.collapsedTranscriptomeWithDecoys.transcripts_quant

gzip $i\_R1_output_forward_paired.fq &

gzip $i\_R2_output_reverse_paired.fq &


done < sample_list
