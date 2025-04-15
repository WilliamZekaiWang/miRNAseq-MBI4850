#!/bin/bash

samples=(
	SRR27929302
SRR27929303
SRR27929304
SRR27929305
SRR27929306
SRR27929307
SRR27929308
SRR27929309
SRR27929310
SRR27929311
	)

for sample in "${samples[@]}"; do
  
  nohup bowtie2 -x ../uni_ind/UniVecInd \
    -1 "${sample}_1_paired.fastq" -2 "${sample}_2_paired.fastq" \
    --un-conc-gz "${sample}_noUni_%.fastq.gz" \
    --al-conc-gz "${sample}_UniVec_%.fastq.gz" \
    --threads 2 \
    -S "${sample}_univec_aligned.sam" > "${sample}_univec_filtering.log" 2>&1 &
  
done

