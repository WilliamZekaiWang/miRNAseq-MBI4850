#!/bin/bash

samples=(
SRR27929302
SRR27929303
SRR27929304
SRR27929305
SRR27929306
	)

# Run Trimmomatic for each sample with nohup
for sample in "${samples[@]}"; do
  
  nohup trimmomatic PE -threads 16 -phred33 \
    "${sample}_1.fastq" "${sample}_2.fastq" \
    "${sample}_1_paired.fastq" "${sample}_1_unpaired.fastq" \
    "${sample}_2_paired.fastq" "${sample}_2_unpaired.fastq" \
    ILLUMINACLIP:adapters.fa:2:30:10:8:true \
    SLIDINGWINDOW:4:20 \
    LEADING:3 \
    TRAILING:3 \
    MINLEN:18 > "${sample}_trimmomatic.log" 2>&1 &
  
done
