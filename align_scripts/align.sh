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

# Define paths
GENOME_DIR="../genome/gen"
THREADS=16

for sample in "${samples[@]}"; do
  
  nohup STAR --runThreadN ${THREADS} \
    --genomeDir ${GENOME_DIR} \
    --readFilesIn "${sample}_clean_1.fastq.gz" "${sample}_clean_2.fastq.gz" \
    --readFilesCommand gunzip -c \
    --outFileNamePrefix "${sample}_" \
    --outSAMtype BAM SortedByCoordinate \
    --outBAMcompression 10 \
    --limitBAMsortRAM 20000000000 > "${sample}_STAR.log" 2>&1 &
  
done
