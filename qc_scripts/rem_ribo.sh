#!/bin/bash

samples=(SRR27929312
SRR27929313
SRR27929314
SRR27929315
SRR27929316
SRR27929317
SRR27929318
SRR27929319
SRR27929320
SRR27929321
	)

for sample in "${samples[@]}"; do
  
  nohup bowtie2 -x ../rna_ind/rRNA_index \
    -1 "${sample}_noUni_1.fastq.gz" -2 "${sample}_noUni_2.fastq.gz" \
    --un-conc-gz "${sample}_clean_%.fastq.gz" \
    --al-conc-gz "${sample}_rRNA_%.fastq.gz" \
    --threads 12 \
    -S "${sample}_rRNA.sam" > "${sample}_rRNA_filtering.log" 2>&1 &
  
done
