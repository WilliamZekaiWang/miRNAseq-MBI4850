#!/bin/bash

# Define the list of sample IDs (without extensions)
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
	SRR27929312
	SRR27929313
	SRR27929314
	SRR27929315
	SRR27929316
	SRR27929317
	SRR27929318
	SRR27929319
	SRR27929320
	SRR27929321
	SRR27929322
	SRR27929323
	SRR27929324
	SRR27929325
	SRR27929326
	SRR27929327
	SRR27929328
	SRR27929329
	SRR27929330
	SRR27929331
	SRR27929332
	SRR27929333
	SRR27929334
	SRR27929335
	SRR27929336
	SRR27929337
	SRR27929338
	SRR27929339
	SRR27929340
	SRR27929341

	)

# Define paths
GTF="hsa.gff3"  # Path to miRNA annotation (GFF/GTF)
THREADS=16

# Run featureCounts for each sample
for sample in "${samples[@]}"; do
  echo "Processing ${sample}..."
  
  nohup featureCounts \
    -T ${THREADS} \
    -a ${GTF} \
    -F gff \
    -o "${sample}_miRNA_counts.txt" \
    -t miRNA \
    -g "Name" \
    -M \
    --fraction \
    --primary \
    -O \
    -s 0 \
    "${sample}_Aligned.sortedByCoord.out.bam" > "${sample}_featureCounts.log" 2>&1 &
  
  echo "Started ${sample} in background (PID $!). Log file: ${sample}_featureCounts.log"
done

echo "All samples submitted for miRNA quantification!"

# After all jobs complete, create a combined counts matrix
echo "To create a combined matrix after all jobs finish, run:"
echo "awk 'BEGIN {OFS=\"\t\"} NR==1 || FNR>1' *_miRNA_counts.txt > combined_miRNA_counts.txt"
