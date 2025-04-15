# Workflow of MBI4850G Project 

This project aimed to reproduce the work of Ibrahim et al., through processing their provided miRNA data
https://doi.org/10.1093/ijnp/pyae013

The associated write-up and more detailed workflow is detailed in `Report.pdf`

## Downloading data

The data is located on the SRA database with id SRP489148

## Cleaning data

The pipeline I immmitated was the exceRpt smallRNA pipeline

### Removing low quality reads (`trim.sh`)

done with trimmomatic

### Remove contaminants (`rem_uni.sh`/`rem_rRNA.sh`)

- Used NCBI UniVec sequences, and ARF (R package) reference sequences for rRNA
- Aligned data with bowtie2 and took un-aligned sequences

### Align and get counts (`align.sh`/`get_counts.sh`)

- Use STAR to align and build index with hg38
- Use subreads to get counts

## Find Differential Expression

### PCA to remove outliers

### DESeq2 to get significantly differentially expressed miRNA
