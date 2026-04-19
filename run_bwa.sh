#!/bin/bash
set -e

REF=~/references/hg38/hg38.fa
R1=~/SRR13187475/SRR13187475_1.fastq
R2=~/SRR13187475/SRR13187475_2.fastq
OUT=~/SRR13187475/SRR13187475_sorted.bam
LOG=~/SRR13187475/bwa_mem.log
RG='@RG\tID:SRR13187475\tSM:SRR13187475\tPL:ILLUMINA\tLB:lib1'

echo "Starting BWA-MEM alignment..."
~/miniconda3/bin/bwa mem -t 8 -R "$RG" "$REF" "$R1" "$R2" 2>"$LOG" \
  | ~/miniconda3/bin/samtools sort -@ 8 -o "$OUT" -

echo "Indexing BAM..."
~/miniconda3/bin/samtools index "$OUT"

echo "Alignment complete: $OUT"
