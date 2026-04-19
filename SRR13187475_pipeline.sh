#!/bin/bash
# =============================================================================
# SRR13187475 Alignment Pipeline
# Covers: SRA download, hg38 reference, BWA-MEM alignment,
#         unique mapping stats, TSS ±1000bp read counting,
#         MACS3 ChIP-seq broad peak calling, Excel export
# =============================================================================
set -euo pipefail

# ── Paths ─────────────────────────────────────────────────────────────────────
CONDA=~/miniconda3/bin
SRA_ID=SRR13187475
WORKDIR=~/${SRA_ID}
REF_DIR=~/references/hg38
REF=${REF_DIR}/hg38.fa
GTF_GZ=${REF_DIR}/gencode.v47.annotation.gtf.gz
THREADS=8
MACS_OUT=${WORKDIR}/macs3_out
EXCEL_OUT=${WORKDIR}/${SRA_ID}_results.xlsx

mkdir -p "${WORKDIR}" "${REF_DIR}" "${MACS_OUT}"

# ── Step 1: Download SRA ──────────────────────────────────────────────────────
echo "[1/11] Downloading SRA: ${SRA_ID}"
${CONDA}/prefetch ${SRA_ID} --output-directory ~/

# ── Step 2: Convert SRA to FASTQ ─────────────────────────────────────────────
echo "[2/11] Converting SRA to paired FASTQ"
${CONDA}/fasterq-dump ${WORKDIR}/${SRA_ID}.sra \
  --outdir ${WORKDIR} \
  --split-files \
  --threads ${THREADS}

R1=${WORKDIR}/${SRA_ID}_1.fastq
R2=${WORKDIR}/${SRA_ID}_2.fastq

# ── Step 3: Download hg38 reference ──────────────────────────────────────────
echo "[3/11] Downloading hg38 reference genome (UCSC)"
if [ ! -f "${REF}" ]; then
  curl -C - -o ${REF_DIR}/hg38.fa.gz \
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
  gzip -d ${REF_DIR}/hg38.fa.gz
else
  echo "  hg38.fa already exists, skipping download"
fi

# ── Step 4: Build BWA index ───────────────────────────────────────────────────
echo "[4/11] Building BWA index (this takes ~1-2 hours)"
if [ ! -f "${REF}.bwt" ]; then
  ${CONDA}/bwa index ${REF}
else
  echo "  BWA index already exists, skipping"
fi

# ── Step 5: Align with BWA-MEM ────────────────────────────────────────────────
echo "[5/11] Aligning reads with BWA-MEM"
BAM=${WORKDIR}/${SRA_ID}_sorted.bam
RG='@RG\tID:'"${SRA_ID}"'\tSM:'"${SRA_ID}"'\tPL:ILLUMINA\tLB:lib1'

${CONDA}/bwa mem \
  -t ${THREADS} \
  -R "${RG}" \
  ${REF} ${R1} ${R2} \
  2> ${WORKDIR}/bwa_mem.log \
  | ${CONDA}/samtools sort -@ ${THREADS} -o ${BAM} -

${CONDA}/samtools index ${BAM}

# ── Step 6: Alignment statistics ─────────────────────────────────────────────
echo "[6/11] Alignment statistics"
echo "--- samtools flagstat ---"
${CONDA}/samtools flagstat ${BAM}

echo ""
echo "--- Unique mapping counts ---"
MAPQ1=$(${CONDA}/samtools view -c -q 1 ${BAM})
MAPQ20=$(${CONDA}/samtools view -c -q 20 ${BAM})
TOTAL=$(${CONDA}/samtools view -c -F 4 ${BAM})
echo "Total mapped reads:           ${TOTAL}"
echo "Uniquely mapped (MAPQ >= 1):  ${MAPQ1}"
echo "Uniquely mapped (MAPQ >= 20): ${MAPQ20}"

# ── Step 7: Download GENCODE GTF and extract TSS windows ─────────────────────
echo "[7/11] Downloading GENCODE v47 GTF"
if [ ! -f "${GTF_GZ}" ]; then
  curl -C - -o ${GTF_GZ} \
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gtf.gz"
else
  echo "  GTF already exists, skipping download"
fi

# Extract strand-aware TSS positions (1bp) from transcript features
echo "  Extracting TSS positions"
TSS_BED=${REF_DIR}/tss.bed
gzip -dc ${GTF_GZ} \
  | awk '$3 == "transcript"' \
  | awk 'BEGIN{OFS="\t"} {
      chrom = $1; strand = $7;
      if (strand == "+") { tss = $4 - 1 }
      else               { tss = $5 - 1 }
      print chrom, tss, tss+1, ".", ".", strand
    }' \
  | sort -k1,1 -k2,2n \
  > ${TSS_BED}

echo "  TSS sites extracted: $(wc -l < ${TSS_BED})"

# ── Step 8: Extend TSS ±1000bp ───────────────────────────────────────────────
echo "[8/11] Extending TSS windows ±1000bp"

${CONDA}/samtools faidx ${REF}
cut -f1,2 ${REF}.fai > ${REF_DIR}/hg38.chrom.sizes

TSS_WINDOW=${REF_DIR}/tss_pm1000.bed
${CONDA}/bedtools slop \
  -i ${TSS_BED} \
  -g ${REF_DIR}/hg38.chrom.sizes \
  -b 1000 \
  > ${TSS_WINDOW}

# ── Step 9: Count reads per TSS window ───────────────────────────────────────
echo "[9/11] Counting reads per TSS ±1000bp window"
COUNTS=${WORKDIR}/tss_pm1000_counts.bed

${CONDA}/bedtools coverage \
  -a ${TSS_WINDOW} \
  -b ${BAM} \
  -counts \
  > ${COUNTS}

echo ""
echo "--- TSS ±1000bp Read Count Summary ---"
awk '{print $7}' ${COUNTS} | sort -n | awk '
  BEGIN { sum=0; n=0; zero=0; gt100=0 }
  { vals[n++]=$1; sum+=$1 }
  $1==0  { zero++ }
  $1>100 { gt100++ }
  END {
    printf "Total TSS windows:        %d\n", n
    printf "Total reads in windows:   %d\n", sum
    printf "Mean reads/window:        %.1f\n", sum/n
    printf "Median reads/window:      %d\n", vals[int(n/2)]
    printf "Max reads/window:         %d\n", vals[n-1]
    printf "Windows with 0 reads:     %d (%.1f%%)\n", zero, zero/n*100
    printf "Windows with >100 reads:  %d (%.1f%%)\n", gt100, gt100/n*100
  }
'

# ── Step 10: MACS3 ChIP-seq peak calling ─────────────────────────────────────
echo "[10/11] Calling broad ChIP-seq peaks with MACS3"
${CONDA}/macs3 callpeak \
  -t ${BAM} \
  -f BAM \
  -g hs \
  --broad \
  --nolambda \
  --outdir ${MACS_OUT} \
  -n ${SRA_ID} \
  2>&1

PEAKS=${MACS_OUT}/${SRA_ID}_peaks.broadPeak
echo "  Peaks called: $(wc -l < ${PEAKS})"

# ── Step 11: Export results to Excel ─────────────────────────────────────────
echo "[11/11] Exporting results to Excel: ${EXCEL_OUT}"
${CONDA}/python3 - <<PYEOF
import pandas as pd

# Sheet 1: TSS ±1000bp read counts
tss_cols = ["chrom", "start", "end", "name", "score", "strand", "read_count"]
tss = pd.read_csv("${COUNTS}", sep="\t", header=None, names=tss_cols)

# Sheet 2: MACS3 broad peaks
peak_cols = ["chrom", "start", "end", "name", "score", "strand",
             "fold_change", "neg_log10_pvalue", "neg_log10_qvalue"]
peaks = pd.read_csv("${PEAKS}", sep="\t", header=None, names=peak_cols)

with pd.ExcelWriter("${EXCEL_OUT}", engine="openpyxl") as writer:
    tss.to_excel(writer, sheet_name="TSS_pm1000_ReadCounts", index=False)
    peaks.to_excel(writer, sheet_name="MACS3_BroadPeaks", index=False)

print(f"  Sheet 1 (TSS_pm1000_ReadCounts): {len(tss):,} rows")
print(f"  Sheet 2 (MACS3_BroadPeaks):      {len(peaks):,} rows")
print(f"  Saved: ${EXCEL_OUT}")
PYEOF

echo ""
echo "Output files:"
echo "  Sorted BAM:          ${BAM}"
echo "  TSS BED:             ${TSS_BED}"
echo "  TSS ±1000bp BED:     ${TSS_WINDOW}"
echo "  Read counts per TSS: ${COUNTS}"
echo "  MACS3 peaks:         ${PEAKS}"
echo "  Excel report:        ${EXCEL_OUT}"
echo "  BWA log:             ${WORKDIR}/bwa_mem.log"
echo ""
echo "Pipeline complete."
