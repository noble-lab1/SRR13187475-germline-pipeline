#!/bin/bash
# =============================================================================
# GATK Germline Variant Calling Pipeline — SRR13187475
# Steps:
#   1. Filter unique concordantly mapped pairs (MAPQ >= 20, proper pairs)
#   2. Mark duplicates
#   3. Download BQSR + VQSR known-sites VCFs (if not present)
#   4. Base Quality Score Recalibration (BQSR)
#   5. HaplotypeCaller (GVCF mode)
#   6. GenotypeGVCFs
#   7. Variant Quality Score Recalibration (VQSR) — SNPs then Indels
#   8. Hard filter fallback (applied to variants that fail VQSR or if VQSR skipped)
# =============================================================================
set -euo pipefail

# ── Paths ─────────────────────────────────────────────────────────────────────
CONDA=~/miniconda3/bin
SRA_ID=SRR13187475
WORKDIR=~/${SRA_ID}
REF_DIR=~/references/hg38
REF=${REF_DIR}/hg38.fa
THREADS=8

BAM_IN=${WORKDIR}/${SRA_ID}_sorted.bam

# BQSR known sites
DBSNP=${REF_DIR}/Homo_sapiens_assembly38.dbsnp138.vcf.gz
MILLS=${REF_DIR}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz

# VQSR training resources (SNPs)
HAPMAP=${REF_DIR}/hapmap_3.3.hg38.vcf.gz
OMNI=${REF_DIR}/1000G_omni2.5.hg38.vcf.gz
G1K_SNP=${REF_DIR}/1000G_phase1.snps.high_confidence.hg38.vcf.gz

# Output files
BAM_FILTERED=${WORKDIR}/${SRA_ID}_filtered.bam
BAM_MARKDUP=${WORKDIR}/${SRA_ID}_markdup.bam
MARKDUP_METRICS=${WORKDIR}/${SRA_ID}_markdup_metrics.txt
RECAL_TABLE=${WORKDIR}/${SRA_ID}_recal.table
BAM_BQSR=${WORKDIR}/${SRA_ID}_bqsr.bam
GVCF=${WORKDIR}/${SRA_ID}.g.vcf.gz
VCF_RAW=${WORKDIR}/${SRA_ID}_genotyped.vcf.gz

# VQSR intermediate files
SNP_RAW=${WORKDIR}/${SRA_ID}_snps_raw.vcf.gz
INDEL_RAW=${WORKDIR}/${SRA_ID}_indels_raw.vcf.gz
SNP_RECAL=${WORKDIR}/${SRA_ID}_snps.recal
SNP_TRANCHES=${WORKDIR}/${SRA_ID}_snps.tranches
SNP_RSCRIPT=${WORKDIR}/${SRA_ID}_snps_plots.R
INDEL_RECAL=${WORKDIR}/${SRA_ID}_indels.recal
INDEL_TRANCHES=${WORKDIR}/${SRA_ID}_indels.tranches
INDEL_RSCRIPT=${WORKDIR}/${SRA_ID}_indels_plots.R

# Final output VCFs
VCF_SNP_VQSR=${WORKDIR}/${SRA_ID}_snps_vqsr.vcf.gz
VCF_INDEL_VQSR=${WORKDIR}/${SRA_ID}_indels_vqsr.vcf.gz

export PATH=${CONDA}:$PATH

# ── Step 1: Filter unique concordantly mapped pairs ───────────────────────────
echo "[1/8] Filtering: unique (MAPQ >= 20) + concordantly mapped pairs"
# -f 0x2   = properly paired (concordant)
# -q 20    = MAPQ >= 20 (unique mapping)
# -F 0x4   = not unmapped
# -F 0x100 = not secondary alignment
# -F 0x800 = not supplementary alignment
${CONDA}/samtools view -@ ${THREADS} \
  -f 0x2 \
  -F 0x4 -F 0x100 -F 0x800 \
  -q 20 \
  -b ${BAM_IN} \
  | ${CONDA}/samtools sort -@ ${THREADS} -o ${BAM_FILTERED} -

${CONDA}/samtools index ${BAM_FILTERED}

TOTAL_IN=$(${CONDA}/samtools view -c ${BAM_IN})
TOTAL_FILT=$(${CONDA}/samtools view -c ${BAM_FILTERED})
echo "  Reads in:      ${TOTAL_IN}"
echo "  Reads kept:    ${TOTAL_FILT}"

# ── Step 2: Mark Duplicates ───────────────────────────────────────────────────
echo "[2/8] Marking duplicates"
${CONDA}/gatk MarkDuplicates \
  -I ${BAM_FILTERED} \
  -O ${BAM_MARKDUP} \
  -M ${MARKDUP_METRICS} \
  --VALIDATION_STRINGENCY SILENT \
  --CREATE_INDEX true

echo "  Duplicate metrics: ${MARKDUP_METRICS}"

# ── Step 3: Create sequence dictionary if missing ────────────────────────────
if [ ! -f "${REF_DIR}/hg38.dict" ]; then
  echo "[3/8] Creating sequence dictionary for hg38..."
  ${CONDA}/gatk CreateSequenceDictionary -R ${REF}
else
  echo "[3/8] Sequence dictionary already present"
fi

# ── Step 4: Base Quality Score Recalibration (BQSR) ──────────────────────────
# BQSR requires known-variant VCFs from the GATK resource bundle (Google Cloud
# Storage — requires a Google account). If they are not present, BQSR is
# skipped and the markdup BAM is used directly, which is an acceptable
# alternative for this exercise.
#
# To download manually:
#   gsutil -u YOUR_PROJECT cp gs://genomics-public-data/resources/broad/hg38/v0/\
#     Homo_sapiens_assembly38.dbsnp138.vcf.gz ~/references/hg38/
#   gsutil -u YOUR_PROJECT cp gs://genomics-public-data/resources/broad/hg38/v0/\
#     Mills_and_1000G_gold_standard.indels.hg38.vcf.gz ~/references/hg38/
#   (also download .tbi index files)

if [ -f "${DBSNP}" ] && [ -f "${MILLS}" ]; then
  echo "[4/8] BQSR — BaseRecalibrator"
  ${CONDA}/gatk BaseRecalibrator \
    -I ${BAM_MARKDUP} \
    -R ${REF} \
    --known-sites ${DBSNP} \
    --known-sites ${MILLS} \
    -O ${RECAL_TABLE}

  echo "[4/8] BQSR — ApplyBQSR"
  ${CONDA}/gatk ApplyBQSR \
    -I ${BAM_MARKDUP} \
    -R ${REF} \
    --bqsr-recal-file ${RECAL_TABLE} \
    -O ${BAM_BQSR}
else
  echo "[4/8] BQSR resource VCFs not found — skipping BQSR"
  echo "  NOTE: Using markdup BAM directly. Base quality scores will not be recalibrated."
  echo "  To enable BQSR, download the GATK resource bundle (requires Google account)."
  echo "  See comments in this script for download instructions."
  BAM_BQSR=${BAM_MARKDUP}
fi

# ── Step 5: HaplotypeCaller (GVCF mode) ──────────────────────────────────────
echo "[5/8] HaplotypeCaller (GVCF mode)"
${CONDA}/gatk HaplotypeCaller \
  -R ${REF} \
  -I ${BAM_BQSR} \
  -O ${GVCF} \
  -ERC GVCF \
  --native-pair-hmm-threads ${THREADS}

# ── Step 6: GenotypeGVCFs ─────────────────────────────────────────────────────
echo "[6/8] GenotypeGVCFs"
${CONDA}/gatk GenotypeGVCFs \
  -R ${REF} \
  -V ${GVCF} \
  -O ${VCF_RAW}

# Split raw VCF into SNPs and Indels for VQSR
${CONDA}/gatk SelectVariants -V ${VCF_RAW} --select-type-to-include SNP   -O ${SNP_RAW}
${CONDA}/gatk SelectVariants -V ${VCF_RAW} --select-type-to-include INDEL -O ${INDEL_RAW}

# ── Step 7: Variant Quality Score Recalibration (VQSR) ───────────────────────
if [ ! -f "${HAPMAP}" ] || [ ! -f "${OMNI}" ] || [ ! -f "${G1K_SNP}" ]; then
  echo "[7/8] VQSR training resources not found — skipping VQSR, using hard filters only"
  cp ${SNP_RAW}   ${VCF_SNP_VQSR}
  cp ${SNP_RAW}.tbi ${VCF_SNP_VQSR}.tbi 2>/dev/null || ${CONDA}/gatk IndexFeatureFile -I ${VCF_SNP_VQSR}
  cp ${INDEL_RAW} ${VCF_INDEL_VQSR}
  cp ${INDEL_RAW}.tbi ${VCF_INDEL_VQSR}.tbi 2>/dev/null || ${CONDA}/gatk IndexFeatureFile -I ${VCF_INDEL_VQSR}
else

echo "[7/8] VQSR — SNPs"
# VariantRecalibrator builds a Gaussian mixture model from known truth/training
# sites to distinguish true variants from sequencing artifacts.
# Resource annotations and their roles:
#   hapmap   — highest confidence truth set (known=false, training=true, truth=true, prior=15)
#   omni     — high-confidence SNP array calls (known=false, training=true, truth=true, prior=12)
#   1000G    — phase 1 SNPs (known=false, training=true, truth=false, prior=10)
#   dbsnp    — known variants but not used for training (known=true, training=false, truth=false, prior=7)
${CONDA}/gatk VariantRecalibrator \
  -V ${SNP_RAW} \
  --trust-all-polymorphic \
  -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 \
  -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 \
  -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 \
  -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP \
  -mode SNP \
  --resource:hapmap,known=false,training=true,truth=true,prior=15   ${HAPMAP}  \
  --resource:omni,known=false,training=true,truth=true,prior=12     ${OMNI}    \
  --resource:1000G,known=false,training=true,truth=false,prior=10   ${G1K_SNP} \
  --resource:dbsnp,known=true,training=false,truth=false,prior=7    ${DBSNP}   \
  -O ${SNP_RECAL} \
  --tranches-file ${SNP_TRANCHES} \
  --rscript-file ${SNP_RSCRIPT}

echo "[7/8] VQSR — ApplyVQSR (SNPs, tranche 99.5%)"
${CONDA}/gatk ApplyVQSR \
  -V ${SNP_RAW} \
  --recal-file ${SNP_RECAL} \
  --tranches-file ${SNP_TRANCHES} \
  --truth-sensitivity-filter-level 99.5 \
  --create-output-variant-index true \
  -mode SNP \
  -O ${VCF_SNP_VQSR}

echo "[7/8] VQSR — VariantRecalibrator (Indels)"
# Indels use Mills + dbSNP only (no SNP array resources)
${CONDA}/gatk VariantRecalibrator \
  -V ${INDEL_RAW} \
  --trust-all-polymorphic \
  -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 \
  -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 \
  -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 \
  -tranche 91.0 -tranche 90.0 \
  -an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP \
  -mode INDEL \
  --max-gaussians 4 \
  --resource:mills,known=false,training=true,truth=true,prior=12  ${MILLS} \
  --resource:dbsnp,known=true,training=false,truth=false,prior=2  ${DBSNP} \
  -O ${INDEL_RECAL} \
  --tranches-file ${INDEL_TRANCHES} \
  --rscript-file ${INDEL_RSCRIPT}

echo "[7/8] VQSR — ApplyVQSR (Indels, tranche 99.0%)"
${CONDA}/gatk ApplyVQSR \
  -V ${INDEL_RAW} \
  --recal-file ${INDEL_RECAL} \
  --tranches-file ${INDEL_TRANCHES} \
  --truth-sensitivity-filter-level 99.0 \
  --create-output-variant-index true \
  -mode INDEL \
  -O ${VCF_INDEL_VQSR}

fi  # end VQSR block

# ── Step 8: Hard filter fallback (variants outside VQSR model) ────────────────
# VQSR labels variants PASS or assigns a tranche filter. Variants that could not
# be modeled (e.g. in low-complexity regions) remain with a VQSR filter tag.
# This step additionally applies hard filters as a safety net.
echo "[8/8] Hard filter fallback on VQSR output"

VCF_SNP_FINAL=${WORKDIR}/${SRA_ID}_snps_final.vcf.gz
VCF_INDEL_FINAL=${WORKDIR}/${SRA_ID}_indels_final.vcf.gz

# Allele Balance (AB): for het calls, alt allele fraction should be 0.2–0.8.
# Uses JEXL to access the AD (allele depth) FORMAT field. Only applied to hets
# — homozygous alt calls are expected to have near-100% alt reads.
# MQ0: number of reads with mapping quality 0 at the site. High MQ0 indicates
# a repetitive/low-complexity region prone to mismapping artifacts.
# MQ0 filter applied as a fraction of total depth (MQ0 / DP > 0.1).

${CONDA}/gatk VariantFiltration \
  -V ${VCF_SNP_VQSR} \
  --filter-expression "QD < 2.0"              --filter-name "QD2"             \
  --filter-expression "FS > 60.0"             --filter-name "FS60"            \
  --filter-expression "MQ < 40.0"             --filter-name "MQ40"            \
  --filter-expression "MQRankSum < -12.5"     --filter-name "MQRankSum-12.5"  \
  --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
  --filter-expression "MQ0 / DP > 0.1"        --filter-name "MQ0Frac"         \
  --filter-expression "vc.isHet() && (vc.getGenotype(0).getAD()[1] * 1.0 / (vc.getGenotype(0).getAD()[0] + vc.getGenotype(0).getAD()[1]) < 0.2 || vc.getGenotype(0).getAD()[1] * 1.0 / (vc.getGenotype(0).getAD()[0] + vc.getGenotype(0).getAD()[1]) > 0.8)" \
    --filter-name "ABfilter" \
  -O ${VCF_SNP_FINAL}

${CONDA}/gatk VariantFiltration \
  -V ${VCF_INDEL_VQSR} \
  --filter-expression "QD < 2.0"               --filter-name "QD2"              \
  --filter-expression "FS > 200.0"             --filter-name "FS200"            \
  --filter-expression "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
  --filter-expression "MQ0 / DP > 0.1"         --filter-name "MQ0Frac"          \
  --filter-expression "vc.isHet() && (vc.getGenotype(0).getAD()[1] * 1.0 / (vc.getGenotype(0).getAD()[0] + vc.getGenotype(0).getAD()[1]) < 0.2 || vc.getGenotype(0).getAD()[1] * 1.0 / (vc.getGenotype(0).getAD()[0] + vc.getGenotype(0).getAD()[1]) > 0.8)" \
    --filter-name "ABfilter" \
  -O ${VCF_INDEL_FINAL}

# ── Summary ───────────────────────────────────────────────────────────────────
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "Pipeline complete."
echo ""
echo "Output files:"
echo "  Filtered BAM (MAPQ>=20, proper pairs): ${BAM_FILTERED}"
echo "  Markdup BAM:                           ${BAM_MARKDUP}"
echo "  Markdup metrics:                       ${MARKDUP_METRICS}"
echo "  BQSR recal table:                      ${RECAL_TABLE}"
echo "  BQSR BAM:                              ${BAM_BQSR}"
echo "  GVCF:                                  ${GVCF}"
echo "  Genotyped VCF (raw):                   ${VCF_RAW}"
echo "  SNPs after VQSR:                       ${VCF_SNP_VQSR}"
echo "  Indels after VQSR:                     ${VCF_INDEL_VQSR}"
echo "  SNPs final (VQSR + hard filter):       ${VCF_SNP_FINAL}"
echo "  Indels final (VQSR + hard filter):     ${VCF_INDEL_FINAL}"
echo "  SNP VQSR R plots:                      ${SNP_RSCRIPT}"
echo "  Indel VQSR R plots:                    ${INDEL_RSCRIPT}"
echo ""

SNP_PASS=$(${CONDA}/gatk SelectVariants -V ${VCF_SNP_FINAL} --exclude-filtered \
  -O /dev/stdout 2>/dev/null | grep -vc "^#" || echo "N/A")
INDEL_PASS=$(${CONDA}/gatk SelectVariants -V ${VCF_INDEL_FINAL} --exclude-filtered \
  -O /dev/stdout 2>/dev/null | grep -vc "^#" || echo "N/A")
echo "  PASS SNPs:   ${SNP_PASS}"
echo "  PASS Indels: ${INDEL_PASS}"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
