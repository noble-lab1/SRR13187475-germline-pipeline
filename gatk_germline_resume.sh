#!/bin/bash
# =============================================================================
# gatk_germline_resume.sh
# Resume GATK pipeline from GenotypeGVCFs onward
# (used after modal_gatk.py returns the GVCF)
# =============================================================================
set -euo pipefail

CONDA=~/miniconda3/bin
SRA_ID=SRR13187475
WORKDIR=~/${SRA_ID}
REF_DIR=~/references/hg38
REF=${REF_DIR}/hg38.fa
THREADS=8

# Use java -jar directly to avoid the gatk wrapper's broken python shebang
JAVA=${CONDA}/java
GATK_JAR=$(ls ${CONDA}/../share/gatk4-*/gatk-package-*-local.jar 2>/dev/null | tail -1)
if [ -z "${GATK_JAR}" ]; then echo "ERROR: GATK jar not found"; exit 1; fi
gatk() { ${JAVA} -jar "${GATK_JAR}" "$@"; }

GVCF=${WORKDIR}/${SRA_ID}.g.vcf.gz
VCF_RAW=${WORKDIR}/${SRA_ID}_genotyped.vcf.gz
SNP_RAW=${WORKDIR}/${SRA_ID}_snps_raw.vcf.gz
INDEL_RAW=${WORKDIR}/${SRA_ID}_indels_raw.vcf.gz
SNP_RECAL=${WORKDIR}/${SRA_ID}_snps.recal
SNP_TRANCHES=${WORKDIR}/${SRA_ID}_snps.tranches
SNP_RSCRIPT=${WORKDIR}/${SRA_ID}_snps_plots.R
INDEL_RECAL=${WORKDIR}/${SRA_ID}_indels.recal
INDEL_TRANCHES=${WORKDIR}/${SRA_ID}_indels.tranches
INDEL_RSCRIPT=${WORKDIR}/${SRA_ID}_indels_plots.R
VCF_SNP_VQSR=${WORKDIR}/${SRA_ID}_snps_vqsr.vcf.gz
VCF_INDEL_VQSR=${WORKDIR}/${SRA_ID}_indels_vqsr.vcf.gz
VCF_SNP_FINAL=${WORKDIR}/${SRA_ID}_snps_final.vcf.gz
VCF_INDEL_FINAL=${WORKDIR}/${SRA_ID}_indels_final.vcf.gz

DBSNP=${REF_DIR}/Homo_sapiens_assembly38.dbsnp138.vcf.gz
MILLS=${REF_DIR}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
HAPMAP=${REF_DIR}/hapmap_3.3.hg38.vcf.gz
OMNI=${REF_DIR}/1000G_omni2.5.hg38.vcf.gz
G1K_SNP=${REF_DIR}/1000G_phase1.snps.high_confidence.hg38.vcf.gz

export PATH=${CONDA}:$PATH

if [ ! -f "${VCF_RAW}" ]; then
  echo "ERROR: Genotyped VCF not found: ${VCF_RAW}"
  echo "Run modal_gatk.py first to generate the genotyped VCF."
  exit 1
fi

# ── SelectVariants ─────────────────────────────────────────────────────────────
echo "[1/4] SelectVariants (SNPs + Indels)"
gatk SelectVariants -V ${VCF_RAW} --select-type-to-include SNP   -O ${SNP_RAW}
gatk SelectVariants -V ${VCF_RAW} --select-type-to-include INDEL -O ${INDEL_RAW}

# ── VQSR ──────────────────────────────────────────────────────────────────────
echo "[2/4] VQSR — SNPs"
gatk VariantRecalibrator \
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
  --rscript-file  ${SNP_RSCRIPT}

gatk ApplyVQSR \
  -V ${SNP_RAW} \
  --recal-file ${SNP_RECAL} \
  --tranches-file ${SNP_TRANCHES} \
  --truth-sensitivity-filter-level 99.5 \
  --create-output-variant-index true \
  -mode SNP \
  -O ${VCF_SNP_VQSR}

echo "[2/4] VQSR — Indels"
gatk VariantRecalibrator \
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
  --rscript-file  ${INDEL_RSCRIPT}

gatk ApplyVQSR \
  -V ${INDEL_RAW} \
  --recal-file ${INDEL_RECAL} \
  --tranches-file ${INDEL_TRANCHES} \
  --truth-sensitivity-filter-level 99.0 \
  --create-output-variant-index true \
  -mode INDEL \
  -O ${VCF_INDEL_VQSR}

# ── Hard filter fallback ──────────────────────────────────────────────────────
echo "[3/4] Hard filter fallback (on top of VQSR)"
gatk VariantFiltration \
  -V ${VCF_SNP_VQSR} \
  --filter-expression "QD < 2.0"              --filter-name "QD2"             \
  --filter-expression "FS > 60.0"             --filter-name "FS60"            \
  --filter-expression "MQ < 40.0"             --filter-name "MQ40"            \
  --filter-expression "MQRankSum < -12.5"     --filter-name "MQRankSum-12.5"  \
  --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
  --filter-expression "MQ0 / DP > 0.1"        --filter-name "MQ0Frac"         \
  -O ${VCF_SNP_FINAL}

gatk VariantFiltration \
  -V ${VCF_INDEL_VQSR} \
  --filter-expression "QD < 2.0"               --filter-name "QD2"              \
  --filter-expression "FS > 200.0"             --filter-name "FS200"            \
  --filter-expression "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
  --filter-expression "MQ0 / DP > 0.1"         --filter-name "MQ0Frac"          \
  -O ${VCF_INDEL_FINAL}

# ── Summary ───────────────────────────────────────────────────────────────────
echo "[4/4] Summary"
SNP_PASS=$(gatk SelectVariants -V ${VCF_SNP_FINAL} --exclude-filtered \
  -O /dev/stdout 2>/dev/null | grep -vc "^#" || echo "N/A")
INDEL_PASS=$(gatk SelectVariants -V ${VCF_INDEL_FINAL} --exclude-filtered \
  -O /dev/stdout 2>/dev/null | grep -vc "^#" || echo "N/A")

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "Pipeline complete."
echo "  Final SNPs VCF:   ${VCF_SNP_FINAL}"
echo "  Final Indels VCF: ${VCF_INDEL_FINAL}"
echo "  PASS SNPs:        ${SNP_PASS}"
echo "  PASS Indels:      ${INDEL_PASS}"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
