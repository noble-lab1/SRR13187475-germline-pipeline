#!/bin/bash
# =============================================================================
# GATK Germline Variant Calling Pipeline — SRR13187475
# Steps:
#   1. Mark duplicates (on the full sorted BAM — uses all reads to anchor mates)
#   2. Filter to unique (MAPQ >= 20) + concordantly mapped pairs
#   3. Scatter intervals + sequence dictionary
#   4. Base Quality Score Recalibration (BQSR), scattered
#   5. HaplotypeCaller (GVCF mode), scattered
#   6. GenotypeGVCFs
#   7. Variant Quality Score Recalibration (VQSR) — SNPs then Indels
#   8. Hard filter fallback (applied to variants that fail VQSR or if VQSR skipped)
# =============================================================================
set -euo pipefail

# All paths, sample ID, and parallelism knobs live in config.sh — override any
# of them by exporting before running (e.g. `SRA_ID=SRR24975641 ./gatk_germline.sh`).
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

# ── Step 1: Mark Duplicates ───────────────────────────────────────────────────
# GATK Best Practices order: duplicate marking runs on the FULL sorted BAM
# before any quality/pair filtering. This is important because markdup uses
# mate positions to anchor duplicates — if a mate fails MAPQ/properly-paired
# filters and is removed first, a real PCR duplicate can be missed.
#
# Using samtools markdup instead of Picard/GATK MarkDuplicates: multithreaded
# (3-5× faster on large WGS BAMs) and produces equivalent flags for downstream
# BQSR / HaplotypeCaller. Pipeline:
#   collate   — groups reads by name (required for mate awareness)
#   fixmate   — adds MC/ms tags used by markdup
#   sort      — restore coordinate order
#   markdup   — flag 0x400 on PCR/optical duplicates; write stats
echo "[1/8] Marking duplicates (samtools markdup) on full sorted BAM"
MARKDUP_TMP=${WORKDIR}/.markdup_tmp
${CONDA}/samtools collate -@ ${THREADS} -O -u ${BAM_IN} ${MARKDUP_TMP} \
  | ${CONDA}/samtools fixmate -@ ${THREADS} -m -u - - \
  | ${CONDA}/samtools sort    -@ ${THREADS} -u - \
  | ${CONDA}/samtools markdup -@ ${THREADS} -s -f ${MARKDUP_METRICS} - ${BAM_MARKDUP}

${CONDA}/samtools index -@ ${THREADS} ${BAM_MARKDUP}
rm -f ${MARKDUP_TMP}.*.bam

echo "  Duplicate metrics: ${MARKDUP_METRICS}"

# ── Step 2: Filter to unique + concordantly mapped pairs ──────────────────────
# Apply after markdup so the 0x400 dup flag is preserved on surviving reads;
# downstream BQSR/HC respect the flag internally.
# -f 0x2   = properly paired (concordant)
# -q 20    = MAPQ >= 20 (unique mapping)
# -F 0x4   = not unmapped
# -F 0x100 = not secondary alignment
# -F 0x800 = not supplementary alignment
echo "[2/8] Filtering: unique (MAPQ >= 20) + concordantly mapped pairs"
${CONDA}/samtools view -@ ${THREADS} \
  -f 0x2 \
  -F 0x4 -F 0x100 -F 0x800 \
  -q 20 \
  -b ${BAM_MARKDUP} \
  -o ${BAM_FILTERED}

${CONDA}/samtools index -@ ${THREADS} ${BAM_FILTERED}

# idxstats reads only the BAM index (mapped + unmapped per ref), much faster
# than re-walking the BAM with `samtools view -c`.
TOTAL_IN=$(${CONDA}/samtools idxstats ${BAM_MARKDUP}   | awk '{s+=$3+$4} END {print s}')
TOTAL_FILT=$(${CONDA}/samtools idxstats ${BAM_FILTERED} | awk '{s+=$3+$4} END {print s}')
echo "  Reads in (markdup):  ${TOTAL_IN}"
echo "  Reads kept:          ${TOTAL_FILT}"

# ── Step 3: Create sequence dictionary + scatter intervals ───────────────────
if [ ! -f "${REF_DIR}/hg38.dict" ]; then
  echo "[3/8] Creating sequence dictionary for hg38..."
  ${CONDA}/gatk CreateSequenceDictionary -R ${REF}
else
  echo "[3/8] Sequence dictionary already present"
fi

# SplitIntervals produces ${SCATTER} evenly-sized .interval_list files used by
# the parallel BQSR (step 4) and HaplotypeCaller (step 5) scatter steps.
# SUBDIVISION_BALANCING_WITHOUT_INTERVAL_SUBDIVISION keeps intervals contiguous
# (important for BQSR and HC — avoids cross-shard edge effects at break points).
if [ ! -d "${INTERVALS_DIR}" ] || [ -z "$(ls -A ${INTERVALS_DIR} 2>/dev/null)" ]; then
  echo "  Scattering intervals (${SCATTER} shards)..."
  mkdir -p ${INTERVALS_DIR}
  ${CONDA}/gatk SplitIntervals \
    -R ${REF} \
    --scatter-count ${SCATTER} \
    --subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION \
    -O ${INTERVALS_DIR}
else
  echo "  Scatter intervals already present (${INTERVALS_DIR})"
fi

# ── Step 4: Base Quality Score Recalibration (BQSR) ──────────────────────────
# BQSR requires known-variant VCFs from the GATK resource bundle. If they are
# not present, BQSR is skipped and the markdup BAM is used directly.
#
# Public (no requester-pays, no billing project needed):
#   BASE=gs://gcp-public-data--broad-references/hg38/v0
#   gcloud storage cp \
#     $BASE/Homo_sapiens_assembly38.dbsnp138.vcf.gz{,.tbi} \
#     $BASE/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz{,.tbi} \
#     $BASE/hapmap_3.3.hg38.vcf.gz{,.tbi} \
#     $BASE/1000G_omni2.5.hg38.vcf.gz{,.tbi} \
#     $BASE/1000G_phase1.snps.high_confidence.hg38.vcf.gz{,.tbi} \
#     ~/references/hg38/

if [ -f "${DBSNP}" ] && [ -f "${MILLS}" ]; then
  mkdir -p ${BQSR_SHARDS_DIR}

  # Exported for the bash -c subshell used by xargs -P below.
  export CONDA REF BAM_FILTERED DBSNP MILLS RECAL_TABLE BAM_BQSR BQSR_SHARDS_DIR

  # Per-shard BaseRecalibrator: each shard writes a small recal table that is
  # later merged with GatherBQSRReports.
  bqsr_base_shard() {
    local interval=$1
    local name; name=$(basename "$interval" .interval_list)
    ${CONDA}/gatk BaseRecalibrator \
      -I ${BAM_FILTERED} \
      -R ${REF} \
      -L "$interval" \
      --known-sites ${DBSNP} \
      --known-sites ${MILLS} \
      -O ${BQSR_SHARDS_DIR}/recal_${name}.table \
      >/dev/null 2>&1 && echo "  [base/$name] done" \
      || { echo "  [base/$name] FAILED"; exit 1; }
  }
  export -f bqsr_base_shard

  # Per-shard ApplyBQSR: each shard writes a piece of the recalibrated BAM,
  # later concatenated with GatherBamFiles.
  bqsr_apply_shard() {
    local interval=$1
    local name; name=$(basename "$interval" .interval_list)
    ${CONDA}/gatk ApplyBQSR \
      -I ${BAM_FILTERED} \
      -R ${REF} \
      -L "$interval" \
      --bqsr-recal-file ${RECAL_TABLE} \
      -O ${BQSR_SHARDS_DIR}/bqsr_${name}.bam \
      >/dev/null 2>&1 && echo "  [apply/$name] done" \
      || { echo "  [apply/$name] FAILED"; exit 1; }
  }
  export -f bqsr_apply_shard

  echo "[4/8] BQSR — BaseRecalibrator (${SCATTER}-way scatter, ${THREADS} parallel)"
  find ${INTERVALS_DIR} -name "*.interval_list" -print0 \
    | xargs -0 -n 1 -P ${THREADS} -I{} bash -c 'bqsr_base_shard "$@"' _ {}

  echo "[4/8] BQSR — GatherBQSRReports"
  ${CONDA}/gatk GatherBQSRReports \
    $(for f in $(ls ${BQSR_SHARDS_DIR}/recal_*.table | sort); do printf -- "-I %s " "$f"; done) \
    -O ${RECAL_TABLE}

  echo "[4/8] BQSR — ApplyBQSR (${SCATTER}-way scatter, ${THREADS} parallel)"
  find ${INTERVALS_DIR} -name "*.interval_list" -print0 \
    | xargs -0 -n 1 -P ${THREADS} -I{} bash -c 'bqsr_apply_shard "$@"' _ {}

  echo "[4/8] BQSR — GatherBamFiles"
  ${CONDA}/gatk GatherBamFiles \
    $(for f in $(ls ${BQSR_SHARDS_DIR}/bqsr_*.bam | sort); do printf -- "-I %s " "$f"; done) \
    -O ${BAM_BQSR} \
    --CREATE_INDEX true
else
  echo "[4/8] BQSR resource VCFs not found — skipping BQSR"
  echo "  NOTE: Using filtered BAM directly. Base quality scores will not be recalibrated."
  echo "  To enable BQSR, download the GATK resource bundle (see comments above)."
  BAM_BQSR=${BAM_FILTERED}
fi

# ── Step 5: HaplotypeCaller (GVCF mode, scattered) ───────────────────────────
# HC is the wall-clock dominator and is effectively single-threaded —
# --native-pair-hmm-threads only accelerates the PairHMM inner loop, not the
# assembly/active-region stages. Scattering across ${SCATTER} intervals and
# running ${THREADS} in parallel gives a near-linear speedup.
# Each shard uses --native-pair-hmm-threads 1 (more shard parallelism > SIMD).
mkdir -p ${HC_SHARDS_DIR}

export CONDA REF BAM_BQSR HC_SHARDS_DIR

hc_shard() {
  local interval=$1
  local name; name=$(basename "$interval" .interval_list)
  ${CONDA}/gatk HaplotypeCaller \
    -R ${REF} \
    -I ${BAM_BQSR} \
    -L "$interval" \
    -O ${HC_SHARDS_DIR}/shard_${name}.g.vcf.gz \
    -ERC GVCF \
    --native-pair-hmm-threads 1 \
    >/dev/null 2>&1 && echo "  [hc/$name] done" \
    || { echo "  [hc/$name] FAILED"; exit 1; }
}
export -f hc_shard

echo "[5/8] HaplotypeCaller (${SCATTER}-way scatter, ${THREADS} parallel)"
find ${INTERVALS_DIR} -name "*.interval_list" -print0 \
  | xargs -0 -n 1 -P ${THREADS} -I{} bash -c 'hc_shard "$@"' _ {}

echo "[5/8] HaplotypeCaller — MergeVcfs (GVCF)"
${CONDA}/gatk MergeVcfs \
  $(for f in $(ls ${HC_SHARDS_DIR}/shard_*.g.vcf.gz | sort); do printf -- "-I %s " "$f"; done) \
  -O ${GVCF}

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

# Allele Balance (AB): for het calls, alt allele fraction should be 0.2–0.8.
# Uses JEXL to access the AD (allele depth) FORMAT field. Only applied to hets
# — homozygous alt calls are expected to have near-100% alt reads.
# MQ0: number of reads with mapping quality 0 at the site. High MQ0 indicates
# a repetitive/low-complexity region prone to mismapping artifacts.
# MQ0 filter applied as a fraction of total depth (MQ0 / DP > 0.1).

# SNP + Indel VariantFiltration are independent — run concurrently.
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
  -O ${VCF_SNP_FINAL} &
SNP_VF_PID=$!

${CONDA}/gatk VariantFiltration \
  -V ${VCF_INDEL_VQSR} \
  --filter-expression "QD < 2.0"               --filter-name "QD2"              \
  --filter-expression "FS > 200.0"             --filter-name "FS200"            \
  --filter-expression "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
  --filter-expression "MQ0 / DP > 0.1"         --filter-name "MQ0Frac"          \
  --filter-expression "vc.isHet() && (vc.getGenotype(0).getAD()[1] * 1.0 / (vc.getGenotype(0).getAD()[0] + vc.getGenotype(0).getAD()[1]) < 0.2 || vc.getGenotype(0).getAD()[1] * 1.0 / (vc.getGenotype(0).getAD()[0] + vc.getGenotype(0).getAD()[1]) > 0.8)" \
    --filter-name "ABfilter" \
  -O ${VCF_INDEL_FINAL} &
INDEL_VF_PID=$!

wait ${SNP_VF_PID} ${INDEL_VF_PID}

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

# bcftools reads the bgzipped VCF directly via htslib (no double-decompression,
# no JVM startup) — ~20× faster than the SelectVariants+grep approach.
SNP_PASS=$(${CONDA}/bcftools view -f PASS -H ${VCF_SNP_FINAL} 2>/dev/null | wc -l | tr -d ' ' || echo "N/A")
INDEL_PASS=$(${CONDA}/bcftools view -f PASS -H ${VCF_INDEL_FINAL} 2>/dev/null | wc -l | tr -d ' ' || echo "N/A")
echo "  PASS SNPs:   ${SNP_PASS}"
echo "  PASS Indels: ${INDEL_PASS}"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
