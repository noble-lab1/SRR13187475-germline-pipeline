#!/bin/bash
# =============================================================================
# Shared configuration for the GATK germline pipeline scripts.
#
# Sourced by: gatk_germline.sh, gatk_germline_resume.sh, modal_gatk.py (via env)
#
# Override any value by exporting it before calling a script, e.g.:
#   SRA_ID=SRR24975641 THREADS=16 SCATTER=32 ./gatk_germline.sh
#
# Every variable uses `:=` so caller-provided env vars win.
# =============================================================================

# ── Toolchain & resources ─────────────────────────────────────────────────────
: "${CONDA:=$HOME/miniconda3/bin}"
: "${THREADS:=8}"
# Scatter count for parallel BQSR + HaplotypeCaller (3-4× THREADS reduces stragglers).
: "${SCATTER:=24}"

# ── Sample + directories ──────────────────────────────────────────────────────
: "${SRA_ID:=SRR13187475}"
: "${WORKDIR:=$HOME/${SRA_ID}}"
: "${REF_DIR:=$HOME/references/hg38}"
: "${REF:=${REF_DIR}/hg38.fa}"

# ── Reference resource files ─────────────────────────────────────────────────
# BQSR known sites
: "${DBSNP:=${REF_DIR}/Homo_sapiens_assembly38.dbsnp138.vcf.gz}"
: "${MILLS:=${REF_DIR}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz}"
# VQSR training resources
: "${HAPMAP:=${REF_DIR}/hapmap_3.3.hg38.vcf.gz}"
: "${OMNI:=${REF_DIR}/1000G_omni2.5.hg38.vcf.gz}"
: "${G1K_SNP:=${REF_DIR}/1000G_phase1.snps.high_confidence.hg38.vcf.gz}"

# ── Input BAM ─────────────────────────────────────────────────────────────────
: "${BAM_IN:=${WORKDIR}/${SRA_ID}_sorted.bam}"

# ── Intermediate BAMs / tables ───────────────────────────────────────────────
: "${BAM_MARKDUP:=${WORKDIR}/${SRA_ID}_markdup.bam}"
: "${BAM_FILTERED:=${WORKDIR}/${SRA_ID}_filtered.bam}"
: "${MARKDUP_METRICS:=${WORKDIR}/${SRA_ID}_markdup_metrics.txt}"
: "${RECAL_TABLE:=${WORKDIR}/${SRA_ID}_recal.table}"
: "${BAM_BQSR:=${WORKDIR}/${SRA_ID}_bqsr.bam}"

# ── Scatter/shard directories ────────────────────────────────────────────────
: "${INTERVALS_DIR:=${WORKDIR}/scatter_intervals}"
: "${BQSR_SHARDS_DIR:=${WORKDIR}/bqsr_shards}"
: "${HC_SHARDS_DIR:=${WORKDIR}/hc_shards}"

# ── GVCF + raw genotyped VCF ─────────────────────────────────────────────────
: "${GVCF:=${WORKDIR}/${SRA_ID}.g.vcf.gz}"
: "${VCF_RAW:=${WORKDIR}/${SRA_ID}_genotyped.vcf.gz}"

# ── VQSR intermediates ───────────────────────────────────────────────────────
: "${SNP_RAW:=${WORKDIR}/${SRA_ID}_snps_raw.vcf.gz}"
: "${INDEL_RAW:=${WORKDIR}/${SRA_ID}_indels_raw.vcf.gz}"
: "${SNP_RECAL:=${WORKDIR}/${SRA_ID}_snps.recal}"
: "${SNP_TRANCHES:=${WORKDIR}/${SRA_ID}_snps.tranches}"
: "${SNP_RSCRIPT:=${WORKDIR}/${SRA_ID}_snps_plots.R}"
: "${INDEL_RECAL:=${WORKDIR}/${SRA_ID}_indels.recal}"
: "${INDEL_TRANCHES:=${WORKDIR}/${SRA_ID}_indels.tranches}"
: "${INDEL_RSCRIPT:=${WORKDIR}/${SRA_ID}_indels_plots.R}"

# ── Output VCFs ──────────────────────────────────────────────────────────────
: "${VCF_SNP_VQSR:=${WORKDIR}/${SRA_ID}_snps_vqsr.vcf.gz}"
: "${VCF_INDEL_VQSR:=${WORKDIR}/${SRA_ID}_indels_vqsr.vcf.gz}"
: "${VCF_SNP_FINAL:=${WORKDIR}/${SRA_ID}_snps_final.vcf.gz}"
: "${VCF_INDEL_FINAL:=${WORKDIR}/${SRA_ID}_indels_final.vcf.gz}"

export PATH=${CONDA}:${PATH}
