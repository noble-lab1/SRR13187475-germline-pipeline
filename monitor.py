#!/usr/bin/env python3
"""
monitor.py — Germline pipeline status dashboard for SRR13187475
Usage:  python3 monitor.py           # one-shot
        python3 monitor.py --watch   # refresh every 30 s
"""

import os
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

WORKDIR  = Path.home() / "SRR13187475"
CONDA    = Path.home() / "miniconda3" / "bin"
MODAL    = Path.home() / "miniconda3" / "bin" / "modal"
SAMTOOLS = Path.home() / "miniconda3" / "bin" / "samtools"

# Modal credentials are auto-loaded by the CLI from ~/.modal.toml
# (run `modal token new` once to populate). No env-var override needed.
VOLUME_NAME        = "gatk-haplotypecaller-data"

# ── Helpers ───────────────────────────────────────────────────────────────────

def gb(path): return path.stat().st_size / 1e9 if path.exists() else None
def mb(path): return path.stat().st_size / 1e6 if path.exists() else None
def age(path):
    if not path.exists(): return "—"
    s = time.time() - path.stat().st_mtime
    if s < 60:   return f"{s:.0f}s ago"
    if s < 3600: return f"{s/60:.0f}m ago"
    return f"{s/3600:.1f}h ago"

def ok(v):  return f"\033[32m✓\033[0m {v}"
def err(v): return f"\033[31m✗\033[0m {v}"
def warn(v):return f"\033[33m⚠\033[0m {v}"
def dim(v): return f"\033[2m{v}\033[0m"
def bold(v):return f"\033[1m{v}\033[0m"

def file_row(label, path, unit="GB"):
    if not path.exists():
        return f"  {label:<30} {err('missing')}"
    size = gb(path) if unit == "GB" else mb(path)
    return f"  {label:<30} {ok(f'{size:.1f} {unit}')}   {dim(age(path))}"

def run(cmd, **kw):
    return subprocess.run(cmd, capture_output=True, text=True, **kw)

# ── Sections ──────────────────────────────────────────────────────────────────

def section(title):
    print(f"\n\033[1;36m{'─'*54}\033[0m")
    print(f"\033[1;36m  {title}\033[0m")
    print(f"\033[1;36m{'─'*54}\033[0m")

def show_local_files():
    section("Local files")
    files = [
        ("SRR13187475_sorted.bam",          "sorted BAM",       "GB"),
        ("SRR13187475_filtered.bam",         "filtered BAM",     "GB"),
        ("SRR13187475_markdup.bam",          "markdup BAM",      "GB"),
        ("SRR13187475_bqsr.bam",             "BQSR BAM",         "GB"),
        ("SRR13187475.g.vcf.gz",             "GVCF (combined)",  "MB"),
        ("SRR13187475_genotyped.vcf.gz",     "genotyped VCF",    "MB"),
        ("SRR13187475_snps_final.vcf.gz",    "final SNPs VCF",   "MB"),
        ("SRR13187475_indels_final.vcf.gz",  "final Indels VCF", "MB"),
    ]
    for fname, label, unit in files:
        path = WORKDIR / fname
        print(file_row(label, path, unit))

def show_pipeline_process():
    section("Running processes")
    procs = run(["ps", "aux"])
    found = []
    keywords = {
        "modal_gatk.py":       "Modal GPU pipeline",
        "gatk_germline.sh":    "gatk_germline.sh",
        "gatk_germline_resume":"gatk_germline_resume.sh",
        "HaplotypeCaller":     "HaplotypeCaller (local)",
        "GenotypeGVCFs":       "GenotypeGVCFs",
        "CombineGVCFs":        "CombineGVCFs",
        "VariantRecalibrator": "VariantRecalibrator (VQSR)",
        "ApplyVQSR":           "ApplyVQSR",
    }
    for line in procs.stdout.splitlines():
        for key, label in keywords.items():
            if key in line and "grep" not in line:
                pid = line.split()[1]
                found.append(f"  {ok(label)}  {dim(f'PID {pid}')}")
                break
    if found:
        for f in found: print(f)
    else:
        print(f"  {dim('No pipeline processes running')}")

def show_modal_log():
    section("Modal pipeline log (last 8 lines)")
    log = WORKDIR / "modal_gatk.log"
    if not log.exists():
        print(f"  {warn('modal_gatk.log not found')}")
        return
    lines = log.read_text().strip().splitlines()
    for line in lines[-8:]:
        if "ERROR" in line or "error" in line.lower():
            print(f"  {err(line)}")
        elif "✓" in line or "done" in line.lower() or "complete" in line.lower():
            print(f"  {ok(line)}")
        else:
            print(f"  {line}")

def show_modal_volume():
    section("Modal volume contents")
    # modal CLI reads ~/.modal.toml directly — no env override.
    dirs = ["/bam", "/ref", "/shards", "/output"]
    for d in dirs:
        r = run([str(MODAL), "volume", "ls", VOLUME_NAME, d])
        lines = [l for l in r.stdout.splitlines() if l.strip()]
        count = len(lines)
        if count:
            # show count and last file
            last = lines[-1].strip()
            print(f"  {d:<12} {ok(f'{count} files')}   {dim(last)}")
        else:
            print(f"  {d:<12} {dim('empty or not found')}")

def show_gatk_log():
    section("GATK germline log (last 6 lines)")
    log = WORKDIR / "gatk_germline.log"
    if not log.exists():
        print(f"  {dim('gatk_germline.log not found')}")
        return
    lines = log.read_text().strip().splitlines()
    for line in lines[-6:]:
        if "ERROR" in line or "Exception" in line:
            print(f"  {err(line)}")
        elif "Progress" in line or "done" in line.lower():
            print(f"  {ok(line)}")
        else:
            print(f"  {line}")

def show_variant_counts():
    section("Variant counts (PASS only)")
    gatk = CONDA / "gatk"
    pairs = [
        ("SNPs",   WORKDIR / "SRR13187475_snps_final.vcf.gz"),
        ("Indels", WORKDIR / "SRR13187475_indels_final.vcf.gz"),
    ]
    for label, vcf in pairs:
        if not vcf.exists():
            print(f"  {label:<10} {dim('not yet generated')}")
            continue
        r = run([str(gatk), "SelectVariants", "-V", str(vcf),
                 "--exclude-filtered", "-O", "/dev/stdout"],
                timeout=120)
        count = sum(1 for l in r.stdout.splitlines() if l and not l.startswith("#"))
        if count:
            print(f"  {label:<10} {ok(f'{count:,} PASS variants')}")
        else:
            print(f"  {label:<10} {warn('0 PASS or count failed')}")

def show_next_step():
    section("Next step")
    gvcf      = WORKDIR / "SRR13187475.g.vcf.gz"
    genotyped = WORKDIR / "SRR13187475_genotyped.vcf.gz"
    snp_final = WORKDIR / "SRR13187475_snps_final.vcf.gz"

    if not gvcf.exists():
        print(f"  Waiting for Modal CombineGVCFs to finish and download GVCF")
        print(f"  Monitor: {dim('tail -f ~/SRR13187475/modal_gatk.log')}")
    elif not genotyped.exists():
        print(f"  {ok('GVCF ready')} — run GenotypeGVCFs:")
        print(f"  {dim('bash ~/SRR13187475/gatk_germline_resume.sh')}")
    elif not snp_final.exists():
        print(f"  {ok('Genotyped VCF ready')} — VQSR/filtering in progress or pending")
        print(f"  {dim('bash ~/SRR13187475/gatk_germline_resume.sh')}")
    else:
        print(f"  {ok('Pipeline complete!')} Final VCFs are ready.")
        print(f"  SNPs:   {snp_final}")
        print(f"  Indels: {WORKDIR / 'SRR13187475_indels_final.vcf.gz'}")

# ── Main ──────────────────────────────────────────────────────────────────────

def run_dashboard():
    os.system("clear")
    print(bold(f"  SRR13187475 Germline Pipeline Monitor  —  {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}"))
    show_local_files()
    show_pipeline_process()
    show_modal_log()
    show_modal_volume()
    show_gatk_log()
    show_variant_counts()
    show_next_step()
    print()

if __name__ == "__main__":
    watch = "--watch" in sys.argv
    interval = 30
    for arg in sys.argv[1:]:
        if arg.startswith("--interval="):
            interval = int(arg.split("=")[1])
    if watch:
        print(f"Watching every {interval}s — Ctrl+C to stop")
        try:
            while True:
                run_dashboard()
                time.sleep(interval)
        except KeyboardInterrupt:
            print("\nStopped.")
    else:
        run_dashboard()
