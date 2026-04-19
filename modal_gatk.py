#!/usr/bin/env python3
"""
modal_gatk.py
=============
Parallelized GATK HaplotypeCaller on Modal cloud infrastructure.

Strategy:
  - Splits BAM by chromosome locally (200-400 MB each) for reliable uploads
  - Reference genome downloaded directly on Modal from UCSC (avoids 3 GB upload)
  - Each chromosome runs as a separate Parabricks GPU job in parallel
  - GVCFs combined and downloaded back locally

Usage:
  python3 modal_gatk.py \\
    --bam  ~/SRR13187475/SRR13187475_bqsr.bam \\
    --ref  ~/references/hg38/hg38.fa          \\
    --out  ~/SRR13187475/SRR13187475.g.vcf.gz
"""

import argparse
import os
import subprocess
import sys
import time
from pathlib import Path

# ── Config ────────────────────────────────────────────────────────────────────
CHROMS      = [f"chr{i}" for i in list(range(1, 23)) + ["X", "Y"]]
VOLUME_NAME = "gatk-haplotypecaller-data"
GATK_VER    = "4.6.2.0"
MODAL_BIN   = os.path.expanduser("~/miniconda3/bin/modal")
SAMTOOLS    = os.path.expanduser("~/miniconda3/bin/samtools")

# Credentials
DEFAULT_MODAL_TOKEN_ID     = "ak-8OHryx1aMl2XEwM0oknuBE"
DEFAULT_MODAL_TOKEN_SECRET = "as-7NggUhXzihTOdJHQc1WD0R"
DEFAULT_NGC_API_KEY        = "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"

# hg38 reference — downloaded on Modal directly from UCSC (avoids uploading 3 GB)
HG38_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"


# ── Args ──────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--modal-token-id",     default=DEFAULT_MODAL_TOKEN_ID)
    p.add_argument("--modal-token-secret", default=DEFAULT_MODAL_TOKEN_SECRET)
    p.add_argument("--ngc-api-key",        default=DEFAULT_NGC_API_KEY)
    p.add_argument("--bam",  required=True, type=Path)
    p.add_argument("--ref",  required=True, type=Path,
                   help="Local reference (used for fai/dict only; Modal downloads its own copy)")
    p.add_argument("--out",  required=True, type=Path)
    p.add_argument("--skip-split",  action="store_true",
                   help="Skip BAM splitting if chr BAMs already exist in /tmp/chr_bams/")
    p.add_argument("--skip-upload", action="store_true",
                   help="Skip upload if chr BAMs already in Modal volume")
    p.add_argument("--skip-shards", action="store_true",
                   help="Skip GPU shards if GVCFs already in Modal volume; go straight to combine")
    p.add_argument("--chroms", nargs="+", default=CHROMS)
    return p.parse_args()


# ── BAM splitting ─────────────────────────────────────────────────────────────

def split_bam_by_chrom(bam: Path, chroms: list, out_dir: Path) -> dict:
    """
    Split BAM into per-chromosome BAMs using samtools.
    Returns {chrom: Path} mapping.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    chr_bams = {}
    for chrom in chroms:
        out_bam = out_dir / f"{bam.stem}_{chrom}.bam"
        if out_bam.exists():
            print(f"  {chrom}: already split, skipping")
            chr_bams[chrom] = out_bam
            continue
        print(f"  Splitting {chrom}...", end=" ", flush=True)
        subprocess.run(
            [SAMTOOLS, "view", "-b", "-o", str(out_bam), str(bam), chrom],
            check=True
        )
        subprocess.run([SAMTOOLS, "index", str(out_bam)], check=True)
        size_mb = out_bam.stat().st_size / 1e6
        print(f"{size_mb:.0f} MB")
        chr_bams[chrom] = out_bam
    return chr_bams


# ── Upload ────────────────────────────────────────────────────────────────────

def file_exists_in_volume(remote_path: str) -> bool:
    """Check if a file already exists in the Modal volume."""
    remote_dir  = str(Path(remote_path).parent)
    remote_name = Path(remote_path).name
    result = subprocess.run(
        [MODAL_BIN, "volume", "ls", VOLUME_NAME, remote_dir],
        capture_output=True, text=True
    )
    return remote_name in result.stdout


def upload_file(local_path: Path, remote_path: str, retries: int = 8):
    """Upload via modal CLI with exponential backoff retries."""
    size_mb = local_path.stat().st_size / 1e6
    print(f"  ↑ {local_path.name}  ({size_mb:.0f} MB)", end=" ", flush=True)
    if file_exists_in_volume(remote_path):
        print("already in volume, skipping")
        return
    for attempt in range(1, retries + 1):
        try:
            subprocess.run(
                [MODAL_BIN, "volume", "put", "--force", VOLUME_NAME,
                 str(local_path), remote_path],
                check=True, timeout=900,
            )
            print("✓")
            return
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as e:
            if attempt < retries:
                wait = min(2 ** attempt, 120)
                print(f"\n    Retry {attempt}/{retries} in {wait}s: {type(e).__name__}", end=" ", flush=True)
                time.sleep(wait)
            else:
                raise RuntimeError(f"Upload failed after {retries} attempts: {e}")


# ── Download ──────────────────────────────────────────────────────────────────

def download_file(remote_path: str, local_path: Path):
    print(f"  ↓ {Path(remote_path).name} → {local_path}")
    subprocess.run(
        [MODAL_BIN, "volume", "get", VOLUME_NAME, remote_path, str(local_path)],
        check=True
    )


# ── Modal app ─────────────────────────────────────────────────────────────────

def run_on_modal(chroms, bam_stem, out_name, ngc_key, skip_shards=False):
    import modal

    volume = modal.Volume.from_name(VOLUME_NAME, create_if_missing=True)

    registry_secret = modal.Secret.from_dict({
        "REGISTRY_USERNAME": "$oauthtoken",
        "REGISTRY_PASSWORD": ngc_key,
    })
    image = (
        modal.Image.from_registry(
            "nvcr.io/nvidia/clara/clara-parabricks:4.3.1-1",
            secret=registry_secret,
            add_python="3.13",
        )
        .apt_install("wget", "samtools", "unzip", "openjdk-17-jre-headless")
        .run_commands(
            f"wget -q -O /tmp/gatk.zip https://github.com/broadinstitute/gatk/releases/download/{GATK_VER}/gatk-{GATK_VER}.zip",
            f"unzip -q /tmp/gatk.zip -d /opt",
            f"ln -sf /opt/gatk-{GATK_VER}/gatk /usr/local/bin/gatk",
            "rm /tmp/gatk.zip",
        )
    )

    app = modal.App("gatk-haplotypecaller")

    # ── Setup function: download reference on Modal once ──────────────────────
    @app.function(
        image=image,
        cpu=4,
        memory=8192,
        volumes={"/data": volume},
        timeout=3600,
        serialized=True,
    )
    def setup_reference():
        import subprocess, os
        ref_path  = "/data/ref/hg38.fa"
        dict_path = "/data/ref/hg38.dict"
        fai_path  = ref_path + ".fai"
        need_commit = False

        if not os.path.exists(ref_path):
            os.makedirs("/data/ref", exist_ok=True)
            print("Downloading hg38 from UCSC...")
            subprocess.run(
                ["wget", "-q", "-O", "/data/ref/hg38.fa.gz", HG38_URL], check=True
            )
            subprocess.run(["gzip", "-d", "/data/ref/hg38.fa.gz"], check=True)
            need_commit = True
        else:
            print("Reference FASTA already present")

        if not os.path.exists(fai_path):
            subprocess.run(["samtools", "faidx", ref_path], check=True)
            need_commit = True

        if not os.path.exists(dict_path):
            print("Creating sequence dictionary...")
            subprocess.run(
                ["samtools", "dict", ref_path, "-o", dict_path], check=True
            )
            need_commit = True

        if need_commit:
            volume.commit()
        print("Reference ready.")
        return ref_path

    # ── Shard: one chromosome ─────────────────────────────────────────────────
    @app.function(
        image=image,
        gpu="A100-80GB",
        cpu=16,
        memory=64 * 1024,
        volumes={"/data": volume},
        timeout=7200,
        retries=1,
        serialized=True,
    )
    def haplotypecaller_shard(chrom: str, bam_stem: str, out_prefix: str) -> str:
        import subprocess, os
        os.makedirs("/data/shards", exist_ok=True)
        bam_path  = f"/data/bam/{bam_stem}_{chrom}.bam"
        ref_path  = "/data/ref/hg38.fa"
        out_gvcf  = f"/data/shards/{out_prefix}_{chrom}.g.vcf.gz"

        cmd = [
            "pbrun", "haplotypecaller",
            "--ref",          ref_path,
            "--in-bam",       bam_path,
            "--out-variants", out_gvcf,
            "--gvcf",
            "--num-gpus",     "1",
        ]
        print(f"[{chrom}] pbrun haplotypecaller started")
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"[{chrom}] failed:\n{result.stderr[-4000:]}")
        print(f"[{chrom}] done → {out_gvcf}")
        volume.commit()
        return out_gvcf

    # ── Combine ───────────────────────────────────────────────────────────────
    @app.function(
        image=image,
        cpu=4,
        memory=16 * 1024,
        volumes={"/data": volume},
        timeout=7200,
        serialized=True,
    )
    def combine_gvcfs(shard_paths: list, out_gvcf: str) -> str:
        import subprocess, os
        os.makedirs("/data/output", exist_ok=True)
        out_path = f"/data/output/{out_gvcf}"
        inputs = []
        for p in sorted(shard_paths):
            inputs += ["-V", p]
        cmd = [
            "gatk", "CombineGVCFs",
            "-R", "/data/ref/hg38.fa",
            *inputs,
            "-O", out_path,
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"CombineGVCFs failed:\n{result.stderr[-4000:]}")
        volume.commit()
        return out_path

    out_prefix = out_name.replace(".g.vcf.gz", "").replace(".vcf.gz", "")

    with app.run():
        print("  Setting up reference genome on Modal...")
        setup_reference.remote()

        if not skip_shards:
            print(f"  Launching {len(chroms)} parallel Parabricks shards...")
            shard_paths = list(
                haplotypecaller_shard.starmap([
                    (chrom, bam_stem, out_prefix) for chrom in chroms
                ])
            )
        else:
            shard_paths = [
                f"/data/shards/{out_prefix}_{chrom}.g.vcf.gz" for chrom in chroms
            ]
            print(f"  Skipping shards (--skip-shards); using {len(shard_paths)} pre-computed GVCFs")

        print(f"\n  Combining {len(shard_paths)} GVCFs...")
        final_remote = combine_gvcfs.remote(shard_paths, out_name)

    return final_remote


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()

    os.environ["MODAL_TOKEN_ID"]     = args.modal_token_id
    os.environ["MODAL_TOKEN_SECRET"] = args.modal_token_secret

    try:
        import modal  # noqa: F401
    except ImportError:
        subprocess.run([sys.executable, "-m", "pip", "install", "modal", "-q"], check=True)

    bam   = args.bam.expanduser().resolve()
    bai   = Path(str(bam) + ".bai")
    if not bai.exists():
        bai = bam.with_suffix(".bai")
    out   = args.out.expanduser().resolve()

    for f in [bam, bai]:
        if not f.exists():
            sys.exit(f"ERROR: File not found: {f}")

    print(f"""
╔══════════════════════════════════════════════════════╗
║    Modal GATK HaplotypeCaller  (Parabricks A100)    ║
╠══════════════════════════════════════════════════════╣
  BAM:         {bam.name}  ({bam.stat().st_size / 1e9:.1f} GB)
  Output:      {out}
  Chromosomes: {len(args.chroms)}
  Strategy:    Split BAM → upload per-chrom → GPU scatter
╚══════════════════════════════════════════════════════╝
""")

    # ── Step 1: Split BAM by chromosome ──────────────────────────────────────
    chr_bam_dir = Path("/tmp/chr_bams")
    if not args.skip_split:
        print("Step 1/4 — Splitting BAM by chromosome...")
        chr_bams = split_bam_by_chrom(bam, args.chroms, chr_bam_dir)
    else:
        print("Step 1/4 — Skipping split (--skip-split)")
        chr_bams = {c: chr_bam_dir / f"{bam.stem}_{c}.bam" for c in args.chroms}

    # ── Step 2: Upload per-chrom BAMs ─────────────────────────────────────────
    if not args.skip_upload:
        print("\nStep 2/4 — Uploading per-chromosome BAMs to Modal volume...")
        print(f"  Volume: {VOLUME_NAME}")
        subprocess.run([MODAL_BIN, "volume", "create", VOLUME_NAME], capture_output=True)
        for chrom, chr_bam in chr_bams.items():
            upload_file(chr_bam,                        f"/bam/{chr_bam.name}")
            upload_file(chr_bam.with_suffix(".bam.bai"), f"/bam/{chr_bam.name}.bai")
        print("All uploads complete.\n")
    else:
        print("Step 2/4 — Skipping upload (--skip-upload)\n")

    # ── Step 3: Run on Modal ──────────────────────────────────────────────────
    print("Step 3/4 — Running Parabricks HaplotypeCaller on Modal A100...")
    remote_gvcf = run_on_modal(
        chroms      = args.chroms,
        bam_stem    = bam.stem,
        out_name    = out.name,
        ngc_key     = args.ngc_api_key,
        skip_shards = args.skip_shards,
    )

    # ── Step 4: Download result ───────────────────────────────────────────────
    print("\nStep 4/4 — Downloading GVCF...")
    download_file(remote_gvcf, out)
    try:
        download_file(remote_gvcf + ".tbi", Path(str(out) + ".tbi"))
    except Exception:
        pass

    if out.exists():
        print(f"""
╔══════════════════════════════════════════════════════╗
║  Done!                                              ║
╠══════════════════════════════════════════════════════╣
  GVCF: {out}

  Next step:
    bash ~/SRR13187475/gatk_germline_resume.sh
╚══════════════════════════════════════════════════════╝
""")
    else:
        sys.exit("ERROR: Output GVCF not found after download.")


if __name__ == "__main__":
    main()
