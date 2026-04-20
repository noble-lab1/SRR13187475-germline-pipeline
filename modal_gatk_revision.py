#!/usr/bin/env python3
"""
modal_gatk.py
=============
Parallelized GATK HaplotypeCaller on Modal cloud infrastructure.

Strategy:
  - Splits BAM by chromosome locally (200-400 MB each) for reliable uploads
  - Reference genome downloaded directly on Modal from UCSC (avoids 3 GB upload)
  - Each chromosome runs as a separate Parabricks GPU job in parallel
  - Per-chromosome VCFs genotyped on Modal, downloaded, merged locally

Credentials:
  Hardcoded defaults are set below for convenience during development.
  Override via CLI flags, environment variables, or `~/.modal_gatk.env`.
  ⚠️  ROTATE THESE BEFORE SHARING THIS FILE OR COMMITTING TO GIT.

Usage:
  python3 modal_gatk.py \\
    --bam  ~/SRR13187475/SRR13187475_bqsr.bam \\
    --ref  ~/references/hg38/hg38.fa          \\
    --out  ~/SRR13187475/SRR13187475.g.vcf.gz
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

# ── Config ────────────────────────────────────────────────────────────────────
CHROMS      = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
VOLUME_NAME = "gatk-haplotypecaller-data"
GATK_VER    = "4.6.2.0"

# hg38 reference — downloaded on Modal directly from UCSC (avoids uploading 3 GB)
HG38_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"

# Manifest for resumability
MANIFEST_PATH = Path.home() / ".modal_gatk_manifest.json"

# Credentials removed — Modal reads ~/.modal.toml, NGC key comes from
# ~/.modal_gatk.env or the NGC_API_KEY env var. See modal_gatk.py for the
# canonical implementation; this revision file is kept for reference only.
DEFAULT_MODAL_TOKEN_ID     = None
DEFAULT_MODAL_TOKEN_SECRET = None
DEFAULT_NGC_API_KEY        = None


# ── Binary resolution ─────────────────────────────────────────────────────────

def find_binary(name: str, fallback_paths: list = None) -> str:
    """Locate an executable, preferring PATH but falling back to common locations."""
    found = shutil.which(name)
    if found:
        return found
    for candidate in (fallback_paths or []):
        p = Path(candidate).expanduser()
        if p.exists() and os.access(p, os.X_OK):
            return str(p)
    sys.exit(f"ERROR: could not find '{name}' on PATH or in fallback locations")


MODAL_BIN = find_binary("modal", ["~/miniconda3/bin/modal", "~/.local/bin/modal"])
SAMTOOLS  = find_binary("samtools", ["~/miniconda3/bin/samtools", "/usr/local/bin/samtools"])


# ── Credential loading ───────────────────────────────────────────────────────

def load_env_file(path: Path) -> dict:
    """Parse a simple KEY=VALUE env file. Returns {} if absent."""
    if not path.exists():
        return {}
    env = {}
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        if "=" in line:
            k, v = line.split("=", 1)
            env[k.strip()] = v.strip().strip('"').strip("'")
    return env


def get_credential(name: str, cli_value: str, env_file: dict,
                   default: str = None, required: bool = True) -> str:
    """Resolve a credential from (1) CLI, (2) env var, (3) env file, (4) hardcoded default."""
    val = cli_value or os.environ.get(name) or env_file.get(name) or default
    if not val and required:
        sys.exit(
            f"ERROR: {name} not set. Provide via --{name.lower().replace('_', '-')}, "
            f"environment variable, or ~/.modal_gatk.env"
        )
    return val


# ── Manifest for resumability ─────────────────────────────────────────────────

def load_manifest() -> dict:
    if MANIFEST_PATH.exists():
        try:
            return json.loads(MANIFEST_PATH.read_text())
        except json.JSONDecodeError:
            print(f"  WARN: manifest at {MANIFEST_PATH} is corrupt, starting fresh")
    return {"split": {}, "uploaded": {}, "shards_done": {}, "genotyped": {}, "downloaded": {}}


def save_manifest(manifest: dict):
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2))


# ── Args ──────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--modal-token-id",     default=None)
    p.add_argument("--modal-token-secret", default=None)
    p.add_argument("--ngc-api-key",        default=None)
    p.add_argument("--bam",  required=True, type=Path)
    p.add_argument("--ref",  required=True, type=Path,
                   help="Local reference (used for fai/dict only; Modal downloads its own copy)")
    p.add_argument("--out",  required=True, type=Path)
    p.add_argument("--chroms", nargs="+", default=CHROMS,
                   help="Chromosomes to process (default: chr1-22,X,Y). "
                        "Use 'auto' to extract from BAM header.")
    p.add_argument("--skip-split",    action="store_true",
                   help="Skip BAM splitting if chr BAMs already exist in /tmp/chr_bams/")
    p.add_argument("--skip-upload",   action="store_true",
                   help="Skip upload if chr BAMs already in Modal volume")
    p.add_argument("--skip-shards",   action="store_true",
                   help="Skip HaplotypeCaller GPU shards if GVCFs already in Modal volume")
    p.add_argument("--skip-modal",    action="store_true",
                   help="Skip all Modal execution (steps 3-4); go straight to local merge")
    p.add_argument("--skip-download", action="store_true",
                   help="Skip download of per-chrom VCFs (assume already local)")
    p.add_argument("--skip-merge",    action="store_true",
                   help="Skip final MergeVcfs step")
    return p.parse_args()


# ── Helpers ───────────────────────────────────────────────────────────────────

def strip_vcf_suffix(name: str) -> str:
    """Strip .g.vcf.gz or .vcf.gz from a filename."""
    for suf in (".g.vcf.gz", ".vcf.gz"):
        if name.endswith(suf):
            return name[: -len(suf)]
    return name


def chroms_from_bam(bam: Path) -> list:
    """Extract chromosome names from BAM header. Filters to canonical chr1-22,X,Y,M."""
    result = subprocess.run(
        [SAMTOOLS, "view", "-H", str(bam)],
        capture_output=True, text=True, check=True,
    )
    chroms = []
    canonical = set(f"chr{i}" for i in range(1, 23)) | {"chrX", "chrY", "chrM"}
    for line in result.stdout.splitlines():
        if line.startswith("@SQ"):
            for field in line.split("\t"):
                if field.startswith("SN:"):
                    name = field[3:]
                    if name in canonical:
                        chroms.append(name)
    return chroms


# ── BAM splitting ─────────────────────────────────────────────────────────────

def split_bam_by_chrom(bam: Path, chroms: list, out_dir: Path, manifest: dict) -> dict:
    """Split BAM into per-chromosome BAMs. Returns {chrom: Path}."""
    out_dir.mkdir(parents=True, exist_ok=True)
    chr_bams = {}
    for chrom in chroms:
        out_bam = out_dir / f"{bam.stem}_{chrom}.bam"
        out_bai = Path(str(out_bam) + ".bai")
        if out_bam.exists() and out_bai.exists():
            print(f"  {chrom}: already split, skipping")
            chr_bams[chrom] = out_bam
            manifest["split"][chrom] = str(out_bam)
            continue
        print(f"  Splitting {chrom}...", end=" ", flush=True)
        subprocess.run(
            [SAMTOOLS, "view", "-b", "-o", str(out_bam), str(bam), chrom],
            check=True,
        )
        subprocess.run([SAMTOOLS, "index", str(out_bam)], check=True)
        size_mb = out_bam.stat().st_size / 1e6
        print(f"{size_mb:.0f} MB")
        chr_bams[chrom] = out_bam
        manifest["split"][chrom] = str(out_bam)
        save_manifest(manifest)
    return chr_bams


# ── Volume inspection ─────────────────────────────────────────────────────────

def list_volume_dir(remote_dir: str) -> set:
    """Return the set of filenames in a volume directory (exact match, not substring)."""
    result = subprocess.run(
        [MODAL_BIN, "volume", "ls", VOLUME_NAME, remote_dir],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        return set()
    names = set()
    for line in result.stdout.splitlines():
        # Modal's table output: strip whitespace and box-drawing chars, take first column
        parts = [p.strip() for p in line.replace("│", "|").split("|") if p.strip()]
        if parts:
            candidate = parts[0]
            # Skip table headers and decorative rows
            if candidate and not candidate.startswith(("-", "━", "═", "Filename", "┏", "┗", "┃")):
                names.add(Path(candidate).name)
    return names


def file_exists_in_volume(remote_path: str, cache: dict = None) -> bool:
    """Check exact filename match in a volume directory, with optional caching."""
    remote_dir  = str(Path(remote_path).parent)
    remote_name = Path(remote_path).name
    if cache is not None:
        if remote_dir not in cache:
            cache[remote_dir] = list_volume_dir(remote_dir)
        return remote_name in cache[remote_dir]
    return remote_name in list_volume_dir(remote_dir)


# ── Upload / Download ─────────────────────────────────────────────────────────

def upload_file(local_path: Path, remote_path: str, vol_cache: dict, retries: int = 8):
    """Upload via modal CLI with exponential backoff retries."""
    size_mb = local_path.stat().st_size / 1e6
    print(f"  ↑ {local_path.name}  ({size_mb:.0f} MB)", end=" ", flush=True)
    if file_exists_in_volume(remote_path, vol_cache):
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
            # Invalidate cache for that dir since we just wrote to it
            vol_cache.pop(str(Path(remote_path).parent), None)
            return
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as e:
            if attempt < retries:
                wait = min(2 ** attempt, 120)
                print(f"\n    Retry {attempt}/{retries} in {wait}s: {type(e).__name__}",
                      end=" ", flush=True)
                time.sleep(wait)
            else:
                raise RuntimeError(f"Upload failed after {retries} attempts: {e}")


def download_file(remote_path: str, local_path: Path, required: bool = True) -> bool:
    """Download a file from volume. Returns True on success."""
    print(f"  ↓ {Path(remote_path).name} → {local_path}")
    result = subprocess.run(
        [MODAL_BIN, "volume", "get", VOLUME_NAME, remote_path, str(local_path)],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        msg = f"download failed for {remote_path}: {result.stderr.strip()}"
        if required:
            raise RuntimeError(msg)
        print(f"    WARN: {msg}")
        return False
    return True


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

    # ── Setup: download reference on Modal once ───────────────────────────────
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
        volume.reload()
        ref_path  = "/data/ref/hg38.fa"
        dict_path = "/data/ref/hg38.dict"
        fai_path  = ref_path + ".fai"
        need_commit = False

        if not os.path.exists(ref_path):
            os.makedirs("/data/ref", exist_ok=True)
            print("Downloading hg38 from UCSC...")
            subprocess.run(
                ["wget", "-q", "-O", "/data/ref/hg38.fa.gz", HG38_URL], check=True,
            )
            subprocess.run(["gzip", "-d", "/data/ref/hg38.fa.gz"], check=True)
            need_commit = True
        else:
            print("Reference FASTA already present")

        if not os.path.exists(fai_path):
            subprocess.run(["samtools", "faidx", ref_path], check=True)
            need_commit = True

        if not os.path.exists(dict_path):
            print("Creating sequence dictionary (GATK CreateSequenceDictionary)...")
            subprocess.run(
                ["gatk", "CreateSequenceDictionary", "-R", ref_path, "-O", dict_path],
                check=True,
            )
            need_commit = True

        if need_commit:
            volume.commit()
        print("Reference ready.")
        return ref_path

    # ── Shard: HaplotypeCaller on one chromosome ──────────────────────────────
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
    def haplotypecaller_shard(chrom: str, bam_stem: str, out_prefix: str) -> dict:
        import subprocess, os
        volume.reload()
        os.makedirs("/data/shards", exist_ok=True)
        bam_path = f"/data/bam/{bam_stem}_{chrom}.bam"
        ref_path = "/data/ref/hg38.fa"
        out_gvcf = f"/data/shards/{out_prefix}_{chrom}.g.vcf.gz"

        if not os.path.exists(bam_path):
            raise FileNotFoundError(f"[{chrom}] BAM missing in volume: {bam_path}")

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
        return {"chrom": chrom, "gvcf": out_gvcf}

    # ── Genotype: one chromosome ──────────────────────────────────────────────
    @app.function(
        image=image,
        cpu=4,
        memory=16 * 1024,
        volumes={"/data": volume},
        timeout=3600,
        retries=1,
        serialized=True,
    )
    def genotype_shard(chrom: str, gvcf_path: str, out_prefix: str) -> dict:
        import subprocess, os
        volume.reload()
        os.makedirs("/data/output", exist_ok=True)
        out_vcf = f"/data/output/{out_prefix}_{chrom}.vcf.gz"

        if not os.path.exists(gvcf_path):
            raise FileNotFoundError(f"[{chrom}] GVCF missing in volume: {gvcf_path}")

        cmd = [
            "gatk", "GenotypeGVCFs",
            "-R", "/data/ref/hg38.fa",
            "-V", gvcf_path,
            "-O", out_vcf,
        ]
        print(f"[{chrom}] GenotypeGVCFs started")
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"[{chrom}] GenotypeGVCFs failed:\n{result.stderr[-4000:]}")
        print(f"[{chrom}] done → {out_vcf}")
        volume.commit()
        return {"chrom": chrom, "vcf": out_vcf}

    out_prefix = strip_vcf_suffix(out_name)

    with app.run():
        print("  Setting up reference genome on Modal...")
        setup_reference.remote()

        if not skip_shards:
            print(f"  Launching {len(chroms)} parallel Parabricks shards...")
            shard_results = list(
                haplotypecaller_shard.starmap(
                    [(chrom, bam_stem, out_prefix) for chrom in chroms]
                )
            )
            done_chroms = {r["chrom"] for r in shard_results}
            missing = set(chroms) - done_chroms
            if missing:
                raise RuntimeError(f"HaplotypeCaller failed for chromosomes: {missing}")
            gvcf_paths = [
                next(r["gvcf"] for r in shard_results if r["chrom"] == c)
                for c in chroms
            ]
        else:
            print(f"  Skipping shards (--skip-shards); verifying GVCFs in volume...")
            existing = list_volume_dir("/shards")
            gvcf_paths = []
            for chrom in chroms:
                name = f"{out_prefix}_{chrom}.g.vcf.gz"
                if name not in existing:
                    raise RuntimeError(
                        f"--skip-shards set but GVCF missing in volume: /shards/{name}"
                    )
                gvcf_paths.append(f"/data/shards/{name}")

        print(f"\n  Launching {len(chroms)} parallel GenotypeGVCFs jobs...")
        geno_results = list(
            genotype_shard.starmap(
                [(chrom, gvcf, out_prefix) for chrom, gvcf in zip(chroms, gvcf_paths)]
            )
        )
        done_chroms = {r["chrom"] for r in geno_results}
        missing = set(chroms) - done_chroms
        if missing:
            raise RuntimeError(f"GenotypeGVCFs failed for chromosomes: {missing}")
        vcf_paths = [
            next(r["vcf"] for r in geno_results if r["chrom"] == c)
            for c in chroms
        ]

    return vcf_paths


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()

    # Credential resolution: CLI > env var > ~/.modal_gatk.env > hardcoded default
    env_file = load_env_file(Path.home() / ".modal_gatk.env")
    modal_id     = get_credential("MODAL_TOKEN_ID",     args.modal_token_id,     env_file, DEFAULT_MODAL_TOKEN_ID)
    modal_secret = get_credential("MODAL_TOKEN_SECRET", args.modal_token_secret, env_file, DEFAULT_MODAL_TOKEN_SECRET)
    ngc_key      = get_credential("NGC_API_KEY",        args.ngc_api_key,        env_file, DEFAULT_NGC_API_KEY)

    os.environ["MODAL_TOKEN_ID"]     = modal_id
    os.environ["MODAL_TOKEN_SECRET"] = modal_secret

    try:
        import modal  # noqa: F401
    except ImportError:
        subprocess.run([sys.executable, "-m", "pip", "install", "modal", "-q"], check=True)

    bam = args.bam.expanduser().resolve()
    bai = Path(str(bam) + ".bai")
    if not bai.exists():
        bai = bam.with_suffix(".bai")
    out = args.out.expanduser().resolve()

    for f in [bam, bai]:
        if not f.exists():
            sys.exit(f"ERROR: File not found: {f}")

    # Chromosome auto-detection
    if args.chroms == ["auto"] or (len(args.chroms) == 1 and args.chroms[0] == "auto"):
        print("Detecting chromosomes from BAM header...")
        chroms = chroms_from_bam(bam)
        print(f"  Found: {chroms}")
    else:
        chroms = args.chroms

    manifest = load_manifest()

    print(f"""
╔══════════════════════════════════════════════════════╗
║    Modal GATK HaplotypeCaller  (Parabricks A100)    ║
╠══════════════════════════════════════════════════════╣
  BAM:         {bam.name}  ({bam.stat().st_size / 1e9:.1f} GB)
  Output:      {out}
  Chromosomes: {len(chroms)}
  Manifest:    {MANIFEST_PATH}
╚══════════════════════════════════════════════════════╝
""")

    # ── Step 1: Split BAM by chromosome ──────────────────────────────────────
    chr_bam_dir = Path("/tmp/chr_bams")
    if not args.skip_split:
        print("Step 1/5 — Splitting BAM by chromosome...")
        chr_bams = split_bam_by_chrom(bam, chroms, chr_bam_dir, manifest)
    else:
        print("Step 1/5 — Skipping split (--skip-split)")
        chr_bams = {c: chr_bam_dir / f"{bam.stem}_{c}.bam" for c in chroms}
        for c, p in chr_bams.items():
            if not p.exists():
                sys.exit(f"ERROR: --skip-split set but {p} does not exist")

    # ── Step 2: Upload per-chrom BAMs ────────────────────────────────────────
    if not args.skip_upload and not args.skip_modal:
        print("\nStep 2/5 — Uploading per-chromosome BAMs to Modal volume...")
        print(f"  Volume: {VOLUME_NAME}")
        subprocess.run(
            [MODAL_BIN, "volume", "create", VOLUME_NAME],
            capture_output=True,
        )
        vol_cache = {}
        for chrom, chr_bam in chr_bams.items():
            bai_path = Path(str(chr_bam) + ".bai")
            upload_file(chr_bam, f"/bam/{chr_bam.name}", vol_cache)
            upload_file(bai_path, f"/bam/{chr_bam.name}.bai", vol_cache)
            manifest["uploaded"][chrom] = True
            save_manifest(manifest)
        print("All uploads complete.\n")
    else:
        reason = "--skip-upload" if args.skip_upload else "--skip-modal"
        print(f"Step 2/5 — Skipping upload ({reason})\n")

    # ── Step 3: Run on Modal ─────────────────────────────────────────────────
    out_prefix = strip_vcf_suffix(out.name)
    if not args.skip_modal:
        label = "Parabricks HaplotypeCaller + GenotypeGVCFs" if not args.skip_shards else "GenotypeGVCFs only"
        print(f"Step 3/5 — Running {label} on Modal...")
        remote_vcf_paths = run_on_modal(
            chroms      = chroms,
            bam_stem    = bam.stem,
            out_name    = out.name,
            ngc_key     = ngc_key,
            skip_shards = args.skip_shards,
        )
        for chrom in chroms:
            manifest["genotyped"][chrom] = True
        save_manifest(manifest)
    else:
        print("Step 3/5 — Skipping Modal run (--skip-modal)")
        remote_vcf_paths = [f"/data/output/{out_prefix}_{c}.vcf.gz" for c in chroms]

    # ── Step 4: Download per-chrom VCFs ──────────────────────────────────────
    vcf_dir = out.parent / "chr_vcfs"
    vcf_dir.mkdir(exist_ok=True)
    local_vcfs = []
    if not args.skip_download and not args.skip_modal:
        print(f"\nStep 4/5 — Downloading {len(remote_vcf_paths)} per-chromosome VCFs...")
        for chrom, remote_path in zip(chroms, remote_vcf_paths):
            local_vcf = vcf_dir / f"{out_prefix}_{chrom}.vcf.gz"
            download_file(remote_path, local_vcf, required=True)
            tbi_remote = remote_path + ".tbi"
            tbi_local  = Path(str(local_vcf) + ".tbi")
            if not download_file(tbi_remote, tbi_local, required=False):
                print(f"    WARN: no .tbi for {chrom}; MergeVcfs may still work")
            local_vcfs.append(local_vcf)
            manifest["downloaded"][chrom] = str(local_vcf)
            save_manifest(manifest)
    else:
        reason = "--skip-download" if args.skip_download else "--skip-modal"
        print(f"Step 4/5 — Skipping download ({reason})")
        local_vcfs = [vcf_dir / f"{out_prefix}_{c}.vcf.gz" for c in chroms]
        for v in local_vcfs:
            if not v.exists():
                sys.exit(f"ERROR: expected local VCF missing: {v}")

    # ── Step 5: Merge per-chrom VCFs ─────────────────────────────────────────
    if args.skip_merge:
        print("Step 5/5 — Skipping merge (--skip-merge)")
        return

    print(f"\nStep 5/5 — Merging {len(local_vcfs)} VCFs with GATK MergeVcfs...")
    gatk = find_binary("gatk", ["~/miniconda3/bin/gatk", "/usr/local/bin/gatk"])
    merge_inputs = []
    for vcf in local_vcfs:
        merge_inputs += ["-I", str(vcf)]
    merged = out.parent / (strip_vcf_suffix(out.name) + "_genotyped.vcf.gz")
    result = subprocess.run(
        [gatk, "MergeVcfs", *merge_inputs, "-O", str(merged)],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        print(f"MergeVcfs stderr:\n{result.stderr[-3000:]}")
        sys.exit("ERROR: MergeVcfs failed.")

    if merged.exists():
        print(f"""
╔══════════════════════════════════════════════════════╗
║  Done!                                              ║
╠══════════════════════════════════════════════════════╣
  Genotyped VCF: {merged}
╚══════════════════════════════════════════════════════╝
""")
    else:
        sys.exit("ERROR: Merged VCF not found after MergeVcfs.")


if __name__ == "__main__":
    main()
