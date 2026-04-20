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
  Modal:  read from `~/.modal.toml` (run `modal token new` once to populate).
  NGC:    read from env var `NGC_API_KEY`, or `~/.modal_gatk.env`, or `--ngc-api-key`.

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
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

# ── Config ────────────────────────────────────────────────────────────────────
CHROMS      = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
VOLUME_NAME = "gatk-haplotypecaller-data"
GATK_VER    = "4.6.2.0"

# hg38 reference — downloaded on Modal directly from UCSC (avoids uploading 3 GB)
HG38_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"

# Manifest for resumability
MANIFEST_PATH = Path.home() / ".modal_gatk_manifest.json"


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
    p.add_argument("--ngc-api-key", default=None,
                   help="NGC API key (overrides env var / ~/.modal_gatk.env)")
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

def split_bam_by_chrom(bam: Path, chroms: list, out_dir: Path, manifest: dict,
                       max_workers: int = 6) -> dict:
    """Split BAM into per-chromosome BAMs in parallel. Returns {chrom: Path}."""
    out_dir.mkdir(parents=True, exist_ok=True)
    chr_bams = {}

    def split_one(chrom: str):
        out_bam = out_dir / f"{bam.stem}_{chrom}.bam"
        out_bai = Path(str(out_bam) + ".bai")
        if out_bam.exists() and out_bai.exists():
            return chrom, out_bam, None
        subprocess.run(
            [SAMTOOLS, "view", "-b", "-o", str(out_bam), str(bam), chrom],
            check=True, capture_output=True,
        )
        subprocess.run([SAMTOOLS, "index", str(out_bam)],
                       check=True, capture_output=True)
        return chrom, out_bam, out_bam.stat().st_size / 1e6

    with ThreadPoolExecutor(max_workers=max_workers) as pool:
        futures = {pool.submit(split_one, c): c for c in chroms}
        for fut in as_completed(futures):
            chrom, out_bam, size_mb = fut.result()
            if size_mb is None:
                print(f"  {chrom}: already split, skipping")
            else:
                print(f"  {chrom}: {size_mb:.0f} MB")
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

def _upload_once(local_path: Path, remote_path: str, retries: int = 8):
    """Upload a single file via modal CLI with exponential backoff. No stdout prints
    so callers can parallelize and format output themselves."""
    for attempt in range(1, retries + 1):
        try:
            subprocess.run(
                [MODAL_BIN, "volume", "put", "--force", VOLUME_NAME,
                 str(local_path), remote_path],
                check=True, timeout=900, capture_output=True,
            )
            return
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as e:
            if attempt < retries:
                time.sleep(min(2 ** attempt, 120))
            else:
                raise RuntimeError(f"Upload failed after {retries} attempts: {e}")


def upload_files_parallel(pairs: list, max_workers: int = 4):
    """Upload (local_path, remote_path) pairs in parallel, skipping files
    already present in the volume."""
    # Prime per-directory caches once (single modal-CLI call per dir).
    dirs = {str(Path(rp).parent) for _, rp in pairs}
    vol_cache = {d: list_volume_dir(d) for d in dirs}

    to_upload = []
    for local, remote in pairs:
        if Path(remote).name in vol_cache.get(str(Path(remote).parent), set()):
            size_mb = local.stat().st_size / 1e6
            print(f"  = {local.name} ({size_mb:.0f} MB) already in volume")
        else:
            to_upload.append((local, remote))

    if not to_upload:
        return

    print(f"  Uploading {len(to_upload)} file(s) with {max_workers} workers...")
    with ThreadPoolExecutor(max_workers=max_workers) as pool:
        futures = {pool.submit(_upload_once, lp, rp): (lp, rp) for lp, rp in to_upload}
        for fut in as_completed(futures):
            lp, rp = futures[fut]
            fut.result()  # raises if upload ultimately failed
            size_mb = lp.stat().st_size / 1e6
            print(f"  ↑ {lp.name} ({size_mb:.0f} MB) ✓")


def vol_path(container_path: str) -> str:
    """Translate an in-container path (/data/foo) to a volume-root CLI path (/foo)."""
    if container_path.startswith("/data/"):
        return container_path[len("/data"):]
    return container_path


def download_file(remote_path: str, local_path: Path, required: bool = True) -> bool:
    """Download a file from volume. Returns True on success."""
    cli_path = vol_path(remote_path)
    print(f"  ↓ {Path(cli_path).name} → {local_path}")
    result = subprocess.run(
        [MODAL_BIN, "volume", "get", VOLUME_NAME, cli_path, str(local_path)],
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

    # ── Shard: HaplotypeCaller + GenotypeGVCFs on one chromosome ──────────────
    # Fused so GenotypeGVCFs runs in the same container as HC — no extra
    # cold-start, no volume commit/reload round-trip between the two stages.
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
    def hc_genotype_shard(chrom: str, bam_stem: str, out_prefix: str) -> dict:
        import subprocess, os
        volume.reload()
        os.makedirs("/data/shards", exist_ok=True)
        os.makedirs("/data/output", exist_ok=True)
        bam_path  = f"/data/bam/{bam_stem}_{chrom}.bam"
        ref_path  = "/data/ref/hg38.fa"
        gvcf_path = f"/data/shards/{out_prefix}_{chrom}.g.vcf.gz"
        vcf_path  = f"/data/output/{out_prefix}_{chrom}.vcf.gz"

        if not os.path.exists(bam_path):
            raise FileNotFoundError(f"[{chrom}] BAM missing in volume: {bam_path}")

        # HaplotypeCaller (GPU). Skip if GVCF already present (resume path).
        if not os.path.exists(gvcf_path):
            print(f"[{chrom}] pbrun haplotypecaller started")
            result = subprocess.run(
                ["pbrun", "haplotypecaller",
                 "--ref",          ref_path,
                 "--in-bam",       bam_path,
                 "--out-variants", gvcf_path,
                 "--gvcf",
                 "--num-gpus",     "1"],
                capture_output=True, text=True,
            )
            if result.returncode != 0:
                raise RuntimeError(f"[{chrom}] haplotypecaller failed:\n{result.stderr[-4000:]}")
            print(f"[{chrom}] HC done → {gvcf_path}")
        else:
            print(f"[{chrom}] GVCF already exists, skipping HC")

        # GenotypeGVCFs (CPU) — runs in the same container, no extra cold-start.
        print(f"[{chrom}] GenotypeGVCFs started")
        result = subprocess.run(
            ["gatk", "GenotypeGVCFs",
             "-R", ref_path,
             "-V", gvcf_path,
             "-O", vcf_path],
            capture_output=True, text=True,
        )
        if result.returncode != 0:
            raise RuntimeError(f"[{chrom}] GenotypeGVCFs failed:\n{result.stderr[-4000:]}")
        print(f"[{chrom}] geno done → {vcf_path}")
        volume.commit()
        return {"chrom": chrom, "gvcf": gvcf_path, "vcf": vcf_path}

    out_prefix = strip_vcf_suffix(out_name)

    with app.run():
        print("  Setting up reference genome on Modal...")
        setup_reference.remote()

        # skip_shards: HC outputs already present in volume; fused function
        # detects existing GVCF and runs GenotypeGVCFs only (still on A100 —
        # speed parity, slight cost trade).
        label = "HC+GenotypeGVCFs" if not skip_shards else "GenotypeGVCFs (HC skipped)"
        print(f"  Launching {len(chroms)} parallel {label} shards...")

        if skip_shards:
            existing = list_volume_dir("/shards")
            missing_gvcfs = [c for c in chroms
                             if f"{out_prefix}_{c}.g.vcf.gz" not in existing]
            if missing_gvcfs:
                raise RuntimeError(
                    f"--skip-shards set but GVCFs missing in volume for: {missing_gvcfs}"
                )

        shard_results = list(
            hc_genotype_shard.starmap(
                [(chrom, bam_stem, out_prefix) for chrom in chroms]
            )
        )
        done_chroms = {r["chrom"] for r in shard_results}
        missing = set(chroms) - done_chroms
        if missing:
            raise RuntimeError(f"Shard failed for chromosomes: {missing}")
        vcf_paths = [
            next(r["vcf"] for r in shard_results if r["chrom"] == c)
            for c in chroms
        ]

    return vcf_paths


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()

    # Modal credentials are auto-loaded by the SDK from ~/.modal.toml
    # (run `modal token new` once to populate). NGC key is still needed at
    # runtime to pull the Parabricks container from nvcr.io.
    env_file = load_env_file(Path.home() / ".modal_gatk.env")
    ngc_key = get_credential("NGC_API_KEY", args.ngc_api_key, env_file)

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
        upload_pairs = []
        for chrom, chr_bam in chr_bams.items():
            bai_path = Path(str(chr_bam) + ".bai")
            upload_pairs.append((chr_bam, f"/bam/{chr_bam.name}"))
            upload_pairs.append((bai_path, f"/bam/{chr_bam.name}.bai"))
        upload_files_parallel(upload_pairs, max_workers=4)
        for chrom in chr_bams:
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
        remote_vcf_paths = [f"/data/output/{out_prefix}_{c}.vcf.gz" for c in chroms]  # translated by vol_path()

    # ── Step 4: Download per-chrom VCFs ──────────────────────────────────────
    vcf_dir = out.parent / "chr_vcfs"
    vcf_dir.mkdir(exist_ok=True)
    local_vcfs = []
    if not args.skip_download:
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
        print("Step 4/5 — Skipping download (--skip-download)")
        local_vcfs = [vcf_dir / f"{out_prefix}_{c}.vcf.gz" for c in chroms]
        for v in local_vcfs:
            if not v.exists():
                sys.exit(f"ERROR: expected local VCF missing: {v}")

    # ── Step 5: Merge per-chrom VCFs ─────────────────────────────────────────
    if args.skip_merge:
        print("Step 5/5 — Skipping merge (--skip-merge)")
        return

    print(f"\nStep 5/5 — Merging {len(local_vcfs)} VCFs with GATK MergeVcfs...")
    java = find_binary("java", ["~/miniconda3/bin/java", "/usr/bin/java"])
    # Locate GATK jar — prefer conda install, fall back to any jar on PATH
    gatk_jar_candidates = list(Path(os.path.expanduser("~/miniconda3/share")).glob("gatk4-*/gatk-package-*-local.jar"))
    if not gatk_jar_candidates:
        sys.exit("ERROR: could not find GATK jar under ~/miniconda3/share/gatk4-*/")
    gatk_jar = str(sorted(gatk_jar_candidates)[-1])
    merge_inputs = []
    for vcf in local_vcfs:
        merge_inputs += ["-I", str(vcf)]
    merged = out.parent / (strip_vcf_suffix(out.name) + "_genotyped.vcf.gz")
    result = subprocess.run(
        [java, "-jar", gatk_jar, "MergeVcfs", *merge_inputs, "-O", str(merged)],
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
