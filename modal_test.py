#!/usr/bin/env python3
"""Quick connectivity test for Modal + NGC Parabricks.

Credentials:
  Modal: auto-loaded from ~/.modal.toml (`modal token new` once)
  NGC:   env var NGC_API_KEY, or ~/.modal_gatk.env
"""
import os, sys
from pathlib import Path

# Load NGC key from ~/.modal_gatk.env if not already in env
NGC_KEY = os.environ.get("NGC_API_KEY")
if not NGC_KEY:
    env_file = Path.home() / ".modal_gatk.env"
    if env_file.exists():
        for line in env_file.read_text().splitlines():
            if line.startswith("NGC_API_KEY="):
                NGC_KEY = line.split("=", 1)[1].strip().strip('"').strip("'")
                break
if not NGC_KEY:
    sys.exit("ERROR: NGC_API_KEY not set (env var or ~/.modal_gatk.env)")

import modal

app = modal.App("gatk-connectivity-test")

# Test 1: basic Modal function (CPU)
basic_image = modal.Image.debian_slim(python_version="3.11")

@app.function(image=basic_image, cpu=1, timeout=60)
def test_modal():
    import platform, subprocess
    arch    = platform.machine()
    cpu_out = subprocess.run(["grep", "-c", "processor", "/proc/cpuinfo"],
                             capture_output=True, text=True).stdout.strip()
    avx     = "yes" if "avx" in open("/proc/cpuinfo").read() else "no"
    return {"arch": arch, "cpus": cpu_out, "avx": avx}

# Test 2: Parabricks NGC image
ngc_secret = modal.Secret.from_dict({
    "REGISTRY_USERNAME": "$oauthtoken",
    "REGISTRY_PASSWORD": NGC_KEY,
})
pbricks_image = (
    modal.Image.from_registry(
        "nvcr.io/nvidia/clara/clara-parabricks:4.3.1-1",
        secret=ngc_secret,
        add_python="3.11",
    )
)

@app.function(
    image=pbricks_image,
    gpu="A100-80GB",
    timeout=120,
)
def test_parabricks():
    import subprocess
    result = subprocess.run(["pbrun", "version"], capture_output=True, text=True)
    gpu_result = subprocess.run(["nvidia-smi", "--query-gpu=name,memory.total",
                                  "--format=csv,noheader"],
                                 capture_output=True, text=True)
    return {
        "pbrun_version": result.stdout.strip() or result.stderr.strip(),
        "gpu":           gpu_result.stdout.strip(),
    }

@app.local_entrypoint()
def main():
    print("\n── Test 1: Modal connectivity ──────────────────────")
    info = test_modal.remote()
    print(f"  ✓ Connected to Modal")
    print(f"  Architecture : {info['arch']}")
    print(f"  CPU cores    : {info['cpus']}")
    print(f"  AVX support  : {info['avx']}")

    print("\n── Test 2: NGC Parabricks image ────────────────────")
    try:
        pb = test_parabricks.remote()
        print(f"  ✓ Parabricks image pulled from NGC")
        print(f"  pbrun version: {pb['pbrun_version']}")
        print(f"  GPU          : {pb['gpu']}")
    except Exception as e:
        print(f"  ✗ Parabricks test failed: {e}")
        sys.exit(1)

    print("\n✓ All tests passed — ready to run modal_gatk.py\n")
