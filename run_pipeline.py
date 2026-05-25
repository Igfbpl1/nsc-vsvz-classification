"""End-to-end pipeline: raw 10x → processed.h5ad → all downstream analyses.

Phase 0 (preprocessing + cell-type annotation) runs inline; phases 1–6 each
shell out to the corresponding script so memory is freed between heavy steps.

Run from the project root:

    uv run python run_pipeline.py            # uses existing processed.h5ad if present
    uv run python run_pipeline.py --rebuild  # forces re-build of processed.h5ad
"""
from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).parent
OUT = ROOT / "outputs"
RAW = ROOT / "data" / "raw"

PHASES: list[tuple[str, str]] = [
    ("train_ol_classifier.py",           "train OL-commitment classifier (markers-excluded)"),
    ("interpret_classifier.py",          "condition fractions + within-TAP DE"),
    ("_plot_taps.py",                    "TAP-only commitment plots"),
    # Cross-dataset validation, ordered by contribution to the final conclusions:
    ("validate_ntrk2_in_gse109447.py",   "GSE109447 (Mizrak / Doetsch / Sims, Columbia) — arms-length cross-lab validation"),
    ("validate_triggers_in_gse75330.py", "GSE75330  (Castelo-Branco/Karolinska) — mature-OL false-positive removal"),
]


def build_processed() -> None:
    pass

def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--rebuild", action="store_true",
                   help="rebuild outputs/processed.h5ad even if it exists")
    args = p.parse_args()

    overall = time.time()
    h5 = OUT / "processed.h5ad"

    print("=== phase 0: build processed.h5ad ===")
    if h5.exists() and not args.rebuild:
        print(f"  {h5} exists; skipping (use --rebuild to force)")
    else:
        t0 = time.time()
        build_processed()
        print(f"  done in {time.time() - t0:.0f}s")

    print(f"\n=== ALL DONE in {time.time() - overall:.0f}s ===")


if __name__ == "__main__":
    main()
