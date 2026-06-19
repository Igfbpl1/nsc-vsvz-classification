"""End-to-end pipeline: raw 10x → processed.h5ad → all downstream analyses.

Phase 0 (preprocessing + cell-type annotation) runs inline; phases 1–6 each
shell out to the corresponding script so memory is freed between heavy steps.

Run from the project root:

    uv run python run_pipeline.py            # uses existing processed.h5ad if present
    uv run python run_pipeline.py --rebuild  # forces re-build of processed.h5ad
"""

from __future__ import annotations

import argparse
import time
import data_io
from pathlib import Path
import scanpy as sc
import preprocess
from compare_tap_fate_methods import run_compare_fate_methods
from markers import MARKERS, OL_LINEAGE
import train_ol_classifier
import tap_analysis
import velocity_build
from rna_velocity_pipeline import run_rna_velocity_pipeline

ROOT = Path(__file__).parent
OUT = ROOT / "outputs"
RAW = ROOT / "data" / "raw"


def build_processed() -> None:
    adata = data_io.load_all(RAW)
    adata = preprocess.qc_filter(adata, output_dir=OUT)
    print("completed quality control")
    adata.to_df().head(100).to_csv(f"{OUT}/raw_barcodes_genes_top100.csv")
    adata = preprocess.normalize_and_embed(adata, output_dir=OUT)
    print("completed normalizing and embeding")
    adata.to_df().head(100).to_csv(f"{OUT}/normalized_barcodes_genes_top100.csv")

    for name, genes in MARKERS.items():
        present = [g for g in genes if g in adata.var_names]
        if not present:
            continue
        sc.tl.score_genes(
            adata,
            gene_list=present,
            score_name=f"score_{name}",
            ctrl_size=min(50, max(10, len(present) * 5)),
            random_state=0,
        )

    score_cols = [f"score_{n}" for n in MARKERS if f"score_{n}" in adata.obs.columns]
    means = adata.obs.groupby("leiden", observed=True)[score_cols].mean()
    # Pure idxmax fails when panels have unequal absolute scales (mature
    # OL out-scores COP just because Mbp/Plp1 are highly expressed). Pure
    # z-score over-corrects (Cd52/Cd69 are weakly expressed in many
    # microglia clusters → Other_Immune wrongly wins). Hybrid: use abs
    # idxmax when the margin is clear; use z-score tie-break only between
    # the top two panels when the abs margin is small.
    z = (means - means.mean(axis=0)) / means.std(axis=0)
    MARGIN = 0.4       # max abs gap that counts as a close call
    RUNNER_FLOOR = 0.5 # runner-up must have this much abs signal to be eligible
    cluster_to_type = {}
    for c in means.index:
        row_abs = means.loc[c].sort_values(ascending=False)
        winner, runner_up = row_abs.index[0], row_abs.index[1]
        margin = row_abs.iloc[0] - row_abs.iloc[1]
        runner_abs = row_abs.iloc[1]
        if margin < MARGIN and runner_abs >= RUNNER_FLOOR:
            # close call AND runner-up has real signal — pick more distinctive
            chosen = z.loc[c, [winner, runner_up]].idxmax()
        else:
            chosen = winner
        cluster_to_type[c] = chosen.replace("score_", "")
    adata.obs["cell_type"] = adata.obs["leiden"].map(cluster_to_type).astype("category")
    adata.uns["cluster_to_type"] = cluster_to_type
    adata.obs["fate_OL_lineage"] = adata.obs["cell_type"].isin(OL_LINEAGE).astype(int)
    sc.settings.figdir = str(OUT)
    sc.settings.dpi_save = 300
    sc.pl.umap(adata, color="cell_type", size=40, alpha=1.0,
               save="leiden_cell_type.png", show=False)

    OUT.mkdir(exist_ok=True)
    h5 = OUT / "processed.h5ad"
    adata.write_h5ad(h5)
    print(f"  wrote {h5}: {adata.n_obs} cells × {adata.n_vars} genes")
    print(adata.obs["cell_type"].value_counts().to_string())

    adata.obs.to_csv(f"{OUT}/barcode_to_cell_type_mapping.csv")


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument(
        "--rebuild",
        action="store_true",
        help="rebuild outputs/processed.h5ad and velocity_combined.h5ad even if they exist",
    )
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
    print("running classifier")
    train_ol_classifier.run_classifier()
    print("running corroboration analysis")
    tap_analysis.run_analysis()
    print("=== velocity build: velocity_combined.h5ad ===")
    vel_h5ad = Path(velocity_build.OUT_H5AD)
    if vel_h5ad.exists() and not args.rebuild:
        print(f"  {vel_h5ad} exists; skipping (use --rebuild to force)")
    else:
        t0 = time.time()
        velocity_build.build_velocity()
        print(f"  done in {time.time() - t0:.0f}s")

    print("running rna velocity pipeline")
    run_rna_velocity_pipeline()
    print("comparing fate methods")
    run_compare_fate_methods()
    print(f"\n=== ALL DONE in {time.time() - overall:.0f}s ===")


if __name__ == "__main__":
    main()
