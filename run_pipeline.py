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
from markers import MARKERS, OL_LINEAGE
import train_ol_classifier
import tap_analysis


ROOT = Path(__file__).parent
OUT = ROOT / "outputs"
RAW = ROOT / "data" / "raw"

def build_processed() -> None:
    adata = data_io.load_all(RAW)
    adata = preprocess.qc_filter(adata)
    print("completed quality control")
    adata.to_df().head(100).to_csv(f"{OUT}/raw_barcodes_genes_top100.csv")
    adata = preprocess.normalize_and_embed(adata)
    print("completed normalizing and embeding")
    adata.to_df().head(100).to_csv(f"{OUT}/normalized_barcodes_genes_top100.csv")

    for name, genes in MARKERS.items():
        present = [g for g in genes if g in adata.var_names]
        if not present:
            continue
        sc.tl.score_genes(
            adata, gene_list=present, score_name=f"score_{name}",
            ctrl_size=min(50, max(10, len(present) * 5)), random_state=0,
        )

    score_cols = [f"score_{n}" for n in MARKERS if f"score_{n}" in adata.obs.columns]
    means = adata.obs.groupby("leiden", observed=True)[score_cols].mean()
    cluster_to_type = {c: row.idxmax().replace("score_", "") for c, row in means.iterrows()}
    adata.obs["cell_type"] = adata.obs["leiden"].map(cluster_to_type).astype("category")
    adata.uns["cluster_to_type"] = cluster_to_type
    adata.obs["fate_OL_lineage"] = adata.obs["cell_type"].isin(OL_LINEAGE).astype(int)

    OUT.mkdir(exist_ok=True)
    h5 = OUT / "processed.h5ad"
    adata.write_h5ad(h5)
    print(f"  wrote {h5}: {adata.n_obs} cells × {adata.n_vars} genes")
    print(adata.obs["cell_type"].value_counts().to_string())

    adata.obs.to_csv(f"{OUT}/barcode_to_cell_type_mapping.csv")




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
    print("running classifier")
    train_ol_classifier.run_classifier()
    print("running corroboration analysis")
    tap_analysis.run_analysis()
    print(f"\n=== ALL DONE in {time.time() - overall:.0f}s ===")


if __name__ == "__main__":
    main()
