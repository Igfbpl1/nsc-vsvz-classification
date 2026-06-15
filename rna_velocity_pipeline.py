"""
RNA Velocity & Trajectory Pipeline — Phase 2–4
Loads kb-python MTX output (spliced + unspliced) and merges with processed.h5ad.

When GSM8253796 download completes, update the Cntl entry in SAMPLE_CONFIGS:
    "h5ad_sample_id": "GSM8253796",
    "kb_dir": "sra_runs/kb_output_GSM8253796/counts_unfiltered",
"""

import os
import pandas as pd
from scipy.io import mmread
import anndata as ad
import scanpy as sc
import scvelo as scv
import matplotlib
matplotlib.use("Agg")

# ── Configuration ──────────────────────────────────────────────────────────────
H5AD_PATH = "outputs/processed.h5ad"
T2G_PATH = "sra_runs/t2g.txt"
OUTPUT_DIR = "outputs/velocity"

# Both paths point to the same folder until GSM8253796 finishes processing.
# To switch to real data: set h5ad_sample_id="GSM8253796" and update kb_dir.
SAMPLE_CONFIGS = [
    {
        "h5ad_sample_id": "GSM8253798",       # placeholder: swap for GSM8253796 when ready
        "velocity_label": "Cntl_placeholder",
        "kb_dir": "sra_runs/kb_output_GSM8253798/counts_unfiltered",
    },
    {
        "h5ad_sample_id": "GSM8253798",
        "velocity_label": "CupRap",
        "kb_dir": "sra_runs/kb_output_GSM8253798/counts_unfiltered",
    },
]

MIN_SHARED_COUNTS = 20
N_TOP_GENES       = 2000
N_PCS             = 30
N_NEIGHBORS       = 30


# ── Gene name mapping (Ensembl ID → gene name) ────────────────────────────────
def load_ensembl_to_name(t2g_path: str) -> dict[str, str]:
    df = pd.read_csv(
        t2g_path, sep="\t", header=None,
        names=["transcript_id", "gene_id", "gene_name",
               "transcript_name", "chrom", "start", "end", "strand"],
    )
    return df.drop_duplicates("gene_id").set_index("gene_id")["gene_name"].to_dict()


# ── Load one MTX layer ─────────────────────────────────────────────────────────
def load_kb_layer(kb_dir: str, layer: str) -> ad.AnnData:
    """Return AnnData (cells × genes) for one kb counts_unfiltered layer."""
    # kb MTX is barcodes × genes (cells × genes)
    matrix = mmread(f"{kb_dir}/{layer}.mtx").tocsr()
    barcodes = pd.read_csv(f"{kb_dir}/{layer}.barcodes.txt", header=None)[0].tolist()
    genes = pd.read_csv(f"{kb_dir}/{layer}.genes.txt", header=None)[0].tolist()
    return ad.AnnData(
        X=matrix,
        obs=pd.DataFrame(index=barcodes),
        var=pd.DataFrame(index=genes),
    )


# ── Build per-sample velocity AnnData ─────────────────────────────────────────
def build_sample_adata(
    config: dict,
    ref: ad.AnnData,
    ensembl_to_name: dict[str, str],
) -> ad.AnnData:
    h5ad_sid = config["h5ad_sample_id"]
    label = config["velocity_label"]
    kb_dir = config["kb_dir"]

    print(f"\n[{label}] Loading kb layers from {kb_dir} ...")
    sp_adata = load_kb_layer(kb_dir, "spliced")
    un_adata = load_kb_layer(kb_dir, "unspliced")
    print(f"  spliced:   {sp_adata.shape[0]:,} barcodes × {sp_adata.shape[1]:,} genes")
    print(f"  unspliced: {un_adata.shape[0]:,} barcodes × {un_adata.shape[1]:,} genes")

    # Subset ref to this sample; strip "GSM#######_" prefix and "-1" suffix
    ref_sub = ref[ref.obs["sample_id"] == h5ad_sid]
    raw_bc = (
        ref_sub.obs_names
        .str.removeprefix(f"{h5ad_sid}_")
        .str.removesuffix("-1")
    )

    # Intersect with barcodes present in both layers
    in_sp = raw_bc.isin(sp_adata.obs_names)
    in_un = raw_bc.isin(un_adata.obs_names)
    keep_mask = in_sp & in_un
    matched_raw_bc = raw_bc[keep_mask]
    matched_ref = ref_sub[keep_mask]
    print(f"  {h5ad_sid} cells in h5ad: {len(raw_bc):,}")
    print(f"  matched in both layers:    {len(matched_raw_bc):,}")

    if len(matched_raw_bc) == 0:
        raise ValueError(
            f"No barcodes matched for {h5ad_sid}. "
            "Check that barcode format strips correctly (prefix/suffix)."
        )

    # Intersect genes between spliced and unspliced
    shared_genes = sp_adata.var_names[sp_adata.var_names.isin(un_adata.var_names)]
    sp_sub = sp_adata[matched_raw_bc.tolist(), shared_genes]
    un_sub = un_adata[matched_raw_bc.tolist(), shared_genes]

    # Build result: raw barcodes as obs_names, Ensembl IDs as var_names
    result = ad.AnnData(
        X=sp_sub.X.copy(),
        obs=matched_ref.obs.copy(),
        var=pd.DataFrame(
            {"gene_name": [ensembl_to_name.get(g, g) for g in shared_genes]},
            index=shared_genes,
        ),
    )
    result.obs_names = matched_raw_bc.tolist()
    result.layers["spliced"] = sp_sub.X.copy()
    result.layers["unspliced"] = un_sub.X.copy()
    result.obs["velocity_label"] = label

    # Transfer UMAP from ref (for visualization)
    result.obsm["X_umap"] = matched_ref.obsm["X_umap"].copy()

    return result


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    print("Loading processed.h5ad ...")
    ref = sc.read_h5ad(H5AD_PATH)

    cluster_to_type = ref.uns.get("cluster_to_type", {})
    ref.obs["cell_type"] = ref.obs["leiden"].map(cluster_to_type).fillna("Unknown")
    print(f"  {ref.shape[0]:,} cells, {ref.shape[1]:,} genes")
    print(f"  cell types: {sorted(ref.obs['cell_type'].unique())}")

    print("\nLoading Ensembl → gene name mapping ...")
    ensembl_to_name = load_ensembl_to_name(T2G_PATH)
    print(f"  {len(ensembl_to_name):,} entries")

    samples = [build_sample_adata(cfg, ref, ensembl_to_name) for cfg in SAMPLE_CONFIGS]

    print("\nConcatenating samples ...")
    keys = [cfg["velocity_label"] for cfg in SAMPLE_CONFIGS]
    combined = ad.concat(samples, keys=keys, label="concat_key", join="inner")
    print(f"  combined: {combined.shape[0]:,} cells × {combined.shape[1]:,} genes")
    print(f"  velocity_label counts:\n{combined.obs['velocity_label'].value_counts()}")

    # ── Preprocessing ──────────────────────────────────────────────────────────
    scv.settings.verbosity = 3
    scv.settings.set_figure_params("scvelo", figsize=(8, 6))

    print("\n[scVelo] Filtering genes ...")
    scv.pp.filter_genes(combined, min_shared_counts=MIN_SHARED_COUNTS)

    print("[scVelo] Normalizing ...")
    sc.pp.normalize_total(combined, target_sum=1e4)
    sc.pp.log1p(combined)

    print(f"[scVelo] Selecting top {N_TOP_GENES} highly variable genes ...")
    sc.pp.highly_variable_genes(combined, n_top_genes=N_TOP_GENES, subset=True)
    print(f"  genes after HVG subset: {combined.n_vars:,}")

    print("[scVelo] Computing PCA and neighbors ...")
    sc.pp.pca(combined, n_comps=N_PCS)
    sc.pp.neighbors(combined, n_pcs=N_PCS, n_neighbors=N_NEIGHBORS)

    print("[scVelo] Computing moments ...")
    scv.pp.moments(combined, n_pcs=None, n_neighbors=None)

    # ── Velocity ───────────────────────────────────────────────────────────────
    print("[scVelo] Recovering dynamics (dynamical model — ~30 min) ...")
    scv.tl.recover_dynamics(combined, n_jobs=4, show_progress_bar=False)

    print("[scVelo] Computing velocity ...")
    scv.tl.velocity(combined, mode="dynamical")

    print("[scVelo] Building velocity graph ...")
    scv.tl.velocity_graph(combined, show_progress_bar=False)

    print("[scVelo] Computing latent time ...")
    scv.tl.latent_time(combined)

    # ── Save ───────────────────────────────────────────────────────────────────
    out_h5ad = f"{OUTPUT_DIR}/velocity_combined.h5ad"
    combined.write_h5ad(out_h5ad)
    print(f"\nSaved → {out_h5ad}")

    # ── Plots ──────────────────────────────────────────────────────────────────
    scv.settings.figdir = OUTPUT_DIR
    color_key = "cell_type"

    print("\nSaving plots ...")

    scv.pl.velocity_embedding_stream(
        combined, basis="umap", color=color_key,
        title="Velocity stream (all cells)",
        save="stream_all.png", show=False,
    )

    scv.pl.scatter(
        combined, color="latent_time", color_map="gnuplot",
        title="Latent time",
        save="latent_time.png", show=False,
    )

    for label in combined.obs["velocity_label"].unique():
        sub = combined[combined.obs["velocity_label"] == label]
        scv.pl.velocity_embedding_stream(
            sub, basis="umap", color=color_key,
            title=f"Velocity — {label}",
            save=f"stream_{label}.png", show=False,
        )

    print("\n[scVelo] Ranking velocity genes ...")
    scv.tl.rank_velocity_genes(combined, groupby=color_key, min_corr=0.3)
    driver_df = pd.DataFrame(combined.uns["rank_velocity_genes"]["names"])
    driver_df.to_csv(f"{OUTPUT_DIR}/velocity_driver_genes.csv", index=False)
    print(f"  Driver genes → {OUTPUT_DIR}/velocity_driver_genes.csv")

    top_drivers = [g for g in pd.unique(driver_df.values.flatten())[:300]
                   if g in combined.var_names]
    scv.pl.heatmap(
        combined,
        var_names=top_drivers,
        sortby="latent_time",
        col_color="cell_type",
        n_convolve=100,
        save="latent_time_heatmap.png",
        show=False,
    )

    print(f"\nDone. Outputs in {OUTPUT_DIR}/")


if __name__ == "__main__":
    main()
