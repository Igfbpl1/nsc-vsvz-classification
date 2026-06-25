"""Build velocity_combined.h5ad from kb-python MTX outputs.

Called from run_pipeline.py with a rebuild guard — only re-run when needed.
"""

from __future__ import annotations

import pandas as pd
from scipy.io import mmread
import anndata as ad
import scanpy as sc
import scvelo as scv
import preprocess

from pathlib import Path

ROOT       = Path(__file__).parent
H5AD_PATH  = str(ROOT / "outputs" / "processed.h5ad")
T2G_PATH   = str(ROOT / "sra_runs" / "t2g.txt")
OUTPUT_DIR = str(ROOT / "outputs" / "velocity")
OUT_H5AD   = str(ROOT / "outputs" / "velocity" / "velocity_combined.h5ad")

SAMPLE_CONFIGS = [
    {"h5ad_sample_id": "GSM8253792", "velocity_label": "CD1_Cntl_0wksRecov",           "kb_dir": "sra_runs/kb_output_GSM8253792/counts_unfiltered"},
    {"h5ad_sample_id": "GSM8253793", "velocity_label": "CD1_Cntl_3wksRecov",      "kb_dir": "sra_runs/kb_output_GSM8253793/counts_unfiltered"},
    {"h5ad_sample_id": "GSM8253794", "velocity_label": "CD1_CupRap_Rep1_3wksRecov",         "kb_dir": "sra_runs/kb_output_GSM8253794/counts_unfiltered"},
    {"h5ad_sample_id": "GSM8253795", "velocity_label": "CD1_CupRap_rep2_3wksRecov",    "kb_dir": "sra_runs/kb_output_GSM8253795/counts_unfiltered"},
    {"h5ad_sample_id": "GSM8253796", "velocity_label": "NesCre_Cntl_Rep1_3wksRecov",               "kb_dir": "sra_runs/kb_output_GSM8253796/counts_unfiltered"},
    {"h5ad_sample_id": "GSM8253797", "velocity_label": "NesCre_Cntl_Rep2_3wksRecov",          "kb_dir": "sra_runs/kb_output_GSM8253797/counts_unfiltered"},
    {"h5ad_sample_id": "GSM8253798", "velocity_label": "NesCre_CR_Rep1_3wksRecov",        "kb_dir": "sra_runs/kb_output_GSM8253798/counts_unfiltered"},
    {"h5ad_sample_id": "GSM8253799", "velocity_label": "NesCre_CR_Rep2_3wksRecov",        "kb_dir": "sra_runs/kb_output_GSM8253799/counts_unfiltered"},
    {"h5ad_sample_id": "GSM8647353", "velocity_label": "CD1_CupRap_Rep2_0wksRecov", "kb_dir": "sra_runs/kb_output_GSM8647353/counts_unfiltered"},
    {"h5ad_sample_id": "GSM8647352", "velocity_label": "CD1_CupRap_Rep1_0wksRecov", "kb_dir": "sra_runs/kb_output_GSM8647352/counts_unfiltered"},
]

MIN_SHARED_COUNTS = 20
N_TOP_GENES       = 2000
N_PCS             = 30
N_NEIGHBORS       = 30


def load_ensembl_to_name(t2g_path: str) -> dict[str, str]:
    df = pd.read_csv(
        t2g_path, sep="\t", header=None,
        names=["transcript_id", "gene_id", "gene_name",
               "transcript_name", "chrom", "start", "end", "strand"],
    )
    return df.drop_duplicates("gene_id").set_index("gene_id")["gene_name"].to_dict()


def load_kb_layer(kb_dir: str, layer: str) -> ad.AnnData:
    matrix   = mmread(f"{kb_dir}/{layer}.mtx").tocsr()
    barcodes = pd.read_csv(f"{kb_dir}/{layer}.barcodes.txt", header=None)[0].tolist()
    genes    = pd.read_csv(f"{kb_dir}/{layer}.genes.txt",    header=None)[0].tolist()
    return ad.AnnData(X=matrix, obs=pd.DataFrame(index=barcodes), var=pd.DataFrame(index=genes))


def build_sample_adata(config: dict, ref: ad.AnnData, ensembl_to_name: dict[str, str]) -> ad.AnnData:
    h5ad_sid = config["h5ad_sample_id"]
    label    = config["velocity_label"]
    kb_dir   = config["kb_dir"]

    print(f"\n[{label}] Loading kb layers from {kb_dir} ...")
    sp_adata = load_kb_layer(kb_dir, "spliced")
    un_adata = load_kb_layer(kb_dir, "unspliced")
    print(f"  spliced:   {sp_adata.shape[0]:,} barcodes × {sp_adata.shape[1]:,} genes")
    print(f"  unspliced: {un_adata.shape[0]:,} barcodes × {un_adata.shape[1]:,} genes")

    ref_sub  = ref[ref.obs["sample_id"] == h5ad_sid]
    raw_bc   = ref_sub.obs_names.str.removeprefix(f"{h5ad_sid}_").str.removesuffix("-1")

    keep_mask    = raw_bc.isin(sp_adata.obs_names) & raw_bc.isin(un_adata.obs_names)
    matched_raw_bc = raw_bc[keep_mask]
    matched_ref    = ref_sub[keep_mask]
    print(f"  {h5ad_sid} cells in h5ad: {len(raw_bc):,}")
    print(f"  matched in both layers:    {len(matched_raw_bc):,}")

    if len(matched_raw_bc) == 0:
        raise ValueError(f"No barcodes matched for {h5ad_sid}.")

    shared_genes = sp_adata.var_names[sp_adata.var_names.isin(un_adata.var_names)]
    sp_sub = sp_adata[matched_raw_bc.tolist(), shared_genes]
    un_sub = un_adata[matched_raw_bc.tolist(), shared_genes]

    result = ad.AnnData(
        X=sp_sub.X.copy(),
        obs=matched_ref.obs.copy(),
        var=pd.DataFrame({"gene_name": [ensembl_to_name.get(g, g) for g in shared_genes]}, index=shared_genes),
    )
    result.obs_names          = matched_raw_bc.tolist()
    result.layers["spliced"]  = sp_sub.X.copy()
    result.layers["unspliced"] = un_sub.X.copy()
    result.obs["velocity_label"] = label
    result.obsm["X_umap_ref"] = matched_ref.obsm["X_umap"].copy()
    return result


def build_velocity() -> None:
    """Preprocess samples, fit dynamics, write velocity_combined.h5ad."""
    import os
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    print("Loading processed.h5ad ...")
    ref = sc.read_h5ad(H5AD_PATH)
    cluster_to_type = ref.uns.get("cluster_to_type", {})
    ref.obs["cell_type"] = ref.obs["leiden"].map(cluster_to_type).fillna("Unknown")
    print(f"  {ref.shape[0]:,} cells, {ref.shape[1]:,} genes")

    print("\nLoading Ensembl → gene name mapping ...")
    ensembl_to_name = load_ensembl_to_name(T2G_PATH)
    print(f"  {len(ensembl_to_name):,} entries")

    samples = [build_sample_adata(cfg, ref, ensembl_to_name) for cfg in SAMPLE_CONFIGS]

    print("\nConcatenating samples ...")
    keys     = [cfg["velocity_label"] for cfg in SAMPLE_CONFIGS]
    combined = ad.concat(samples, keys=keys, label="concat_key", join="inner")
    print(f"  combined: {combined.shape[0]:,} cells × {combined.shape[1]:,} genes")
    print(f"  velocity_label counts:\n{combined.obs['velocity_label'].value_counts()}")

    scv.settings.verbosity = 3
    scv.settings.set_figure_params("scvelo", figsize=(8, 6))

    print("\n[scVelo] Filtering genes with low shared counts ...")
    scv.pp.filter_genes(combined, min_shared_counts=MIN_SHARED_COUNTS)

    print("[scVelo] Normalizing, selecting HVGs, computing neighbors ...")
    combined = preprocess.normalize_embed_velocity(
        combined, n_top_genes=N_TOP_GENES, n_pcs=N_PCS, n_neighbors=N_NEIGHBORS,
    )
    combined.obsm["X_umap"] = combined.obsm["X_umap_ref"]
    print(f"  {combined.n_obs:,} velocity cells; using ref UMAP for plots")

    print("[scVelo] Computing moments ...")
    scv.pp.moments(combined, n_pcs=None, n_neighbors=None)

    print("[scVelo] Recovering dynamics (dynamical model — ~30 min) ...")
    scv.tl.recover_dynamics(combined, n_jobs=8, show_progress_bar=False)

    print("[scVelo] Computing velocity ...")
    scv.tl.velocity(combined, mode="dynamical")

    print("[scVelo] Building velocity graph ...")
    scv.tl.velocity_graph(combined, show_progress_bar=False)

    combined.write_h5ad(OUT_H5AD)
    print(f"\nSaved → {OUT_H5AD}")
