"""
RNA Velocity & Trajectory Pipeline — Phase 2–4
Loads kb-python MTX output (spliced + unspliced) and merges with processed.h5ad.
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

SAMPLE_CONFIGS = [
    {
        "h5ad_sample_id": "GSM8253792",
        "velocity_label": "CD1_Cntl",
        "kb_dir": "sra_runs/kb_output_GSM8253792/counts_unfiltered",
    },
    {
        "h5ad_sample_id": "GSM8253793",
        "velocity_label": "CD1_Cntl_3wks",
        "kb_dir": "sra_runs/kb_output_GSM8253793/counts_unfiltered",
    },
    {
        "h5ad_sample_id": "GSM8253794",
        "velocity_label": "CD1_CupRap",
        "kb_dir": "sra_runs/kb_output_GSM8253794/counts_unfiltered",
    },
    {
        "h5ad_sample_id": "GSM8253796",
        "velocity_label": "Cntl",
        "kb_dir": "sra_runs/kb_output_GSM8253796/counts_unfiltered",
    },
    {
        "h5ad_sample_id": "GSM8253798",
        "velocity_label": "CupRap_Rep1",
        "kb_dir": "sra_runs/kb_output_GSM8253798/counts_unfiltered",
    },
    {
        "h5ad_sample_id": "GSM8253799",
        "velocity_label": "CupRap_Rep2",
        "kb_dir": "sra_runs/kb_output_GSM8253799/counts_unfiltered",
    },
    {
        "h5ad_sample_id": "GSM8647353",
        "velocity_label": "CD1_CupRap_0wks_Rep2",
        "kb_dir": "sra_runs/kb_output_GSM8647353/counts_unfiltered",
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

    # Keep ref UMAP for reference; fresh UMAP computed after PCA on velocity cells
    result.obsm["X_umap_ref"] = matched_ref.obsm["X_umap"].copy()

    return result


# ── Main ──────────────────────────────────────────────────────────────────────
def run_rna_velocity_pipeline():
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

    # ── Preprocessing (mirrors preprocess.normalize_and_embed) ────────────────
    scv.settings.verbosity = 3
    scv.settings.set_figure_params("scvelo", figsize=(8, 6))

    print("\n[scVelo] Filtering genes with low shared counts ...")
    scv.pp.filter_genes(combined, min_shared_counts=MIN_SHARED_COUNTS)

    # Normalize and log-transform X (layers stay as raw counts for scVelo)
    print("[scVelo] Normalizing ...")
    combined.layers["counts"] = combined.X.copy()
    sc.pp.normalize_total(combined, target_sum=1e4)
    sc.pp.log1p(combined)

    # Batch-aware HVG selection (same as run_pipeline: seurat flavor, per sample_id)
    print(f"[scVelo] Selecting top {N_TOP_GENES} HVGs (batch-aware) ...")
    sc.pp.highly_variable_genes(
        combined, n_top_genes=N_TOP_GENES, flavor="seurat", batch_key="sample_id"
    )
    print(f"  HVGs: {int(combined.var['highly_variable'].sum())}")

    # Scale a copy for PCA/UMAP/neighbors — keeps combined.X unscaled for velocity
    # (same pattern as normalize_and_embed: adata_hvg is the scaled working copy)
    adata_hvg = combined[:, combined.var["highly_variable"]].copy()
    sc.pp.scale(adata_hvg, max_value=10)
    sc.tl.pca(adata_hvg, n_comps=N_PCS)
    sc.pp.neighbors(adata_hvg, n_pcs=N_PCS, n_neighbors=N_NEIGHBORS)
    sc.tl.umap(adata_hvg)

    # Copy PCA, neighbors, and graph back into combined (used for velocity computation)
    combined.obsm["X_pca"]          = adata_hvg.obsm["X_pca"]
    combined.obsp["connectivities"] = adata_hvg.obsp["connectivities"]
    combined.obsp["distances"]      = adata_hvg.obsp["distances"]
    combined.uns["neighbors"]       = adata_hvg.uns["neighbors"]
    del adata_hvg

    # Use ref UMAP for plotting — it has well-separated clusters from all 10 samples.
    # Velocity is computed on the local neighborhood graph above; visualization
    # on the established embedding keeps cell type boundaries interpretable.
    combined.obsm["X_umap"] = combined.obsm["X_umap_ref"]
    print(f"  neighbors computed on {combined.n_obs:,} velocity cells; using ref UMAP for plots")

    # Subset combined to HVGs so scVelo moments use the same gene space
    combined = combined[:, combined.var["highly_variable"]].copy()

    print("[scVelo] Computing moments ...")
    scv.pp.moments(combined, n_pcs=None, n_neighbors=None)  # uses pre-computed neighbors

    # ── Velocity ───────────────────────────────────────────────────────────────
    print("[scVelo] Recovering dynamics (dynamical model — ~30 min) ...")
    scv.tl.recover_dynamics(combined, n_jobs=4, show_progress_bar=False)

    print("[scVelo] Computing velocity ...")
    scv.tl.velocity(combined, mode="dynamical")

    print("[scVelo] Building velocity graph ...")
    scv.tl.velocity_graph(combined, show_progress_bar=False)

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

    for label in combined.obs["velocity_label"].unique():
        sub = combined[combined.obs["velocity_label"] == label]
        scv.pl.velocity_embedding_stream(
            sub, basis="umap", color=color_key,
            title=f"Velocity — {label}",
            save=f"stream_{label}.png", show=False,
        )

    # load ensembl → gene name mapping
    ensembl_to_name = load_ensembl_to_name(T2G_PATH)

    # SHAP lookup (gene name → rank + value)
    shap_df = pd.read_csv("outputs/trigger_genes.csv")
    shap_lookup = shap_df.set_index("gene")[["shap_mean_abs"]].copy()
    shap_lookup["shap_rank"] = range(1, len(shap_lookup) + 1)

    # canonical marker gene names (to flag known markers)
    from markers import MARKERS
    canonical_genes = {g for genes in MARKERS.values() for g in genes}

    def write_driver_csvs(adata, label_suffix=""):
        scv.tl.rank_velocity_genes(adata, groupby=color_key, min_corr=0.3)
        rvg = adata.uns["rank_velocity_genes"]
        driver_df = pd.DataFrame(rvg["names"])
        scores_df = pd.DataFrame(rvg["scores"])

        var_df = adata.var[
            [c for c in ["fit_likelihood", "spearmans_score"] if c in adata.var.columns]
        ]

        lineage_cols = [c for c in ["NSC", "TAP", "Neuroblast", "OPC", "COP", "OL"]
                        if c in driver_df.columns]

        for col in lineage_cols:
            ensembl_ids = driver_df[col]
            gene_names  = ensembl_ids.map(lambda x: ensembl_to_name.get(x, x))

            cell_df = pd.DataFrame()
            cell_df["rank"]       = range(1, len(ensembl_ids) + 1)
            cell_df["ensembl_id"] = ensembl_ids.values
            cell_df["gene_name"]  = gene_names.values
            cell_df["corr"]       = scores_df[col].values
            cell_df["likelihood"] = ensembl_ids.map(
                lambda x: var_df["fit_likelihood"].get(x, None) if "fit_likelihood" in var_df.columns else None
            ).values
            cell_df["spearmans"]  = ensembl_ids.map(
                lambda x: var_df["spearmans_score"].get(x, None) if "spearmans_score" in var_df.columns else None
            ).values
            cell_df["shap_rank"]  = gene_names.map(
                lambda g: shap_lookup.loc[g, "shap_rank"] if g in shap_lookup.index else None
            ).values
            cell_df["shap_value"] = gene_names.map(
                lambda g: shap_lookup.loc[g, "shap_mean_abs"] if g in shap_lookup.index else None
            ).values
            cell_df["canonical"]  = gene_names.map(
                lambda g: 1 if g in canonical_genes else 0
            ).values

            suffix = f"_{label_suffix}" if label_suffix else ""
            out_path = f"{OUTPUT_DIR}/velocity_drivers_{col}{suffix}.csv"
            cell_df.to_csv(out_path, index=False)
            print(f"  {col}{suffix} drivers → {out_path}")

    print("\n[scVelo] Ranking velocity genes (pooled) ...")
    write_driver_csvs(combined)
    # capture pooled rankings before per-condition runs can overwrite combined.uns
    driver_df_pooled = pd.DataFrame(combined.uns["rank_velocity_genes"]["names"])

    print("\n[scVelo] Ranking velocity genes per condition ...")
    for condition in combined.obs["treatment"].unique():
        subset = combined[combined.obs["treatment"] == condition].copy()
        print(f"  condition={condition}: {subset.n_obs:,} cells")
        write_driver_csvs(subset, label_suffix=condition)

    print(f"\nDone. Outputs in {OUTPUT_DIR}/")
