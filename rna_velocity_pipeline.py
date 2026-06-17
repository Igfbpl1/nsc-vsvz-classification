"""RNA velocity analysis: plots and driver gene rankings from velocity_combined.h5ad."""

from __future__ import annotations

import pandas as pd
import scanpy as sc
import scvelo as scv
import matplotlib
matplotlib.use("Agg")

from velocity_build import load_ensembl_to_name, OUT_H5AD, OUTPUT_DIR, T2G_PATH
from markers import MARKERS


def run_rna_velocity_pipeline() -> None:
    combined = sc.read_h5ad(OUT_H5AD)

    scv.settings.figdir    = OUTPUT_DIR
    scv.settings.verbosity = 3
    scv.settings.set_figure_params("scvelo", figsize=(8, 6))
    color_key = "cell_type"

    # ── Stream plots ───────────────────────────────────────────────────────────
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

    # ── Driver gene rankings ───────────────────────────────────────────────────
    ensembl_to_name = load_ensembl_to_name(T2G_PATH)

    shap_df    = pd.read_csv("outputs/trigger_genes.csv")
    shap_lookup = shap_df.set_index("gene")[["shap_mean_abs"]].copy()
    shap_lookup["shap_rank"] = range(1, len(shap_lookup) + 1)

    canonical_genes = {g for genes in MARKERS.values() for g in genes}

    def write_driver_csvs(adata, label_suffix=""):
        scv.tl.rank_velocity_genes(adata, groupby=color_key, min_corr=0.3)
        rvg       = adata.uns["rank_velocity_genes"]
        driver_df = pd.DataFrame(rvg["names"])
        scores_df = pd.DataFrame(rvg["scores"])
        var_df    = adata.var[
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

            suffix   = f"_{label_suffix}" if label_suffix else ""
            out_path = f"{OUTPUT_DIR}/velocity_drivers_{col}{suffix}.csv"
            cell_df.to_csv(out_path, index=False)
            print(f"  {col}{suffix} drivers → {out_path}")

    print("\n[scVelo] Ranking velocity genes (pooled) ...")
    write_driver_csvs(combined)
    driver_df_pooled = pd.DataFrame(combined.uns["rank_velocity_genes"]["names"])

    print("\n[scVelo] Ranking velocity genes per condition ...")
    for condition in combined.obs["treatment"].unique():
        subset = combined[combined.obs["treatment"] == condition].copy()
        print(f"  condition={condition}: {subset.n_obs:,} cells")
        write_driver_csvs(subset, label_suffix=condition)

    # ── Phase portraits ────────────────────────────────────────────────────────
    combined_named = combined.copy()
    combined_named.var_names = pd.Index(
        [ensembl_to_name.get(g, g) for g in combined.var_names]
    )
    combined_named.var_names_make_unique()

    lineage_cols_plot = [c for c in ["NSC", "TAP", "Neuroblast", "OPC", "COP", "OL"]
                         if c in driver_df_pooled.columns]

    print("\nSaving phase portraits (top 4 driver genes per cell type) ...")
    for col in lineage_cols_plot:
        top_genes = [
            ensembl_to_name.get(g, g)
            for g in driver_df_pooled[col].iloc[:4].tolist()
        ]
        top_genes = [g for g in top_genes if g in combined_named.var_names]
        if not top_genes:
            continue
        scv.pl.velocity(
            combined_named, var_names=top_genes, basis="umap",
            color=color_key,
            save=f"phase_{col}.png", show=False,
        )
        print(f"  {col}: {top_genes}")

    print(f"\nDone. Outputs in {OUTPUT_DIR}/")
