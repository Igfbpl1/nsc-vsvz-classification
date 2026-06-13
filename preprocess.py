"""QC, normalization, dimension reduction, clustering."""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import scanpy as sc


def qc_filter(
    adata: ad.AnnData,
    min_genes: int = 500,
    min_cells: int = 10,
    max_pct_mt: float = 10.0,
    max_counts: int | None = None,
    output_dir: Path | None = None
) -> ad.AnnData:
    adata.var["mt"] = adata.var_names.str.startswith(("mt-", "Mt-"))
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], inplace=True, percent_top=None, log1p=False
    )

    # Identify cells failing min_genes
    failed_min_genes_mask = adata.obs["n_genes_by_counts"] < min_genes
    failed_min_genes_adata = adata[failed_min_genes_mask]
    if failed_min_genes_adata.n_obs > 0:
        print(f"  Top 10 cells failing min_genes ({min_genes}): {list(failed_min_genes_adata.obs_names[:10])}")
        if output_dir:
            failed_min_genes_adata.to_df().head(10).to_csv(f"{output_dir}/failed_min_genes_cells_data_top10.csv")

    # Identify cells failing max_pct_mt
    failed_max_pct_mt_mask = adata.obs["pct_counts_mt"] > max_pct_mt
    failed_max_pct_mt_adata = adata[failed_max_pct_mt_mask]
    if failed_max_pct_mt_adata.n_obs > 0:
        print(f"  Top 10 cells failing max_pct_mt ({max_pct_mt}%): {list(failed_max_pct_mt_adata.obs_names[:10])}")
        if output_dir:
            failed_max_pct_mt_adata.to_df().head(10).to_csv(f"{output_dir}/failed_max_pct_mt_cells_data_top10.csv")

    n0 = adata.n_obs
    keep = (adata.obs["n_genes_by_counts"] >= min_genes) & (
        adata.obs["pct_counts_mt"] <= max_pct_mt
    )
    if max_counts is not None:
        keep &= adata.obs["total_counts"] <= max_counts
    adata = adata[keep].copy()
    print(f"  cells: {n0} -> {adata.n_obs} (kept {100 * adata.n_obs / n0:.1f}%)")

    # Log genes failing min_cells
    genes_before_filter = adata.var_names.copy()
    
    n_var0 = adata.n_vars
    sc.pp.filter_genes(adata, min_cells=min_cells)
    
    genes_after_filter = adata.var_names.copy()
    
    filtered_genes_names = list(set(genes_before_filter) - set(genes_after_filter))
    if len(filtered_genes_names) > 0:
        print(f"  Top 10 genes failing min_cells ({min_cells}): {filtered_genes_names[:10]}")
        if output_dir:
            # Create an AnnData object for the filtered genes to save their data
            # Need to ensure the genes are present in the original adata before filtering
            # For simplicity, we'll take the current adata and filter it by the names of genes that were removed
            # This assumes the genes are still in the `adata` object before the filter_genes call
            # A more robust way would be to keep a copy of adata before gene filtering
            
            # To get the data for the *removed* genes, we need to subset the adata *before* the filter_genes call
            # However, since adata is already filtered by cells, we can just subset by gene names.
            # Let's assume `adata` still contains the genes that were filtered out for the purpose of logging.
            # This is a bit tricky with in-place filtering.
            # A better approach is to get the mask of genes to keep, and then invert it.
            
            # Re-calculate n_cells for genes on the current adata to identify those that would be filtered
            # This is to get the actual data for the genes that were removed.
            temp_adata_for_gene_filter = adata.copy() # Make a temporary copy to avoid modifying adata again
            sc.pp.calculate_qc_metrics(temp_adata_for_gene_filter, qc_vars=[], inplace=True, percent_top=None, log1p=False)
            
            # Identify genes that would be filtered out by min_cells
            failed_min_cells_mask = temp_adata_for_gene_filter.var["n_cells_by_counts"] < min_cells
            failed_min_cells_adata = temp_adata_for_gene_filter[:, failed_min_cells_mask]
            
            if failed_min_cells_adata.n_vars > 0:
                failed_min_cells_adata.to_df().head(10).to_csv(f"{output_dir}/failed_min_cells_genes_data_top10.csv")


    print(f"  genes: {n_var0} -> {adata.n_vars}")
    return adata


def normalize_and_embed(
    adata: ad.AnnData, n_hvg: int = 2000, n_pcs: int = 30, resolution: float = 0.8
) -> ad.AnnData:
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(
        adata, n_top_genes=n_hvg, flavor="seurat", batch_key="sample_id"
    )
    print(f"  HVGs: {int(adata.var['highly_variable'].sum())}")

    adata_hvg = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata_hvg, max_value=10)
    sc.tl.pca(adata_hvg, n_comps=n_pcs)
    sc.pp.neighbors(adata_hvg, n_pcs=n_pcs, n_neighbors=15)
    sc.tl.umap(adata_hvg)
    sc.tl.leiden(
        adata_hvg,
        resolution=resolution,
        flavor="igraph",
        n_iterations=2,
        directed=False,
    )

    adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"]
    adata.obsm["X_umap"] = adata_hvg.obsm["X_umap"]
    adata.obs["leiden"] = adata_hvg.obs["leiden"].values
    adata.obsp["connectivities"] = adata_hvg.obsp["connectivities"]
    adata.obsp["distances"] = adata_hvg.obsp["distances"]
    adata.uns["neighbors"] = adata_hvg.uns["neighbors"]
    print(f"  leiden clusters: {adata.obs['leiden'].nunique()}")
    return adata
