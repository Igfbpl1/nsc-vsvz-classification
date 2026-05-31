"""QC, normalization, dimension reduction, clustering."""
from __future__ import annotations

import anndata as ad
import scanpy as sc


def qc_filter(adata: ad.AnnData,
              min_genes: int = 500,
              min_cells: int = 10,
              max_pct_mt: float = 10.0,
              max_counts: int | None = None) -> ad.AnnData:
    adata.var["mt"] = adata.var_names.str.startswith(("mt-", "Mt-"))
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None, log1p=False)

    n0 = adata.n_obs
    keep = (adata.obs["n_genes_by_counts"] >= min_genes) & \
           (adata.obs["pct_counts_mt"] <= max_pct_mt)
    if max_counts is not None:
        keep &= (adata.obs["total_counts"] <= max_counts)
    adata = adata[keep].copy()
    print(f"  cells: {n0} -> {adata.n_obs} (kept {100*adata.n_obs/n0:.1f}%)")

    n_var0 = adata.n_vars
    sc.pp.filter_genes(adata, min_cells=min_cells)
    print(f"  genes: {n_var0} -> {adata.n_vars}")
    return adata


def normalize_and_embed(adata: ad.AnnData,
                        n_hvg: int = 2000,
                        n_pcs: int = 30,
                        resolution: float = 0.8) -> ad.AnnData:
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg, flavor="seurat", batch_key="sample_id")
    print(f"  HVGs: {int(adata.var['highly_variable'].sum())}")

    adata_hvg = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata_hvg, max_value=10)
    sc.tl.pca(adata_hvg, n_comps=n_pcs)
    sc.pp.neighbors(adata_hvg, n_pcs=n_pcs, n_neighbors=15)
    sc.tl.umap(adata_hvg)
    sc.tl.leiden(adata_hvg, resolution=resolution, flavor="igraph",
                 n_iterations=2, directed=False)

    adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"]
    adata.obsm["X_umap"] = adata_hvg.obsm["X_umap"]
    adata.obs["leiden"] = adata_hvg.obs["leiden"].values
    adata.obsp["connectivities"] = adata_hvg.obsp["connectivities"]
    adata.obsp["distances"] = adata_hvg.obsp["distances"]
    adata.uns["neighbors"] = adata_hvg.uns["neighbors"]
    print(f"  leiden clusters: {adata.obs['leiden'].nunique()}")
    return adata
