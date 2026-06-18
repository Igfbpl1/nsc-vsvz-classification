"""Marker gene dotplot: validates cell type labels against canonical marker expression."""

from __future__ import annotations

import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")

from markers import MARKERS
from velocity_build import load_ensembl_to_name, OUT_H5AD, OUTPUT_DIR, T2G_PATH

CT_ORDER = [
    "NSC", "TAP", "Neuroblast",
    "OPC", "COP", "OL",
    "Astrocyte", "Microglia", "Ependymal", "Endothelial", "Mural",
]


def _build_name_to_eid(adata: sc.AnnData, t2g_path: str) -> dict[str, str]:
    """Map gene symbol → versioned Ensembl ID present in adata, via base-ID matching."""
    ensembl_to_name = load_ensembl_to_name(t2g_path)
    adata_base = {v.split(".")[0]: v for v in adata.var_names}
    name_to_eid: dict[str, str] = {}
    for versioned_eid, sym in ensembl_to_name.items():
        base = versioned_eid.split(".")[0]
        if base in adata_base:
            name_to_eid[sym] = adata_base[base]
    return name_to_eid


def _subset_to_markers(adata: sc.AnnData, name_to_eid: dict[str, str]) -> sc.AnnData:
    """Return adata subset to marker genes with gene-symbol var_names.

    Genes shared across multiple cell types are assigned to their first
    appearance in CT_ORDER to avoid duplicate columns in the dotplot.
    """
    seen: set[str] = set()
    eids, syms = [], []
    for ct in CT_ORDER:
        for gene in MARKERS.get(ct, []):
            eid = name_to_eid.get(gene)
            if eid and eid not in seen:
                seen.add(eid)
                eids.append(eid)
                syms.append(gene)

    sub = adata[:, eids].copy()
    sub.var_names = pd.Index(syms)
    return sub


def run_marker_dotplot() -> None:
    adata = sc.read_h5ad(OUT_H5AD)

    name_to_eid = _build_name_to_eid(adata, T2G_PATH)
    adata_sub = _subset_to_markers(adata, name_to_eid)

    ct_order = [ct for ct in CT_ORDER if ct in adata_sub.obs["cell_type"].values]
    adata_sub.obs["cell_type"] = pd.Categorical(
        adata_sub.obs["cell_type"], categories=ct_order, ordered=True
    )

    var_names = {
        ct: [g for g in MARKERS[ct] if g in adata_sub.var_names]
        for ct in ct_order
    }

    sc.settings.figdir = OUTPUT_DIR
    sc.pl.dotplot(
        adata_sub,
        var_names=var_names,
        groupby="cell_type",
        standard_scale="var",
        colorbar_title="Mean expr\n(scaled)",
        size_title="Fraction\nexpressing",
        figsize=(24, 8),
        save="_marker_dotplot.png",
        show=False,
    )

    found = sum(len(v) for v in var_names.values())
    total = sum(len(v) for v in MARKERS.values())
    print(f"Saved marker dotplot ({found}/{total} markers found) → {OUTPUT_DIR}/dotplot__marker_dotplot.png")


if __name__ == "__main__":
    run_marker_dotplot()
