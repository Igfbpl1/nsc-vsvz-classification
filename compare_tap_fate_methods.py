"""
Compare CellRank and ML fate predictions for TAP cells.

  1. ML model   → XGBoost-trained P(OL) per TAP (outputs/ol_commitment.csv)
  2. CellRank   → Velocity-based fate probability P(OL) per TAP

Reports binary agreement and mean P(OL) per condition.

Output:
  outputs/tap_fate_comparison.csv     — per-TAP cellrank_P_OL and ml_P_OL
  outputs/tap_fate_per_condition.csv  — mean P(OL) per condition per method
"""

from __future__ import annotations
from pathlib import Path

import cellrank as cr
import pandas as pd
import scanpy as sc
import scvelo as scv
from cellrank.kernels import ConnectivityKernel, VelocityKernel

from markers import OL_LINEAGE

ROOT = Path(__file__).parent
OUT  = ROOT / "outputs"


def run_compare_fate_methods():
    print("loading velocity_combined.h5ad ...")
    adata = sc.read_h5ad(OUT / "velocity" / "velocity_combined.h5ad")
    adata.obs_names_make_unique()
    print(f"  cells: {adata.n_obs:,}  genes: {adata.n_vars:,}")

    if "velocity_graph" not in adata.uns:
        print("recomputing velocity graph ...")
        scv.tl.velocity_graph(adata, show_progress_bar=False, n_jobs=1)

    # ── CellRank fate probabilities ───────────────────────────────────────────
    print("\n[CellRank] building transition matrix from velocity ...")
    vk = VelocityKernel(adata).compute_transition_matrix(n_jobs=4)
    ck = ConnectivityKernel(adata).compute_transition_matrix()
    combined_kernel = 0.8 * vk + 0.2 * ck

    # Terminal cells defined from cell_type labels. The OL/NB cluster
    # assignments are well-supported by canonical markers in
    # outputs/velocity/dotplot__marker_dotplot.png (Mbp/Mog/Plp1 tight in OL,
    # Dcx/Tubb3 tight in Neuroblast). Using cell_type avoids the cost of
    # GPCCA's Schur decomposition (~30 min dense on 36k cells without SLEPc).
    print("[CellRank] setting terminal states from cell_type labels ...")
    g = cr.estimators.GPCCA(combined_kernel)
    # Match the ML positive class (OL_LINEAGE = {OPC, COP, OL}) so CellRank's
    # P(OL_lineage) is apples-to-apples with the XGBoost classifier. OPC/COP
    # aren't true absorbing states velocity-wise, but bundling them into the
    # same terminal group preserves their forward flux to OL while letting us
    # report one combined commitment probability.
    ol_cells = adata.obs.index[adata.obs["cell_type"].isin(OL_LINEAGE)]
    nb_cells = adata.obs.index[adata.obs["cell_type"] == "Neuroblast"]
    g.set_terminal_states({"OL": ol_cells, "Neuroblast": nb_cells})

    print("[CellRank] computing fate probabilities ...")
    g.compute_fate_probabilities()

    adata.obs["cellrank_P_OL"] = g.fate_probabilities[:, "OL"].X.squeeze()
    adata.obs["cellrank_P_NB"] = g.fate_probabilities[:, "Neuroblast"].X.squeeze()
    print("  CellRank P(OL) computed.")

    # ── Persist CellRank results ──────────────────────────────────────────────
    # Write fate probs back into the h5ad so downstream scripts can load them
    # without recomputing the transition matrix (~35s).
    vk.write_to_adata()
    adata.write(OUT / "velocity" / "velocity_combined.h5ad", compression="gzip")
    print(f"  → wrote cellrank_P_OL/P_NB + transition matrix to {OUT}/velocity/velocity_combined.h5ad")

    # Also save a lightweight per-cell CSV for all cells (not just TAPs)
    cellrank_all = adata.obs[["sample_id", "cell_type", "velocity_label",
                              "cellrank_P_OL", "cellrank_P_NB"]].copy()
    cellrank_all.index.name = "barcode"
    cellrank_all.to_csv(OUT / "cellrank_fate_probabilities.csv")
    print(f"  → wrote {OUT}/cellrank_fate_probabilities.csv ({len(cellrank_all):,} cells)")

    # ── Align TAPs with ML predictions ───────────────────────────────────────
    print("\n[ML] loading XGBoost predictions ...")
    ml = pd.read_csv(OUT / "ol_commitment.csv").set_index("cell_id")

    adata.obs["original_cell_id"] = (
        adata.obs["sample_id"].astype(str) + "_" + adata.obs_names + "-1"
    )
    tap_mask = adata.obs["cell_type"] == "TAP"
    tap_df = adata.obs[tap_mask][
        ["original_cell_id", "velocity_label", "cellrank_P_OL"]
    ].copy()
    tap_df["ml_P_OL"] = tap_df["original_cell_id"].map(ml["prob_OL"])
    tap_df["bias"]    = tap_df["original_cell_id"].map(ml["bias"])

    before = len(tap_df)
    tap_df = tap_df.dropna(subset=["cellrank_P_OL", "ml_P_OL", "bias"])
    print(f"  TAPs with both predictions + bias: {len(tap_df)}/{before}")

    # Binary calls — P(OL) > 0.5 is natural since OL + NB are the only
    # terminal states, so P(OL) + P(NB) = 1
    tap_df["cellrank_ol"] = (tap_df["cellrank_P_OL"] > 0.5).astype(int)
    tap_df["ml_ol"]       = (tap_df["ml_P_OL"] > 0.5).astype(int)
    tap_df["bias_ol"]     = (tap_df["bias"] > 0).astype(int)

    tap_df.to_csv(OUT / "tap_fate_comparison.csv", index=False)
    print(f"  → {OUT}/tap_fate_comparison.csv")

    # ── Agreements ──────────────────────────────────────────────
    print(f"\nMethod Comparison (n={len(tap_df)} TAPs):")

    ml_ol_mask = tap_df["ml_ol"] == 1
    cr_ol_mask = tap_df["cellrank_ol"] == 1
    bias_ol_mask = tap_df["bias_ol"] == 1

    ml_cr_intersect = (ml_ol_mask & cr_ol_mask).sum()
    ml_bias_intersect = (ml_ol_mask & bias_ol_mask).sum()

    print(f"\nOL-Specific Agreement (Focusing only on OL predictions):")
    print(f"  ML predicts OL:       {ml_ol_mask.sum()} TAPs")
    print(f"  CellRank predicts OL: {cr_ol_mask.sum()} TAPs")
    print(f"  Bias predicts OL:     {bias_ol_mask.sum()} TAPs")
    print(f"  Of the TAPs ML predicts as OL, CellRank also predicts OL: {(ml_cr_intersect / ml_ol_mask.sum() * 100):.2f}%")
    print(f"  Of the TAPs ML predicts as OL, Bias also predicts OL:     {(ml_bias_intersect / ml_ol_mask.sum() * 100):.2f}%")

    # ── Per-condition OL and NB counts and means ─────────────────────────────
    grp = tap_df.groupby("velocity_label")
    n_taps             = grp["ml_P_OL"].count()
    cellrank_n_ol      = grp["cellrank_ol"].sum()
    ml_n_ol            = grp["ml_ol"].sum()
    bias_n_ol          = grp["bias_ol"].sum()
    cellrank_mean_p_ol = grp["cellrank_P_OL"].mean()
    ml_mean_p_ol       = grp["ml_P_OL"].mean()
    bias_mean          = grp["bias"].mean()

    per_cond = pd.DataFrame({
        "n_TAPs":             n_taps,
        "cellrank_n_OL":      cellrank_n_ol,
        "cellrank_mean_P_OL": cellrank_mean_p_ol.round(3),
        "cellrank_n_NB":      n_taps - cellrank_n_ol,
        "cellrank_mean_P_NB": (1 - cellrank_mean_p_ol).round(3),
        "ml_n_OL":            ml_n_ol,
        "ml_mean_P_OL":       ml_mean_p_ol.round(3),
        "ml_n_NB":            n_taps - ml_n_ol,
        "ml_mean_P_NB":       (1 - ml_mean_p_ol).round(3),
        "bias_n_OL":          bias_n_ol,
        "bias_mean":          bias_mean.round(3),
    }).sort_index()

    print("\nPer-condition OL/NB counts (threshold > 0.5 for probs, > 0 for bias) and mean scores:")
    print(per_cond.to_string())

    per_cond.to_csv(OUT / "tap_fate_per_condition.csv")
    print(f"\n  → {OUT}/tap_fate_per_condition.csv")


if __name__ == "__main__":
    run_compare_fate_methods()
