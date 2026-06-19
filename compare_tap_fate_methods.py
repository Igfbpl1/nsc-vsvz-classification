"""
Compare CellRank and ML fate predictions for TAP cells.

  1. ML model   → XGBoost-trained P(OL) per TAP (outputs/ol_commitment.csv)
  2. CellRank   → Velocity-based fate probability P(OL) per TAP

Reports Pearson/Spearman correlation and mean P(OL) per condition.

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
from cellrank.kernels import ConnectivityKernel, PrecomputedKernel
from scipy.stats import pearsonr, spearmanr

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
    # Build velocity transition matrix via scvelo, then wrap in PrecomputedKernel.
    # This avoids cellrank 2.0.7's VelocityKernel.__init__ shape-check that calls
    # np.testing.assert_array_equal(x=, y=) — kwargs renamed to actual/desired in
    # numpy 2 (fixed in cellrank >=2.1.0, which conflicts with scipy>=1.17).
    print("\n[CellRank] building transition matrix from velocity ...")
    T_vel = scv.utils.get_transition_matrix(adata)
    vk = PrecomputedKernel(T_vel, adata=adata)
    ck = ConnectivityKernel(adata).compute_transition_matrix()
    combined_kernel = 0.8 * vk + 0.2 * ck

    # Terminal cells defined from cell_type labels. The OL/NB cluster
    # assignments are well-supported by canonical markers in
    # outputs/velocity/dotplot__marker_dotplot.png (Mbp/Mog/Plp1 tight in OL,
    # Dcx/Tubb3 tight in Neuroblast). Using cell_type avoids the cost of
    # GPCCA's Schur decomposition (~30 min dense on 36k cells without SLEPc).
    print("[CellRank] setting terminal states from cell_type labels ...")
    g = cr.estimators.GPCCA(combined_kernel)
    ol_cells = adata.obs.index[adata.obs["cell_type"] == "OL"]
    nb_cells = adata.obs.index[adata.obs["cell_type"] == "Neuroblast"]
    g.set_terminal_states({"OL": ol_cells, "Neuroblast": nb_cells})

    print("[CellRank] computing fate probabilities ...")
    # cellrank 2.0.7's row-sum tolerance is hardcoded at rtol=1e-3 (the
    # `check_sum_tol` kwarg only exists in >=2.1.0). Tighten the solver tol
    # so probabilities converge closely enough to satisfy that check.
    g.compute_fate_probabilities(tol=1e-10)

    adata.obs["cellrank_P_OL"] = g.fate_probabilities[:, "OL"].X.squeeze()
    print("  CellRank P(OL) computed.")

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

    before = len(tap_df)
    tap_df = tap_df.dropna(subset=["cellrank_P_OL", "ml_P_OL"])
    print(f"  TAPs with both predictions: {len(tap_df)}/{before}")

    # Binary calls — P(OL) > 0.5 is natural since OL + NB are the only
    # terminal states, so P(OL) + P(NB) = 1
    tap_df["cellrank_ol"] = (tap_df["cellrank_P_OL"] > 0.5).astype(int)
    tap_df["ml_ol"]       = (tap_df["ml_P_OL"] > 0.5).astype(int)

    tap_df.to_csv(OUT / "tap_fate_comparison.csv", index=False)
    print(f"  → {OUT}/tap_fate_comparison.csv")

    # ── Correlation ───────────────────────────────────────────────────────────
    pearson_r, pearson_p   = pearsonr(tap_df["cellrank_P_OL"], tap_df["ml_P_OL"])
    spearman_r, spearman_p = spearmanr(tap_df["cellrank_P_OL"], tap_df["ml_P_OL"])
    print(f"\nCellRank vs ML (n={len(tap_df)} TAPs):")
    print(f"  Pearson r  = {pearson_r:.3f}  (p={pearson_p:.2e})")
    print(f"  Spearman r = {spearman_r:.3f}  (p={spearman_p:.2e})")

    # ── Per-condition OL and NB counts and means ─────────────────────────────
    # Since OL + NB are the only terminal states, P(NB) = 1 - P(OL)
    grp = tap_df.groupby("velocity_label")
    n_taps             = grp["ml_P_OL"].count()
    cellrank_n_ol      = grp["cellrank_ol"].sum()
    ml_n_ol            = grp["ml_ol"].sum()
    cellrank_mean_p_ol = grp["cellrank_P_OL"].mean()
    ml_mean_p_ol       = grp["ml_P_OL"].mean()

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
    }).sort_index()

    print("\nPer-condition OL/NB counts (threshold > 0.5) and mean P:")
    print(per_cond.to_string())

    per_cond.to_csv(OUT / "tap_fate_per_condition.csv")
    print(f"\n  → {OUT}/tap_fate_per_condition.csv")


if __name__ == "__main__":
    run_compare_fate_methods()
