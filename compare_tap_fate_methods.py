"""
Compare three methods for predicting TAP cell fate:

  1. ML model      → XGBoost-trained P(OL) per cell (outputs/ol_commitment.csv)
  2. Bias score    → Canonical marker panel scoring (existing in ol_commitment.csv)
  3. CellRank      → Velocity-based fate probabilities (most principled for intermediate cells)

For each TAP, compute three scores. Report pairwise correlations, agreement,
Cohen's kappa, and per-condition OL-leaning fractions.

Output:
  outputs/tap_fate_comparison.csv          — per-TAP scores from all 3 methods
  outputs/tap_fate_method_summary.csv      — pairwise agreement statistics
"""

from __future__ import annotations
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
from scipy.stats import pearsonr, spearmanr
import cellrank as cr
from cellrank.kernels import VelocityKernel, ConnectivityKernel


def cohen_kappa(a, b):
    """Cohen's kappa for two binary classifications."""
    p_observed = (a == b).mean()
    p_a = a.mean()
    p_b = b.mean()
    p_chance = p_a * p_b + (1 - p_a) * (1 - p_b)
    return (p_observed - p_chance) / (1 - p_chance)


def main():
    print("loading velocity_combined.h5ad ...")
    adata = sc.read_h5ad("outputs/velocity/velocity_combined.h5ad")
    adata.obs_names_make_unique()
    print(f"  cells: {adata.n_obs:,}  genes: {adata.n_vars:,}")

    if "velocity_graph" not in adata.uns:
        print("recomputing velocity graph ...")
        scv.tl.velocity_graph(adata, show_progress_bar=False, n_jobs=1)

    # ── CellRank fate probabilities ───────────────────────────────────────────
    print("\n[CellRank] building transition matrix from velocity ...")
    vk = VelocityKernel(adata).compute_transition_matrix(n_jobs=1)
    ck = ConnectivityKernel(adata).compute_transition_matrix()
    combined_kernel = 0.8 * vk + 0.2 * ck

    print("[CellRank] computing fate probabilities ...")
    g = cr.estimators.GPCCA(combined_kernel)
    ol_cells = adata.obs.index[adata.obs["cell_type"] == "OL"]
    nb_cells = adata.obs.index[adata.obs["cell_type"] == "Neuroblast"]
    g.set_terminal_states({"OL": ol_cells, "Neuroblast": nb_cells})
    g.compute_fate_probabilities(check_sum_tol=1e-1, tol=1e-8)

    fate_df = g.fate_probabilities.X.copy()
    fate_states = list(g.fate_probabilities.names)
    ol_idx = fate_states.index("OL")
    adata.obs["cellrank_P_OL"] = fate_df[:, ol_idx]
    print(f"  CellRank P(OL) computed.")

    # ── ML and bias score from existing predictions ──────────────────────────
    print("\n[ML/Bias] loading XGBoost predictions ...")
    ml = pd.read_csv("outputs/ol_commitment.csv").set_index("cell_id")

    # ── Align on common TAP cells ────────────────────────────────────────────
    adata.obs["original_cell_id"] = (
        adata.obs["sample_id"].astype(str) + "_" + adata.obs_names + "-1"
    )
    tap_mask = adata.obs["cell_type"] == "TAP"
    tap_df = adata.obs[tap_mask][[
        "original_cell_id", "sample_id", "velocity_label", "cellrank_P_OL"
    ]].copy()
    tap_df["ml_P_OL"] = tap_df["original_cell_id"].map(ml["prob_OL"])
    tap_df["bias_score"] = tap_df["original_cell_id"].map(ml["bias"])

    before = len(tap_df)
    tap_df = tap_df.dropna(subset=["cellrank_P_OL", "ml_P_OL", "bias_score"])
    print(f"  TAPs with all 3 predictions: {len(tap_df)}/{before}")

    # binary calls
    tap_df["cellrank_is_ol"] = (tap_df["cellrank_P_OL"] > 0.5).astype(int)
    tap_df["ml_is_ol"]       = (tap_df["ml_P_OL"] > 0.5).astype(int)
    tap_df["bias_is_ol"]     = (tap_df["bias_score"] > 0).astype(int)

    tap_df.to_csv("outputs/tap_fate_comparison.csv", index=False)
    print(f"  → outputs/tap_fate_comparison.csv  (n={len(tap_df)})")

    # ── Comparison stats ─────────────────────────────────────────────────────
    print(f"\n{'='*80}")
    print("DISTRIBUTION (n = {} TAPs)".format(len(tap_df)))
    print(f"{'='*80}")
    print(tap_df[["cellrank_P_OL", "ml_P_OL", "bias_score"]].describe().round(3).to_string())

    # Pairwise stats — build a summary CSV
    methods = [
        ("ml_P_OL", "ml_is_ol", "ML"),
        ("cellrank_P_OL", "cellrank_is_ol", "CellRank"),
        ("bias_score", "bias_is_ol", "Bias"),
    ]
    rows = []
    for i, (s1, b1, n1) in enumerate(methods):
        for j, (s2, b2, n2) in enumerate(methods):
            if i >= j:
                continue
            pearson_r, _ = pearsonr(tap_df[s1], tap_df[s2])
            spearman_r, _ = spearmanr(tap_df[s1], tap_df[s2])
            agreement = (tap_df[b1] == tap_df[b2]).mean()
            kappa = cohen_kappa(tap_df[b1], tap_df[b2])
            both_ol = ((tap_df[b1] == 1) & (tap_df[b2] == 1)).sum()
            m1_ol = tap_df[b1].sum()
            m2_ol = tap_df[b2].sum()
            pos_agreement_m1 = both_ol / m1_ol if m1_ol > 0 else 0
            pos_agreement_m2 = both_ol / m2_ol if m2_ol > 0 else 0

            rows.append({
                "method_1": n1,
                "method_2": n2,
                "n_TAPs": len(tap_df),
                "pearson_r": round(pearson_r, 3),
                "spearman_r": round(spearman_r, 3),
                "overall_agreement": round(agreement * 100, 1),
                "cohen_kappa": round(kappa, 3),
                f"pct_{n1}_OL_confirmed_by_{n2}": round(pos_agreement_m1 * 100, 1),
                f"pct_{n2}_OL_confirmed_by_{n1}": round(pos_agreement_m2 * 100, 1),
            })

    print(f"\n{'='*80}")
    print("PAIRWISE METHOD COMPARISON")
    print(f"{'='*80}")
    for r in rows:
        print(f"\n{r['method_1']} vs {r['method_2']}:")
        print(f"  Pearson r:          {r['pearson_r']}")
        print(f"  Spearman r:         {r['spearman_r']}")
        print(f"  Overall agreement:  {r['overall_agreement']}%")
        print(f"  Cohen's kappa:      {r['cohen_kappa']}  ({'good' if r['cohen_kappa']>0.6 else 'moderate' if r['cohen_kappa']>0.4 else 'fair' if r['cohen_kappa']>0.2 else 'slight'} agreement)")
        for k, v in r.items():
            if k.startswith("pct_") and "_confirmed_" in k:
                print(f"  {k}: {v}%")

    summary_df = pd.DataFrame([{
        "method_1": r["method_1"],
        "method_2": r["method_2"],
        "n_TAPs": r["n_TAPs"],
        "pearson_r": r["pearson_r"],
        "spearman_r": r["spearman_r"],
        "overall_agreement_pct": r["overall_agreement"],
        "cohen_kappa": r["cohen_kappa"],
    } for r in rows])
    summary_df.to_csv("outputs/tap_fate_method_summary.csv", index=False)
    print(f"\n  → outputs/tap_fate_method_summary.csv")

    # ── Per-condition OL-leaning fractions ───────────────────────────────────
    print(f"\n{'='*80}")
    print("OL-LEANING FRACTION PER CONDITION")
    print(f"{'='*80}")
    pivot_rows = []
    for method, col in [("CellRank", "cellrank_is_ol"),
                        ("ML", "ml_is_ol"),
                        ("Bias score", "bias_is_ol")]:
        row = {"method": method}
        for cond in sorted(tap_df["velocity_label"].unique()):
            n = (tap_df["velocity_label"] == cond).sum()
            n_ol = ((tap_df["velocity_label"] == cond) & (tap_df[col] == 1)).sum()
            row[cond] = f"{n_ol/n*100:.1f}% ({n_ol}/{n})"
        pivot_rows.append(row)

    pivot = pd.DataFrame(pivot_rows)
    print(pivot.to_string(index=False))
    pivot.to_csv("outputs/tap_fate_per_condition.csv", index=False)
    print(f"\n  → outputs/tap_fate_per_condition.csv")


if __name__ == "__main__":
    main()
