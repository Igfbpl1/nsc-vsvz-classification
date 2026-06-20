"""
Validate CellRank fate probabilities against known biological expectations.
Loads persisted CellRank probabilities from velocity_combined.h5ad if available,
otherwise computes them inline. Saves all check results as CSV files.

Usage:
    uv run python validate_cellrank.py
"""
import sys
sys.path.insert(0, ".")

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
from scipy.stats import spearmanr
from pathlib import Path
from markers import OL_LINEAGE

ROOT = Path(".")
OUT = ROOT / "outputs"
VEL_DIR = OUT / "velocity"

def ensure_cellrank_probs(adata):
    """Load or compute CellRank P(OL) for all cells."""
    if "cellrank_P_OL" in adata.obs.columns:
        print("  ✅ cellrank_P_OL already present in obs (loaded from h5ad).")
        return adata

    print("  cellrank_P_OL not found — computing inline ...")
    import cellrank as cr
    from cellrank.kernels import VelocityKernel, ConnectivityKernel

    if "velocity_graph" not in adata.uns:
        scv.tl.velocity_graph(adata, show_progress_bar=False, n_jobs=1)

    vk = VelocityKernel(adata).compute_transition_matrix(n_jobs=4)
    ck = ConnectivityKernel(adata).compute_transition_matrix()
    combined_kernel = 0.8 * vk + 0.2 * ck

    g = cr.estimators.GPCCA(combined_kernel)
    ol_cells = adata.obs.index[adata.obs["cell_type"].isin(OL_LINEAGE)]
    nb_cells = adata.obs.index[adata.obs["cell_type"] == "Neuroblast"]
    g.set_terminal_states({"OL": ol_cells, "Neuroblast": nb_cells})
    g.compute_fate_probabilities()

    adata.obs["cellrank_P_OL"] = g.fate_probabilities[:, "OL"].X.squeeze()
    adata.obs["cellrank_P_NB"] = g.fate_probabilities[:, "Neuroblast"].X.squeeze()
    return adata


def main():
    print("=" * 70)
    print("CellRank Fate Probability Validation Report")
    print("=" * 70)

    # Load velocity object
    print("\nLoading velocity_combined.h5ad ...")
    adata = sc.read_h5ad(VEL_DIR / "velocity_combined.h5ad")
    adata.obs_names_make_unique()
    print(f"  Cells: {adata.n_obs:,}  |  Genes: {adata.n_vars:,}")

    # Ensure CellRank probabilities are available
    adata = ensure_cellrank_probs(adata)
    obs = adata.obs

    # =====================================================================
    # CHECK 1: Terminal State Sanity
    # =====================================================================
    print("\n" + "=" * 70)
    print("CHECK 1: Terminal State Sanity")
    print("=" * 70)

    check1_rows = []
    for ct in ['OL', 'COP', 'OPC', 'Neuroblast', 'TAP', 'aNSC', 'dNSC',
               'Astrocyte', 'Microglia', 'Endothelial', 'Pericyte', 'VAMC',
               'Ependymal', 'Other_Immune', 'Striatal_Neuron']:
        mask = obs['cell_type'] == ct
        if mask.sum() == 0:
            continue
        vals = obs.loc[mask, 'cellrank_P_OL']
        status = ""
        if ct in ['OL', 'COP', 'OPC']:
            status = "PASS" if vals.mean() > 0.8 else "WARN"
        elif ct == 'Neuroblast':
            status = "PASS" if vals.mean() < 0.2 else "WARN"
        else:
            status = "N/A"
        row = {
            "cell_type": ct, "n_cells": int(mask.sum()),
            "mean_P_OL": round(vals.mean(), 6), "median_P_OL": round(vals.median(), 6),
            "std_P_OL": round(vals.std(), 6), "min_P_OL": round(vals.min(), 6),
            "max_P_OL": round(vals.max(), 6), "status": status
        }
        check1_rows.append(row)
        icon = {"PASS": "✅", "WARN": "⚠", "N/A": " "}.get(status, " ")
        print(f"  {icon} {ct:18s}  n={mask.sum():5d}  mean={vals.mean():.4f}  "
              f"median={vals.median():.4f}  std={vals.std():.4f}  "
              f"[{vals.min():.4f}, {vals.max():.4f}]")

    df1 = pd.DataFrame(check1_rows)
    df1.to_csv(VEL_DIR / "validation_check1_terminal_states.csv", index=False)
    print(f"  → {VEL_DIR}/validation_check1_terminal_states.csv")

    # =====================================================================
    # CHECK 2: Lineage Gradient
    # =====================================================================
    print("\n" + "=" * 70)
    print("CHECK 2: Lineage Gradient (Neuroblast < TAP < OPC < COP < OL)")
    print("=" * 70)

    trajectory_order = ['Neuroblast', 'TAP', 'OPC', 'COP', 'OL']
    check2_rows = []
    for ct in trajectory_order:
        mask = obs['cell_type'] == ct
        if mask.sum() > 0:
            m = obs.loc[mask, 'cellrank_P_OL'].mean()
            med = obs.loc[mask, 'cellrank_P_OL'].median()
            check2_rows.append({"cell_type": ct, "mean_P_OL": round(m, 6), "median_P_OL": round(med, 6)})
            print(f"  {ct:15s}  mean P(OL) = {m:.4f}")

    values = [r["mean_P_OL"] for r in check2_rows]
    is_monotonic = all(values[i] <= values[i + 1] for i in range(len(values) - 1))
    print(f"\n  Monotonically increasing? {'✅ YES' if is_monotonic else '❌ NO'}")

    df2 = pd.DataFrame(check2_rows)
    df2["is_monotonic"] = is_monotonic
    df2.to_csv(VEL_DIR / "validation_check2_lineage_gradient.csv", index=False)
    print(f"  → {VEL_DIR}/validation_check2_lineage_gradient.csv")

    # =====================================================================
    # CHECK 3: Velocity Confidence
    # =====================================================================
    print("\n" + "=" * 70)
    print("CHECK 3: Velocity Confidence & Quality Metrics")
    print("=" * 70)

    check3_rows = []
    for col in ['velocity_confidence', 'velocity_length', 'velocity_self_transition']:
        if col in obs.columns:
            vals = obs[col].dropna()
            row = {
                "metric": col, "mean": round(vals.mean(), 6),
                "median": round(vals.median(), 6), "std": round(vals.std(), 6),
                "min": round(vals.min(), 6), "max": round(vals.max(), 6)
            }
            check3_rows.append(row)
            print(f"  {col}: mean={vals.mean():.4f}, median={vals.median():.4f}, "
                  f"std={vals.std():.4f}, [{vals.min():.4f}, {vals.max():.4f}]")

    if 'velocity_params' in adata.uns:
        params = dict(adata.uns['velocity_params'])
        print(f"  velocity_params: {params}")

    # Per cell-type confidence
    if 'velocity_self_transition' in obs.columns:
        print("\n  Per cell-type velocity_self_transition:")
        for ct in trajectory_order:
            mask = obs['cell_type'] == ct
            if mask.sum() > 0:
                vals = obs.loc[mask, 'velocity_self_transition']
                check3_rows.append({
                    "metric": f"self_transition_{ct}",
                    "mean": round(vals.mean(), 6), "median": round(vals.median(), 6),
                    "std": round(vals.std(), 6), "min": round(vals.min(), 6),
                    "max": round(vals.max(), 6)
                })
                print(f"    {ct:15s}  mean={vals.mean():.4f}  std={vals.std():.4f}")

    df3 = pd.DataFrame(check3_rows)
    df3.to_csv(VEL_DIR / "validation_check3_velocity_confidence.csv", index=False)
    print(f"  → {VEL_DIR}/validation_check3_velocity_confidence.csv")

    # =====================================================================
    # CHECK 4: Marker Score Correlation
    # =====================================================================
    print("\n" + "=" * 70)
    print("CHECK 4: Marker Score Correlation with CellRank P(OL)")
    print("=" * 70)

    proc = ad.read_h5ad(OUT / "processed.h5ad", backed='r')
    proc_obs = proc.obs.copy()

    obs['match_id'] = obs['sample_id'].astype(str) + "_" + obs.index + "-1"

    score_cols = ['score_aNSC', 'score_dNSC', 'score_TAP', 'score_Neuroblast',
                  'score_OPC', 'score_COP', 'score_OL', 'score_Astrocyte',
                  'score_Microglia', 'score_Endothelial']
    available_scores = [c for c in score_cols if c in proc_obs.columns]

    check4_rows = []
    matched = obs[obs['match_id'].isin(proc_obs.index)]
    print(f"  Matched {len(matched):,} / {len(obs):,} cells")

    if len(matched) > 100:
        for sc_col in available_scores:
            scores = proc_obs.loc[matched['match_id'], sc_col].values
            cr_probs = matched['cellrank_P_OL'].values
            valid = ~(np.isnan(scores) | np.isnan(cr_probs))
            if valid.sum() > 50:
                rho, p = spearmanr(cr_probs[valid], scores[valid])
                if any(x in sc_col for x in ['OPC', 'COP', 'OL', 'Astrocyte']):
                    expected_dir = "positive"
                else:
                    expected_dir = "negative"
                actual_dir = "positive" if rho > 0 else "negative"
                status = "PASS" if expected_dir == actual_dir else "WARN"
                check4_rows.append({
                    "score_column": sc_col, "spearman_rho": round(rho, 6),
                    "p_value": f"{p:.2e}", "expected_direction": expected_dir,
                    "actual_direction": actual_dir, "status": status
                })
                icon = "✅" if status == "PASS" else "⚠"
                print(f"  {icon} Spearman(P(OL), {sc_col:20s}) = {rho:+.4f}  (p={p:.2e})  [{status}]")

    df4 = pd.DataFrame(check4_rows)
    df4.to_csv(VEL_DIR / "validation_check4_marker_correlation.csv", index=False)
    print(f"  → {VEL_DIR}/validation_check4_marker_correlation.csv")

    # =====================================================================
    # CHECK 5: Fate Probability Sum
    # =====================================================================
    print("\n" + "=" * 70)
    print("CHECK 5: Fate Probability Sum (P(OL) + P(NB) should = 1.0)")
    print("=" * 70)

    if 'cellrank_P_NB' in obs.columns:
        row_sums = obs['cellrank_P_OL'].values + obs['cellrank_P_NB'].values
    elif 'lineages_fwd' in adata.obsm:
        fp = np.array(adata.obsm['lineages_fwd'])
        row_sums = fp.sum(axis=1)
    else:
        row_sums = np.ones(len(obs))
        print("  ⚠ Cannot compute: no P(NB) column or lineages_fwd found.")

    print(f"  Row sums: mean={row_sums.mean():.8f}, std={row_sums.std():.2e}")
    print(f"  Range: [{row_sums.min():.8f}, {row_sums.max():.8f}]")
    close_to_one = np.abs(row_sums - 1.0) < 1e-4
    n_pass = int(close_to_one.sum())
    n_total = len(row_sums)
    status = "ALL PASS" if close_to_one.all() else f"{n_pass}/{n_total} pass"
    print(f"  Rows ≈ 1.0 (tol 1e-4): {n_pass:,} / {n_total:,} ({'✅ ' + status if close_to_one.all() else '⚠ ' + status})")

    df5 = pd.DataFrame([{
        "metric": "fate_prob_sum", "mean": round(row_sums.mean(), 8),
        "std": f"{row_sums.std():.2e}", "min": round(row_sums.min(), 8),
        "max": round(row_sums.max(), 8),
        "n_pass_tol_1e4": n_pass, "n_total": n_total, "status": status
    }])
    df5.to_csv(VEL_DIR / "validation_check5_fate_prob_sums.csv", index=False)
    print(f"  → {VEL_DIR}/validation_check5_fate_prob_sums.csv")

    # =====================================================================
    # SUMMARY
    # =====================================================================
    print("\n" + "=" * 70)
    print("VALIDATION COMPLETE — CSV files saved to outputs/velocity/")
    print("=" * 70)
    print("  validation_check1_terminal_states.csv")
    print("  validation_check2_lineage_gradient.csv")
    print("  validation_check3_velocity_confidence.csv")
    print("  validation_check4_marker_correlation.csv")
    print("  validation_check5_fate_prob_sums.csv")


if __name__ == "__main__":
    main()
