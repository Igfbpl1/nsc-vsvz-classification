from __future__ import annotations

from pathlib import Path

import pandas as pd
import scanpy as sc

from markers import MARKERS

ROOT = Path(__file__).parent
OUT = ROOT / "outputs"

ALL_MARKER_GENES = {g for genes in MARKERS.values() for g in genes}

# GSM accession -> condition label (matches velocity_build.py velocity_label naming)
CONDITION = {
    "GSM8253792": "CD1_Cntl_0wksRecov",
    "GSM8253793": "CD1_Cntl_3wksRecov",
    "GSM8253794": "CD1_CupRap_Rep1_3wksRecov",
    "GSM8253795": "CD1_CupRap_rep2_3wksRecov",
    "GSM8253796": "NesCre_Cntl_Rep1_3wksRecov",
    "GSM8253797": "NesCre_Cntl_Rep2_3wksRecov",
    "GSM8253798": "NesCre_CR_Rep1_3wksRecov",
    "GSM8253799": "NesCre_CR_Rep2_3wksRecov",
    "GSM8647352": "CD1_CupRap_Rep1_0wksRecov",
    "GSM8647353": "CD1_CupRap_Rep2_0wksRecov",
}

# Order conditions appear in the per-condition agreement view.
CONDITION_ORDER = [
    "CD1_Cntl_0wksRecov", "CD1_Cntl_3wksRecov", "CD1_CupRap_Rep1_0wksRecov", 
    "CD1_CupRap_Rep2_0wksRecov", "CD1_CupRap_Rep1_3wksRecov", "CD1_CupRap_rep2_3wksRecov", 
    "NesCre_Cntl_Rep1_3wksRecov", "NesCre_Cntl_Rep2_3wksRecov", 
    "NesCre_CR_Rep1_3wksRecov", "NesCre_CR_Rep2_3wksRecov",
]


def write_agreement_view(tap_frame: pd.DataFrame) -> None:
    """Per-condition bias-vs-prediction agreement breakdown (2x2 quadrants)."""
    b, p = tap_frame["bias_is_ol"], tap_frame["prediction_is_ol"]
    df = tap_frame.assign(
        _nn=(b == 0) & (p == 0),   # both NB
        _oo=(b == 1) & (p == 1),   # both OL
        _bo=(b == 1) & (p == 0),   # bias-OL only
        _po=(b == 0) & (p == 1),   # pred-OL only
    )
    g = df.groupby("condition", observed=True)
    view = pd.DataFrame({
        "n": g.size(),
        "both NB (0/0)": g["_nn"].sum(),
        "both OL (1/1)": g["_oo"].sum(),
        "bias-OL only (1/0)": g["_bo"].sum(),
        "pred-OL only (0/1)": g["_po"].sum(),
    }).reindex([c for c in CONDITION_ORDER if c in g.groups])
    view["agree %"] = (
        (view["both NB (0/0)"] + view["both OL (1/1)"]) / view["n"] * 100
    ).round(1)

    total = view.drop(columns="agree %").sum()
    total["agree %"] = round(
        (view["both NB (0/0)"].sum() + view["both OL (1/1)"].sum())
        / view["n"].sum() * 100, 1
    )
    view = pd.concat([view, pd.DataFrame([total], index=["Grand Total"])])
    view.index.name = "condition"
    view.to_csv(f"{OUT}/test2_agreement_view.csv")


def write_ol_rate_comparison(tap_frame: pd.DataFrame) -> None:
    """Per-condition OL-leaning *rate* for each method (marginal proportions).

    Compares how many TAPs each method calls OL-leaning, as a rate, with the
    percentage-point difference between them. Deliberately omits a single
    "agreement %" because OL-leaning TAPs are rare (~8-10%), so 100 - diff
    always reads near-perfect and oversells the match; the two rates plus their
    pp difference are the honest comparison.
    """
    g = tap_frame.groupby("condition", observed=True)
    v = pd.DataFrame({
        "n": g.size(),
        "bias_OL": g["bias_is_ol"].sum(),
        "pred_OL": g["prediction_is_ol"].sum(),
    }).reindex([c for c in CONDITION_ORDER if c in g.groups])

    total = v.sum()
    total.name = "Grand Total"
    v = pd.concat([v, total.to_frame().T])

    v["bias_OL_%"] = (v["bias_OL"] / v["n"] * 100).round(1)
    v["pred_OL_%"] = (v["pred_OL"] / v["n"] * 100).round(1)
    v["diff_pp"] = (v["pred_OL_%"] - v["bias_OL_%"]).abs().round(1)
    v = v[["n", "bias_OL", "bias_OL_%", "pred_OL", "pred_OL_%", "diff_pp"]]
    v.index.name = "condition"
    v.to_csv(f"{OUT}/test2_ol_rate_comparison.csv")


def run_analysis() -> None:
    print("loading data")
    adata = sc.read_h5ad(OUT / "processed.h5ad")
    predictions = pd.read_csv(OUT / "ol_commitment.csv").set_index("cell_id")
    adata.obs["pred_p_OL"] = predictions.loc[adata.obs_names, "prob_OL"].values
    adata.obs["bias"] = predictions.loc[adata.obs_names, "bias"].values

    tap = adata[adata.obs["cell_type"] == "TAP"].copy()
    print(f"TAPs: {tap.n_obs}")
    tap_frame = tap.obs[["sample_id", "pred_p_OL", "bias"]].copy()
    tap_frame.insert(1, "condition", tap_frame["sample_id"].map(CONDITION))
    tap_frame["bias_is_ol"] = tap_frame["bias"].apply(lambda p: 1 if p > 0 else 0)
    tap_frame["prediction_is_ol"] = tap_frame["pred_p_OL"].apply(
        lambda p: 1 if p > 0.5 else 0
    )
    tap_frame.to_csv(f"{OUT}/out_of_sample_tap_comparison_test2.csv")
    write_agreement_view(tap_frame)
    write_ol_rate_comparison(tap_frame)
