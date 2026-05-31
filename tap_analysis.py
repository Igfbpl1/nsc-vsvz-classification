from __future__ import annotations

from pathlib import Path

import pandas as pd
import scanpy as sc

from markers import MARKERS

ROOT = Path(__file__).parent
OUT = ROOT / "outputs"

ALL_MARKER_GENES = {g for genes in MARKERS.values() for g in genes}


def run_analysis() -> None:
    print("loading data")
    adata = sc.read_h5ad(OUT / "processed.h5ad")
    predictions = pd.read_csv(OUT / "ol_commitment.csv").set_index("cell_id")
    adata.obs["pred_p_OL"] = predictions.loc[adata.obs_names, "prob_OL"].values
    adata.obs["bias"] = predictions.loc[adata.obs_names, "bias"].values

    tap = adata[adata.obs["cell_type"] == "TAP"].copy()
    print(f"TAPs: {tap.n_obs}")
    tap_frame = tap.obs[["sample_id", "pred_p_OL", "bias"]]
    tap_frame["bias_is_ol"] = tap_frame["bias"].apply(lambda p: 1 if p > 0 else 0)
    tap_frame["prediction_is_ol"] = tap_frame["pred_p_OL"].apply(
        lambda p: 1 if p > 0.5 else 0
    )
    tap_frame.to_csv(f"{OUT}/out_of_sample_tap_comparison_test2.csv")
