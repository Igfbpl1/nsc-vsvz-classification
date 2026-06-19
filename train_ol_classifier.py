"""Train OL-lineage vs Neuroblast similarity classifier; predict OL commitment on TAPs.

Goal: given a cell's expression, score how committed it looks toward
the oligodendrocyte repair lineage. The features of the model are the candidate trigger genes.

The model is trained on HVGs **minus the lineage markers**. The exclusion set is
built from two sources in markers.py:
  * MARKERS[POS_TYPES + NEG_TYPES] — subtype-specific markers (Pdgfra, Cspg4,
    Mog, Mag, Opalin, Mbp, Plp1, Dlx1, Dlx2, …)
  * LINEAGE_GENES (OL_lineage + Neuroblast_lineage) — pan-lineage TFs and
    structural genes that identify lineage commitment but not subtype
    (Olig1, Olig2, Sox10, Cnp, Lhfpl3, Mobp, Tubb3, Cd24a).
Forcing the model to find non-marker signal makes the SHAP importance list a
useful trigger-gene candidate set rather than a re-discovery of canonical
oligodendrocyte / neuroblast genes.

Train on cells with confident terminal labels:
  positive class: OPC + COP + OL
  negative class: Neuroblast
Predicted on: all TAPs (the cells whose fate we actually care about because they are undecided).
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import shap
import xgboost as xgb
from scipy.stats import spearmanr
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import train_test_split

from markers import MARKERS, OL_LINEAGE, LINEAGE_GENES

ROOT = Path(__file__).parent
OUT = ROOT / "outputs"
MODELS = ROOT / "models"
for d in (OUT, MODELS):
    d.mkdir(parents=True, exist_ok=True)

HOLDOUT_SAMPLES = [
    "GSM8253796",
    "GSM8253797",
    "GSM8253798",
    "GSM8253799",
]  # NesCre CR Rep2/Rep1

POS_TYPES = list(OL_LINEAGE)  # OPC, COP, OL
NEG_TYPES = ["Neuroblast"]
TAP_TYPE = "TAP"

EXCLUDED_MARKERS = sorted(
    {g for ct in POS_TYPES + NEG_TYPES for g in MARKERS[ct]}
    | set(LINEAGE_GENES["OL_lineage"])
    | set(LINEAGE_GENES["Neuroblast_lineage"])
)


def _to_dense(X) -> np.ndarray:
    return X.toarray() if hasattr(X, "toarray") else np.asarray(X)


def labels_for(adata: ad.AnnData) -> pd.Series:
    """1 = OL-lineage, 0 = Neuroblast, NaN = neither (excluded from train/evaluation)."""
    y = pd.Series(np.nan, index=adata.obs_names, dtype=float)
    y[adata.obs["cell_type"].isin(POS_TYPES)] = 1.0
    y[adata.obs["cell_type"].isin(NEG_TYPES)] = 0.0
    return y


def run_classifier() -> None:
    print("loading processed.h5ad")
    adata = sc.read_h5ad(OUT / "processed.h5ad")
    print(f"  {adata.n_obs} cells x {adata.n_vars} genes")

    counts = adata.obs["cell_type"].value_counts()
    print("  cell_type counts:")
    for k, v in counts.items():
        print(f"    {k}: {v}")

    if not all(s in adata.obs["sample_id"].unique().tolist() for s in HOLDOUT_SAMPLES):
        raise SystemExit(
            f"held-out samples {', '.join(HOLDOUT_SAMPLES)} not found in obs['sample_id']"
        )

    bias = (
        adata.obs["score_OPC"] + adata.obs["score_COP"] + adata.obs["score_OL"]
    ) / 3 - adata.obs["score_Neuroblast"]
    adata.obs["bias"] = bias.values

    # ------------------------------------------------------------------
    # Build feature set: HVGs minus 29 known lineage markers
    # ------------------------------------------------------------------
    hvg_mask = adata.var["highly_variable"].values
    feat = adata.var_names[hvg_mask].tolist()
    before = len(feat)
    feat = [g for g in feat if g not in EXCLUDED_MARKERS]
    print(
        f"\n=== features: {before} HVGs - {before - len(feat)} markers = {len(feat)} ==="
    )

    # ------------------------------------------------------------------
    # Model Training / held-out split
    # ------------------------------------------------------------------
    y_all = labels_for(adata)
    is_eligible = y_all.notna()
    held_out_mask = adata.obs["sample_id"].isin(HOLDOUT_SAMPLES).values
    train_mask = is_eligible.values & ~held_out_mask
    test_mask = is_eligible.values & held_out_mask

    print(
        f"  train cells: {train_mask.sum()} (POS={int(y_all[train_mask].sum())}, "
        f"NEG={int((y_all[train_mask] == 0).sum())})"
    )
    print(
        f"  held-out cells: {test_mask.sum()} (POS={int(y_all[test_mask].sum())}, "
        f"NEG={int((y_all[test_mask] == 0).sum())})"
    )

    X_full = adata[:, feat].X
    X_train = _to_dense(X_full[train_mask])
    X_test = _to_dense(X_full[test_mask])
    y_train = y_all[train_mask].values.astype(int)
    y_test = y_all[test_mask].values.astype(int)

    X_tr, X_val, y_tr, y_val = train_test_split(
        X_train, y_train, test_size=0.2, stratify=y_train, random_state=42
    )

    pos = int(y_tr.sum())
    neg = int(len(y_tr) - pos)
    spw = neg / max(pos, 1)
    print(f"  scale_pos_weight = {spw:.3f}")

    # ------------------------------------------------------------------
    # Model Fitting
    # ------------------------------------------------------------------
    clf = xgb.XGBClassifier(
        n_estimators=500,
        max_depth=6,
        learning_rate=0.05,
        subsample=0.8,
        colsample_bytree=0.6,
        reg_lambda=1.0,
        scale_pos_weight=spw,
        objective="binary:logistic",
        eval_metric="auc",
        early_stopping_rounds=20,
        n_jobs=-1,
        random_state=42,
        tree_method="hist",
    )
    clf.fit(X_tr, y_tr, eval_set=[(X_val, y_val)], verbose=False)
    print(f"  best_iteration: {clf.best_iteration}, val AUC: {clf.best_score:.4f}")

    p_test = clf.predict_proba(X_test)[:, 1]
    auc = roc_auc_score(y_test, p_test)
    ap = average_precision_score(y_test, p_test)
    print(f"  held-out ({', '.join(HOLDOUT_SAMPLES)}) AUC: {auc:.4f}, AP: {ap:.4f}")
    y_test_frame = pd.DataFrame(y_test, columns=["fate_OL_lineage"])
    p_test_frame = pd.DataFrame(p_test, columns=["p_ol_nb"])
    out_of_sample_comparison = y_test_frame.join(
        p_test_frame,
    )
    out_of_sample_comparison["prediction_is_ol"] = out_of_sample_comparison[
        "p_ol_nb"
    ].apply(lambda p: 1 if p > 0.5 else 0)
    out_of_sample_comparison.to_csv(f"{OUT}/out_of_sample_comparison.csv")

    # ------------------------------------------------------------------
    # SHAP feature importance
    # ------------------------------------------------------------------
    print("  computing SHAP on held-out test set")
    rng = np.random.default_rng(42)
    n_shap = min(2000, X_test.shape[0])
    idx = rng.choice(X_test.shape[0], size=n_shap, replace=False)
    explainer = shap.TreeExplainer(clf)
    sv = explainer.shap_values(X_test[idx])
    shap_abs_mean = np.abs(sv).mean(axis=0)
    importance = (
        pd.DataFrame(
            {
                "gene": feat,
                "shap_mean_abs": shap_abs_mean,
                "gain": [
                    clf.get_booster()
                    .get_score(importance_type="gain")
                    .get(f"f{i}", 0.0)
                    for i in range(len(feat))
                ],
            }
        )
        .sort_values("shap_mean_abs", ascending=False)
        .reset_index(drop=True)
    )

    top = importance.head(50).copy()

    ol_mask = adata.obs["cell_type"].isin(POS_TYPES)
    nb_mask = adata.obs["cell_type"].isin(NEG_TYPES)
    directions = []
    for gene in top["gene"]:
        if gene not in adata.var_names:
            directions.append("NOT_FOUND")
            continue
        import scipy.sparse as _sp
        x = adata[:, gene].X
        if _sp.issparse(x):
            x = x.toarray().flatten()
        ol_mean = float(x[ol_mask].mean())
        nb_mean = float(x[nb_mask].mean())
        ratio = ol_mean / nb_mean if nb_mean > 0 else float("inf")
        if ratio > 1.5:
            directions.append("POSITIVE_OL")
        elif nb_mean > 0 and (1 / ratio) > 1.5:
            directions.append("NEGATIVE_OL")
        else:
            directions.append("AMBIGUOUS")
    top["direction"] = directions

    top.to_csv(OUT / "trigger_genes.csv", index=False)
    print("  top-50 features -> trigger_genes.csv")
    print(importance.head(15).to_string(index=False))

    # ------------------------------------------------------------------
    # Per-cell predictions on the full dataset
    # ------------------------------------------------------------------
    print("  predicting on full dataset (per-cell OL commitment)")
    X_dense_all = _to_dense(X_full)
    p_all = clf.predict_proba(X_dense_all)[:, 1]

    df = pd.DataFrame(
        {
            "cell_id": adata.obs_names,
            "sample_id": adata.obs["sample_id"].values,
            "cell_type": adata.obs["cell_type"].values,
            "bias": adata.obs["bias"].values,
            "prob_OL": p_all,
        }
    )
    df.to_csv(OUT / "ol_commitment.csv", index=False)
    print("  per-cell predictions -> ol_commitment.csv")

    tap_mask = adata.obs["cell_type"].values == TAP_TYPE
    rho, p_rho = spearmanr(p_all[tap_mask], adata.obs["bias"].values[tap_mask])
    print(f"  TAP sanity: Spearman(predicted_prob, bias) = {rho:.3f} (p={p_rho:.2e})")

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    summary = pd.DataFrame(
        [
            {
                "auc": auc,
                "ap": ap,
                "spw": spw,
                "best_iter": int(clf.best_iteration),
                "val_auc": float(clf.best_score),
                "n_features": len(feat),
                "n_top_features": len(importance),
            }
        ]
    )
    print("\nclassifier summary:")
    print(summary.to_string(index=False))
    print("\ndone")
