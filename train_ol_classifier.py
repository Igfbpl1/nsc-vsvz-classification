"""Train OL-lineage vs Neuroblast similarity classifier; predict OL commitment on TAPs.

Goal (per CLAUDE.md): given a cell's expression, score how committed it looks toward
the oligodendrocyte repair lineage. Feature importances = candidate trigger genes.

The model is trained on HVGs **minus the 29 known lineage markers** (Pdgfra, Cspg4,
Olig1, Olig2, Sox10, Cnp, Lhfpl3, Mbp, Mog, Plp1, Mobp, …). Forcing the model to find
non-marker signal makes the SHAP importance list a useful trigger-gene candidate set.

Train on cells with confident terminal labels:
  positive class: OPC + COP + OL
  negative class: Neuroblast
Held out for honest generalization AUC: one entire sample (GSM8253799).
Predicted on: all TAPs (the cells whose fate we actually care about).
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

from markers import MARKERS, OL_LINEAGE

ROOT = Path(__file__).parent
OUT = ROOT / "outputs"
MODELS = ROOT / "models"
for d in (OUT, MODELS):
    d.mkdir(parents=True, exist_ok=True)

HOLDOUT_SAMPLE = "GSM8253799"  # NesCre CR Rep2

POS_TYPES = list(OL_LINEAGE)         # OPC, COP, OL
NEG_TYPES = ["Neuroblast"]
TAP_TYPE = "TAP"

EXCLUDED_MARKERS = sorted({g for ct in POS_TYPES + NEG_TYPES for g in MARKERS[ct]})


def _to_dense(X) -> np.ndarray:
    return X.toarray() if hasattr(X, "toarray") else np.asarray(X)


def labels_for(adata: ad.AnnData) -> pd.Series:
    """1 = OL-lineage, 0 = Neuroblast, NaN = neither (excluded from train/eval)."""
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

    if HOLDOUT_SAMPLE not in adata.obs["sample_id"].unique().tolist():
        raise SystemExit(f"held-out sample {HOLDOUT_SAMPLE} not found in obs['sample_id']")

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
    print(f"\n=== features: {before} HVGs - {before - len(feat)} markers = {len(feat)} ===")

    # ------------------------------------------------------------------
    # Train / held-out split
    # ------------------------------------------------------------------
    y_all = labels_for(adata)
    is_eligible = y_all.notna()
    held_out_mask = (adata.obs["sample_id"] == HOLDOUT_SAMPLE).values
    train_mask = is_eligible.values & ~held_out_mask
    test_mask = is_eligible.values & held_out_mask

    print(f"  train cells: {train_mask.sum()} (POS={int(y_all[train_mask].sum())}, "
          f"NEG={int((y_all[train_mask] == 0).sum())})")
    print(f"  held-out cells: {test_mask.sum()} (POS={int(y_all[test_mask].sum())}, "
          f"NEG={int((y_all[test_mask] == 0).sum())})")

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
    # Fit
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
    print(f"  held-out ({HOLDOUT_SAMPLE}) AUC: {auc:.4f}, AP: {ap:.4f}")

    model_path = MODELS / "ol_classifier.json"
    clf.save_model(model_path)
    print(f"  model -> {model_path}")

    # ------------------------------------------------------------------
    # SHAP feature importance
    # ------------------------------------------------------------------
    print("  computing SHAP on a 2000-cell sample of the training set")
    rng = np.random.default_rng(42)
    n_shap = min(2000, X_train.shape[0])
    idx = rng.choice(X_train.shape[0], size=n_shap, replace=False)
    explainer = shap.TreeExplainer(clf)
    sv = explainer.shap_values(X_train[idx])
    shap_abs_mean = np.abs(sv).mean(axis=0)
    importance = pd.DataFrame({
        "gene": feat,
        "shap_mean_abs": shap_abs_mean,
        "gain": [clf.get_booster().get_score(importance_type="gain").get(f"f{i}", 0.0)
                 for i in range(len(feat))],
    }).sort_values("shap_mean_abs", ascending=False).reset_index(drop=True)

    importance.head(50).to_csv(OUT / "trigger_genes.csv", index=False)
    print(f"  top-50 features -> trigger_genes.csv")
    print(importance.head(15).to_string(index=False))

    # ------------------------------------------------------------------
    # Per-cell predictions on the full dataset
    # ------------------------------------------------------------------
    print("  predicting on full dataset (per-cell OL commitment)")
    X_dense_all = _to_dense(X_full)
    p_all = clf.predict_proba(X_dense_all)[:, 1]

    df = pd.DataFrame({
        "cell_id": adata.obs_names,
        "sample_id": adata.obs["sample_id"].values,
        "cell_type": adata.obs["cell_type"].values,
        "bias": adata.obs["bias"].values,
        "prob_OL": p_all,
    })
    df.to_csv(OUT / "ol_commitment.csv", index=False)
    print(f"  per-cell predictions -> ol_commitment.csv")

    tap_mask = adata.obs["cell_type"].values == TAP_TYPE
    rho, p_rho = spearmanr(p_all[tap_mask], adata.obs["bias"].values[tap_mask])
    print(f"  TAP sanity: Spearman(predicted_prob, bias) = {rho:.3f} (p={p_rho:.2e})")

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    summary = pd.DataFrame([{
        "auc": auc, "ap": ap, "spw": spw,
        "best_iter": int(clf.best_iteration), "val_auc": float(clf.best_score),
        "n_features": len(feat), "n_top_features": len(importance),
    }])
    summary.to_csv(OUT / "classifier_summary.csv", index=False)
    print("\nclassifier summary:")
    print(summary.to_string(index=False))
    print("\ndone")
