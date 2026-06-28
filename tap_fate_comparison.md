# TAP Fate Prediction — Method Comparison

Comparison of two methods for predicting OL vs NB fate in TAP cells from the V-SVZ.

---

## 1. Motivation

The project's central deliverable is a per-TAP fate prediction: "is this TAP heading toward OL or NB lineage?" Two methods produce candidate answers:

| Method | Approach | Input |
|---|---|---|
| **ML model (XGBoost)** | Supervised classifier trained on terminal cells | ~1,970 non-canonical HVGs (2000 HVGs − 31 excluded lineage genes from `MARKERS[POS+NEG]` ∪ `LINEAGE_GENES`) |
| **Hand-rolled bias score** | Arithmetic marker expression score difference | `(score_OPC + score_COP + score_OL)/3 - score_Neuroblast` calculated via Scanpy `score_genes` |

These methods serve as mutual references to evaluate how robustly we can capture lineage commitment. The hand-rolled bias score acts as an expression-based canonical baseline.

---

## 2. The Question Each Method Answers

These methods are different in the *signal* they read from the data:

| Method | What it asks | What it captures |
|---|---|---|
| **ML model** | "Does this TAP's overall pattern look like cells that become OL?" | Trajectory — learned boundary from real terminal cells, uses positive + negative signals + non-linear gene interactions across ~1,970 HVGs |
| **Bias Score** | "Does this TAP express more canonical OL than NB markers?" | Canonical marker enrichment — basic relative expression levels of 29 hand-picked marker genes |

ML uses **steady-state gene expression**. The bias score uses **canonical gene expression levels**. The signals themselves are biologically distinct, so the methods *can* disagree on individual cells. However, both methods anchor on the same upstream data structure (see Section 4), so they are not statistically independent.

---

## 3. Results — 10-Sample Run (3,597 TAPs)

### Per-cell agreement (OL Specific)

| Method | Total TAPs called OL | Percentage of total TAPs |
|---|:---:|:---:|
| **ML** | 366 | 10.2% |
| **Bias Score** | 298 | 8.3% |

### Per-condition counts
* ML probs are thresholded at `P(OL) > 0.5` (since OL and NB are the only absorbing states).
* Bias score is thresholded at `bias > 0` (positive leans OL, negative leans NB).

**OL-leaning TAPs (leaning OL / Total TAPs)**

| Method | CD1_Cntl_0wksRecov | CD1_Cntl_3wksRecov | CD1_CupRap_Rep1_3wksRecov | CD1_CupRap_Rep1_0wksRecov | CD1_CupRap_Rep2_0wksRecov | CD1_CupRap_rep2_3wksRecov | NesCre_Cntl_Rep1_3wksRecov | NesCre_Cntl_Rep2_3wksRecov | NesCre_CR_Rep1_3wksRecov | NesCre_CR_Rep2_3wksRecov |
|---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| **ML** | 10.2% (58/568) | 10.5% (78/746) | 11.2% (36/320) | 11.1% (10/90) | **37.7% (29/77)** | 9.8% (27/275) | 8.3% (45/541) | 8.6% (25/290) | 10.3% (31/302) | 7.0% (27/388) |
| **Bias** | 7.9% (45/568) | 6.7% (50/746) | 6.2% (20/320) | 6.7% (6/90) | **35.1% (27/77)** | 10.9% (30/275) | 6.1% (33/541) | 10.0% (29/290) | 13.6% (41/302) | 4.4% (17/388) |

### Key Observations:
1. **The 0-week Outlier**: Both expression-based methods independently flag the acute injury sample `CD1_CupRap_Rep2_0wksRecov` (no recovery) as a major outlier. The proportion of OL-leaning TAPs in this sample is roughly **3x to 5x higher** than in the control/recovery samples across both methods.
2. **CupRap vs. Control at 3 weeks**: For the recovery samples, Cup-Rap TAPs are **not** more OL-leaning than Control TAPs. Across both methods, the 3-week recovery conditions sit in a low baseline band. This reinforces the finding that Cup-Rap recovery does not skew the TAP population toward OL differentiation compared to controls.

---

## 4. Methodology Notes — Structural Sharing

Before reading these numbers as "independent validations," it's worth being explicit about what these methods share upstream:

1. **The `cell_type` labels**: ML trains on `cell_type` definitions; the Bias score uses marker panels that defined those cell types in the first place.
2. **The expression matrix**: Same QC filters, normalization, and Harmony batch integration.

**What is NOT shared**:
* **XGBoost (ML)**: Learns a non-linear decision boundary across 1,970 non-canonical highly variable genes.
* **Bias Score**: Uses the raw arithmetic difference of average expression scores across 29 canonical genes.

---

## 5. Output Files

- `compare_tap_fate_methods.py` — comparison script
- `outputs/tap_fate_comparison.csv` — per-cell columns (`ml_P_OL`, `bias`, `ml_ol`, `bias_ol`)
- `outputs/tap_fate_per_condition.csv` — per-condition summarized statistics

---

## 6. Reproducibility

```bash
uv sync  # Installs scipy<1.17, numpy>=2
uv run python compare_tap_fate_methods.py
```
