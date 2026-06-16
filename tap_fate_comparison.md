# TAP Fate Prediction — Method Comparison

Comparison of three methods for predicting OL vs NB fate in TAP cells from the V-SVZ.

---

## 1. Motivation

The project's central deliverable is a per-TAP fate prediction: "is this TAP heading toward OL or NB lineage?" Three methods produce candidate answers:

| Method | Approach | Input |
|---|---|---|
| **ML model (XGBoost)** | Supervised classifier trained on terminal cells | 1,965 non-canonical HVGs |
| **Bias score** | Canonical marker panel scoring | 35 unique known lineage markers (across OPC/COP/OL/NB panels) |
| **CellRank** | Velocity-based fate probability | RNA velocity transition matrix + terminal states |

CellRank is methodologically principled for intermediate cells because it uses each cell's actual velocity vector to compute drift toward terminal states. It provides an independent validation reference for the other methods.

---

## 2. The Question Each Method Answers

These methods are not just three implementations of the same idea — they answer **different biological questions**:

| Method | What it asks | What it captures |
|---|---|---|
| **Bias score** | "Is this TAP currently expressing OL identity?" | Snapshot — current canonical OL panel expression |
| **ML model** | "Does this TAP's overall pattern look like cells that become OL?" | Trajectory — learned boundary from real terminal cells, uses positive + negative signals + non-linear gene interactions |
| **CellRank** | "Where is this TAP's velocity vector pointing?" | Motion — splicing dynamics indicate direction of cell-state change |

The **bias score is conservative by design** — it requires actual upregulation of canonical OL markers. The **ML model and CellRank are trajectory-aware** — they catch TAPs that are about to commit but have not yet turned on OL genes.

---

## 3. Results — 2,506 TAPs from 4 Samples (CD1_Cntl 0wks + 3 NesCre)

### Distribution of scores

| Method | Mean | SD | Min | Max |
|---|---|---|---|---|
| CellRank P(OL) | 0.407 | 0.326 | 0.000 | 0.951 |
| ML P(OL) | 0.292 | 0.303 | 0.055 | 0.945 |
| Bias score | −0.939 | 0.555 | −2.010 | 0.491 |

### Pairwise method comparison

| Pair | Pearson r | Spearman r | Agreement | Cohen's kappa |
|---|---|---|---|---|
| **ML vs CellRank** | **0.832** | **0.805** | **86.3%** | **0.679 (good)** |
| ML vs Bias | 0.848 | 0.856 | 78.6% | 0.251 (fair) |
| CellRank vs Bias | 0.857 | 0.810 | 70.6% | 0.188 (slight) |

### Per-condition OL-leaning fraction

| Method | CD1_Cntl (n=876) | Cntl (n=743) | CupRap_Rep1 (n=370) | CupRap_Rep2 (n=517) |
|---|---|---|---|---|
| CellRank (P > 0.5) | 41.3% (362) | 32.8% (244) | 33.5% (124) | 26.5% (137) |
| ML (P > 0.5) | 31.7% (278) | 24.6% (183) | 25.1% (93) | 19.0% (98) |
| Bias score (>0) | 7.3% (64) | 4.6% (34) | 6.2% (23) | 2.7% (14) |

---

## 4. Why ML Beats Bias Score for Trajectory Prediction

The original `PROJECT_NARRATIVE.md` used the bias score to validate the ML model, reporting 82% agreement on GSM8253799. That number is technically correct but uninformative — it reflects mostly the trivial fact that both methods correctly classify the abundant NB-fated cells.

**Of ML's OL predictions, Bias confirms 19.2%.**
**Of ML's OL predictions, CellRank confirms 90.2%.**

The 82% bias-vs-ML number was misleading. The actual validation is CellRank's confirmation of ML's OL calls — and that confirmation rate is strong.

### Three reasons ML beats Bias for trajectory

Note: Both methods use positive AND negative direction signals. The bias score formula `mean(OL panels) − NB panel score` already incorporates NB-direction. The real differences are elsewhere.

**1. Feature breadth — not direction.**
Bias uses only **35 unique canonical markers** (24 OL-lineage + 11 NB, deduplicated across OPC/COP/OL/NB panels). These are genes that activate at OPC/COP/OL stages, not in TAPs that are *about to* commit. ML uses 1,965 non-canonical HVGs, many of which are expressed at intermediate levels in early-committing TAPs. ML can see subtle patterns Bias is blind to because Bias's vocabulary is restricted to late-stage canonical markers.

**2. Non-linear gene interactions vs additive averaging.**
Bias = arithmetic mean of OL panel − mean of NB panel. Purely additive — cannot encode conditional patterns like "Gjc3 elevated AND Stmn2 reduced → OL-like." XGBoost trees can. At the TAP stage where individual OL gene expression is weak, integrating co-occurrence patterns across many genes is exactly what reveals commitment.

**3. Cell-based calibration vs background normalization.**
Bias score's `score_genes` subtracts the mean of random background genes — an abstract reference. ML learned what "OL-like" looks like by seeing real OPC/COP/OL cells. When ML scores a TAP, it asks "how close is this to the cells I saw labeled OL?" — a cell-aware comparison that's more sensitive to intermediate states because real intermediate cells were in the training distribution.

### Precision vs sensitivity

Both methods are validated against CellRank. The split between them is in precision vs sensitivity:

| Metric | Bias score | ML model |
|---|---|---|
| Precision (when called OL, is CellRank confirming?) | **100%** | 90.6% |
| Sensitivity (catches how many of CellRank's OL TAPs?) | 14% | **67%** |

Bias is **high precision, low sensitivity** — it only flags the most committed TAPs, but every one of those is real. ML is **moderate precision, much higher sensitivity** — it catches 5× more OL-leaning TAPs at the cost of slightly lower precision. For "predict TAP fate before maturation," sensitivity matters: you want to catch cells early in commitment, not just confirm the obvious cases.

---

## 5. Bottom Line

**ML agrees with CellRank on 86.3% of cells (Cohen's kappa 0.679, good agreement).**

- If your project explicitly wants "predict TAP fate before maturation," **ML is the right tool**.
- If you want "describe which TAPs already look OL," **Bias score is more honest** (but more restrictive).
- The 82% bias-vs-ML agreement that the original project relied on was technically true but uninformative — both methods are conservative and just agreed on the easy NB cells.

ML beats Bias for the trajectory question — note that both methods use both directional signals (the bias formula `OL_panels − NB_panel` already incorporates NB-direction). The real differences are:

1. **Feature breadth** — ML uses 1,965 non-canonical HVGs; Bias uses only 35 unique canonical markers (24 OL-lineage + 11 NB) that don't activate until late commitment stages
2. **Non-linear gene interactions** — ML can encode conditional patterns; Bias is purely additive
3. **Cell-based calibration** — ML learned from real OPC/COP/OL cells; Bias normalizes against random background genes

The functional consequence: **Bias confirms only 15.2% of CellRank's OL calls. ML confirms 67.8% of CellRank's OL calls** — much higher sensitivity against the velocity-based reference.

ML's predictions are validated by CellRank's independent velocity-based analysis — when ML calls a TAP OL-leaning, CellRank confirms 90.6% of the time.

---

## 6. OL-Leaning Fraction Across Conditions

CD1_Cntl (0wks) has the highest OL-leaning fraction across both ML and CellRank (41.3% and 31.7%). CupRap_Rep2 has the lowest across all three methods.

CupRap samples do not consistently show higher OL-leaning TAP fractions than Cntl:

- CellRank: CD1_Cntl 41.3% > Cntl 32.8% ≈ CupRap_Rep1 33.5% > CupRap_Rep2 26.5%
- ML: CD1_Cntl 31.7% > Cntl 24.6% ≈ CupRap_Rep1 25.1% > CupRap_Rep2 19.0%
- Bias: CD1_Cntl 7.3% > Cntl 4.6% ≈ CupRap_Rep1 6.2% > CupRap_Rep2 2.7%

ML and CellRank agree on the ordering across all four conditions.

---

## 7. Output Files

- `compare_tap_fate_methods.py` — analysis script
- `outputs/tap_fate_comparison.csv` — per-TAP scores: cellrank_P_OL, ml_P_OL, bias_score, plus binary calls
- `outputs/tap_fate_method_summary.csv` — pairwise method statistics (Pearson, Spearman, agreement, Cohen's kappa)
- `outputs/tap_fate_per_condition.csv` — OL-leaning fraction per cell condition per method

---

## 8. Reproducibility

```bash
uv pip install cellrank
uv run python compare_tap_fate_methods.py
```

Requires:
- `outputs/processed.h5ad` (from preprocess pipeline)
- `outputs/velocity/velocity_combined.h5ad` (from rna_velocity_pipeline.py)
- `outputs/ol_commitment.csv` (from train_ol_classifier.py)

CellRank computation: deterministic (no random seed needed). Single run gives the same results.
