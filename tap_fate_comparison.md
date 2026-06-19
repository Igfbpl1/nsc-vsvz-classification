# TAP Fate Prediction — Method Comparison

Comparison of two methods for predicting OL vs NB fate in TAP cells from the V-SVZ.

---

## 1. Motivation

The project's central deliverable is a per-TAP fate prediction: "is this TAP heading toward OL or NB lineage?" Two methods produce candidate answers:

| Method | Approach | Input |
|---|---|---|
| **ML model (XGBoost)** | Supervised classifier trained on terminal cells | 1,965 non-canonical HVGs |
| **CellRank** | Velocity-based fate probability | RNA velocity transition matrix + terminal states |

CellRank is methodologically principled for intermediate cells because it uses each cell's actual velocity vector to compute drift toward terminal states. It provides an independent validation reference for the XGBoost classifier.

A third method — canonical marker bias score — was used in earlier iterations of the project but is excluded here because (a) it is conservative by design and just confirms obvious cases, and (b) its low Cohen's kappa against both other methods was uninformative. The bias score still exists in `outputs/ol_commitment.csv` for completeness.

---

## 2. The Question Each Method Answers

These methods are not just two implementations of the same idea — they answer **different biological questions** with different signals:

| Method | What it asks | What it captures |
|---|---|---|
| **ML model** | "Does this TAP's overall pattern look like cells that become OL?" | Trajectory — learned boundary from real terminal cells, uses positive + negative signals + non-linear gene interactions across 1,965 HVGs |
| **CellRank** | "Where is this TAP's velocity vector pointing?" | Motion — splicing dynamics indicate direction of cell-state change |

ML uses **gene expression patterns**. CellRank uses **splicing kinetics**. The two signals are biologically independent — agreement between them is therefore non-trivial.

---

## 3. Results — 7-Sample Run (4,268 TAPs)

### Per-cell agreement

| Statistic | Value |
|---|---|
| n TAPs with both predictions | 4,268 |
| Pearson r | **0.835** (p ≈ 0) |
| Spearman r | **0.855** (p ≈ 0) |

Two methods using independent biological signals correlate at r > 0.8 per cell. This is the strongest single statistic supporting the ML classifier's per-cell scoring.

### Per-condition counts (per-cell threshold P > 0.5)

A cell is OL-leaning if its individual P(OL) > 0.5, NB-leaning otherwise. Since OL + NB are the only terminal states, every TAP falls into one of the two.

**OL-leaning TAPs**

| Method | CD1_Cntl (876) | CD1_Cntl_3wks (1129) | CD1_CupRap (463) | CD1_CupRap_0wks_Rep2 (172) | Cntl (742) | CupRap_Rep1 (370) | CupRap_Rep2 (516) |
|---|---:|---:|---:|---:|---:|---:|---:|
| CellRank | 36.3% (318) | 31.6% (357) | 27.2% (126) | **67.4% (116)** | 29.9% (222) | 28.9% (107) | 20.3% (105) |
| ML | 31.7% (278) | 33.0% (373) | 23.5% (109) | **66.9% (115)** | 24.7% (183) | 25.1% (93) | 18.8% (97) |

**NB-leaning TAPs**

| Method | CD1_Cntl (876) | CD1_Cntl_3wks (1129) | CD1_CupRap (463) | CD1_CupRap_0wks_Rep2 (172) | Cntl (742) | CupRap_Rep1 (370) | CupRap_Rep2 (516) |
|---|---:|---:|---:|---:|---:|---:|---:|
| CellRank | 63.7% (558) | 68.4% (772) | 72.8% (337) | **32.6% (56)** | 70.1% (520) | 71.1% (263) | 79.7% (411) |
| ML | 68.3% (598) | 67.0% (756) | 76.5% (354) | **33.1% (57)** | 75.3% (559) | 74.9% (277) | 81.2% (419) |

`CD1_CupRap_0wks_Rep2` (n=172) is the only no-recovery sample and a large outlier on both methods. At this timepoint, rapamycin has blocked OPC→OL maturation, so OPCs are accumulating (paper Fig 3). TAPs in this sample may be projecting toward the stalled OPC pool — but n=172 is small; interpret with caution.

### Key direction

CupRap TAPs are **not** more OL-leaning than Control TAPs. By both CellRank and ML, the share of OL-leaning TAPs in CupRap conditions is at or below corresponding Controls:

- CellRank ordering (excluding outlier): CD1_Cntl (36.3%) > CD1_Cntl_3wks (31.6%) > Cntl (29.9%) ≈ CupRap_Rep1 (28.9%) ≈ CD1_CupRap (27.2%) > CupRap_Rep2 (20.3%)
- ML ordering: CD1_Cntl_3wks (33.0%) ≈ CD1_Cntl (31.7%) > CupRap_Rep1 (25.1%) ≈ Cntl (24.7%) > CD1_CupRap (23.5%) > CupRap_Rep2 (18.8%)

ML and CellRank agree on the ordering for most conditions. CupRap_Rep2 (NesCre, 3wks) is consistently lowest by both methods.

---

## 4. Methodology Notes on CellRank

CellRank's terminal states are defined from `cell_type` labels to match the ML classifier's positive class. The `OL` terminal group is the full OL lineage — cells labeled `OL`, `OPC`, or `COP` (i.e. `OL_LINEAGE` from `markers.py`); the `Neuroblast` group is cells labeled `Neuroblast`; every other cell (TAPs, NSCs, Astrocyte, etc.) is transient. For each transient cell, `cellrank_P_OL` is the probability that a random walk following the velocity transition matrix is absorbed into the OL-lineage group.

Bundling OPC/COP into the OL terminal group keeps the CellRank vs ML comparison apples-to-apples: both methods now report `P(OL lineage)` rather than ML reporting "OL-or-OPC-or-COP" and CellRank reporting strict "OL only". An OL-only terminal definition was tested and gives near-identical results in this dataset (Pearson r = 0.835 either way; mean P shifts by ≤0.006 per condition, OL counts shift by 1–2 cells per condition), because the OPC→COP→OL velocity flux is functional enough that strict-OL absorbing already captures the lineage mass via forward flow.

This approach inherits whatever errors exist in the Leiden clustering. We accept this because:

1. The canonical marker dotplot (`outputs/velocity/dotplot__marker_dotplot.png`) shows that the OL/OPC/COP and Neuroblast cluster labels are well-supported — `Mbp`/`Mog`/`Plp1` tight in OL, `Pdgfra`/`Cspg4` tight in OPC, `Dcx`/`Tubb3` tight in Neuroblast.
2. An alternative using GPCCA's `compute_macrostates()` + `predict_terminal_states()` to discover terminal cells from velocity topology was tested. It requires dense Schur decomposition on a 36k × 36k matrix (~30 min CPU on M2 Pro without SLEPc), and given the clean clustering, the methodological gain is small.

The CellRank transition matrix combines the velocity kernel (80%) and connectivity kernel (20%):

```
T_vel = scv.utils.get_transition_matrix(adata)
combined_kernel = 0.8 * PrecomputedKernel(T_vel) + 0.2 * ConnectivityKernel
```

(`PrecomputedKernel` is used in place of `VelocityKernel(adata)` to side-step a numpy-2 incompatibility in cellrank 2.0.7; the transition matrix it wraps is what `VelocityKernel.compute_transition_matrix` would have produced.)

---

## 5. Bottom Line

**ML and CellRank agree on per-cell P(OL lineage) at Pearson r = 0.835 / Spearman r = 0.855 across 4,268 TAPs.**

Two methods built on independent biological signals (gene expression vs splicing kinetics) converge on the same per-cell ranking. This is strong cross-validation for the XGBoost classifier's per-TAP fate score.

By both methods, **CupRap does not enrich OL fate at the TAP stage**. The share of OL-leaning TAPs in CupRap conditions is at or below corresponding Controls. This is consistent with the published claim (Willis et al. Figure 2I, r = 0.995 for TAP transcriptome Cntl vs CupRap) that TAP transcriptional identity is unchanged under CupRap, and adds a velocity-based confirmation that the velocity vectors at the TAP stage are not pushing toward OL fate either.

The OL commitment signal — both transcriptionally and via velocity dynamics — appears post-TAP, at the OPC/COP transition (see velocity driver rankings in `outputs/velocity/velocity_drivers_*.csv`).

---

## 6. Output Files

- `compare_tap_fate_methods.py` — analysis script
- `outputs/tap_fate_comparison.csv` — per-TAP `cellrank_P_OL`, `ml_P_OL`, and binary OL calls (threshold > 0.5)
- `outputs/tap_fate_per_condition.csv` — per-condition OL/NB counts for both methods (CSV includes continuous mean P(OL) / P(NB) columns as well)

---

## 7. Reproducibility

```bash
uv pip install cellrank
uv run python compare_tap_fate_methods.py
```

Requires:
- `outputs/velocity/velocity_combined.h5ad` (from `rna_velocity_pipeline.py`)
- `outputs/ol_commitment.csv` (from `train_ol_classifier.py`)

CellRank computation: deterministic for a given input. Single run gives the same results.
