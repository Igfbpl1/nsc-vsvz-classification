# TAP Fate Prediction — Method Comparison

Comparison of two methods for predicting OL vs NB fate in TAP cells from the V-SVZ.

---

## 1. Motivation

The project's central deliverable is a per-TAP fate prediction: "is this TAP heading toward OL or NB lineage?" Two methods produce candidate answers:

| Method | Approach | Input |
|---|---|---|
| **ML model (XGBoost)** | Supervised classifier trained on terminal cells | ~1,970 non-canonical HVGs (2000 HVGs − 31 excluded lineage genes from `MARKERS[POS+NEG]` ∪ `LINEAGE_GENES`) |
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

## 3. Results — 7-Sample Run (2,944 TAPs)

### Per-cell agreement

| Statistic | Value |
|---|---|
| n TAPs with both predictions | 2,944 |
| Pearson r | **0.828** (p ≈ 0) |
| Spearman r | **0.838** (p ≈ 0) |

Two methods using independent biological signals correlate at r > 0.8 per cell. This is the strongest single statistic supporting the ML classifier's per-cell scoring.

### Per-condition counts (per-cell threshold P > 0.5)

A cell is OL-leaning if its individual P(OL) > 0.5, NB-leaning otherwise. Since OL + NB are the only terminal states, every TAP falls into one of the two.

**OL-leaning TAPs**

| Method | CD1_Cntl (568) | CD1_Cntl_3wks (746) | CD1_CupRap (320) | CD1_CupRap_0wks_Rep2 (77) | Cntl (542) | CupRap_Rep1 (302) | CupRap_Rep2 (389) |
|---|---:|---:|---:|---:|---:|---:|---:|
| CellRank | 22.7% (129) | 18.1% (135) | 21.3% (68) | **53.2% (41)** | 18.6% (101) | 24.5% (74) | 16.5% (64) |
| ML | 10.2% (58) | 10.5% (78) | 11.3% (36) | **37.7% (29)** | 8.3% (45) | 10.3% (31) | 6.9% (27) |

**NB-leaning TAPs**

| Method | CD1_Cntl (568) | CD1_Cntl_3wks (746) | CD1_CupRap (320) | CD1_CupRap_0wks_Rep2 (77) | Cntl (542) | CupRap_Rep1 (302) | CupRap_Rep2 (389) |
|---|---:|---:|---:|---:|---:|---:|---:|
| CellRank | 77.3% (439) | 81.9% (611) | 78.8% (252) | **46.8% (36)** | 81.4% (441) | 75.5% (228) | 83.5% (325) |
| ML | 89.8% (510) | 89.5% (668) | 88.8% (284) | **62.3% (48)** | 91.7% (497) | 89.7% (271) | 93.1% (362) |

`CD1_CupRap_0wks_Rep2` (n=77) is the only no-recovery sample and a large outlier on both methods (the only condition where OL-leaning approaches ~50%). At this timepoint, rapamycin has blocked OPC→OL maturation, so OPCs are accumulating (paper Fig 3). TAPs in this sample may be projecting toward the stalled OPC pool — but n=77 is small; interpret with caution.

### Key direction

CupRap TAPs are **not** more OL-leaning than Control TAPs. Across both methods, every 3wks condition (Control and CupRap alike) sits in a tight low band — CellRank ~16–25% OL, ML ~7–11% OL — with CupRap *not* above Control:

- CellRank OL-leaning (excluding outlier): CupRap_Rep1 (24.5%) ≈ CD1_Cntl (22.7%) ≈ CD1_CupRap (21.3%) ≈ Cntl (18.6%) ≈ CD1_Cntl_3wks (18.1%) ≈ CupRap_Rep2 (16.5%)
- ML OL-leaning: CD1_CupRap (11.3%) ≈ CD1_Cntl_3wks (10.5%) ≈ CupRap_Rep1 (10.3%) ≈ CD1_Cntl (10.2%) ≈ Cntl (8.3%) ≈ CupRap_Rep2 (6.9%)

The conditions are essentially indistinguishable at the TAP stage, which is the point: CupRap does not shift TAP fate toward OL. This matches the paper's r = 0.995 TAP-transcriptome result. Only the acute no-recovery sample (`CD1_CupRap_0wks_Rep2`) stands apart, where committed OPCs are accumulating downstream.

---

## 4. Methodology Notes on CellRank

CellRank's terminal states are defined from `cell_type` labels to match the ML classifier's positive class. The `OL` terminal group is the full OL lineage — cells labeled `OL`, `OPC`, or `COP` (i.e. `OL_LINEAGE` from `markers.py`); the `Neuroblast` group is cells labeled `Neuroblast`; every other cell (TAPs, **dNSCs, aNSCs**, Astrocyte, etc.) is transient. For each transient cell, `cellrank_P_OL` is the probability that a random walk following the velocity transition matrix is absorbed into the OL-lineage group.

> **Note on the velocity object's cell-type labels.** The `velocity_combined.h5ad` consumed by this script still carries the pre-refactor labels (`NSC`, `Mural`) rather than the post-refactor split (`aNSC`/`dNSC`, `Pericyte`/`VAMC`). The terminal-state bundle (`OL`, `OPC`, `COP`, `Neuroblast`) is unaffected because those keys exist in both vocabularies. Aggregate per-condition numbers are therefore unchanged. To get aNSC/dNSC-level transient-state interpretation, rebuild the velocity object: `python velocity_build.py`.

Bundling OPC/COP into the OL terminal group keeps the CellRank vs ML comparison apples-to-apples: both methods now report `P(OL lineage)` rather than ML reporting "OL-or-OPC-or-COP" and CellRank reporting strict "OL only". An OL-only terminal definition was tested and gives near-identical results in this dataset (Pearson r = 0.835 either way; mean P shifts by ≤0.006 per condition, OL counts shift by 1–2 cells per condition), because the OPC→COP→OL velocity flux is functional enough that strict-OL absorbing already captures the lineage mass via forward flow.

This approach inherits whatever errors exist in the Leiden clustering. We accept this because:

1. The canonical marker dotplot (`outputs/velocity/dotplot__marker_dotplot.png`) shows that the OL/OPC/COP and Neuroblast cluster labels are well-supported — `Mbp`/`Mog`/`Plp1` tight in OL, `Pdgfra`/`Cspg4` tight in OPC, `Dcx`/`Tubb3` tight in Neuroblast.
2. An alternative using GPCCA's `compute_macrostates()` + `predict_terminal_states()` to discover terminal cells from velocity topology was tested. It requires dense Schur decomposition on a 36k × 36k matrix (~30 min CPU on M2 Pro without SLEPc), and given the clean clustering, the methodological gain is small.

The CellRank transition matrix combines the velocity kernel (80%) and connectivity kernel (20%):

```python
vk = VelocityKernel(adata).compute_transition_matrix(n_jobs=4)
ck = ConnectivityKernel(adata).compute_transition_matrix()
combined_kernel = 0.8 * vk + 0.2 * ck
```

---

## 5. Bottom Line

**ML and CellRank agree on per-cell P(OL lineage) at Pearson r = 0.828 / Spearman r = 0.838 across 2,944 TAPs.**

Two methods built on independent biological signals (gene expression vs splicing kinetics) converge on the same per-cell ranking. This is strong cross-validation for the XGBoost classifier's per-TAP fate score.

By both methods, **CupRap does not enrich OL fate at the TAP stage**. The share of OL-leaning TAPs in CupRap conditions is at or below corresponding Controls. This is consistent with the published claim (Willis et al. Figure 2I, r = 0.995 for *average* TAP transcriptome Cntl vs CupRap): TAP steady-state expression identity does not shift toward OL under CupRap.

Velocity adds a finer-grained observation that is **not** in tension with the r = 0.995 result. Velocity measures **splicing kinetics**, not average expression — two TAP populations can have identical average expression vectors (r = 0.995) while differing in how quickly specific genes are being upregulated. At the TAP stage, the only condition-dependent velocity signal we detect points *away from* OL: Meis2 climbs from rank 10 (Cntl) to rank 2 (CupRap) in TAP velocity drivers (corr 71.99, likelihood 0.79, Spearman 0.93) and is independently SHAP rank 4 NEGATIVE_OL — the model uses its *absence* as an OL signal. Its velocity-driven rise in CupRap TAPs is consistent with TAPs being more NB-committed under CupRap, not less. No positive OL-leaning gene appears in the top 100 TAP velocity drivers in either condition. Net: CupRap leaves TAP average identity unchanged (paper's r = 0.995) and reinforces, rather than reverses, the NB-default trajectory in velocity.

The OL-positive commitment signal — both transcriptionally and via velocity dynamics — first appears **post-TAP**, at the OPC/COP transition. Fa2h is rank 1 in the COP velocity drivers (corr 138.32) AND SHAP rank 9 POSITIVE_OL; Gjc3 is rank 11 in COP velocity drivers AND SHAP rank 8 POSITIVE_OL — same cross-method convergence that's absent at TAP. Driver tables: `outputs/velocity/velocity_drivers_*.csv`.

---

## 6. Output Files

- `compare_tap_fate_methods.py` — analysis script
- `outputs/tap_fate_comparison.csv` — per-TAP `cellrank_P_OL`, `ml_P_OL`, and binary OL calls (threshold > 0.5)
- `outputs/tap_fate_per_condition.csv` — per-condition OL/NB counts for both methods (CSV includes continuous mean P(OL) / P(NB) columns as well)

---

## 7. Reproducibility

```bash
uv sync  # picks up cellrank>=2.3.1, scipy<1.17, numpy>=2
uv run python compare_tap_fate_methods.py
```

Requires:
- `outputs/velocity/velocity_combined.h5ad` (from `rna_velocity_pipeline.py`)
- `outputs/ol_commitment.csv` (from `train_ol_classifier.py`)

CellRank computation: deterministic for a given input. Single run gives the same results.

Dependency notes (current pyproject.toml):
- `cellrank>=2.3.1` — clean upstream, no monkey-patches in `compare_tap_fate_methods.py`
- `scipy>=1.13,<1.17` — capped to satisfy `pygam` (transitive of cellrank)
- The previous numpy-2 monkey-patch and `tol=1e-10` workaround were removed; outputs are unchanged
