# RNA Velocity & Trajectory Analysis

---

## 1. Background & Motivation

The project uses a static snapshot of gene expression to classify cell types and an XGBoost model (with SHAP scores) to identify trigger genes that predict whether a TAP will commit to the Neuroblast or OL lineage.

> Cell-type annotation now uses paper-aligned marker panels (Willis et al. 2025 STAR Methods) with NSC split into **aNSC** (Egfr, Ascl1) and **dNSC** (Meg3, Sparc, Fbxo2, Id3), and Mural split into **Pericyte** and **VAMC**. Cluster assignment is `idxmax` over absolute panel scores with a z-score tie-break when the top two panels are within margin 0.4 and the runner-up has absolute signal ≥ 0.5. See `markers.py` and `run_pipeline.py:41-69`.

While SHAP scores tell us *which* genes distinguish the lineages, static data cannot tell us *when* those genes act — are they early drivers or late consequence markers?

**RNA velocity answers this.** By modeling the ratio of unspliced (newly transcribed) to spliced (mature) mRNA, scVelo acts as a molecular clock. It infers directional flow, allowing genes to be ranked by how tightly their splicing dynamics align with cell movement along the trajectory — separating early drivers from late state markers.

---

## 2. Biological Questions

By comparing CD1 Cntl (0wks baseline) + NesCre Cntl (3wks) vs NesCre CupRap (3wks × 2 replicates):

1. **Directionality** — Does TAP→OL transit accelerate under CupRap treatment?
2. **Early Gene Identification** — Which genes initiate the OL transcriptional cascade before the cell looks like an OL?
3. **Perturbation Response** — Do velocity drivers differ between Cntl and CupRap conditions?
4. **Trigger Gene Validation** — Do the top SHAP trigger genes appear in velocity driver rankings?

---

## 3. How RNA Velocity Works

Each kb count output contains a cells × genes matrix with two layers:

| Layer | Source | Meaning |
|---|---|---|
| `spliced` | cDNA reads (exonic) | Mature mRNA — already processed |
| `unspliced` | intronic reads | Pre-mRNA — newly transcribed, not yet spliced |

scVelo models splicing kinetics per gene using this ODE system:

```
du/dt = α - βu       (unspliced: produced at rate α, consumed at rate β)
ds/dt = βu - γs      (spliced: produced from unspliced, degraded at rate γ)
```

By fitting α (transcription), β (splicing), γ (degradation) per gene, scVelo predicts where spliced RNA is heading — that is the **velocity vector** for each cell.

---

## 4. Data Processing Steps

### Step 1 — Build kb Reference Index
```bash
cd sra_runs
nohup kb ref \
  -i index.idx -g t2g.txt -f1 cdna.fa -f2 intron.fa \
  -c1 cdna_t2c.txt -c2 intron_t2c.txt \
  --workflow lamanno \
  Mus_musculus.GRCm39.dna.primary_assembly.fa.gz \
  Mus_musculus.GRCm39.110.gtf.gz > ref.log 2>&1 &
```

### Step 2 — Download FASTQs and Run kb count

**Important — two FASTQ layouts exist depending on SRA submission:**

| Submission type | Files produced | Use for kb count |
|---|---|---|
| Dual-indexed (e.g. SRR314436xx) | `_1` I1, `_2` I2, `_3` R1 (28bp), `_4` R2 (91bp) | `_3` + `_4` |
| Single-indexed (e.g. SRR289128xx) | `_1` I1, `_2` R1 (28bp), `_3` R2 (91bp) | `_2` + `_3` |

Check layout by reading first sequence length:
```bash
awk 'NR==2{print length($0); exit}' SRR_ID_1.fastq   # 8bp = I1
awk 'NR==2{print length($0); exit}' SRR_ID_2.fastq   # 28bp = R1 or 91bp = R2
awk 'NR==2{print length($0); exit}' SRR_ID_3.fastq
```

`process_sample.sh` auto-detects the layout (checks for `_4.fastq`; falls back to `_2`+`_3`).

### Step 3 — Run Velocity Pipeline
```bash
uv run python rna_velocity_pipeline.py
uv run python compare_tap_fate_methods.py
```

Output: `outputs/velocity/`

> **Dependency update (cellrank 2.0.7 → 2.3.1, scipy 1.17.1 → 1.16.3)**
> `compare_tap_fate_methods.py` previously included a 16-line monkey-patch on `numpy.testing.assert_array_equal` to work around a kwarg-renaming incompatibility in cellrank 2.0.7's `VelocityKernel.__init__`. With cellrank ≥2.3.1 the upstream fix is in place. The patch is removed and `compute_fate_probabilities()` runs at defaults (no `tol=1e-10` workaround). The per-condition `tap_fate_per_condition.csv` numbers are byte-identical before and after this upgrade, confirming pure cleanup with no behavior change.

---

## 5. Output Files

| File | What it shows |
|---|---|
| `scvelo_stream_all.png` | Global velocity arrows on UMAP — direction of differentiation for all cells |
| `scvelo_stream_{label}.png` | Per-condition velocity stream (CD1_Cntl, Cntl, CupRap_Rep1, CupRap_Rep2) |
| `velocity_drivers_{CellType}.csv` | Pooled velocity driver genes per cell type, ranked by correlation |
| `velocity_drivers_{CellType}_Cntl.csv` | Driver genes for Cntl cells of that type |
| `velocity_drivers_{CellType}_CupRap.csv` | Driver genes for CupRap cells of that type |
| `velocity_combined.h5ad` | Full AnnData with spliced/unspliced layers + velocity vectors |

Each driver CSV contains: `rank`, `ensembl_id`, `gene_name`, `corr`, `likelihood`, `spearmans`, `shap_rank`, `shap_value`, `canonical`.

---

## 6. Interpreting Velocity Driver Genes

### What Each Column Represents

Each driver CSV corresponds to **where in the trajectory the gene's velocity is active**:

```
dNSC → aNSC → TAP → Neuroblast
                  ↘ OPC → COP → OL

Early ←——————————————————————————→ Late
(dNSC/aNSC/TAP)                 (OL/NB)
```

| Column | When the gene is active | What it means |
|---|---|---|
| **dNSC** | In dormant NSCs (Meg3/Sparc/Fbxo2/Id3 high) | Quiescent stem cell state; pre-activation |
| **aNSC** | In activated NSCs (Egfr/Ascl1 high) | Stem cell that just entered cycle; fate-uncommitted |
| **TAP** | In TAPs, before commitment | **Early fate-decision gene** — most valuable for this project |
| **OPC** | After OL commitment, driving OPC→COP | Mid-lineage marker |
| **COP** | Late, driving COP→OL maturation | Late commitment marker |
| **OL** | In mature OLs | Terminal differentiation marker |
| **Neuroblast** | In NB arm after commitment | NB terminal marker |

> **Note on the velocity driver CSVs**: `velocity_drivers_TAP.csv`, `velocity_drivers_OPC.csv`, etc. were generated *before* the markers.py refactor that split NSC → {dNSC, aNSC}. The current `outputs/velocity/` therefore has no `velocity_drivers_aNSC.csv` / `velocity_drivers_dNSC.csv`, and TAP velocity drivers include ~1591 cells that are now correctly classified as aNSC (not TAP). Rebuilding `velocity_combined.h5ad` from the refactored `processed.h5ad` (via `python velocity_build.py`) will produce per-aNSC/dNSC driver CSVs and a cleaner TAP-only driver list.

### How to Identify a True Early Marker

A gene qualifies as an early marker if it satisfies:

| Criterion | What to look for |
|---|---|
| Active before commitment | Appears in NSC or TAP driver CSV at high rank |
| Dynamically regulated | High Spearman correlation + likelihood |
| Discriminates fate | Appears in SHAP top 20 (trigger_genes.csv) |

**The sweet spot is low velocity rank + low SHAP rank** — dynamically active during TAP transition AND predictive of which fate the cell takes.

---

## 7. Results — 7-Sample Run

Velocity computed on 7 samples (36,728 cells; 2,949 TAPs):
- **GSM8253792** — CD1, Cntl, 0wks
- **GSM8253793** — CD1, Cntl, 3wks
- **GSM8253794** — CD1, CupRap, 3wks
- **GSM8253796** — NesCre, Cntl, 3wks
- **GSM8253798** — NesCre, CupRap, 3wks, Rep1
- **GSM8253799** — NesCre, CupRap, 3wks, Rep2
- **GSM8647353** — CD1, CupRap, 0wks, Rep2

### Stream Plot

The stream plot shows the same trajectory topology as earlier (4-sample) runs:
- **dNSC → aNSC → TAP**: clear flow along the precursor cascade
- **TAP → Neuroblast**: strong dominant arm
- **TAP bifurcation**: rightward arm toward OPC → COP → OL visible alongside dominant Neuroblast arm
- Non-lineage cells (Microglia, Astrocyte, Endothelial, Pericyte, VAMC, Ependymal) self-contained

Adding samples did not change the trajectory structure. The current stream PNGs were generated against the pre-refactor labels (NSC, Mural); cell-type *colors* in those PNGs will differ once the velocity object is rebuilt against the refactored `processed.h5ad`, but trajectory direction is independent of label and unchanged.

### TAP Velocity Drivers — Srrm4 Rank 1, Meis2 Rank 2 in CupRap

| Rank | Pooled | Cntl | CupRap |
|---|---|---|---|
| 1 | Srrm4 | Srrm4 | **Srrm4** |
| 2 | Meis2 | Draxin | **Meis2** |
| 3 | Plxna2 | Dpysl3 | Plxna2 |
| 4 | Dpysl3 | Igfbpl1 | Negr1 |
| 5 | Nfib | Nsg2 | Nfib |

**Meis2** is rank 2 in CupRap TAPs (corr 71.99, likelihood 0.79, Spearman 0.93) and rank 10 in Cntl — it climbs specifically in the treatment condition. It is also **SHAP rank 4** — two independent methods converge on the same gene. However, Meis2 is **NB-leaning** (NEGATIVE_OL in the SHAP direction analysis); the model uses its *absence* as an OL signal.

**No positive OL-leaning signal exists at the TAP stage.** All 16 SHAP-confirmed positive OL markers are absent from the top 100 TAP drivers in both conditions. Olig2 and Ncam1 appear but at low ranks with weak correlations.

### OPC/COP Velocity Drivers — Key Patterns

**Pattern 1: ML model and velocity converge on the same OL genes at the COP stage**

The ML model (trained on non-canonical HVGs, no prior gene list) identified positive OL markers. Velocity independently identifies two of these as strong COP drivers:

| Gene | SHAP rank | OPC rank | COP rank |
|---|---|---|---|
| Fa2h | 9 | 37 | **1** |
| Gjc3 | 8 | 22 | 11 |

Neither method used the other's output. Both arrived at the same genes independently. (Dock10, flagged in an earlier run, dropped out of both the SHAP top-50 and the OPC/COP driver lists in this run.)

**Pattern 2: Myelin lipid biosynthesis program across OPC→COP**

```
OPC stage:
  Ugt8a (rank 30) — galactosylceramide synthase
  Fa2h  (rank 38) — fatty acid 2-hydroxylase

COP stage:
  Fa2h   (rank  1) — dominant COP driver (corr 138.32)
  Cldn11 (rank 22) — myelin tight junction protein
```

**Pattern 3: Ncam1 is a continuous driver across all three OL lineage stages**

| Stage | Rank | Corr | Likelihood |
|---|---|---|---|
| OPC | 12 | 35.3 | 0.57 |
| COP | 34 | 16.5 | 0.57 |
| OL | 4 | 123.6 | 0.57 |

Active across all three stages. Likelihood 0.57 is among the highest in the OPC list.

**Pattern 4: Ntn1 (Netrin-1) is the top OPC driver (rank 1, corr 86.93)**

The strongest velocity signal at the OPC stage, ahead of all myelin genes.

---

## 8. SHAP Direction Analysis — Critical Methodological Finding

SHAP measures discrimination *magnitude*, not direction. A high SHAP score means a gene strongly distinguishes OL from NB cells — but the model can use that gene as a **positive** OL predictor (high gene → OL fate) OR as a **negative** OL predictor (high gene → NB fate, low gene → OL fate). Direction is determined by comparing mean expression in OL-lineage (OPC+COP+OL) vs Neuroblast cells. Full rankings with direction in `outputs/trigger_genes.csv`.

Top 10 SHAP genes (SHAP computed on held-out NesCre test set): Stmn2, Pllp, Tmsb10, Meis2, Dlx6os1, Igfbpl1, Sox11, Gjc3, Fa2h, Tubb2b. Of these, Pllp (rank 2), Gjc3 (rank 8), and Fa2h (rank 9) are POSITIVE_OL; the remaining 7 are NEGATIVE_OL (NB markers used in reverse).

**Summary across top 50:** 22 POSITIVE_OL, 28 NEGATIVE_OL.

### True Positive OL-Leaning Markers (SHAP top 50)

_SHAP ranks below are from the current 7-sample run (Pllp 2, Gjc3 8, Fa2h 9). The OL% / NB% / ratio columns are carried over from the prior run and are pending re-verification against the current `outputs/trigger_genes.csv`; treat them as approximate. Dock10 dropped out of the SHAP top-50 in this run._

| SHAP rank | Gene | OL% | NB% | OL/NB ratio | Notes |
|---|---|---|---|---|---|
| 2 | Pllp | 90% | 2% | 114× | Plasmolipin — myelin lipid raft |
| 8 | Gjc3 | 89% | 1% | 137× | Connexin 30.2 — gap junction |
| 9 | Fa2h | 75% | 1% | 116× | Myelin sphingolipid synthesis — top COP velocity driver |
| 13 | Cryab | 90% | 9% | 30× | αB-crystallin — myelin maintenance |
| 17 | Tspan2 | 86% | 4% | 57× | Early OL surface marker |
| 18 | Cnp | 94% | 18% | 17× | 2',3'-CNPase — myelin marker |
| 29 | Myrf | 72% | 1% | 244× | Master myelin TF |
| 31 | Trf | 76% | 8% | 26× | Transferrin — iron transport for myelination |
| 35 | Gatm | 89% | 4% | 53× | Creatine synthesis |
| 36 | Mal | 68% | 5% | 45× | Myelin and lymphocyte protein |
| 40 | Gsn | 78% | 10% | 15× | Gelsolin — cytoskeletal |
| 41 | Neat1 | 71% | 7% | 20× | lncRNA expressed in OL |
| 42 | Ndrg1 | 63% | 2% | 54× | Myelin maintenance |
| 44 | Car2 | 69% | 9% | 23× | Carbonic anhydrase II — myelin |
| 45 | S100a1 | 65% | 3% | 34× | S100 calcium-binding protein |

**Fa2h** is cross-validated: SHAP rank 9 (POSITIVE_OL) AND top COP velocity driver (rank 1) in both conditions.

---

## 9. Current Sample Status

| GSM                                  | Strain | Treatment | Timepoint | Velocity status |
|--------------------------------------|---|---|---|-----------------|
| GSM8253792                           | CD1 | Cntl | 0wks | ✓ included      |
| GSM8253793                           | CD1 | Cntl | 3wks | ✓ included      |
| GSM8253794                           | CD1 | CupRap | 3wks | ✓ included      |
| GSM8253796                           | NesCre | Cntl | 3wks | ✓ included      |
| GSM8253798                           | NesCre | CupRap | 3wks | ✓ included      |
| GSM8253799                           | NesCre | CupRap | 3wks | ✓ included      |
| GSM8647353                           | CD1 | CupRap | 0wks | ✓ included      |
| GSM8253795 (SRR28912877)             | CD1 | CupRap | 3wks | x discluded (counted, not in current velocity build) |
| GSM8253797 (SRR28912867 SRR28912868) | NesCre | Cntl | 3wks | x discluded     |
| GSM8647352  (SRR31443698 SRR31443699)                       | CD1 | CupRap | 0wks | pending         |

---

## 10. What 0wks vs 3wks Tells Us About OL Biology

### Timepoint Asymmetry

- **0wks** = sampled immediately after treatment ends (acute injury/response phase)
- **3wks** = sampled 3 weeks after treatment ends (recovery/remyelination phase)

The CD1_Cntl 0wks stream plot showing the quietest OL lineage flow is consistent with this: at baseline (no injury, no recovery), the OL commitment machinery is minimally active. The CupRap 3wks samples showing the strongest OL commitment dynamics reflects active remyelination.

### Cell-Type Composition Shift (CD1 CupRap only, strain-controlled)

| Cell type | CupRap 0wks (n=3,928) | CupRap 3wks (n=3,550) |
|---|---|---|
| TAP | 2.0% | 9.0% |
| Neuroblast | 1.9% | 28.7% |
| Microglia | 34.6% | 30.3% |
| OPC | 6.0% | 1.5% |
| COP | 11.4% | 3.0% |
| OL | 20.4% | 13.4% |
| **OL-lineage total** | **37.8%** | **17.9%** |

(CD1 CupRap only: 0wks = GSM8647353; 3wks = GSM8253794. Strain-matched so the timepoint effect is not confounded by strain.)

The OL-lineage proportion (OPC+COP+OL) is much higher at 0wks (37.8%) than 3wks (17.9%). OPC and COP fractions collapse by 3wks (OPC 6.0%→1.5%, COP 11.4%→3.0%), while Neuroblast surges (1.9%→28.7%) and TAP rises (2.0%→9.0%). This is the acute-injury → recovery transition: at 0wks the niche is at peak OL-precursor generation with neurogenesis suppressed; by 3wks the OL-precursor wave has receded and the neurogenic program has rebounded. (These are *proportions* of a fixed pie — a population's share moving does not by itself prove its absolute count changed.)

The current velocity analysis includes GSM8647353 (CD1_CupRap_0wks_Rep2). GSM8647352 (CD1_CupRap_0wks_Rep1) is still pending.

### aNSC fate-bias observation (independent of velocity)

After the markers.py refactor (split NSC → aNSC/dNSC), the per-cell `bias` score (`(score_OPC + score_COP + score_OL)/3 − score_Neuroblast`) reveals a clean acute-injury signal at the aNSC stage. Using all CD1 samples:

| Condition | n aNSCs | bias median | fraction OL-leaning (bias > 0) |
|---|---|---|---|
| Control (0+3wks) | 1029 | -0.07 | 43% |
| CupRap 3wkRecov | 428 | -0.12 | 41% |
| CupRap **NoRecov** | 134 | **+0.15** | **67%** |

Mann-Whitney U test, aNSC bias Control vs CupRap_NoRecov: **p = 2.9 × 10⁻¹⁴**. The 3-week recovery aNSCs are statistically indistinguishable from controls (p = 0.77), so the OL-commitment push is **transient** — acute demyelination tilts the activated stem cell pool toward oligodendrogenesis, and by 3 weeks the distribution returns to baseline. This matches the paper's central claim about injury-driven NSC activation toward repair and is recovered here without telling the model the conditions.

Caveat: the ML model's `prob_OL` on aNSCs is unreliable (~0.45 across conditions) because aNSCs are not in the training set (POS = OL+COP+OPC, NEG = Neuroblast). For aNSC inference, use the `bias` score, not `prob_OL`. See `ol_commitment.csv` and `train_ol_classifier.py:56-61`.
