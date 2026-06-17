# RNA Velocity & Trajectory Analysis

---

## 1. Background & Motivation

The project uses a static snapshot of gene expression to classify cell types and an XGBoost model (with SHAP scores) to identify trigger genes that predict whether a TAP will commit to the Neuroblast or OL lineage.

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

---

## 5. Output Files

| File | What it shows |
|---|---|
| `scvelo_stream_all.png` | Global velocity arrows on UMAP — direction of differentiation for all cells |
| `scvelo_stream_{label}.png` | Per-condition velocity stream (CD1_Cntl, Cntl, CupRap_Rep1, CupRap_Rep2) |
| `scvelo_phase_{CellType}.png` | Phase portraits for top 4 velocity driver genes per cell type |
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
NSC → TAP → Neuroblast
          ↘ OPC → COP → OL

Early ←————————————————————→ Late
(NSC/TAP)                  (OL/NB)
```

| Column | When the gene is active | What it means |
|---|---|---|
| **NSC** | In NSCs, upstream of TAPs | May mark cells about to become TAPs |
| **TAP** | In TAPs, before commitment | **Early fate-decision gene** — most valuable for this project |
| **OPC** | After OL commitment, driving OPC→COP | Mid-lineage marker |
| **COP** | Late, driving COP→OL maturation | Late commitment marker |
| **OL** | In mature OLs | Terminal differentiation marker |
| **Neuroblast** | In NB arm after commitment | NB terminal marker |

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

Velocity computed on 7 samples:
- **GSM8253792** — CD1, Cntl, 0wks
- **GSM8253793** — CD1, Cntl, 3wks
- **GSM8253794** — CD1, CupRap, 3wks
- **GSM8253796** — NesCre, Cntl, 3wks
- **GSM8253798** — NesCre, CupRap, 3wks, Rep1
- **GSM8253799** — NesCre, CupRap, 3wks, Rep2
- **GSM8647353** — CD1, CupRap, 0wks, Rep2

### Stream Plot

The 7-sample stream plot shows the same trajectory topology as the 4-sample run:
- **NSC → TAP**: clear flow
- **TAP → Neuroblast**: strong dominant arm
- **TAP bifurcation**: rightward arm toward OPC → COP → OL visible alongside dominant Neuroblast arm
- Non-lineage cells (Microglia, Astrocyte, Endothelial, Mural, Ependymal) self-contained

Adding 3 samples did not change the trajectory structure.

### TAP Velocity Drivers — Meis2 is Rank 1 in CupRap

| Rank | Pooled | Cntl | CupRap |
|---|---|---|---|
| 1 | Dpysl3 | Dpysl3 | **Meis2** |
| 2 | Srrm4 | Nfib | Srrm4 |
| 3 | Nfib | Elavl2 | Plxna2 |
| 4 | Plxna2 | Bcl11a | Dpysl3 |
| 5 | Fut9 | Rnd3 | Nfib |

**Meis2** is rank 1 in CupRap TAPs (corr 79.77, likelihood 0.79, Spearman 0.93) and rank 17 in Cntl. It is also **SHAP rank 5** — two independent methods converge on the same gene specifically in the treatment condition. However, Meis2 is **NB-leaning** (OL mean 0.275 vs NB mean 2.934 — NEGATIVE_OL in the SHAP direction analysis). The model uses its *absence* as an OL signal.

**No positive OL-leaning signal exists at the TAP stage.** All 16 SHAP-confirmed positive OL markers are absent from the top 100 TAP drivers in both conditions. Olig2 and Ncam1 appear but at low ranks with weak correlations.

### OPC/COP Velocity Drivers — Key Patterns

**Pattern 1: ML model and velocity converge on the same OL genes at the COP stage**

The ML model (trained on non-canonical HVGs, no prior gene list) identified 16 positive OL markers. Velocity independently identifies 3 of these as COP drivers:

| Gene | SHAP rank | OPC rank | COP rank |
|---|---|---|---|
| Fa2h | 15 | 38 | **1** |
| Gjc3 | 10 | 18 | 9 |
| Dock10 | 12 | — | 49 |

Neither method used the other's output. Both arrived at the same 3 genes independently.

**Pattern 2: Myelin lipid biosynthesis program across OPC→COP**

```
OPC stage:
  Ugt8a (rank 30) — galactosylceramide synthase
  Fa2h  (rank 38) — fatty acid 2-hydroxylase

COP stage:
  Fa2h   (rank  1) — top driver (corr 60.74, up from 30.91 in 4-sample run)
  Cldn11 (rank 22) — myelin tight junction protein
```

**Pattern 3: Ncam1 is a continuous driver across all three OL lineage stages**

| Stage | Rank | Corr | Likelihood |
|---|---|---|---|
| OPC | 11 | 39.0 | 0.57 |
| COP | 41 | 12.6 | 0.57 |
| OL | 4 | 123.6 | 0.57 |

Active across all three stages. Likelihood 0.57 is among the highest in the OPC list.

**Pattern 4: Ntn1 (Netrin-1) is the top OPC driver (rank 1, corr 96.87)**

The strongest velocity signal at OPC stage, ahead of all myelin genes. Corr increased from 62.85 (4-sample) to 96.87 (7-sample).

---

## 8. SHAP Direction Analysis — Critical Methodological Finding

SHAP measures discrimination *magnitude*, not direction. A high SHAP score means a gene strongly distinguishes OL from NB cells — but the model can use that gene as a **positive** OL predictor (high gene → OL fate) OR as a **negative** OL predictor (high gene → NB fate, low gene → OL fate). Direction is determined by comparing mean expression in OL-lineage (OPC+COP+OL) vs Neuroblast cells. Full rankings with direction in `outputs/trigger_genes.csv`.

Top 10 SHAP genes (SHAP computed on held-out NesCre test set): Stmn2, Igfbpl1, Pllp, Dlx6os1, Meis2, Bcl11a, Celf4, Sox11, Tmsb10, Gjc3. Of these, Pllp (rank 3) and Gjc3 (rank 10) are POSITIVE_OL; the remaining 8 are NEGATIVE_OL (NB markers used in reverse).

**Summary across top 50:** 24 POSITIVE_OL, 26 NEGATIVE_OL.

### True Positive OL-Leaning Markers (SHAP top 50)

| SHAP rank | Gene | OL% | NB% | OL/NB ratio | Notes |
|---|---|---|---|---|---|
| 3 | Pllp | 90% | 2% | 114× | Plasmolipin — myelin lipid raft |
| 10 | Gjc3 | 89% | 1% | 137× | Connexin 30.2 — gap junction |
| 12 | Dock10 | 79% | 3% | 53× | Cytoskeletal regulator |
| 13 | Cryab | 90% | 9% | 30× | αB-crystallin — myelin maintenance |
| 15 | Fa2h | 75% | 1% | 116× | Myelin sphingolipid synthesis — top COP velocity driver |
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

**Fa2h** is cross-validated: SHAP rank 15 (POSITIVE_OL) AND top COP velocity driver in both conditions.

---

## 9. Current Sample Status

| GSM | Strain | Treatment | Timepoint | Velocity status |
|---|---|---|---|---|
| GSM8253792 | CD1 | Cntl | 0wks | ✓ included |
| GSM8253793 | CD1 | Cntl | 3wks | ✓ included |
| GSM8253794 | CD1 | CupRap | 3wks | ✓ included |
| GSM8253796 | NesCre | Cntl | 3wks | ✓ included |
| GSM8253798 | NesCre | CupRap | 3wks | ✓ included |
| GSM8253799 | NesCre | CupRap | 3wks | ✓ included |
| GSM8647353 | CD1 | CupRap | 0wks | ✓ included |
| GSM8253795 | CD1 | CupRap | 3wks | pending |
| GSM8253797 | NesCre | Cntl | 3wks | pending |
| GSM8647352 | CD1 | CupRap | 0wks | pending |

---

## 10. What 0wks vs 3wks Tells Us About OL Biology

### Timepoint Asymmetry

- **0wks** = sampled immediately after treatment ends (acute injury/response phase)
- **3wks** = sampled 3 weeks after treatment ends (recovery/remyelination phase)

The CD1_Cntl 0wks stream plot showing the quietest OL lineage flow is consistent with this: at baseline (no injury, no recovery), the OL commitment machinery is minimally active. The CupRap 3wks samples showing the strongest OL commitment dynamics reflects active remyelination.

### Cell-Type Composition Shift (CD1 CupRap only, strain-controlled)

| Cell type | CupRap 0wks (n=8,062) | CupRap 3wks (n=5,943) |
|---|---|---|
| TAP | 3.8% | 13.9% |
| Neuroblast | 3.8% | 31.0% |
| Microglia | 53.6% | 29.9% |
| OPC | 4.1% | 1.7% |
| COP | 7.3% | 4.1% |
| OL | 12.6% | 11.4% |
| **OL-lineage total** | **24.0%** | **17.2%** |

The OL-lineage proportion (OPC+COP+OL) is higher at 0wks (24.0%) than 3wks (17.2%). OPC and COP fractions both decrease by 3wks while OL fraction is similar. TAP and Neuroblast fractions increase substantially at 3wks. Microglia decrease from 53.6% to 29.9%.

The current velocity analysis uses 3wks samples. The CD1 0wks samples (GSM8647352, GSM8647353) are not yet included.
