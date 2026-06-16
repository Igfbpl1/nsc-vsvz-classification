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
NSC → TAP → OPC → COP → OL
              ↘ Neuroblast

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

## 7. Results — 4-Sample Run

Velocity computed on 4 samples:
- **GSM8253792** — CD1, Cntl, 0wks (baseline anchor)
- **GSM8253796** — NesCre, Cntl, 3wks
- **GSM8253798** — NesCre, CupRap, 3wks, Rep1
- **GSM8253799** — NesCre, CupRap, 3wks, Rep2

### Stream Plot — Biologically Coherent Velocity

The global stream plot (`scvelo_stream_all.png`) shows clear, biologically sensible directional flow:

- **NSC → TAP**: confirmed flow from NSC into TAP cluster
- **TAP → Neuroblast**: strong well-directed flow (dominant arm)
- **TAP → OPC → COP → OL**: rightward flow through the OL lineage is visible — the bifurcation the ML model is built around
- **Non-lineage cells** (Microglia, Astrocyte, Endothelial, Mural, Ependymal): mostly self-contained, not transitioning into lineage — correct

The TAP bifurcation toward Neuroblast and OPC is structurally visible.

### Per-Condition Stream Plots — CupRap Amplifies OL Commitment

| Condition | OL lineage flow | Notes |
|---|---|---|
| CD1_Cntl (0wks) | Sparse, weak | Quietest condition — pre-recovery baseline |
| Cntl (NesCre, 3wks) | Moderate | Some OPC/COP flow but relatively quiet |
| CupRap Rep1 | Stronger, more directed | Arrows through OPC→COP→OL visibly denser |
| CupRap Rep2 | Consistent with Rep1 | Replicate agreement confirms the signal |

The OL commitment arm is clearly more active in CupRap vs both Cntl conditions. The two CupRap replicates are consistent with each other, ruling out a single-sample artifact. CD1_Cntl at 0wks being the quietest establishes a true pre-recovery baseline that makes the CupRap signal sharper by comparison.

### TAP Velocity Drivers — No Positive OL Signal; Meis2 NB Convergence

| Rank | Pooled | Cntl | CupRap |
|---|---|---|---|
| 1 | Bcl11a | Bcl11a | Rnf165 |
| 2 | Rnf165 | Dpysl3 | Srrm4 |
| 3 | Srrm4 | Nfib | Bcl11a |
| 4 | Dpysl3 | Mn1 | **Meis2** |
| 5 | Nfib | Rnf165 | Plxna2 |

**Meis2** ranks 4th in CupRap TAPs (corr 64.85) and is absent from the Cntl top 15. It is also **SHAP rank 5** — two independent methods converge on the same gene specifically in the treatment condition. However, Meis2 is **NB-leaning** (OL mean 0.275 vs NB mean 2.934 — NEGATIVE_OL in the SHAP direction analysis). The model uses its *absence* as an OL signal. Its stronger velocity in CupRap TAPs likely reflects more TAPs actively committing to NB fate under treatment, with a separate subset defaulting to OL by failing to engage this program.

**No positive OL-leaning signal exists at the TAP stage.** The full top 100 TAP drivers were searched for the 7 SHAP-confirmed positive OL markers (Pllp, Gjc3, Dock10, Cryab, Fa2h, Cnp, Tspan2) — none appear. The only OL-associated genes present are Olig2 (rank 73 CupRap, corr 21; rank 50 Cntl, corr 28) and Ncam1 (rank ~63-78, corr ~22) — both ~3× weaker than Meis2 and not CupRap-enriched. Olig2 is actually stronger in Cntl than CupRap.

This confirms that OL commitment is not encoded by a positive transcriptional switch at the TAP stage. The commitment decision is defined by failure to engage the NB program (Bcl11a, Meis2, Srrm4, Nfib), with OL identity only becoming velocity-active at the OPC/COP stage downstream.

### OPC Velocity Drivers — Condition-Specific Signals

Both conditions share the same core OPC program (Ntn1, Gpr37l1, Dpysl3, Cthrc1 dominant in both). Condition-specific differences:

| Gene | Cntl rank | CupRap rank | Function |
|---|---|---|---|
| **Ntrk2** (TrkB) | not in top 10 | 6 | BDNF receptor — neurotrophin signalling activated under treatment |
| **Sgk1** | 6 | not in top 10 | Stress kinase — baseline OPC stress response |
| Fut9 | 9 | 3 | Jumps to rank 3 in CupRap |

Ntrk2 emerging in CupRap OPCs suggests active neurotrophin signalling as part of the injury-driven commitment response.

### COP Velocity Drivers — Strongest Treatment Effect

**Fa2h (fatty acid 2-hydroxylase)** is the top driver in both conditions — a direct myelin sphingolipid synthesis gene, confirming the OL commitment signal is real at the COP stage.

Overall CupRap COP correlation scores are ~2x higher than Cntl (Fa2h: 26.7 vs 15.4), indicating dramatically more active OL commitment dynamics under treatment.

CupRap-specific COP drivers not present in Cntl top 10:

| Gene | CupRap rank | Function |
|---|---|---|
| **Hspa4l** | 5 | Heat shock protein — active myelination machinery |
| **Kif26a** | 7 | Kinesin — microtubule motor for myelin membrane extension |
| **Ptprb** | 9 | Receptor tyrosine phosphatase (canonical marker) |

### Summary: Where the CupRap Effect Lives

```
Stage    Cntl vs CupRap                         CupRap-specific signal
--------------------------------------------------------------------------
TAP      Mostly shared NB program;              Meis2 rank 4 (NB-leaning, SHAP #5)
         no positive OL signal in either        No OL marker in top 100
OPC      ~90% shared core program               Ntrk2 (TrkB/BDNF), Fut9 amplified
COP      ~2x stronger scores in CupRap          Hspa4l, Kif26a, Fa2h amplified
OL       Same program, stronger dynamics        Consistent with prior runs
```

The CupRap effect accumulates downstream of commitment. At the TAP stage, both conditions run the same NB-default program — CupRap intensifies it (Meis2 rises) but introduces no positive OL signal. OL identity only becomes velocity-active at OPC/COP stage.

---

## 8. SHAP Direction Analysis — Critical Methodological Finding

### The Problem

SHAP measures discrimination *magnitude*, not direction. A high SHAP score means a gene strongly distinguishes OL from NB cells — but the model can use that gene as a **positive** OL predictor (high gene → OL fate) OR as a **negative** OL predictor (high gene → NB fate, low gene → OL fate).

### Results

```
shap_rank,gene,shap_value,ol_mean,nb_mean,direction
1,Stmn2,0.711,0.071,2.965,NEGATIVE_OL
2,Igfbpl1,0.440,0.054,2.199,NEGATIVE_OL
3,Pllp,0.258,1.836,0.016,POSITIVE_OL
4,Dlx6os1,0.209,0.050,2.294,NEGATIVE_OL
5,Meis2,0.187,0.275,2.934,NEGATIVE_OL
6,Bcl11a,0.139,0.035,1.521,NEGATIVE_OL
7,Sox11,0.100,0.218,2.999,NEGATIVE_OL
8,Celf4,0.100,0.062,1.613,NEGATIVE_OL
9,Tmsb10,0.094,0.802,3.754,NEGATIVE_OL
10,Gjc3,0.070,1.955,0.014,POSITIVE_OL
11,Meg3,0.067,0.214,2.244,NEGATIVE_OL
12,Arx,0.060,0.035,1.487,NEGATIVE_OL
13,Dock10,0.042,1.231,0.023,POSITIVE_OL
14,Cryab,0.036,2.707,0.091,POSITIVE_OL
15,Grin2b,0.029,0.017,0.602,NEGATIVE_OL
16,Fa2h,0.028,1.207,0.010,POSITIVE_OL
17,Nfib,0.024,0.851,2.888,NEGATIVE_OL
18,Cnp,0.023,3.171,0.186,POSITIVE_OL
19,Tspan2,0.021,2.211,0.039,POSITIVE_OL
20,Hmgn2,0.021,0.731,1.948,NEGATIVE_OL
```

**Summary:** 7 positive OL markers, 13 negative OL markers (NB markers used in reverse).

### The 7 True Positive OL-Leaning Markers

| SHAP rank | Gene | OL% | NB% | OL/NB ratio | Notes |
|---|---|---|---|---|---|
| 3 | Pllp | 90% | 2% | 114× | Plasmolipin — myelin lipid raft |
| 10 | Gjc3 | 89% | 1% | 137× | Connexin 30.2 — gap junction |
| 13 | Dock10 | 79% | 3% | 53× | Cytoskeletal regulator |
| 14 | Cryab | 90% | 9% | 30× | αB-crystallin — myelin maintenance |
| 16 | Fa2h | 75% | 1% | 116× | Myelin sphingolipid synthesis — top COP velocity driver |
| 18 | Cnp | 94% | 18% | 17× | 2',3'-CNPase — myelin marker |
| 19 | Tspan2 | 86% | 4% | 57× | Early OL surface marker |

All 7 are biologically coherent myelin/membrane structural genes. **Fa2h is the key cross-validated gene**: SHAP rank 16 (positive OL marker) AND top COP velocity driver in both Cntl and CupRap conditions (~2x stronger under treatment).

---

## 9. Current Sample Status

| GSM | Strain | Treatment | Timepoint | Velocity status |
|---|---|---|---|---|
| GSM8253792 | CD1 | Cntl | 0wks | ✓ included (baseline anchor) |
| GSM8253796 | NesCre | Cntl | 3wks | ✓ included |
| GSM8253798 | NesCre | CupRap | 3wks | ✓ included |
| GSM8253799 | NesCre | CupRap | 3wks | ✓ included |
| GSM8253793 | CD1 | Cntl | 3wks | downloading |
| GSM8253794 | CD1 | CupRap | 3wks | downloading |
| GSM8253795 | CD1 | CupRap | 3wks | pending |
| GSM8253797 | NesCre | Cntl | 3wks | pending |
| GSM8647352 | CD1 | CupRap | 0wks | pending |
| GSM8647353 | CD1 | CupRap | 0wks | pending |

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

The OL-lineage proportion is higher at 0wks — the system expands precursors immediately in response to injury. By 3wks, OPCs and COPs have matured into OLs and the pool empties. TAPs and Neuroblasts increase dramatically at 3wks as normal V-SVZ neurogenesis resumes.

The current velocity analysis (3wks samples dominant) captures the recovery phase. The acute commitment moment — when microglial IGF1/OSM signaling first hits NSC/TAP populations — occurs at 0wks and will be testable once the remaining CD1 0wks samples (GSM8647352, GSM8647353) are added.
