# RNA Velocity & Trajectory Analysis

---

## 1. Background & Motivation

The project uses a static snapshot of gene expression to classify cell types and an XGBoost model (with SHAP scores) to identify trigger genes that predict whether a TAP will commit to the Neuroblast or OL lineage.

While SHAP scores tell us *which* genes distinguish the lineages, static data cannot tell us *when* those genes act — are they early drivers or late consequence markers?

**RNA velocity answers this.** By modeling the ratio of unspliced (newly transcribed) to spliced (mature) mRNA, scVelo acts as a molecular clock. It infers a continuous latent time and directional flow, allowing genes to be sorted by when their transcription spikes along the trajectory — clearly separating early drivers from late state markers.

---

## 2. Biological Questions

By comparing NesCre Cntl (GSM8253796) vs NesCre CupRap (GSM8253798):

1. **Directionality** — Does TAP→OL transit accelerate or stall under CupRap treatment?
2. **Early Gene Identification** — Which genes initiate the OL transcriptional cascade before the cell looks like an OL?
3. **Perturbation Response** — Do early genes activate differently in CupRap vs Cntl?
4. **Trigger Gene Validation** — Do the top SHAP trigger genes precede the lineage split in developmental time?

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

**Pipeline:**

```
GSM8253796 MTX (Cntl)  ──┐
                           ├── merge on barcodes ──> combined AnnData
GSM8253798 MTX (CupRap)──┘                          spliced + unspliced layers
                                                      +
processed.h5ad ───────────────────────────────────> UMAP, cell_type labels

scVelo fits kinetics per gene → velocity vector per cell → directional arrows on UMAP
```

Both samples stay in one AnnData with a `velocity_label` column. Plots are stratified by condition to compare Cntl vs CupRap side by side.

---

## 4. Data Processing Steps

### Step 1 — Download Reference Files
```bash
cd sra_runs
wget ftp://ftp.ensembl.org/pub/release-110/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
wget ftp://ftp.ensembl.org/pub/release-110/gtf/mus_musculus/Mus_musculus.GRCm39.110.gtf.gz
```

### Step 2 — Build kb Reference Index
```bash
nohup kb ref \
  -i index.idx -g t2g.txt -f1 cdna.fa -f2 intron.fa \
  -c1 cdna_t2c.txt -c2 intron_t2c.txt \
  --workflow lamanno \
  Mus_musculus.GRCm39.dna.primary_assembly.fa.gz \
  Mus_musculus.GRCm39.110.gtf.gz > ref.log 2>&1 &
```

### Step 3 — Download FASTQs and Run kb count

**Local Mac (small samples):**
```bash
# download with --include-technical to get barcode reads (_3 files)
fasterq-dump --split-files --include-technical SRR_ID1 SRR_ID2 --outdir .

# run kb count (script handles per-sample FASTQ lists)
bash run_kb_count.sh GSM8253796
bash run_kb_count.sh GSM8253798
```

**EC2 (large samples — download in parallel):**
```bash
# parallel prefetch + fasterq-dump per SRR to maximize throughput
nohup fasterq-dump --split-files --include-technical SRR_ID1 --outdir . > fastq_1.log 2>&1 &
nohup fasterq-dump --split-files --include-technical SRR_ID2 --outdir . > fastq_2.log 2>&1 &
nohup fasterq-dump --split-files --include-technical SRR_ID3 --outdir . > fastq_3.log 2>&1 &

# run kb count directly (without run_kb_count.sh) so binary paths use system defaults
nohup kb count \
  -i index.idx \
  -g t2g.txt \
  -x 10xv3 \
  -o kb_output_GSM8253799 \
  -c1 cdna_t2c.txt \
  -c2 intron_t2c.txt \
  --workflow lamanno \
  SRR28912869_3.fastq SRR28912869_4.fastq \
  SRR28912870_3.fastq SRR28912870_4.fastq \
  SRR28912871_3.fastq SRR28912871_4.fastq \
  > kb_count.log 2>&1 &

tail -f kb_count.log
```

**Transfer counts_unfiltered back to local:**
```bash
# on EC2
tar czf gsm799_counts.tgz -C kb_output_GSM8253799 counts_unfiltered

# on Mac
scp -i ~/.ssh/vsvz-key.pem \
  ubuntu@<ec2-dns>:~/project/sra_runs/gsm799_counts.tgz \
  sra_runs/

mkdir -p sra_runs/kb_output_GSM8253799
tar xzf sra_runs/gsm799_counts.tgz -C sra_runs/kb_output_GSM8253799
```

**Notes:**
- FASTQ files have 4 reads per spot. Use `_3` (R1, 28bp barcode+UMI) and `_4` (R2, 90bp cDNA) only.
- `--include-technical` is required — without it fasterq-dump omits the barcode reads (`_3` files).
- After kb count completes, delete FASTQs + BUS files to free space; keep only `counts_unfiltered/`.

### Step 4 — Run Velocity Pipeline
```bash
uv run python rna_velocity_pipeline.py
```

Output: `outputs/velocity/`

---

## 5. Output Files

| File | What it shows |
|---|---|
| `scvelo_stream_all.png` | Global velocity arrows on UMAP — direction of differentiation for all cells |
| `scvelo_stream_Cntl.png` | Velocity arrows for Cntl condition only (GSM8253796) |
| `scvelo_stream_CupRap.png` | Velocity arrows for CupRap condition only (GSM8253798) |
| `scvelo_latent_time.png` | Each cell coloured by latent time (0=NSC/earliest, 1=committed/latest) |
| `scvelo_heatmap_latent_time_heatmap.png` | Driver genes ordered by when they peak along latent time — left=early, right=late |
| `scvelo_zoom_NSC_TAP.png` | Zoomed velocity at the NSC→TAP boundary |
| `scvelo_zoom_TAP_Neuroblast.png` | Zoomed velocity along the NB arm |
| `scvelo_zoom_TAP_OL.png` | Zoomed velocity along the OL arm |
| `scvelo_zoom_TAP_bifurcation.png` | The bifurcation point where TAPs split between NB and OL fate |
| `scvelo_phase_NSC_TAP.png` | Phase portraits for NSC→TAP markers: Sox2, Egfr, Mki67, Ascl1, Ccnd2 |
| `scvelo_phase_TAP_Neuroblast.png` | Phase portraits for TAP→NB markers: Dcx, Dlx1, Tubb3, Dlx2 |
| `scvelo_phase_TAP_OL.png` | Phase portraits for TAP→OL markers: Olig2, Pdgfra, Mbp, Mog, Gpr17 |
| `velocity_drivers_NSC.csv` | Velocity driver genes for NSC, ranked by correlation |
| `velocity_drivers_TAP.csv` | Velocity driver genes for TAP, ranked by correlation |
| `velocity_drivers_Neuroblast.csv` | Velocity driver genes for Neuroblast, ranked by correlation |
| `velocity_drivers_OPC.csv` | Velocity driver genes for OPC, ranked by correlation |
| `velocity_drivers_COP.csv` | Velocity driver genes for COP, ranked by correlation |
| `velocity_drivers_OL.csv` | Velocity driver genes for OL, ranked by correlation |

A gene with a clean **almond shape** in a phase portrait has dynamic spliced/unspliced kinetics consistent with velocity. A cloud shape means noisy signal.

Each driver CSV contains: `rank`, `ensembl_id`, `gene_name`, `corr`, `likelihood`, `spearmans`, `shap_rank`, `shap_value`, `canonical`.

---

## 6. Interpreting Velocity Driver Genes

### What Each Column Represents

Each driver CSV corresponds to **where in the trajectory the gene's velocity is active**, not what the gene does:

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

### What Velocity Rank Means

Rank 1 means that gene's splicing dynamics **most tightly and consistently align with the velocity direction** across all cells of that type.

- **Rank 1 (Bcl11a in TAP)** — every time a TAP moves toward NB fate, Bcl11a's unspliced RNA accumulates in lockstep. Strongest correlation of any gene in TAPs.
- **Rank 10 (Epha5 in TAP)** — also correlated but weaker; may only track velocity in a subset of TAPs, or signal is noisier, or it activates slightly later.

High rank alone does not tell you **why** the gene correlates. A rank 1 gene could be:
1. **A driver** — causes the fate decision
2. **An early responder** — reacts quickly to an upstream signal
3. **A strong bystander** — highly expressed gene whose kinetics happen to correlate

### How to Identify a True Early Marker

A gene qualifies as an early marker if it satisfies all three:

| Criterion | What to look for |
|---|---|
| Active before commitment | Appears in NSC or TAP driver CSV at high rank |
| Dynamically regulated | Clean almond shape in phase portrait |
| Discriminates fate | Appears in SHAP top 20 (trigger_genes.csv) |

**The sweet spot is low velocity rank + low SHAP rank** — the gene is both kinetically active during the TAP transition AND predictive of which fate the cell takes.

| Gene | TAP vel rank | SHAP rank | Interpretation |
|---|---|---|---|
| Bcl11a | 1 | 6 | Strong driver candidate — dynamic in TAPs AND discriminates fate |
| Nfib | 3 | 17 | Same — both analyses agree |
| Epha5 | 10 | none | Active in TAPs but doesn't discriminate OL vs NB — likely general progenitor gene |
| Stmn2 | none | 1 | Discriminates fate strongly but not dynamic in TAPs — state marker, not a driver |

**A gene in TAP + SHAP top 20 is the gold standard.** This applies to Bcl11a, Nfib, and Meis2 for the NB arm.

---

## 7. Results (Cntl GSM8253796 vs CupRap GSM8253798)

### SHAP Top 20 × Velocity Driver Cross-Reference

```
gene,shap_rank,shap_value,velocity_column,velocity_rank,interpretation
Bcl11a,6,0.139,TAP,1,strongest early NB marker — TAP velocity rank 1
Bcl11a,6,0.139,Neuroblast,43,confirmed in committed NB too
Nfib,17,0.024,TAP,3,early NB — present in TAP and NSC
Nfib,17,0.024,NSC,58,detectable even at NSC stage
Meis2,5,0.187,TAP,5,strong NB fate discriminator
Meis2,5,0.187,Neuroblast,52,confirmed downstream
Grin2b,15,0.029,Neuroblast,10,NB commitment marker
Grin2b,15,0.029,TAP,41,weak TAP signal
Fa2h,16,0.028,COP,1,strongest OL-lineage velocity driver
Fa2h,16,0.028,OL,8,confirmed in mature OL
Tspan2,19,0.021,OL,3,early OL surface marker
Gjc3,10,0.070,OL,22,OL marker but weaker in combined Cntl+CupRap analysis
Gjc3,10,0.070,OPC,21,present across OL lineage
Igfbpl1,2,0.440,TAP,69,weak TAP velocity signal — not a robust cell-intrinsic OL driver
Stmn2,1,0.711,none,none,stable NB state marker — not dynamically induced
```

### Cntl vs CupRap: What Changed in Velocity Rankings

Comparing CupRap-only placeholder run to real Cntl + CupRap:

```
gene,previous_velocity_rank,new_velocity_rank,column,interpretation
Bcl11a,6,1,TAP,strengthened — constitutive early NB driver in both conditions
Nfib,2,3,TAP,stable — robust early NB marker in both conditions
Fa2h,19,1,COP,strengthened — now top COP driver across both conditions
Tspan2,16,3,OL,strengthened — more prominent OL marker with Cntl added
Gjc3,3,22,OL,weakened — was CupRap-specific; less prominent with Cntl included
Gjc3,24,66,COP,weakened — same pattern confirms CupRap-dependence
Igfbpl1,99,69,TAP,marginal — not a robust velocity driver in either run
```

**Gjc3 dropping from OL rank 3 to rank 22** when Cntl is added suggests its OL dynamics are partially CupRap-dependent — stronger when microglia are activated (IGF1/OSM signaling). This directly connects to the project's central question about microglial signaling driving TAP→OL commitment.

**Fa2h and Tspan2 strengthening** with Cntl added confirms they are constitutive OL commitment markers, present regardless of treatment condition.

**Bcl11a at TAP rank 1** in both conditions confirms it as the dominant early NB fate marker, constitutive and condition-independent.

### Final Summary: Best Cross-Validated Early Markers

```
gene,fate,shap_rank,shap_value,best_velocity_rank,velocity_column,condition_dependence,conclusion
Bcl11a,NB,6,0.139,1,TAP,constitutive,strongest early NB marker — both analyses agree
Nfib,NB,17,0.024,3,TAP,constitutive,early NB — detectable at NSC stage
Meis2,NB,5,0.187,5,TAP,constitutive,strong NB fate discriminator
Fa2h,OL,16,0.028,1,COP,constitutive,strongest OL-lineage velocity driver
Tspan2,OL,19,0.021,3,OL,constitutive,early OL surface marker
Gjc3,OL,10,0.070,22,OL,CupRap-enriched,OL marker but condition-dependent — weaker in Cntl
Stmn2,NB,1,0.711,none,none,constitutive,stable NB state marker — not an early driver
Igfbpl1,OL,2,0.440,69,TAP,none,SHAP high but velocity weak — not a cell-intrinsic OL driver
```

### OL-Fate Velocity Drivers Across the Lineage

Checking Ptprz1 and OL lineage genes across TAP, OPC, COP, OL columns:

**Ptprz1 appears only in TAP (rank 10) and nowhere downstream.** Its splicing dynamics are active during the TAP→OPC transition and have already stabilized by the time the cell becomes an OPC. This is the definition of an early commitment marker — being switched on in TAPs heading toward OL fate before any canonical OPC marker is visible.

OPC top drivers are dominated by axon guidance and adhesion genes (Ntn1, Spry4, Ntm, F3) rather than canonical OPC markers like Pdgfra or Olig1. The OPC velocity field captures structural remodeling dynamics.

COP is where OL commitment first becomes clearly visible in velocity: Fa2h rank 1 (SHAP 16), Olig2 rank 4, Cldn11 rank 19.

OL has the strongest and cleanest signal of any cell type: Olig2 corr=98.45, Mag rank 5, Fa2h rank 8, Tspan2 rank 3 (SHAP 19).

OL-associated velocity drivers per cell type (independent ranked lists — ranks are NOT comparable across columns):

| Cell type | Gene | Rank within cell type | Note |
|---|---|---|---|
| TAP | Ptprz1 | 10 | Only OL-associated gene with strong TAP velocity signal |
| TAP | Smpd3 | 16 | Myelin sphingolipid metabolism |
| OPC | Enpp2 | 16 | OL-associated; top OPC drivers are mostly axon guidance genes |
| COP | Fa2h | 1 | Strongest COP driver; SHAP rank 16 |
| COP | Olig2 | 4 | Canonical OL TF |
| COP | Cldn11 | 19 | Myelin tight junction |
| OL | Olig2 | 1 | Corr=98, dominant OL driver |
| OL | Tspan2 | 3 | SHAP rank 19 |
| OL | Mag | 5 | Myelin-associated glycoprotein |
| OL | Fa2h | 8 | Also top COP driver — spans both stages |

These genes are active in their respective cell types independently. Ptprz1 and Enpp2 are unrelated genes — they do not form a cascade or handoff.

### Key Biological Finding

**Ptprz1 is the only gene showing velocity dynamics specifically at the TAP→OPC transition.** Everything else only becomes prominent at COP or OL stage. This makes Ptprz1 the best candidate for an early OL-fate marker in TAPs — even though SHAP did not rank it, because it does not discriminate OL vs NB at steady-state expression (it marks OPC-committed cells regardless of how they got there).

**The NB arm has one dominant early driver (Bcl11a, TAP rank 1). The OL arm has a staged cascade** — Ptprz1 at TAP stage, Fa2h/Olig2 at COP stage, Olig2/Tspan2/Mag at OL stage. OL fate commitment appears to be a gradual multi-step process rather than a single transcriptional switch.

```
gene,fate,tap_rank,cop_rank,ol_rank,shap_rank,notes
Ptprz1,OL,10,none,none,none,only TAP-stage OL velocity signal — early commitment marker
Smpd3,OL,16,none,none,none,myelin sphingolipid metabolism — early TAP signal
Fa2h,OL,none,1,8,16,constitutive — strongest cross-validated OL driver
Tspan2,OL,none,none,3,19,constitutive OL surface marker
Olig2,OL,none,4,1,none,canonical OL TF — strongest OL velocity signal (corr=98)
Mag,OL,none,none,5,none,late OL maturation marker
Bcl11a,NB,1,none,none,6,dominant early NB driver — TAP rank 1
Nfib,NB,3,none,none,17,early NB — present in NSC too
Meis2,NB,5,none,none,5,strong NB fate discriminator
```

---

## 8. SHAP Direction Analysis — Critical Methodological Finding

### The Problem

SHAP measures discrimination *magnitude*, not direction. A high SHAP score means a gene strongly distinguishes OL from NB cells — but the model can use that gene as a **positive** OL predictor (high gene → OL fate) OR as a **negative** OL predictor (high gene → NB fate, low gene → OL fate).

The project's goal is to find **positive OL-leaning markers** — genes upregulated when a TAP commits to OL fate. Negative markers (NB markers used in reverse) are biologically interesting but do not meet the deliverable.

### Method

For each of the top 20 SHAP genes, mean expression was compared between OL lineage (OPC + COP + OL) and Neuroblasts:
- **POSITIVE_OL**: ol_mean > 1.5 × nb_mean
- **NEGATIVE_OL (NB marker)**: nb_mean > 1.5 × ol_mean
- **AMBIGUOUS**: similar expression

### Results

```
shap_rank,gene,shap_value,ol_mean,nb_mean,tap_mean,ol_pct,nb_pct,ol_to_nb_ratio,direction
1,Stmn2,0.711,0.071,2.965,0.897,5.1,96.1,0.02,NEGATIVE_OL
2,Igfbpl1,0.440,0.054,2.199,1.119,4.8,93.2,0.02,NEGATIVE_OL
3,Pllp,0.258,1.836,0.016,0.018,90.2,1.7,114.01,POSITIVE_OL
4,Dlx6os1,0.209,0.050,2.294,0.840,4.0,93.3,0.02,NEGATIVE_OL
5,Meis2,0.187,0.275,2.934,1.668,27.5,98.2,0.09,NEGATIVE_OL
6,Bcl11a,0.139,0.035,1.521,0.972,3.5,84.5,0.02,NEGATIVE_OL
7,Sox11,0.100,0.218,2.999,2.763,17.8,97.2,0.07,NEGATIVE_OL
8,Celf4,0.100,0.062,1.613,0.570,6.7,86.9,0.04,NEGATIVE_OL
9,Tmsb10,0.094,0.802,3.754,2.780,53.1,97.9,0.21,NEGATIVE_OL
10,Gjc3,0.070,1.955,0.014,0.010,89.1,1.4,137.08,POSITIVE_OL
11,Meg3,0.067,0.214,2.244,0.232,16.0,72.5,0.10,NEGATIVE_OL
12,Arx,0.060,0.035,1.487,1.180,3.6,77.8,0.02,NEGATIVE_OL
13,Dock10,0.042,1.231,0.023,0.016,79.4,2.7,52.63,POSITIVE_OL
14,Cryab,0.036,2.707,0.091,0.072,90.1,8.5,29.87,POSITIVE_OL
15,Grin2b,0.029,0.017,0.602,0.241,2.1,42.3,0.03,NEGATIVE_OL
16,Fa2h,0.028,1.207,0.010,0.008,75.1,1.1,116.07,POSITIVE_OL
17,Nfib,0.024,0.851,2.888,2.416,69.7,96.6,0.29,NEGATIVE_OL
18,Cnp,0.023,3.171,0.186,0.244,93.5,17.9,17.02,POSITIVE_OL
19,Tspan2,0.021,2.211,0.039,0.058,86.3,4.1,56.93,POSITIVE_OL
20,Hmgn2,0.021,0.731,1.948,2.962,65.3,88.1,0.38,NEGATIVE_OL
```

**Summary:** 7 positive OL markers, 13 negative OL markers (NB markers used in reverse), 0 ambiguous.

### The 7 True Positive OL-Leaning Markers

These are the genes that genuinely mark OL commitment by their **presence** (not absence):

| SHAP rank | Gene | OL pct | NB pct | OL/NB ratio | Notes |
|---|---|---|---|---|---|
| 3 | Pllp | 90% | 2% | 114× | Plasmolipin — myelin lipid raft component |
| 10 | Gjc3 | 89% | 1% | 137× | Connexin 30.2 — gap junction protein; CupRap-enriched in velocity |
| 13 | Dock10 | 79% | 3% | 53× | Cytoskeletal regulator |
| 14 | Cryab | 90% | 9% | 30× | αB-crystallin — stress response, myelin maintenance |
| 16 | Fa2h | 75% | 1% | 116× | Fatty acid 2-hydroxylase — myelin sphingolipid synthesis |
| 18 | Cnp | 94% | 18% | 17× | 2',3'-cyclic nucleotide phosphodiesterase — myelin marker |
| 19 | Tspan2 | 86% | 4% | 57× | Tetraspanin 2 — early oligodendrocyte surface marker |

All 7 are biologically coherent **myelin/membrane structural genes** — consistent with the metabolic/structural priming hypothesis identified earlier (vs the single dominant TF that NB fate uses with Bcl11a).

### The 13 Negative OL Markers (NB Markers in Reverse)

These genes have high SHAP rank because they are **strong Neuroblast markers** — the model uses their *absence* as an OL signal. Their presence indicates a TAP is heading to NB fate, not OL fate.

Includes the highest SHAP genes: Stmn2 (#1), Igfbpl1 (#2), Dlx6os1 (#4), Meis2 (#5), Bcl11a (#6), Sox11 (#7).

### Conclusion: The Refined Panel

The project's deliverable should be the **7 positive OL-leaning markers**:
1. Pllp
2. Gjc3
3. Dock10
4. Cryab
5. Fa2h
6. Cnp
7. Tspan2

The other 13 are valid lineage discriminators but do not answer the project's question.

### Decision: Keep the Symmetric Classifier

Retraining with NB markers excluded was considered and rejected. The asymmetric filter (removing only NB-biased features while keeping all OL-biased features) imposes a direction-imposing constraint that changes what's being measured — from "OL vs NB discrimination" to "OL identity expression." This is methodologically different from the canonical marker removal (which symmetrically excludes textbook markers from both classes).

The symmetric classifier is retained as the source of truth. The 7 positive OL-leaning markers above are derived by post-hoc filtering of the SHAP top 20 by direction — a description of what the panel contains, not a new model.

---

## 9. Plan: Expand to All 10 Samples

### Why

The current velocity analysis on 2 samples (10,246 cells) shows a clear early NB driver (Bcl11a at TAP velocity rank 1) but **no comparable early OL driver in TAPs**. The only TAP-stage OL signal is Ptprz1 at rank 10 — moderate at best.

Three possibilities for the absent early OL signal:

1. **Sample size** — Only ~5,200 TAPs total in 2 samples. Of those, only a fraction are actively committing to OL fate (most go to the default NB fate). The signal may be present but underpowered.
2. **Biology** — OL commitment may genuinely lack a dominant early TF switch, instead relying on gradual structural/metabolic priming visible only at OPC/COP stage.
3. **Methodology** — The 2-sample velocity neighborhood graph differs from the 10-sample UMAP embedding, introducing noise.

Adding all 10 samples addresses (1) and (3) simultaneously:
- Doubles TAP cell count (more statistical power)
- Includes 5 CupRap samples (more cells actively transitioning to OL)
- Velocity neighborhood graph matches the 10-sample UMAP exactly

If a positive OL-leaning TAP marker exists, the 10-sample run is the most likely place to find it. If none emerges, that supports interpretation (2) as a real biological finding — OL commitment is structurally different from NB commitment.

### Sample Inventory

| GSM | Strain | Treatment | Timepoint | Cells | Priority |
|---|---|---|---|---|---|
| GSM8253796 | NesCre | Cntl | 3wks | 6,035 | ✓ done |
| GSM8253798 | NesCre | CupRap | 3wks | 4,211 | ✓ done |
| GSM8253799 | NesCre | CupRap | 3wks | 10,567 | running on EC2 |
| GSM8253797 | NesCre | Cntl | 3wks | 4,710 | next |
| GSM8253792 | CD1 | Cntl | 0wks | 4,376 | cross-strain |
| GSM8253793 | CD1 | Cntl | 3wks | 4,061 | cross-strain |
| GSM8253794 | CD1 | CupRap | 3wks | 3,550 | cross-strain |
| GSM8253795 | CD1 | CupRap | 3wks | 2,393 | cross-strain |
| GSM8647352 | CD1 | CupRap | 0wks | 4,134 | early timepoint |
| GSM8647353 | CD1 | CupRap | 0wks | 3,928 | early timepoint |

Total: 47,965 cells across 10 samples. Approximately 50,000 cells will provide ~5x the TAP signal of the current 2-sample run.

### Expected Storage and Compute

- FASTQ data per sample: ~80–200 GB (kept only during kb count, deleted after)
- kb count output (`counts_unfiltered/`): ~150 MB per sample
- Velocity pipeline RAM: ~30 GB peak with all 10 samples
- Velocity pipeline time: ~5 hours (dominated by `recover_dynamics`)
- AWS cost estimate: ~$10–15 total on r6i.2xlarge

### What to Look For After the 10-Sample Run

Re-run the `velocity_drivers_TAP.csv` analysis and look specifically for:

1. **Positive OL markers appearing in TAP velocity** — any of Pllp, Gjc3, Dock10, Cryab, Fa2h, Cnp, Tspan2 showing up in the TAP driver list at high rank would be the deliverable: a confirmed early OL-fate marker.
2. **Ptprz1 rank stability** — if it stays high or strengthens, it's confirmed as the dominant early OL marker. If it drops, it was sample-specific noise.
3. **Cntl vs CupRap differential in TAP** — genes that appear as TAP drivers only in CupRap samples are condition-dependent OL commitment drivers, directly linked to microglial IGF1/OSM signaling.

If no positive OL markers emerge in TAP velocity even with 10 samples and 5x the TAP cell count, the finding is biological: **OL commitment lacks a dominant early transcriptional driver** — a publishable mechanistic insight in itself.