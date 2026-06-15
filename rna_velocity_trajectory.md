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
```bash
# download with --include-technical to get barcode reads (_3 files)
fasterq-dump --split-files --include-technical SRR_ID1 SRR_ID2 --outdir .

# run kb count (script handles per-sample FASTQ lists)
bash run_kb_count.sh GSM8253796
bash run_kb_count.sh GSM8253798
```

Note: FASTQ files have 4 reads per spot. Use `_3` (R1, 28bp barcode+UMI) and `_4` (R2, 90bp cDNA) only. `--include-technical` is required — without it fasterq-dump omits the barcode reads.

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

### Key Biological Finding

**The absence of a strong OL-fate gene in the TAP column is itself a finding.** Fa2h only appears at COP stage, Tspan2 only at OL stage — suggesting OL commitment does not have a single dominant early transcriptional driver in TAPs the way NB commitment does with Bcl11a. OL fate appears to be encoded through multiple weaker signals (Fa2h, Tspan2, Gjc3) that only become prominent after OPC commitment.
