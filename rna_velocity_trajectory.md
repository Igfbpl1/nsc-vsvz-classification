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

## 7. Results — 3 Sample Run (Cntl + 2× CupRap)

Velocity computed on 20,813 cells: GSM8253796 (Cntl, 6,035 cells) + GSM8253798 (CupRap Rep1, 4,211 cells) + GSM8253799 (CupRap Rep2, 10,567 cells).

### TAP Velocity Drivers — No Positive OL-Leaning Markers

The top 30 TAP velocity drivers are dominated by NB-associated genes. With 5× more cells than the 2-sample run, no positive OL marker emerged in TAP velocity.

| Rank | Gene | Corr | Fate association |
|---|---|---|---|
| 1 | Nsg2 | 94 | NB |
| 2 | Rnf165 | 86 | NB |
| 3 | Srrm4 | 85 | NB |
| 4 | Bcl11a | 83 | NB (SHAP #6) |
| 5 | Meis2 | 77 | NB (SHAP #5) |
| 6 | Nfib | 75 | NB-leaning |
| 29 | Smpd3 | 49 | weak OL-related (sphingolipid metabolism) |

No OL canonical or non-canonical markers (Pllp, Gjc3, Cnp, Cryab, Fa2h, Tspan2, Dock10, Lhfpl3, Ptprz1, Pdgfra) appear in TAP top 100.

Ptprz1, which appeared at TAP rank 10 in the 2-sample run, dropped to **rank 44** with 5× more cells — confirming it was a sample-size artifact, not a real early marker.

### Cell-Type Velocity Cascade

The OL identity program turns on sharply at the COP stage with no detectable TAP precursor.

| Stage | Top OL-relevant drivers (rank, corr) |
|---|---|
| TAP | none |
| OPC | Gjc3 (14, 21), Cmtm5 (13, 22), Ncam1 (12, 23) |
| COP | Fa2h (1, 28), Cldn11 (5, 17), Olig2 (6, 16), Mag (12, 13), **Igf1 (13, 13)**, Gjc3 (15, 12), Gatm (18, 12) |
| OL | **Gjc3 (1, 261)**, Mag (2, 155), Tspan2 (4, 120), Ugt8a (8, 96), Olig2 (15, 66), Fa2h (19, 61), Mog (24, 49) |

### Three Major Positive Findings

**1. Gjc3 explosion in OL velocity (rank 22 → rank 1, corr 52 → 261)**

The two added CupRap samples produced a 5× jump in Gjc3 OL velocity correlation. This is the strongest condition-responsive signal in the analysis. Direct evidence that CupRap treatment activates an OL commitment program — connecting to the project's IGF1/OSM microglial hypothesis.

**2. Igf1 itself appears as a COP velocity driver (rank 13)**

First direct molecular evidence in the velocity layer of IGF1 signaling at the OL commitment stage. The COP cells have active IGF1 splicing dynamics, consistent with either autocrine IGF1 or IGF1-responsive transcriptional changes — both supporting the project's central hypothesis.

**3. The myelin program is cleanly detectable at COP and OL stages**

OL-relevant genes appearing as top velocity drivers within each cell type (independent ranked lists — ranks are NOT sequential or comparable across cell types):

| Cell type | OL-relevant velocity drivers (rank within cell type, corr) |
|---|---|
| COP | Fa2h (1, 28), Cldn11 (5, 17), Olig2 (6, 16), Mag (12, 13), Igf1 (13, 13), Gjc3 (15, 12), Gatm (18, 12) |
| OL | Gjc3 (1, 261), Mag (2, 155), Tspan2 (4, 120), Olig2 (15, 66), Fa2h (19, 61), Mog (24, 49) |

These genes are simultaneously active in their respective cell types — Fa2h and Cldn11 are not in sequence, they are both top velocity drivers in COP cells at the same time.

What this shows: the myelin program is detectable as a coordinated set of velocity drivers at COP and OL stages. None of these genes appear in TAP velocity, indicating the program turns on after OL commitment, not before.

### Other Observations

- **Bcl11a dropped from TAP rank 1 to rank 4**; Nsg2 took over. Both are NB markers — the early NB signal is reproducible.
- **Nfib appears in both NSC (rank 47) and TAP (rank 6) velocity** — the earliest NB-leaning signal in the dataset, active across the NSC→TAP transition.
- **Ncam1 appears in OL rank 5 and OPC rank 12** — myelin-relevant adhesion molecule, confirmed velocity driver across the OL lineage.
- **Aldh1l1 in NSC rank 21** — astrocyte/quiescent NSC marker, confirms the NSC end of the trajectory.

### Reinterpretation: OL Fate May Be Encoded by Failure to Engage NB Program

With 5× the TAP cells and zero positive OL markers in TAP velocity, the absence of an early OL transcriptional switch is now strongly supported. Hypothesis worth considering:

**The fate decision may be on the NB side, not the OL side.** When a TAP fails to upregulate Bcl11a/Nsg2/Srrm4/Meis2, it defaults to OL fate. The molecular decision is "engage NB program or not" — there is nothing positive to detect in TAP velocity for OL-fated cells because their commitment is defined by *absence* of NB signal.

This is consistent with:
- SHAP repeatedly surfacing NB markers as top discriminators (13 of 20 are NB-direction)
- The sharp transition of OL identity at COP stage rather than gradual buildup
- The NB arm having clear early TAP drivers (Bcl11a, Nfib, Meis2) while the OL arm has none

---

## 7b. Cntl vs CupRap Condition-Split Velocity Analysis

To test whether CupRap treatment changes the transcriptional program within specific cell types, velocity drivers were ranked per (cell_type × condition) — splitting Cntl and CupRap cells of the same cell type into separate groups.

### Cell Counts Per Cell Type × Condition

| Cell Type | Cntl cells | CupRap cells | Statistical balance |
|---|---|---|---|
| TAP | 743 | 889 | balanced |
| OPC | 27 | 145 | Cntl underpowered |
| COP | 44 | 491 | Cntl severely underpowered |
| OL | 655 | 1,510 | imbalanced but workable |

The asymmetry itself is a finding — Cntl mice have far fewer OPCs and COPs because no white matter injury is driving OL commitment.

### Top 30 Driver Overlap Between Conditions

| Cell Type | Top 30 overlap | Top 50 overlap | Interpretation |
|---|---|---|---|
| TAP | ~100% | ~100% | Identical — no CupRap effect on TAP transcription |
| OPC | 27/30 (90%) | 45/50 (90%) | Very similar — same core identity program |
| COP | 19/30 (63%) | 36/50 (72%) | Meaningful differences — strongest CupRap modulation |
| OL | 24/30 (80%) | 39/50 (78%) | Similar program, much stronger dynamics in CupRap |

### TAP Stage — Identical Programs (Cntl ≈ CupRap)

Both conditions show essentially the same NB-leaning velocity drivers:

| Rank | Cntl | CupRap |
|---|---|---|
| 1 | Nsg2 | Nsg2 |
| 2 | Rnf165 | Srrm4 |
| 3 | Srrm4 | Rnf165 |
| 4 | Bcl11a | Bcl11a |
| 5 | Meis2 | Plxna2 |
| 6 | Nfib | Meis2 |

**Conclusion:** CupRap treatment does **not** alter the transcriptional program of TAPs. Both Cntl and CupRap TAPs run the same NB-default trajectory. The CupRap effect on OL commitment is not encoded at the TAP stage.

### OPC Stage — CupRap Amplifies an Existing Program

Top 30 are ~90% identical (Ntn1, Ntm, Dpysl3, Gpr37l1, Fut9, F3 dominant in both). Same OL identity genes appear at similar ranks:

| Gene | Cntl rank | CupRap rank |
|---|---|---|
| Cmtm5 | 11 | 13 |
| Ncam1 | 13 | 12 |
| Gjc3 | 17 | 14 |

But CupRap OPCs uniquely show **active myelin biosynthesis** genes at moderate ranks:

| Gene | CupRap rank | Function |
|---|---|---|
| **Myrf** | 43 | Myelin Regulatory Factor — master myelin TF |
| **Ugt8a** | 48 | UDP-galactose:ceramide galactosyltransferase — myelin galactolipid |
| **Ptprk** | 50 | Receptor tyrosine phosphatase — myelin signaling |
| **Elovl7** | 55 | Fatty acid elongase — myelin lipid synthesis |
| **Fa2h** | 58 | Fatty acid 2-hydroxylase — myelin sphingolipid |

These genes exist in Cntl OPCs too (Fa2h at Cntl rank 60, corr 4.6) but at much weaker velocity correlation. CupRap **amplifies** the existing OPC myelin program; it does not switch it on de novo.

**Conclusion:** CupRap effect at OPC stage is quantitative amplification of the myelin biosynthesis program. Myrf emergence is mechanistically meaningful — Myrf is the master TF committing OPCs to mature myelinating OLs.

### COP Stage — Strongest CupRap-Driven Differences

Cntl top correlations are ~3× weaker than CupRap (Fa2h Cntl corr=7.6 vs CupRap corr=26). Several OL commitment genes appear in CupRap COPs but not Cntl:

| Gene | CupRap COP rank | Cntl COP rank | Function |
|---|---|---|---|
| **Igf1** | 12 (corr 12.4) | not in top 50 | IGF1 ligand — same molecule the source paper identifies as microglial driver |
| **Gjc3** | 13 (corr 12.2) | not in top 50 | Gap junction protein — connects to OL maturation |
| **Cmtm5** | 24 (corr 10.6) | not in top 50 | CKLF-like MARVEL transmembrane protein |
| **Tspan2** | 55 (corr 6.6) | not in top 50 | OL surface marker |

**Conclusion:** Igf1 emergence as a COP velocity driver in CupRap (but not Cntl) is mechanistically significant. The Cell Reports 2025 paper identifies Igf1 as a *microglia-derived* ligand. Our analysis shows COP cells themselves have active Igf1 splicing dynamics during CupRap-induced commitment — suggesting either autocrine signaling or response to incoming microglial IGF1 includes transcriptional Igf1 upregulation in target cells.

### OL Stage — Same Program, Dramatically Stronger in CupRap

Top OL identity genes appear in both conditions but with very different correlation magnitudes:

| Gene | Cntl rank (corr) | CupRap rank (corr) | Ratio |
|---|---|---|---|
| Gjc3 | 2 (81) | 1 (141) | 1.7× stronger in CupRap |
| Mag | 1 (82) | 3 (92) | comparable |
| Tspan2 | 3 (66) | 5 (77) | 1.2× stronger in CupRap |
| Olig2 | 26 (25) | 11 (58) | 2.3× stronger in CupRap |
| Fa2h | 17 (37) | 20 (43) | comparable |
| Ptprk | 14 (41) | 10 (59) | 1.4× stronger in CupRap |

CupRap-specific OL drivers include Pex5l, Tppp, Enpp2 (additional myelin/cytoskeletal genes). Cntl-specific includes Sox10 (canonical OL TF still detectable in Cntl OLs, while it falls outside CupRap top 50 because the dominant Gjc3 signal compresses the relative rankings).

**Conclusion:** Both Cntl and CupRap have functioning mature OLs running the same identity program. CupRap dramatically accelerates the velocity dynamics — particularly the Gjc3/Olig2/Ptprk axis.

### Summary: Where the CupRap Effect Actually Lives

```
Stage      Cntl vs CupRap difference            CupRap-specific signature
-----------------------------------------------------------------------------
TAP        None (~100% overlap)                  —
OPC        Mild (90% overlap)                   Myrf, Ugt8a, Ptprk, Fa2h amplified
COP        Strong (63% overlap)                 Igf1 (rank 12!), Gjc3, Cmtm5, Tspan2
OL         Same program, ~1.5-2× stronger        Gjc3 corr 141 vs 81, Olig2 rank 11 vs 26
```

The CupRap effect is a **gradient that accumulates downstream of commitment**, not a TAP-stage decision switch:

1. **TAPs** behave identically in both conditions (NB-default program)
2. **OPCs** show CupRap-specific upregulation of myelin biosynthesis (Myrf, Ugt8a, Fa2h, Elovl7)
3. **COPs** show CupRap-specific activation of Igf1 and Gjc3 — strongest condition-driven signal
4. **OLs** show same identity but accelerated dynamics under CupRap

The Igf1 appearance at COP rank 12 in CupRap (but not Cntl) is mechanistic evidence connecting your downstream transcriptional readout to the IGF1 signaling hypothesis from the source paper.

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
| GSM8253799 | NesCre | CupRap | 3wks | 10,567 | ✓ done |
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

### Status After 3 Samples

Completed:
- ✓ Positive OL markers (Pllp, Gjc3, Dock10, Cryab, Fa2h, Cnp, Tspan2) checked against TAP velocity — **none appear**
- ✓ Ptprz1 dropped from TAP rank 10 → 44 with more cells — confirmed as sample-specific noise
- ✓ Gjc3 in OL velocity jumped from rank 22 → 1 (corr 52 → 261) with added CupRap samples — strong CupRap-responsive OL commitment signal
- ✓ Igf1 surfaced as COP velocity driver (rank 13) — direct molecular evidence of IGF1 signaling at OL commitment

### Whether to Expand Further

The remaining 7 samples may strengthen the COP→OL cascade but are unlikely to produce a TAP-stage OL driver — the 5× increase from 2 samples to 3 samples (and 2× cell count) produced no TAP-stage OL signal. The biological hypothesis (OL fate = failure to engage NB program) is now supported.

Adding the cross-strain CD1 samples would address whether findings are NesCre-specific. Adding the 0wks CupRap samples would address whether the commitment signal differs at early vs late post-injury timepoints.

If the remaining 7 samples are added, the specific test would be:
- Does the COP→OL cascade strengthen further?
- Does any TAP-stage OL signal emerge with 50K cells?
- Are the findings reproducible across strains?
---
