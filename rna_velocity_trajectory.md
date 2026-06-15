# RNA Velocity & Trajectory Analysis Plan

## Context
Currently, the project uses a static snapshot of gene expression to classify cell types and uses an XGBoost model (with SHAP scores) to identify "trigger genes" that predict whether a Transit Amplifying Progenitor (TAP) will commit to the Neuroblast or Oligodendrocyte (OL) lineage.

While SHAP scores are excellent at telling us *which* genes distinguish the lineages, static data makes it difficult to definitively say *when* those genes act (i.e., are they early drivers or late consequence markers?). 

## Can RNA Velocity Identify "Early Genes"?
**Yes. This is exactly where RNA velocity excels.** 
By modeling the ratio of unspliced (newly transcribed) to spliced (mature) mRNA, RNA velocity acts as a "molecular clock." It doesn't just show you clusters; it infers a continuous **latent time** and directional flow. You can explicitly sort genes by when their transcription spikes along this latent time axis, allowing you to clearly separate "early driver genes" from "late mature genes."

## Questions We Can Answer (Comparing NesCre Cntl vs. NesCre CR)
By processing the NesCre control and treated (Cuprizone/CR) samples, we can answer:
1. **Directionality:** Does the transit from TAP to OL actively accelerate or stall under Cuprizone treatment?
2. **Early Gene Identification:** Which genes initiate the transcriptional cascade towards the OL lineage before the cell even looks like an OL?
3. **Perturbation Response:** Do the "early genes" activate differently or at different times in the CR sample compared to the control?
4. **Trigger Gene Validation:** Do the top SHAP "trigger genes" actually precede the lineage split in developmental time?

---

## The Action Plan

### Phase 1: Raw Data Processing (The Bottleneck)
*   **Action:** Download FASTQ files from the Sequence Read Archive (SRA) for `GSM8253796` (Cntl) and `GSM8253798` (CR).
*   **Processing:** Use `kb-python` (kallisto|bustools) instead of CellRanger to align reads to a reference mouse genome and quantify spliced vs. unspliced transcripts. 
    *   **Step 1: Download Reference Files (run from `sra_runs/`):**
        ```sh
        cd sra_runs
        wget ftp://ftp.ensembl.org/pub/release-110/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
        wget ftp://ftp.ensembl.org/pub/release-110/gtf/mus_musculus/Mus_musculus.GRCm39.110.gtf.gz
        gunzip *.gz
        ```
    *   **Step 2: Build the Reference Index (run from `sra_runs/`):**
        ```sh
        # All reference files and output live in sra_runs/
        # Script: run_kb_ref.sh
        kb ref \
          -i index.idx \
          -g t2g.txt \
          -f1 cdna.fa \
          -f2 intron.fa \
          -c1 cdna_t2c.txt \
          -c2 intron_t2c.txt \
          --workflow lamanno \
          --kallisto /opt/homebrew/bin/kallisto \
          --bustools /opt/homebrew/bin/bustools \
          Mus_musculus.GRCm39.dna.primary_assembly.fa \
          Mus_musculus.GRCm39.110.gtf
        ```
    *   **Step 3: Run the Count Command (Pooling 3 runs for a single biological sample):**
        *   Note: FASTQ files have 4 reads per spot (`_1`=I1 10bp, `_2`=I2 10bp, `_3`=R1 28bp barcode+UMI, `_4`=R2 90bp cDNA). Use `_3` and `_4` only.
        ```sh
        # Script: run_kb_count.sh (run from repo root; script cd's into sra_runs/)
        cd sra_runs
        kb count \
          -i index.idx \
          -g t2g.txt \
          -x 10xv3 \
          -o kb_output_GSM8253796 \
          -c1 cdna_t2c.txt \
          -c2 intron_t2c.txt \
          --workflow lamanno \
          --loom \
          --kallisto /opt/homebrew/bin/kallisto \
          --bustools /opt/homebrew/bin/bustools \
          SRR28912872_3.fastq SRR28912872_4.fastq \
          SRR28912873_3.fastq SRR28912873_4.fastq \
          SRR28912874_3.fastq SRR28912874_4.fastq
        ```
*   **Output:** One `.loom` file per sample containing `spliced` and `unspliced` as matrix layers, written to `sra_runs/kb_output_<GSM>/counts_unfiltered/adata.loom`.
*   **Validation:** Verify unspliced fraction is 30–60% of aligned reads (expected range for single-nucleus brain data; confirmed ~41% for GSM8253796: spliced=99M reads, unspliced=68M reads).

---

## How the Two Loom Files Produce Velocity

Each loom file contains a cells × genes matrix with two layers:

| Layer | Source | Meaning |
|-------|--------|---------|
| `spliced` | cDNA reads (exonic) | Mature mRNA — already processed |
| `unspliced` | intronic reads | Pre-mRNA — newly transcribed, not yet spliced |

scVelo models splicing kinetics per gene using this ODE system:

```
du/dt = α - βu       (unspliced: produced at rate α, consumed at rate β)
ds/dt = βu - γs      (spliced: produced from unspliced, degraded at rate γ)
```

By fitting α (transcription), β (splicing), γ (degradation) for each gene, scVelo predicts where `s` is heading — that is the **velocity vector** for each cell.

**Merge + velocity pipeline:**

```
GSM8253796 loom (Cntl)  ──┐
                           ├── merge on barcodes ──> combined AnnData
GSM8253798 loom (CR)    ──┘                          spliced + unspliced layers
                                                      +
existing h5ad (subset to 2 samples) ──────────────> UMAP, cell_type labels

scVelo fits kinetics per gene → velocity vector per cell
Project vectors onto UMAP → directional arrows
```

Both samples stay in one AnnData with a `sample_id` column. scVelo fits one model across all cells; plots are stratified by `sample_id` to compare Cntl vs CR side by side.

**What the comparison answers:**
- Does CR **accelerate or stall** TAP→OL commitment? (arrow magnitude/direction change)
- Do SHAP trigger genes show velocity earlier in latent time under CR vs Cntl?

---

### Phase 2: Velocity Inference & Integration
*   **Action:** In Python (using the `scvelo` library), load the `.loom` files from `sra_runs/kb_output_<GSM>/counts_unfiltered/adata.loom` and merge them with your existing, fully-annotated `processed.h5ad` file.
*   **Processing:** Calculate RNA velocity using scVelo's dynamical model.
*   **Output:** 
    *   Velocity stream plots mapped onto your existing UMAP. (You will visually see arrows pointing from NSCs to TAPs to OLs/Neuroblasts).
    *   A calculation of "Latent Time" for every cell.
*   **Validation:** The velocity arrows should broadly agree with known biology (e.g., arrows should not point from mature OLs back to NSCs).

### Phase 3: Identifying Early Drivers (The Core Goal)
*   **Action:** Use scVelo to identify "dynamical driver genes."
*   **Processing:** scVelo tests which genes' expression patterns best correlate with the inferred directional vector. We will then plot these genes along the continuous Latent Time axis.
*   **Output:** 
    *   A ranked list of velocity driver genes.
    *   Expression heatmaps cascaded by time: Visually showing Phase 1 (early activation genes), Phase 2 (mid-transit genes), and Phase 3 (mature terminal genes).
*   **Validation:** Cross-reference the identified "early/driver" genes with known literature on OL and Neuroblast development to ensure biological sanity.

### Phase 4: Integration with XGBoost SHAP Scores
*   **Action:** Merge the findings from the classifier pipeline with the velocity pipeline.
*   **Processing:** 
    *   Map your current `trigger_genes.csv` (ranked by SHAP) against the scVelo Latent Time cascade.
    *   Determine if the highest SHAP-scoring genes are primarily "early" (drivers) or "late" (state markers).
*   **Output:** A finalized, highly-confident list of "Early Committment Trigger Genes" that are supported by BOTH predictive machine learning (SHAP) AND physical transcriptional kinetics (RNA velocity).
*   **Validation:** Do the dynamics of these specific trigger genes change drastically between the Cntl and CR (treated) samples? If a gene's velocity goes up in CR, it validates it as a biologically active participant in the treatment response.