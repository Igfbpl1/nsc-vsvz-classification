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
*   **Output:** Two `.loom` files containing `spliced` and `unspliced` count matrices.
*   **Validation:** Verify that the unspliced fraction is reasonably high (typically 15-25% of total reads in single-nucleus/single-cell brain data).

### Phase 2: Velocity Inference & Integration
*   **Action:** In Python (using the `scvelo` library), load the `.loom` files and merge them with your existing, fully-annotated `processed.h5ad` file.
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