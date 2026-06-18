# Next Steps

---

## What Has Already Been Built

Before extending the work — what exists is genuinely non-trivial. Most undergraduates don't touch RNA velocity. The current project has:

- A full scRNA-seq preprocessing pipeline (QC, normalization, HVG selection, Leiden clustering)
- A binary fate classifier with 99.85% accuracy (XGBoost + SHAP explainability)
- RNA velocity on spliced/unspliced counts built from raw FASTQs via kb-python
- A cross-validated biological finding: CupRap acts post-commitment on OPCs, not in TAPs
- A complete scientific narrative (PROJECT_NARRATIVE.md, rna_velocity_trajectory.md, GSE266687_summary.md)

---

## Immediate Next Steps (Inside This Project)

These can be done in the next few weeks and directly strengthen the existing work.

### 1. Complete the Velocity Run With All 10 Samples

Three GSMs are still pending:

| GSM | Strain | Treatment | Timepoint |
|---|---|---|---|
| GSM8253795 (SRR28912877) | CD1 | CupRap | 3wks |
| GSM8253797 (SRR28912867, SRR28912868) | NesCre | Cntl | 3wks |
| GSM8647352 (SRR31443698, SRR31443699) | CD1 | CupRap | 0wks |

Adding them increases statistical power and makes the CupRap vs Control comparison more robust. This is the lowest-hanging fruit in the project.

### 2. Pseudotime Analysis (Monocle3 or PAGA)

RNA velocity gives direction. Pseudotime gives a scalar — how far along each cell is in the trajectory. Running Monocle3 or PAGA on the TAP→OPC→COP→OL path would let you ask:

> At what pseudotime point does the CupRap-driven transcriptional change first appear?

This directly tests the post-commitment claim computationally with a third independent method. Libraries: `scvelo.tl.paga`, `cellrank`, or `pyroe` / Monocle3 via R.

### 3. Differential Expression in OPCs: CupRap vs Control

The paper shows TAPs are transcriptionally unchanged by CupRap (r = 0.995, Figure 2I). What has not been done is a formal differential expression analysis on **OPCs** stratified by condition. Running `sc.tl.rank_genes_groups` on the OPC population with `groupby='treatment'` would produce a gene list for what CupRap is doing to committed OPCs at the transcriptional level. This is one function call on data already in memory and could produce a new finding.

### 4. CellRank — Fate Probability From Velocity Dynamics

CellRank uses velocity vectors + Markov chains to compute the probability of each cell reaching each terminal state. Applied to TAPs, it gives a per-cell P(OL fate) derived purely from dynamics — completely independent of gene expression. Comparing CellRank P(OL) to XGBoost P(OL) per TAP is a strong cross-validation. If they agree, it closes the loop on the central argument. Library: `cellrank`.

---

## Medium-Term: Apply These Skills to a New Dataset

Taking everything built here and applying it to a related disease dataset makes the work comparative and more publishable.

### Multiple Sclerosis Datasets on GEO

Search GEO for `scRNA-seq "multiple sclerosis" oligodendrocyte` — several large human datasets exist. The OL vs NB classifier, retrained or applied zero-shot, could ask: do MS OPCs show the same velocity signature as CupRap OPCs? This connects a mouse model finding to human disease.

### Human V-SVZ

Postmortem human V-SVZ scRNA-seq datasets exist on GEO. Applying the trigger genes (Pllp, Gjc3, Fa2h, Tspan2, Dock10) to a human dataset tests cross-species conservation of the OL commitment signature — a real biological question that has not been answered.

### Spinal Cord Injury

Same OL remyelination biology, different injury context. The entire analytical pipeline applies directly, and the comparison would test whether the post-commitment CupRap effect is injury-model-specific or general.

---

## Deeper Techniques Worth Learning

Given what is already working, these are the natural next tools in roughly increasing order of difficulty:

| Technique | What It Adds | Library |
|---|---|---|
| CellRank | Fate probability from velocity — third independent fate method | `cellrank` |
| pySCENIC | Transcription factor regulon analysis — which TFs drive the OL program? | `pyscenic` |
| scVI / scANVI | Deep learning for batch correction + latent representation | `scvi-tools` |
| NicheNet / LIANA | Ligand-receptor communication modeling (extends the IGF1/OSM story) | `liana-py` |
| Squidpy | Spatial analysis — the Xenium data for this paper is on GEO | `squidpy` |

**CellRank** is the most immediately useful given scVelo is already running.

**pySCENIC** is probably the highest-impact for discovering something new: it would identify which transcription factors are regulatorily active in OPCs under CupRap but not in TAPs — a mechanistic answer to why the commitment line is where it is, and directly actionable for experimental follow-up.

---

## Where This Goes

### Science Fairs — This Work Is Competition-Ready

Regeneron ISEF and its regional feeders accept computational biology. This project has a clear question, a novel method, a validated result, and a biological interpretation. The CupRap post-commitment finding is a defensible, original claim. RNA velocity is not common in high school ISEF entries. The SHAP trigger gene list + velocity convergence on the same genes (Fa2h, Gjc3, Dock10) independently is a strong scientific argument.

### Preprint on bioRxiv

A preprint can be posted without a PI or institutional affiliation. The bar is: does it say something true and non-obvious? The trigger gene list, the OL vs NB classifier, and the velocity-based CupRap post-commitment finding meet that bar. A short computational resource paper is a realistic target, especially after completing the OPC differential expression and CellRank analyses.

### Summer Research Programs

Programs like RSI (Research Science Institute), PRIMES at MIT, and university HSSP programs look for students who have already done independent research and want mentorship or lab access to extend it. The GitHub repository and this project serve as a portfolio demonstrating real analytical capability.

### Contact the Miller Lab

Freda Miller (freda.miller@msl.ubc.ca) is at UBC. The paper is open access. The velocity analysis performed here found something the paper does not report — specifically the post-commitment specificity of the CupRap effect demonstrated through velocity driver analysis. A short, respectful email summarizing the finding and the methods is not unusual; computational re-analyses of published data are welcomed when they add something. Worst case: no reply. Best case: collaboration or acknowledgment in follow-up work.

---

## The Single Highest-Priority Action

Run the differential expression analysis on OPCs stratified by CupRap vs Control. It is one function call on data already loaded, it uses no new tools, and it would produce a concrete gene list that directly tests the post-commitment claim with a third independent method — complementing the velocity drivers and the paper's Figure 2I. That closes the loop on the central argument and makes the project ready to write up.

```python
# In preprocess.py or a new script, after loading processed.h5ad:
opc_cells = adata[adata.obs['cell_type'] == 'OPC'].copy()
sc.tl.rank_genes_groups(opc_cells, groupby='treatment', groups=['CupRap'], reference='Cntl', method='wilcoxon')
sc.pl.rank_genes_groups(opc_cells, n_genes=25, sharey=False)
```
