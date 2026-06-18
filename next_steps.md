# Next Steps

---

## What Has Already Been Built

The current project has:

- A full scRNA-seq preprocessing pipeline (QC, normalization, HVG selection, Leiden clustering)
- A binary fate classifier (XGBoost + SHAP explainability) → `ol_commitment.csv`, `trigger_genes.csv`
- RNA velocity on spliced/unspliced counts built from raw FASTQs via kb-python → per-condition stream plots, velocity driver gene rankings per cell type and condition
- CellRank fate probability analysis (velocity-based P(OL) per TAP cell) → `tap_fate_comparison.csv`
- Three-method TAP fate comparison (XGBoost, CellRank, bias score) with pairwise Pearson/Spearman/Cohen's kappa statistics → `tap_fate_method_summary.csv`, `tap_fate_per_condition.csv`
- Canonical marker dotplot validating cell type labels across all 11 cell types → `dotplot__marker_dotplot.png`
- Three positive computational findings that are not in the published paper:
  1. **Trigger gene list** — SHAP identified Fa2h, Gjc3, Dock10, Pllp, Tspan2 as the genes most discriminating of OL fate, derived independently from gene expression alone
  2. **Fa2h is the #1 COP velocity driver in both conditions** — the strongest transcriptional momentum signal toward OL maturation localizes to a single gene at a specific lineage stage
  3. **Meis2 is the rank-1 CupRap TAP velocity driver** (corr 79.77, Spearman 0.93) — a neuroblast gene, meaning CupRap TAPs are being pushed toward NB fate, not OL; this is unexpected and not reported by Willis et al.
  - Context: CellRank P(OL) and XGBoost P(OL) are both lower in CupRap TAPs than control TAPs across all conditions, consistent with the Meis2 finding
- A complete scientific narrative (PROJECT_NARRATIVE.md, rna_velocity_trajectory.md, GSE266687_summary.md, tap_fate_comparison.md)

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
| pySCENIC | Transcription factor regulon analysis — which TFs drive the OL program? | `pyscenic` |
| scVI / scANVI | Deep learning for batch correction + latent representation | `scvi-tools` |
| NicheNet / LIANA | Ligand-receptor communication modeling (extends the IGF1/OSM story) | `liana-py` |
| Squidpy | Spatial analysis — the Xenium data for this paper is on GEO | `squidpy` |

**pySCENIC** is probably the highest-impact for discovering something new: it would identify which transcription factors are regulatorily active in OPCs under CupRap but not in TAPs — a mechanistic answer to why the commitment line is where it is, and directly actionable for experimental follow-up.

---

## Where This Goes

### Science Fairs — This Work Is Competition-Ready

Regeneron ISEF and its regional feeders accept computational biology. This project has a clear question, a novel method, and positive original findings. RNA velocity is not common in high school ISEF entries. The strongest argument is the convergence: SHAP (gene expression) and velocity (splicing kinetics) independently identify the same genes — Fa2h, Gjc3, Dock10 — as marking OL commitment. Two methods, different biological signals, same answer. The Meis2 finding (CupRap TAPs pushed toward NB fate) is a specific, unexpected, and testable claim not reported in the paper.

### Preprint on bioRxiv

A preprint can be posted without a PI or institutional affiliation. The bar is: does it say something true and non-obvious? The trigger gene list (Fa2h, Gjc3, Dock10, Pllp, Tspan2), the SHAP-velocity convergence on the same genes, and the Meis2 finding in CupRap TAPs all meet that bar. 
### Summer Research Programs

Programs like RSI (Research Science Institute), PRIMES at MIT, and university HSSP programs look for students who have already done independent research and want mentorship or lab access to extend it. The GitHub repository and this project serve as a portfolio demonstrating real analytical capability.

### Contact the Miller Lab

Freda Miller (freda.miller@msl.ubc.ca) is at UBC. The paper is open access. The velocity analysis performed here found something the paper does not report — specifically the post-commitment specificity of the CupRap effect demonstrated through velocity driver analysis. A short, respectful email summarizing the finding and the methods is not unusual; computational re-analyses of published data are welcomed when they add something. Worst case: no reply. Best case: collaboration or acknowledgment in follow-up work.

---
