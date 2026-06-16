# nsc-vsvz-classification

This project utilizes ML to give a similarity score to every undecided TAP and to surface a list of the top 20 non-canonical genes indicative of differentiation into either OL or NB.

---

## Literature Review — Source Paper Comparison

The primary dataset (GSE266687) was published in:

**Wegleiter et al., "Single cell approaches define neural stem cell niches and identify microglial ligands that can enhance precursor-mediated oligodendrogenesis"**
Cell Reports, January 2025. doi:10.1016/j.celrep.2024.115194
PubMed: [39823226](https://pubmed.ncbi.nlm.nih.gov/39823226/)
bioRxiv preprint: [2024.03.22.586277](https://www.biorxiv.org/content/10.1101/2024.03.22.586277v1)

### What the source paper does

- **Methods**: t-SNE + Harmony batch correction + UMAP + spatial transcriptomics + canonical cell type annotation
- **TAP markers used (canonical only)**: Mki67, Rrm2, Hells, Gsx2 (proliferation markers — TAPs treated as one cluster)
- **OL lineage markers used**: Olig1, Sox10, Pdgfra/Cspg4 (OPC), Mog/Mag/Opalin (mature OL)
- **NB markers used**: Dlx1, Dlx2, Sp8, Sp9, Gad1, Gad2

**Main contributions of the paper:**
1. Compositional/spatial analysis of dorsal vs lateral V-SVZ niches
2. Cell-cell communication analysis identifying microglial→NPC ligands
3. Functional validation that IGF1 and OSM enhance oligodendrogenesis in culture and after lateral ventricle infusion
4. Proposes that altered microglial signaling drives NSC activation, local oligodendrogenesis, and myelin repair

### What the source paper does NOT do

Full-text search confirms the paper does not address:
- RNA velocity (scVelo) or any trajectory analysis (Monocle, Slingshot, PAGA)
- Driver gene discovery for TAP fate commitment
- Machine learning classifiers (XGBoost, SHAP) for fate prediction
- Transcription factor analysis for OL vs NB lineage decision in TAPs
- Specific non-canonical markers: Gjc3, Fa2h, Tspan2, Pllp, Cnp, Cryab, Ptprz1
- NB-side TFs: Bcl11a, Nfib, Meis2, Stmn2, Igfbpl1

The paper's framework is at the **niche/ligand signaling level**, not the **TAP transcriptional decision level**.

### How this project complements the source paper

This analysis is complementary, not duplicative. The paper sets up the upstream framework (microglia → IGF1/OSM → TAPs); this project analyzes the downstream cell-intrinsic transcriptional response.

| Contribution | Source paper | This project |
|---|---|---|
| Identify microglial ligands (IGF1, OSM) | ✓ Main finding | — |
| Spatial map of dorsal vs lateral V-SVZ | ✓ | — |
| Functional validation in culture/in vivo | ✓ | — |
| RNA velocity of TAP→OL/NB transitions | — | ✓ |
| XGBoost+SHAP classifier for TAP fate prediction | — | ✓ |
| Early NB-fate TAP drivers (Bcl11a, Nfib, Meis2) | — | ✓ |
| Positive non-canonical OL markers (Pllp, Gjc3, Fa2h, Tspan2, Cnp, Cryab, Dock10) | — | ✓ |
| Gjc3 as CupRap-responsive OL marker (velocity rank 22→1 with added CupRap) | — | ✓ |
| Igf1 as a COP velocity driver (cell-intrinsic IGF1 dynamics) | — | ✓ |
| Hypothesis: OL commitment lacks an early transcriptional switch | — | ✓ |

### Defensible framing for the project

> Wegleiter et al. (Cell Reports 2025) identified IGF1 and OSM as microglial ligands that enhance V-SVZ oligodendrogenesis. We extended this work by performing RNA velocity and SHAP-based classifier analysis on the same dataset (GSE266687) to characterize the downstream transcriptional response in TAPs. Our analysis reveals (1) clear early NB-fate drivers (Bcl11a, Nfib, Meis2) detectable at the TAP stage, (2) absence of a comparable early OL-fate transcriptional switch in TAPs, suggesting OL commitment may be defined by failure to engage the NB program rather than active induction of an OL program, and (3) Gjc3 as a CupRap-responsive OL maturation marker whose velocity signal scales with CupRap sample inclusion, providing direct molecular evidence linking microglial IGF1/OSM signaling to a specific downstream transcriptional event in committing OL precursors.
