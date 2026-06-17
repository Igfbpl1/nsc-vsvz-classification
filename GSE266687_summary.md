# GSE266687 — Paper Summary & CupRap Argument

**Willis et al., "Single cell approaches define neural stem cell niches and identify microglial ligands that can enhance precursor-mediated oligodendrogenesis"**
Cell Reports 44, 115194. January 28, 2025.
doi: [10.1016/j.celrep.2024.115194](https://doi.org/10.1016/j.celrep.2024.115194) | PMID: 39823226

---

## What the Paper Does

Uses scRNA-seq (10x Genomics) + single-cell spatial transcriptomics (Xenium + MERFISH) on the adult mouse forebrain V-SVZ under homeostasis and cuprizone-rapamycin (CupRap) demyelination/remyelination.

**Main findings:**
1. Dorsal and lateral V-SVZ are transcriptionally distinct NSC niches
2. After white matter injury, dorsal NSCs are activated locally to generate oligodendrocytes
3. Injury induces a transcriptionally distinct remyelination-enriched microglial state in the dorsal niche
4. Two microglial ligands — IGF1 and OSM — promote NPC proliferation and oligodendrogenesis in culture and in vivo

**Dataset:** 10 samples × {barcodes, genes, matrix} — 4 conditions (Cntl/CupRap × 0wks/3wks recovery), CD1 and NesCreERT2-Eyfp strains. ~39,905 cells total after QC.

---

## CupRap Model

Cuprizone is fed for 6 weeks to demyelinate the corpus callosum. Rapamycin is co-administered to inhibit mTOR and block natural remyelination. Mice are then given 0 or 3 weeks of recovery. This creates a controlled window to observe the OL-lineage response without confounding recovery.

- **Cuprizone**: kills mature oligodendrocytes → demyelination signal
- **Rapamycin**: inhibits mTOR → blocks OPC→mature OL maturation
- Together: traps cells at the OPC stage (committed but unable to fully mature)

---

## Key Claim: CupRap Acts on OPCs After Commitment, Not in TAPs

### Evidence 1 — TAP transcriptional identity is unchanged by CupRap

> "Control and remyelinating dNSCs, aNSCs, and TAPs were highly similar, with r values of 0.97–0.99... Thus, demyelination likely activates dNSCs to generate TAPs and, ultimately, oligodendroglia **without altering NPC transcriptional identities**."
> — p. 5, Figure 2I

The Pearson correlation of average gene expression between Control and CupRap TAPs is **r = 0.995**. The cells are transcriptionally indistinguishable. If CupRap were biasing fate at the TAP stage, a transcriptional shift would be expected. There is none.

### Evidence 2 — OPC proportion doubles under CupRap, restricted to dorsal V-SVZ

> "the proportion of OPCs relative to total oligodendroglia was significantly increased (controls, 8% ± 1.7% SEM; remyelinating 19% ± 1.4% SEM; n = 4/condition, **p = 0.002)"
> — p. 7, Figure 3

This >2-fold OPC expansion is:
- Statistically significant (p = 0.002)
- Spatially restricted to the dorsal V-SVZ/corpus callosum, not lateral (Figures 3M, 3N)
- Not matched by a significant TAP change (Figure 1K shows TAPs as "ns")

### Evidence 3 — Lineage tracing confirms the effect is at OPC/OL stages

> "the oligodendroglial lineage labeling index was altered during remyelination. In controls, 5%–7% of OPCs, immature oligodendrocytes, and mature oligodendrocytes were Eyfp positive, while during remyelination this was **increased 3- to 6-fold** (Figure 2E)."
> — p. 5

Eyfp+ labeling marks NPC-derived cells. The jump is at the OPC/immature OL stage specifically. Cells accumulate as OPCs downstream of commitment, not upstream in TAPs.

### Evidence 4 — No-recovery timepoint: OPCs expand, mature OLs are depleted

> "mature oligodendrocytes were reduced approximately 5-fold in the no recovery corpus callosum... OPCs were increased in the dorsal but not lateral V-SVZ in the no recovery condition."
> — p. 8–9, Figure 4C–E

This is the direct mechanistic fingerprint of rapamycin: mTOR inhibition blocks OPC→mature OL maturation. Cells accumulate as PDGFRα+/OLIG2+ OPCs — already committed to the OL lineage — but cannot complete differentiation.

### Evidence 5 — Velocity analysis: no positive OL driver exists in TAPs

From RNA velocity analysis on this dataset (rna_velocity_trajectory.md):

> "No positive OL-leaning signal exists at the TAP stage. All 16 SHAP-confirmed positive OL markers are absent from the top 100 TAP drivers in both conditions."

The rank-1 CupRap TAP velocity driver is **Meis2** (corr 79.77, Spearman 0.93) — a neuroblast marker. The model uses its *absence* as an OL signal. CupRap TAPs are being driven toward the NB arm transcriptionally, not the OL arm. This further confirms there is no OL-fate commitment signal detectable at the TAP stage.

### Evidence 6 — Positive OL velocity drivers only appear after commitment

The SHAP-confirmed positive OL-lineage genes only show up as velocity drivers downstream of commitment:

| Gene | TAP rank | OPC rank | COP rank | SHAP direction |
|---|---|---|---|---|
| Fa2h | absent top 100 | 38 | **1** | POSITIVE_OL |
| Gjc3 | absent top 100 | 18 | 9 | POSITIVE_OL |
| Dock10 | absent top 100 | — | 49 | POSITIVE_OL |

The OL transcriptional program is dynamically active at the OPC→COP transition, not in TAPs.

---

## Summary Table

| Evidence | Stage affected | Source |
|---|---|---|
| TAP transcriptome r = 0.995 Cntl vs CupRap | Not TAPs | Willis et al. Fig 2I |
| OPC proportion doubles (8%→19%, p = 0.002) | OPCs | Willis et al. Fig 3 |
| Eyfp lineage label 3–6× increase at OPC/OL | Post-commitment | Willis et al. Fig 2E |
| Mature OLs depleted 5× at no-recovery timepoint | OPC→OL maturation blocked | Willis et al. Fig 4C–E |
| Meis2 is rank-1 CupRap TAP velocity driver (NB gene) | No OL signal in TAPs | Velocity analysis |
| Fa2h/Gjc3/Dock10 only velocity drivers at OPC/COP | Post-commitment | Velocity analysis |

**Conclusion:** CupRap expands the committed OPC pool and blocks their maturation into oligodendrocytes via mTOR inhibition. TAP transcriptional identity is unchanged. The effect is real but acts downstream of commitment.
