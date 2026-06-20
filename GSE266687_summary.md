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

## Experimental Design — Why These 10 Samples

The 10 samples vary along **three axes**: strain (CD1 vs NesCre) × treatment (Cntl vs CupRap) × timepoint (0wks vs 3wks). The design is deliberately unbalanced — each axis answers a different question.

**CD1 block (wild-type) — a clean 2×2:**

| | 0wks | 3wks |
|---|---|---|
| **Cntl** | GSM8253792 | GSM8253793 |
| **CupRap** | GSM8647352, GSM8647353 | GSM8253794, GSM8253795 |

**NesCre block (NesCreERT2-Eyfp lineage tracing) — 3wks only:**

| | 3wks |
|---|---|
| **Cntl** | GSM8253796, GSM8253797 |
| **CupRap** | GSM8253798, GSM8253799 |

**What each axis proves:**

- **Treatment (Cntl vs CupRap)** — the core perturbation: does demyelination push the V-SVZ toward making oligodendrocytes, and where in the lineage? CupRap is replicated (×2) because it is the condition of interest. *Proves:* OPC proportion roughly doubles, but TAP identity is unchanged (r = 0.995) → the response is downstream of TAPs.
- **Timepoint (0wks vs 3wks, CD1 only)** — the kinetics: acute injury (0wks, harvested at the end of the 6-week treatment = peak generation) vs recovery window (3wks). Studied in wild-type because counting cell-type proportions needs no genetic label. *Proves:* OL-precursor generation peaks acutely and recedes by 3wks.
- **Strain (CD1 vs NesCre)** — characterization vs causal attribution. CD1 gives clean reference biology (atlas, composition, velocity); NesCre adds lineage tracing to prove new OLs are NSC-derived. Only needed at 3wks (the remyelination window). *Proves:* EYFP+ labeling jumps 3–6× at the OPC/OL stage → new oligodendroglia genuinely come from V-SVZ stem cells.

**Important — the growth factors (IGF1, OSM) are NOT in these sequenced samples.** Every GSM sample is one cell of the {strain × treatment × timepoint} grid; none received recombinant IGF1 or OSM. Those ligands were (1) *identified computationally* from this atlas via cell-cell communication analysis, then (2) *validated in separate functional experiments* (cultured NPCs + lateral-ventricle infusion) that are not part of the sequencing dataset. The 10 samples are the characterization + lineage-attribution engine; the growth-factor mechanism is a downstream chapter built on what they localized.

---

## CupRap Model

Cuprizone is fed for 6 weeks to demyelinate the corpus callosum. Rapamycin is co-administered to inhibit mTOR and block natural remyelination. Mice are then given 0 or 3 weeks of recovery. This creates a controlled window to observe the OL-lineage response without confounding recovery.

- **Cuprizone**: kills mature oligodendrocytes → demyelination signal
- **Rapamycin**: inhibits mTOR → blocks OPC→mature OL maturation
- Together: traps cells at the OPC stage (committed but unable to fully mature)

**Why "remyelination blocked" but OPCs still form — the lineage is a pipeline, and rapamycin only blocks the end:**

```
NSC → TAP → OPC → COP → immature OL → mature myelinating OL → wraps axon (REMYELINATION)
            └──────── generation ────────┘   └──── maturation + myelination ────┘
                   (still happens —                (this is what rapamycin blocks;
                    driven harder by injury)        mTOR is required for myelin synthesis)
```

mTOR is required for the terminal maturation/myelination step, not for generating or proliferating OPCs upstream. So under CupRap, cuprizone (upstream) drives the niche to *generate* OPCs while rapamycin (terminal) prevents them maturing — they accumulate as a **traffic jam**. "OL-lineage cells forming" (generation) and "remyelination" (functional myelin wrapping) are different events; only the latter is blocked.

> **Reading the proportions correctly:** cell-type fractions in the scRNA-seq are *compositional* (they sum to 100%), so one population's share moving does not mean its absolute count changed — when microglia crash and neuroblasts surge at 3wks, every other fraction is mechanically diluted. Also, the mature OLs that cuprizone kills live in the **corpus callosum (white matter)**, a different location from the **V-SVZ niche** being sequenced; the OL-lineage cells captured here are mostly newly-generated OPCs/COPs, not the dying white-matter OLs.

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

**What EYFP-positive means (NesCreERT2-Eyfp system).** EYFP is a genetic lineage label, not a natural gene — present only in the NesCre strain (CD1 is wild-type, unlabeled). The strain name decodes the mechanism:

- **Nes** (Nestin promoter) — active in neural stem cells → controls *where* the label can switch on.
- **CreERT2** — a tamoxifen-inducible Cre recombinase → controls *when* (only fires when tamoxifen is given).
- **Eyfp** — a reporter behind a "stop" cassette; once Cre fires it is permanently deleted and EYFP turns on.

```
give tamoxifen → Cre activates only in Nestin+ NSCs → permanently deletes the stop
              → EYFP on, and inherited by ALL descendants (it is a DNA edit)
```

So **EYFP+ = a cell that was an NSC at tamoxifen time, or any of its progeny** — a permanent stamp on the entire V-SVZ NSC lineage. It decouples "was generated from a stem cell" (the label) from "completed myelination" (the function), which is exactly why it can show new OL-lineage cells are being *generated* even while rapamycin blocks their maturation. Without it you could not tell whether the OPC expansion was newly-made cells or rearrangement of pre-existing ones.

### Evidence 4 — No-recovery timepoint: OPCs expand, mature OLs are depleted

> "mature oligodendrocytes were reduced approximately 5-fold in the no recovery corpus callosum... OPCs were increased in the dorsal but not lateral V-SVZ in the no recovery condition."
> — p. 8–9, Figure 4C–E

This is the direct mechanistic fingerprint of rapamycin: mTOR inhibition blocks OPC→mature OL maturation. Cells accumulate as PDGFRα+/OLIG2+ OPCs — already committed to the OL lineage — but cannot complete differentiation.

### Evidence 5 — Velocity analysis: no positive OL driver exists in TAPs

From RNA velocity analysis on this dataset (rna_velocity_trajectory.md):

> "No positive OL-leaning signal exists at the TAP stage. All 16 SHAP-confirmed positive OL markers are absent from the top 100 TAP drivers in both conditions."

The top CupRap TAP velocity drivers are **Srrm4** (rank 1, corr 73.07) and **Meis2** (rank 2, corr 71.99, Spearman 0.93) — both neuroblast/neuronal genes, not OL genes. Meis2 climbs from rank 10 in Cntl to rank 2 in CupRap; the model uses its *absence* as an OL signal. CupRap TAPs are being driven toward the NB arm transcriptionally, not the OL arm. This further confirms there is no OL-fate commitment signal detectable at the TAP stage.

### Evidence 6 — Positive OL velocity drivers only appear after commitment

The SHAP-confirmed positive OL-lineage genes only show up as velocity drivers downstream of commitment:

| Gene | TAP rank | OPC rank | COP rank | OL rank | SHAP direction |
|---|---|---|---|---|---|
| Fa2h | absent top 100 | 37 | **1** | — | POSITIVE_OL |
| Gjc3 | absent top 100 | 22 | 11 | **1** | POSITIVE_OL |

The OL transcriptional program is dynamically active at the OPC→COP transition, not in TAPs.

---

## Summary Table

| Evidence | Stage affected | Source |
|---|---|---|
| TAP transcriptome r = 0.995 Cntl vs CupRap | Not TAPs | Willis et al. Fig 2I |
| OPC proportion doubles (8%→19%, p = 0.002) | OPCs | Willis et al. Fig 3 |
| Eyfp lineage label 3–6× increase at OPC/OL | Post-commitment | Willis et al. Fig 2E |
| Mature OLs depleted 5× at no-recovery timepoint | OPC→OL maturation blocked | Willis et al. Fig 4C–E |
| Srrm4/Meis2 top CupRap TAP velocity drivers (NB genes) | No OL signal in TAPs | Velocity analysis |
| Fa2h (COP rank 1), Gjc3 (OL rank 1) drive only post-TAP | Post-commitment | Velocity analysis |

**Conclusion:** CupRap expands the committed OPC pool and blocks their maturation into oligodendrocytes via mTOR inhibition. TAP transcriptional identity is unchanged. The effect is real but acts downstream of commitment.
