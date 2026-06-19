"""Cell-type marker panels for V-SVZ scRNA-seq annotation.

Markers are taken verbatim from the STAR Methods of:
  Willis et al. (2025) "Single cell approaches define neural stem cell niches
  and identify microglial ligands that can enhance precursor-mediated
  oligodendrogenesis." Cell Reports 44, 115194.

The relevant passage (p. 23):
  "We used gene expression overlays on t-SNE plots (FeaturePlots, Seurat) to
  annotate cell types and state based on the following well-defined marker
  genes, including: For all NSCs, Rcn3, Nes, Veph1, Notum, Tspan18; for
  activated NSCs, Egfr, Ascl1; for dormant NSCs, Meg3, Sparc, Fbxo2, Id3;
  for TAPs, Mki67, Rrm2, Hells, Gsx2; for astrocytes, Aqp4, Agt, Hsd11b1,
  Lcat; for all oligolineage cells: Olig1, Sox10; for immature
  oligodendrocytes, Enpp4, Gpr17, Fyn, Nkx2-2; for OPCs, Pdgfra, Cspg4;
  for mature oligodendrocytes, Mog, Mag, Opalin; for all neuroblasts,
  Dlx1, Dlx2, Sp8, Sp9, Gad1, Gad2; for all microglia, Aif1, Trem2, Cx3cr1,
  Tmem119; for other immune cells including T cells and neutrophils,
  Cd52, Cd69; for striatal neurons, Calb1, Gad1, Gad2, Bcl11b; for
  endothelial cells, Pecam1, Esam, Plvap; for pericytes, Carmn, Cspg4,
  Ano1; for VAMCs, Pdgfrb, Myh11, Mylk; and for ependymal cells, Foxj1,
  Pvalb."

Differences from the prior local panel:
  * NSC split into aNSC and dNSC (paper's whole NSC story depends on this).
  * Mural split into Pericyte and VAMC.
  * Added Striatal_Neuron and Other_Immune categories that the paper uses.
  * TAP panel now matches the paper (was previously a mix of aNSC + Neuroblast
    markers that would have caused mis-assignment).
  * Microglia panel uses the paper's 4 markers; `Iba1` (a protein name, not
    an MGI gene symbol) was removed in favor of its gene `Aif1`.
  * Astrocyte panel uses the paper's V-SVZ-specific markers (Hsd11b1, Lcat
    were missing); pan-astrocyte markers that also light up NSCs (Gfap,
    Sox9, Slc1a3, Aldoc) were removed.
  * Cell-type key names kept singular (Neuroblast, OPC, etc.) to preserve
    compatibility with downstream pipeline scripts.
"""

MARKERS = {
    # ---- NSC lineage ----------------------------------------------------
    "aNSC": ["Egfr", "Ascl1"],
    "dNSC": ["Meg3", "Sparc", "Fbxo2", "Id3"],
    "TAP":  ["Mki67", "Rrm2", "Hells", "Gsx2"],
    "Neuroblast": ["Dlx1", "Dlx2", "Sp8", "Sp9", "Gad1", "Gad2", "Dcx"],

    # ---- Oligo lineage --------------------------------------------------
    # Paper markers in *first* position; subtype-specific extensions
    # follow. Excluded any pan-oligo TFs (Sox10, Olig1, Olig2) that would
    # cause idxmax ties between OPC / COP / OL — those live in
    # LINEAGE_GENES below for ML-feature exclusion.
    "OPC": ["Pdgfra", "Cspg4"],
    "COP": ["Enpp4", "Gpr17", "Fyn", "Nkx2-2",
            "Bcas1", "Tcf7l2", "Bmp4"],              # immature OLGs
    "OL":  ["Mog", "Mag", "Opalin",
            "Mbp", "Plp1", "Cldn11", "Ermn"],        # mature OLGs

    # ---- Glia + ependymal ----------------------------------------------
    "Astrocyte": ["Aqp4", "Agt", "Hsd11b1", "Lcat"],
    "Ependymal": ["Foxj1", "Pvalb"],

    # ---- Immune ---------------------------------------------------------
    "Microglia":   ["Aif1", "Trem2", "Cx3cr1", "Tmem119"],
    "Other_Immune": ["Cd52", "Cd69"],

    # ---- Vasculature ----------------------------------------------------
    "Endothelial": ["Pecam1", "Esam", "Plvap"],
    "Pericyte":    ["Carmn", "Cspg4", "Ano1"],
    "VAMC":        ["Pdgfrb", "Myh11", "Mylk"],

    # ---- Neurons --------------------------------------------------------
    "Striatal_Neuron": ["Calb1", "Gad1", "Gad2", "Bcl11b"],
}

# Lineage-level genes that identify *which lineage* a cell is committed to
# (oligodendrocyte vs neuronal) but not the specific subtype. These are
# DELIBERATELY kept out of MARKERS because:
#   * Adding pan-oligo TFs (Sox10, Olig1/2) to OPC / COP / OL panels would
#     break idxmax annotation — all three subtypes would score equally.
#   * Tubb3 is expressed in mature neurons too, so adding it to Neuroblast
#     would confuse it with Striatal_Neuron.
#
# Used by train_ol_classifier.py to *exclude* these genes from the feature
# matrix when training the OL-vs-Neuroblast classifier, so SHAP importance
# surfaces non-canonical trigger candidates rather than these well-known
# lineage genes.
LINEAGE_GENES = {
    "OL_lineage":        ["Olig1", "Olig2", "Sox10", "Cnp", "Lhfpl3", "Mobp"],
    "Neuroblast_lineage": ["Tubb3", "Cd24a"],
}

OL_LINEAGE = {"OPC", "COP", "OL"}
