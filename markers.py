"""V-SVZ + OL-lineage marker gene panels from the literature."""

MARKERS = {
    "qNSC":        ["Gfap", "Slc1a3", "Aldh1l1", "Vim", "Sox2", "Sox9", "Hes5", "Aldoc", "Apoe", "Id3"],
    "aNSC":        ["Egfr", "Ascl1", "Mki67", "Top2a", "Sox2", "Hes5", "Nes"],
    "TAP":         ["Egfr", "Mki67", "Ccnd2", "Dlx1", "Dlx2", "Ascl1", "Top2a"],
    "Neuroblast":  ["Dcx", "Sox4", "Sox11", "Stmn2", "Stmn3", "Nrxn3", "Tubb3", "Dlx1", "Dlx2"],
    "OPC":         ["Pdgfra", "Cspg4", "Olig1", "Olig2", "Sox10", "Cnp", "Lhfpl3"],
    "COP":         ["Gpr17", "Sox10", "Olig2", "Itpr2", "Bmp4", "Bcas1", "Tcf7l2"],
    "OL":          ["Mbp", "Mog", "Mag", "Plp1", "Mobp", "Cldn11", "Mal", "Trf", "Cnp"],
    "Astrocyte":   ["Aqp4", "Slc1a2", "Aldh1l1", "S100b", "Mfge8", "Ntsr2", "Slc1a3"],
    "Microglia":   ["Cx3cr1", "Tmem119", "P2ry12", "C1qa", "Csf1r", "Trem2", "Hexb", "Aif1"],
    "Ependymal":   ["Foxj1", "Ccdc153", "Tmem212", "Rsph1", "Sox2"],
    "Endothelial": ["Cldn5", "Pecam1", "Flt1", "Cdh5"],
    "Mural":       ["Pdgfrb", "Acta2", "Rgs5"],
}

OL_LINEAGE = {"OPC", "COP", "OL"}
