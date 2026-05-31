# Marker Gene Panels — PanglaoDB-derived (Mouse)

Source: **PanglaoDB** (Franzén et al. 2019, *Database*), file `PanglaoDB_markers_27_Mar_2020.tsv`.
Filtered to mouse (`species` contains `Mm`), sorted by `canonical marker` → `specificity_mouse` → `sensitivity_mouse`, top 25 per type.
Gene symbols converted to mouse Title Case (e.g., `GFAP` → `Gfap`).

TAP and COP are not separate cell types in PanglaoDB and were added manually (see notes).

| Panel | Marker genes |
|---|---|
| **qNSC** | `Apoe`, `Ccnd2`, `Dcx`, `Id4`, `Aldoc`, `Igf1r`, `Dlx1`, `Neurod1`, `Hes5`, `Dlx2`, `Eomes`, `Prom1`, `Hopx`, `Arx`, `Gfap`, `Dynlt1c`, `Ascl2`, `Pou3f4`, `Prss56`, `Sox1`, `Sp8`, `Sox3`, `Slc1a3`, `Rnd3`, `Sp9` |
| **aNSC** | `Egfr`, `Ascl1`, `Mki67`, `Mcm2`, `Ccnd2`, `Apoe`, `Dcx`, `Id4`, `Aldoc`, `Igf1r`, `Dlx1`, `Neurod1`, `Hes5`, `Dlx2`, `Eomes`, `Prom1`, `Hopx`, `Arx`, `Gfap`, `Dynlt1c`, `Ascl2`, `Pou3f4`, `Prss56`, `Sox1`, `Sp8` |
| **TAP** | `Ascl1`, `Egfr`, `Dlx1`, `Dlx2`, `Mki67`, `Top2a`, `Hmgb2`, `Cenpf` |
| **Neuroblast** | `Ncam1`, `Ncan`, `Pbx1`, `Igfbpl1`, `Ezh2`, `Dcx`, `Epha4`, `Cnr1`, `Itga6`, `Dab1`, `Mark2`, `Nes`, `Neurod1`, `Grm5`, `Dlx2`, `Lrp8`, `Isl1`, `Cux2`, `Pros1`, `Eomes`, `Neurog2`, `Ntng1`, `Egfr`, `Erbb4`, `Draxin` |
| **OPC** | `Nnat`, `Cnp`, `Cspg5`, `Aldoc`, `Fyn`, `Olig1`, `Epn2`, `Gpr37l1`, `Olig2`, `Cdo1`, `Pdgfra`, `Nkx6-2`, `Lad1`, `Etv5`, `Gpr17`, `Cspg4`, `Matn4`, `Ascl1`, `Lhfpl3`, `C1ql1`, `Neu4`, `Kcnip3`, `Pcdh15`, `Tmem100`, `Vcan` |
| **OL** | `Plp1`, `Mbp`, `Enpp2`, `Fasn`, `Olig1`, `Ptgds`, `Hdac11`, `Mcam`, `Olig2`, `Hepacam`, `Gamt`, `Pnpla2`, `Prkcz`, `Efnb3`, `Cldn11`, `Fgfr2`, `Pllp`, `Creb5`, `Bace1`, `Ddc`, `Pou3f1`, `Pde8a`, `Il33`, `Klhl2`, `Anln` |
| **Astrocyte** | `Apoe`, `Fabp7`, `Ets1`, `Gja1`, `Aldoc`, `Gsta4`, `Htra1`, `Cmtm5`, `Lcn2`, `Pla2g7`, `Luzp2`, `Gpr37l1`, `Nkain4`, `Hspa2`, `Cldn10`, `Bysl`, `Acsbg1`, `Entpd2`, `Aqp4`, `Ccr7`, `Cd40`, `Fgfr3`, `Fzd2`, `Agt`, `Hmg20a` |
| **Microglia** | `Fos`, `Egr1`, `Cd53`, `Ctss`, `Fcgr3`, `C1qb`, `Aif1`, `Apbb1ip`, `Csf1r`, `Cyth4`, `Mafb`, `Olfml3`, `Itgb5`, `Pf4`, `Ccl4`, `Cbr2`, `Colec12`, `C5ar1`, `Clec3b`, `Ccl3`, `Cebpa`, `Ptgs1`, `Ltc4s`, `Fcrls`, `Lpcat2` |
| **Ependymal** | `Efnb3`, `Pcp4l1`, `Aqp1`, `Calml4`, `Myb`, `Aqp4`, `Celsr2`, `Hspa2`, `Hdc`, `Krt15`, `Angptl2`, `Gfap`, `Enkur`, `Pltp`, `Pvalb`, `Ak8`, `Foxj1`, `Ak7`, `Cfap44`, `Odf3b`, `Ccdc153`, `Crocc`, `Trim71`, `Six3`, `Usp18` |
| **Endothelial** | `Cd9`, `Gpm6a`, `Id3`, `Ifi27`, `Clic4`, `Plac8`, `Igfbp7`, `Gng11`, `Cd82`, `Cdkn1c`, `Id1`, `Ly6a`, `Ets1`, `Pdpn`, `Mgp`, `Ednrb`, `Emp1`, `Ptgds`, `Plec`, `Mgll`, `Mylk`, `Ltc4s`, `Cd34`, `Epas1`, `Lifr` |
| **Pericyte** | `Mfge8`, `Col1a1`, `Pecam1`, `Ifitm1`, `Dlk1`, `Gnb4`, `Myo1b`, `Acta2`, `Mcam`, `Pth1r`, `Des`, `Pdgfrb`, `Anpep`, `Cd248`, `Msx1`, `Angpt1`, `Ctgf`, `Higd1b`, `Alpl`, `Cspg4`, `Cog7`, `Notch3`, `Abcc9`, `Kcnj8`, `Inpp4b` |
| **VSMC** | `Bgn`, `Myl9`, `Pcp4l1`, `Itga1`, `Nrp2`, `Mylk`, `Ehd2`, `Fabp4`, `Acta2`, `Ogn`, `Fhl2`, `Gjc1`, `Lmo7`, `Fbln5`, `Crabp1`, `Mfap5`, `Notch3`, `Des`, `Lamb2`, `Hexim1`, `Gja4`, `Lox`, `Aaed1`, `Pde4dip`, `Akt2` |

## Notes
- **qNSC vs aNSC**: PanglaoDB lumps both under `Neural stem/precursor cells`. aNSC list above adds canonical activation markers (`Egfr`, `Ascl1`, `Mki67`, `Mcm2`, `Ccnd2`) on top of the PanglaoDB pull.
- **TAP**: Not in PanglaoDB. Added manually from Cebrián-Silla et al. 2021 *eLife*.
- Full per-gene table (with specificity/sensitivity scores) saved to `panglao_panel.csv`.

---

## Criteria used to build the panel

### 1. Filter — which rows to keep

```python
df_mouse = df[df["species"].str.contains("Mm", na=False)]
```

Keeps any row where `species` contains `Mm` (so both pure-mouse rows and `Mm Hs` shared rows are included). Rows that are human-only (`Hs`) are dropped.

### 2. Cell-type matching — exact string match on `cell type` column

| Panel label | PanglaoDB `cell type` value matched |
|---|---|
| qNSC | `Neural stem/precursor cells` |
| aNSC | `Neural stem/precursor cells` (same bucket — PanglaoDB doesn't separate q/a) |
| Neuroblast | `Neuroblasts` |
| OPC | `Oligodendrocyte progenitor cells` |
| OL | `Oligodendrocytes` |
| Astrocyte | `Astrocytes` |
| Microglia | `Microglia` |
| Ependymal | `Ependymal cells` |
| Endothelial | `Endothelial cells` |
| Pericyte | `Pericytes` |
| VSMC | `Smooth muscle cells` |
| TAP | (not in PanglaoDB — manual) |
| COP | (not in PanglaoDB — manual) |

### 3. Ranking — sort order within each cell type

```python
sub.sort_values(["canonical marker", "specificity_mouse", "sensitivity_mouse"], ascending=False)
```

Three-tier sort, in this priority:
1. `canonical marker` (1 vs 0) — textbook markers float to the top
2. `specificity_mouse` — how unique the gene is to that cell type (high = more specific)
3. `sensitivity_mouse` — how frequently the gene is detected in that cell type (high = more reliable)

`NaN` values in `canonical marker` are treated as `0`. Missing specificity/sensitivity values sort to the bottom.

### 4. Cutoff — top N

`head(25)` per cell type. No specificity threshold applied — the ranking does the filtering. For a stricter list, add `sub[sub["specificity_mouse"] > 0.5]` before the head().

### 5. Symbol formatting

```python
def to_mouse_symbol(s):
    return s[:1].upper() + s[1:].lower()
```

PanglaoDB stores symbols in human-style ALL CAPS (`GFAP`, `PDGFRA`). Converted to mouse Title Case (`Gfap`, `Pdgfra`).

### 6. Manual additions

| Panel | Added | Why | Source |
|---|---|---|---|
| **TAP** | `Ascl1, Egfr, Dlx1, Dlx2, Mki67, Top2a, Hmgb2, Cenpf` | Not a cell type in PanglaoDB at all | Cebrián-Silla 2021 *eLife* |
| **COP** | `Gpr17, Tcf7l2, Bmp4, Sox6, Neu4` (then top-10 OPC genes appended) | Not separated from OPC in PanglaoDB | Marques 2016 *Science* |
| **aNSC** | `Egfr, Ascl1, Mki67, Mcm2, Ccnd2` prepended | qNSC and aNSC share PanglaoDB's NSC bucket — these activation markers force differentiation | Codega 2014; Llorens-Bobadilla 2015 |

### One-sentence version for your methods

> "Mouse markers were filtered from PanglaoDB (Franzén et al. 2019), ranked by `canonical marker` then `specificity_mouse` then `sensitivity_mouse`, and the top 25 genes per cell type were taken. Cell types not present in PanglaoDB (transit-amplifying progenitors and committed oligodendrocyte progenitors) were supplemented with markers from Cebrián-Silla et al. 2021 and Marques et al. 2016 respectively."