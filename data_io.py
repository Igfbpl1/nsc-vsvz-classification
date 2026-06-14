from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import anndata as ad
import pandas as pd
import scipy.io
import scipy.sparse


@dataclass
class Sample:
    gsm: str
    label_in_filename: str
    strain: str
    treatment: str
    timepoint: str
    replicate: int


SAMPLES = [
    Sample("GSM8253792", "CD1_Cntl_0wksRecov", "CD1", "Cntl", "0wks", 1),
    Sample("GSM8253793", "CD1_Cntl_3wksRecov", "CD1", "Cntl", "3wks", 1),
    Sample("GSM8253794", "CD1_CR_Rep1", "CD1", "CupRap", "3wks", 1),
    Sample("GSM8253795", "CD1_CR_Rep2", "CD1", "CupRap", "3wks", 2),
    Sample("GSM8253796", "NesCre_Cntl_Rep1", "NesCre", "Cntl", "3wks", 1),
    Sample("GSM8253797", "NesCre_Cntl_Rep2", "NesCre", "Cntl", "3wks", 2),
    Sample("GSM8253798", "NesCre_CR_Rep1", "NesCre", "CupRap", "3wks", 1),
    Sample("GSM8253799", "NesCre_CR_Rep2", "NesCre", "CupRap", "3wks", 2),
    Sample("GSM8647352", "CD1_CR1_NoRecov1", "CD1", "CupRap", "0wks", 1),
    Sample("GSM8647353", "CD1_CR1_NoRecov2", "CD1", "CupRap", "0wks", 2),
]


def load_sample(raw_dir: Path, gsm: str, label: str) -> ad.AnnData:
    prefix = raw_dir / f"{gsm}_{label}"
    barcodes = pd.read_csv(f"{prefix}_barcodes.tsv", header=None, sep="\t")[0].values
    genes = pd.read_csv(
        f"{prefix}_genes.tsv",
        header=None,
        sep="\t",
        names=["gene_id", "gene_symbol", "feature_type"],
    )
    mtx = scipy.io.mmread(f"{prefix}_matrix.mtx").tocsr()
    if mtx.shape[0] != len(genes):
        mtx = mtx.T.tocsr()
    X = mtx.T.tocsr()  # cells × genes

    var = genes.copy()
    var.index = pd.Index(var["gene_symbol"].astype(str), name=None)
    var = var.drop(columns=["gene_symbol"])
    obs = pd.DataFrame(index=pd.Index(barcodes.astype(str), name="barcode"))
    return ad.AnnData(X=X, obs=obs, var=var)


def load_all(raw_dir: Path) -> ad.AnnData:
    parts = []
    for sample in SAMPLES:
        a = load_sample(raw_dir, sample.gsm, sample.label_in_filename)
        a.obs["sample_id"] = sample.gsm
        a.obs["sample_label"] = sample.label_in_filename
        a.obs["strain"] = sample.strain
        a.obs["treatment"] = sample.treatment
        a.obs["timepoint"] = sample.timepoint
        a.obs["replicate"] = sample.replicate
        a.obs_names = [f"{sample.gsm}_{b}" for b in a.obs_names]
        a.var_names_make_unique()
        parts.append(a)
        print(
            f"  {sample.gsm} {sample.label_in_filename}: {a.n_obs} cells × {a.n_vars} genes"
        )

    combined = ad.concat(parts, axis=0, join="inner", merge="same", index_unique=None)
    combined.var_names_make_unique()
    return combined
