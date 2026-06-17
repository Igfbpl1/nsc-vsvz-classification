import argparse
import sys


def validate_matrix(kb_output_dir: str) -> bool:
    from scipy.io import mmread
    import pandas as pd

    cu = f"{kb_output_dir}/counts_unfiltered"

    s = mmread(f"{cu}/spliced.mtx").tocsr()
    u = mmread(f"{cu}/unspliced.mtx").tocsr()
    genes = pd.read_csv(f"{cu}/spliced.genes.txt", header=None)[0]
    u_genes = pd.read_csv(f"{cu}/unspliced.genes.txt", header=None)[0]

    total = s.sum() + u.sum()
    unspliced_pct = u.sum() / total * 100 if total else 0

    print(f"KB output    : {kb_output_dir}")
    print()
    print(f"{'':20} {'Spliced':>12} {'Unspliced':>12}")
    print("-" * 46)
    print(f"{'Barcodes':<20} {s.shape[0]:>12,} {u.shape[0]:>12,}")
    print(f"{'Genes':<20} {s.shape[1]:>12,} {u.shape[1]:>12,}")
    print(f"{'Non-zero entries':<20} {s.nnz:>12,} {u.nnz:>12,}")
    print(f"{'Total counts':<20} {int(s.sum()):>12,} {int(u.sum()):>12,}")
    print()
    print(f"Unspliced fraction : {unspliced_pct:.1f}%  (expected 30–60% for single-nucleus brain data)")
    print(f"Sample gene names  : {list(genes[:3])}")
    print()

    passes = True
    if s.shape[1] != u_genes.shape[0]:
        print("FAIL — spliced and unspliced gene counts differ.")
        passes = False
    if s.nnz == 0 or u.nnz == 0:
        print("FAIL — matrix has zero non-zero entries.")
        passes = False
    if not (30 <= unspliced_pct <= 65):
        print(f"WARN — unspliced fraction {unspliced_pct:.1f}% outside expected 30–60% range.")
    if passes:
        print("PASS — matrices look valid.")
    return passes


def validate(raw_barcodes_path: str, kb_output_dir: str) -> bool:
    raw = set(line.strip().replace("-1", "") for line in open(raw_barcodes_path) if line.strip())

    results = {}
    for layer in ("spliced", "unspliced"):
        kb_path = f"{kb_output_dir}/counts_unfiltered/{layer}.barcodes.txt"
        kb = set(line.strip() for line in open(kb_path) if line.strip())
        overlap = raw & kb
        pct = len(overlap) / len(raw) * 100 if raw else 0
        results[layer] = (len(raw), len(kb), len(overlap), pct)

    print(f"Ground truth : {raw_barcodes_path}  ({results['spliced'][0]:,} cells)")
    print(f"KB output    : {kb_output_dir}")
    print()
    print(
        f"{'Layer':<12} {'Raw cells':>10} {'KB barcodes':>12} {'Overlap':>10} {'Match %':>8}"
    )
    print("-" * 56)
    for layer, (n_raw, n_kb, n_overlap, pct) in results.items():
        print(f"{layer:<12} {n_raw:>10,} {n_kb:>12,} {n_overlap:>10,} {pct:>7.1f}%")

    print()
    spliced_pct = results["spliced"][3]
    if spliced_pct >= 85:
        print(
            f"PASS — {spliced_pct:.1f}% of ground-truth barcodes found in kb spliced output."
        )
        return True
    else:
        print(
            f"FAIL — only {spliced_pct:.1f}% overlap. Likely wrong sample or wrong SRR accessions."
        )
        return False


SAMPLES = {
    "GSM8253792": {
        "raw": "data/raw/GSM8253792_CD1_Cntl_0wksRecov_barcodes.tsv",
        "kb_output": "sra_runs/kb_output_GSM8253792",
    },
    "GSM8253793": {
        "raw": "data/raw/GSM8253793_CD1_Cntl_3wksRecov_barcodes.tsv",
        "kb_output": "sra_runs/kb_output_GSM8253793",
    },
    "GSM8253794": {
        "raw": "data/raw/GSM8253794_CD1_CR_Rep1_barcodes.tsv",
        "kb_output": "sra_runs/kb_output_GSM8253794",
    },
    "GSM8253795": {
        "raw": "data/raw/GSM8253795_CD1_CR_Rep2_barcodes.tsv",
        "kb_output": "sra_runs/kb_output_GSM8253795",
    },
    "GSM8253796": {
        "raw": "data/raw/GSM8253796_NesCre_Cntl_Rep1_barcodes.tsv",
        "kb_output": "sra_runs/kb_output_GSM8253796",
    },
    "GSM8253797": {
        "raw": "data/raw/GSM8253797_NesCre_Cntl_Rep2_barcodes.tsv",
        "kb_output": "sra_runs/kb_output_GSM8253797",
    },
    "GSM8253798": {
        "raw": "data/raw/GSM8253798_NesCre_CR_Rep1_barcodes.tsv",
        "kb_output": "sra_runs/kb_output_GSM8253798",
    },
    "GSM8253799": {
        "raw": "data/raw/GSM8253799_NesCre_CR_Rep2_barcodes.tsv",
        "kb_output": "sra_runs/kb_output_GSM8253799",
    },
    "GSM8647352": {
        "raw": "data/raw/GSM8647352_CD1_CR1_NoRecov1_barcodes.tsv",
        "kb_output": "sra_runs/kb_output_GSM8647352",
    },
    "GSM8647353": {
        "raw": "data/raw/GSM8647353_CD1_CR1_NoRecov2_barcodes.tsv",
        "kb_output": "sra_runs/kb_output_GSM8647353",
    },
}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Validate kb-python output against GEO ground truth."
    )
    parser.add_argument("--gsm", choices=list(SAMPLES), help="GSM accession to validate")
    parser.add_argument("--raw", help="Path to ground-truth barcodes.tsv")
    parser.add_argument("--kb-output", help="Path to kb count output directory")
    parser.add_argument(
        "--matrix-only", action="store_true",
        help="Skip barcode overlap check; only validate matrix dimensions and counts"
    )
    args = parser.parse_args()

    if args.gsm:
        s = SAMPLES[args.gsm]
        raw_path, kb_dir = s["raw"], s["kb_output"]
    elif args.kb_output and (args.raw or args.matrix_only):
        raw_path, kb_dir = args.raw, args.kb_output
    else:
        print("Usage:")
        print("  python3 validate_barcodes.py --gsm GSM8253792")
        print(f"  valid --gsm values: {', '.join(SAMPLES)}")
        print("  python3 validate_barcodes.py --gsm GSM8253796 --matrix-only")
        print("  python3 validate_barcodes.py --raw <path> --kb-output <dir>")
        sys.exit(1)

    ok = True
    if not args.matrix_only:
        ok = validate(raw_path, kb_dir)
        print()
    ok = validate_matrix(kb_dir) and ok
    sys.exit(0 if ok else 1)
