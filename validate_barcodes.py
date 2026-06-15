import argparse
import sys


def validate(raw_barcodes_path: str, kb_output_dir: str) -> bool:
    raw = set(l.strip().replace("-1", "") for l in open(raw_barcodes_path) if l.strip())

    results = {}
    for layer in ("spliced", "unspliced"):
        kb_path = f"{kb_output_dir}/counts_unfiltered/{layer}.barcodes.txt"
        kb = set(l.strip() for l in open(kb_path) if l.strip())
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
    "GSM8253796": {
        "raw": "data/raw/GSM8253796_NesCre_Cntl_Rep1_barcodes.tsv",
        "kb_output": "sra_runs/kb_output_GSM8253798",
    },
    "GSM8253798": {
        "raw": "data/raw/GSM8253798_NesCre_CR_Rep1_barcodes.tsv",
        "kb_output": "sra_runs/kb_output_GSM8253798",
    },
}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Validate kb-python barcodes against GEO ground truth."
    )
    parser.add_argument(
        "--gsm", choices=list(SAMPLES), help="GSM accession to validate"
    )
    parser.add_argument("--raw", help="Path to ground-truth barcodes.tsv")
    parser.add_argument("--kb-output", help="Path to kb count output directory")
    args = parser.parse_args()

    if args.gsm:
        s = SAMPLES[args.gsm]
        raw_path, kb_dir = s["raw"], s["kb_output"]
    elif args.raw and args.kb_output:
        raw_path, kb_dir = args.raw, args.kb_output
    else:
        print("Usage:")
        print("  python validate_barcodes.py --gsm GSM8253796")
        print("  python validate_barcodes.py --gsm GSM8253798")
        print("  python validate_barcodes.py --raw <path> --kb-output <dir>")
        sys.exit(1)

    ok = validate(raw_path, kb_dir)
    sys.exit(0 if ok else 1)
