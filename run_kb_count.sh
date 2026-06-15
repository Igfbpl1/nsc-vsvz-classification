#!/bin/bash
# Run kb count for a single GSM sample.
# Usage: bash run_kb_count.sh <GSM_ID>
# All reference files and FASTQs must be in sra_runs/.

set -euo pipefail

GSM="${1:-}"
if [[ -z "$GSM" ]]; then
    echo "Usage: bash run_kb_count.sh <GSM_ID>"
    echo "  bash run_kb_count.sh GSM8253796"
    echo "  bash run_kb_count.sh GSM8253798"
    exit 1
fi

cd "$(dirname "$0")/sra_runs"

case "$GSM" in
    GSM8253796)
        # NesCre Cntl — 4 runs
        FASTQS=(
            SRR28912875_3.fastq SRR28912875_4.fastq
            SRR28912876_3.fastq SRR28912876_4.fastq
            SRR28912884_3.fastq SRR28912884_4.fastq
            SRR28912887_3.fastq SRR28912887_4.fastq
        )
        ;;
    GSM8253798)
        # NesCre CR — 3 runs
        FASTQS=(
            SRR28912872_3.fastq SRR28912872_4.fastq
            SRR28912873_3.fastq SRR28912873_4.fastq
            SRR28912874_3.fastq SRR28912874_4.fastq
        )
        ;;
    GSM8253799)
        # NesCre CupRap Rep2 — 3 runs (SRR28912869: 32GB, SRR28912870: 32GB, SRR28912871: 13GB)
        FASTQS=(
            SRR28912869_3.fastq SRR28912869_4.fastq
            SRR28912870_3.fastq SRR28912870_4.fastq
            SRR28912871_3.fastq SRR28912871_4.fastq
        )
        ;;
    *)
        echo "Unknown GSM: $GSM. Add it to run_kb_count.sh."
        exit 1
        ;;
esac

echo "Running kb count for $GSM with ${#FASTQS[@]} files..."

kb count \
    -i index.idx \
    -g t2g.txt \
    -x 10xv3 \
    -o "kb_output_${GSM}" \
    -c1 cdna_t2c.txt \
    -c2 intron_t2c.txt \
    --workflow lamanno \
    --loom \
    --kallisto /opt/homebrew/bin/kallisto \
    --bustools /opt/homebrew/bin/bustools \
    "${FASTQS[@]}"

echo "Done. Output: sra_runs/kb_output_${GSM}/"
