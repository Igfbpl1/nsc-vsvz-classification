#!/bin/bash
# process_sample.sh — download & process one sample for velocity analysis.
#
# Usage:    bash process_sample.sh <SAMPLE_NAME> <SRR_LIST>
# Example:  bash process_sample.sh "CupRap_Rep2" "SRR28912885 SRR28912886"
#
# Can also be set via environment variable:
#   SRR_LIST="SRR..." bash process_sample.sh <SAMPLE_NAME>

set -euo pipefail

# --- Argument Parsing ---
SAMPLE_NAME="${1:-}"
# If SRR_LIST is passed as the second argument, use it. Otherwise, check env var.
SRR_LIST_ARG="${2:-}"
SRR_LIST="${SRR_LIST_ARG:-$SRR_LIST}" # Use arg, fallback to env var

if [[ -z "$SAMPLE_NAME" ]] || [[ -z "$SRR_LIST" ]]; then
    echo "Usage: bash process_sample.sh <SAMPLE_NAME> \"<SRR_LIST>\""
    echo "Example: bash process_sample.sh \"CupRap_Rep2\" \"SRR28912885 SRR28912886\""
    exit 1
fi

cd ~/project/sra_runs

# ── Skip if already done ──────────────────────────────────────────────────
if [[ -f "${SAMPLE_NAME}_counts.tgz" ]]; then
    echo "✓ Already packaged: ${SAMPLE_NAME}_counts.tgz"
    ls -lh "${SAMPLE_NAME}_counts.tgz"
    exit 0
fi
if [[ -d "kb_output_${SAMPLE_NAME}/counts_unfiltered" ]] && \
   [[ -f "kb_output_${SAMPLE_NAME}/counts_unfiltered/spliced.mtx" ]]; then
    echo "✓ counts_unfiltered exists, just packaging..."
    tar czf "${SAMPLE_NAME}_counts.tgz" -C "kb_output_${SAMPLE_NAME}" counts_unfiltered
    ls -lh "${SAMPLE_NAME}_counts.tgz"
    exit 0
fi

# ── Set SRR IDs ───────────────────────────────────────────────────────────
echo "[1/5] Using provided SRR IDs for ${SAMPLE_NAME} ..."
echo "  SRRs: $SRR_LIST"

# ── Download FASTQs (sequentially to avoid disk spikes) ──────────────────
echo "[2/5] Downloading FASTQs ..."
for SRR in $SRR_LIST; do
    if { [[ -f "${SRR}_3.fastq" ]] && [[ -f "${SRR}_4.fastq" ]]; } || \
       { [[ -f "${SRR}_2.fastq" ]] && [[ -f "${SRR}_3.fastq" ]]; }; then
        echo "  ${SRR}: already present"
        continue
    fi
    echo "  prefetching ${SRR} ..."
    prefetch "$SRR"
    echo "  fasterq-dumping ${SRR} ..."
    fasterq-dump --split-files --include-technical "$SRR" --outdir .
done

# ── Build FASTQ list and run kb count ────────────────────────────────────
echo "[3/5] Running kb count for ${SAMPLE_NAME} ..."
FASTQS=()
for SRR in $SRR_LIST; do
    if [[ -f "${SRR}_4.fastq" ]]; then
        # layout: _3=R1 (barcode+UMI), _4=R2 (cDNA)
        FASTQS+=("${SRR}_3.fastq" "${SRR}_4.fastq")
    else
        # layout: _1=I1 (index), _2=R1 (barcode+UMI), _3=R2 (cDNA)
        FASTQS+=("${SRR}_2.fastq" "${SRR}_3.fastq")
    fi
done

kb count \
    -i index.idx \
    -g t2g.txt \
    -x 10xv3 \
    -o "kb_output_${SAMPLE_NAME}" \
    -c1 cdna_t2c.txt \
    -c2 intron_t2c.txt \
    --workflow lamanno \
    "${FASTQS[@]}"

# ── Package counts_unfiltered ────────────────────────────────────────────
echo "[4/5] Packaging counts_unfiltered ..."
tar czf "${SAMPLE_NAME}_counts.tgz" -C "kb_output_${SAMPLE_NAME}" counts_unfiltered

# ── Cleanup to free disk ─────────────────────────────────────────────────
echo "[5/5] Cleaning up ..."
for SRR in $SRR_LIST; do
    rm -f "${SRR}"_*.fastq
    rm -rf "$SRR"
done

# remove BUS files and other intermediates from kb_output
cd "kb_output_${SAMPLE_NAME}"
rm -f *.bus matrix.ec transcripts.txt inspect*.json 10x_version3_whitelist.txt
rm -rf tmp
cd ..

# Optional: also delete kb_output_<SAMPLE_NAME>/ since the tarball has counts_unfiltered
# Uncomment if you want maximum disk space recovery:
# rm -rf "kb_output_${SAMPLE_NAME}"

echo ""
echo "✓ Done with ${SAMPLE_NAME}"
echo "  Output: $(pwd)/${SAMPLE_NAME}_counts.tgz"
ls -lh "${SAMPLE_NAME}_counts.tgz"
echo ""
echo "Disk usage:"
df -h . | tail -1
