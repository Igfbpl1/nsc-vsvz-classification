#!/bin/bash
# process_all_remaining.sh — process all 7 remaining samples sequentially.
#
# Usage: nohup bash process_all_remaining.sh > process_all.log 2>&1 &
# Monitor: tail -f process_all.log
#
# Runs process_sample.sh for each remaining GSM, one at a time.
# Each sample is fully processed + packaged + cleaned up before the next starts,
# so peak disk usage is bounded by the largest single sample (~150GB FASTQs).
#
# Expected runtime: ~3-4 hours per sample × 7 = ~25 hours on r6i.2xlarge.

set -euo pipefail

cd ~/project

SAMPLES=(
    GSM8253792   # CD1 Cntl 0wks
    GSM8253793   # CD1 Cntl 3wks
    GSM8253794   # CD1 CupRap 3wks
    GSM8253795   # CD1 CupRap 3wks
    GSM8253797   # NesCre Cntl 3wks
    GSM8647352   # CD1 CupRap 0wks
    GSM8647353   # CD1 CupRap 0wks
)

for GSM in "${SAMPLES[@]}"; do
    echo ""
    echo "==============================================================="
    echo "  Processing $GSM at $(date)"
    echo "==============================================================="
    bash process_sample.sh "$GSM"
    echo "  Finished $GSM at $(date)"
done

echo ""
echo "==============================================================="
echo "  All ${#SAMPLES[@]} samples done at $(date)"
echo "==============================================================="
ls -lh ~/project/sra_runs/GSM*_counts.tgz
