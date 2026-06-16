#!/bin/bash
# process_all_remaining.sh — process all remaining samples sequentially.
#
# This script now defines samples as groups of SRR runs and processes each
# group by calling the updated process_sample.sh script.

set -euo pipefail
cd "$(dirname "$0")"

# --- Define Sample Groups ---
# Each key is the sample name, and the value is a space-separated list of SRR IDs.

declare -A SAMPLES
SAMPLES["GSM8253793"]="SRR28912882 SRR28912883"
SAMPLES["GSM8253794"]="SRR28912878 SRR28912879 SRR28912880 SRR28912881"
SAMPLES["GSM8253797"]="SRR28912867 SRR28912868"
SAMPLES["GSM8647352"]="SRR31443698 SRR31443699"
SAMPLES["GSM8647353"]="SRR31443695 SRR31443696 SRR31443697 SRR31443700"

# --- Process Samples ---

for SAMPLE_NAME in "${!SAMPLES[@]}"; do
    SRR_LIST="${SAMPLES[$SAMPLE_NAME]}"
    echo ""
    echo "==============================================================="
    echo "  Processing $SAMPLE_NAME at $(date)"
    echo "  SRRs: $SRR_LIST"
    echo "==============================================================="

    # Call the updated processing script, passing sample name and SRR list
    bash process_sample.sh "$SAMPLE_NAME" "$SRR_LIST"

    echo "  Finished $SAMPLE_NAME at $(date)"
done

echo ""
echo "==============================================================="
echo "  All ${#SAMPLES[@]} samples done at $(date)"
echo "==============================================================="
ls -lh sra_runs/*_counts.tgz
