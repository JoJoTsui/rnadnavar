#!/usr/bin/env bash
# Run the deepest available SEQC2 WGS tumor-normal pair (WGS_IL).
# The WGS_LL pair is not present in the mounted dataset.
# Uses the shared SEQC2 config and its ukb.pad50 capture BED.

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

python3 "$HERE/scripts/run_pipeline.py" \
    --config "$HERE/config/runner.yaml" \
    --input "$HERE/csv/seqc2_wgs_il.csv" \
    --outdir "$HERE/output/seqc2.wgs.il" \
    "$@"
