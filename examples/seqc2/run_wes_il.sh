#!/usr/bin/env bash
# Run the deeper SEQC2 WES_IL tumor-normal pair (BAM input, DNA-only).
#
# The original WES_LL launcher remains available as run_wes_ll.sh.
#
# Usage:
#   bash run_wes_il.sh             (resume automatically)
#   bash run_wes_il.sh --dry-run   (print the Nextflow command)

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

python3 "$HERE/scripts/run_pipeline.py" \
    --config "$HERE/config/runner.yaml" \
    --input "$HERE/csv/seqc2_wes_il.csv" \
    --outdir "$HERE/output/seqc2.wes.il" \
    "$@"
