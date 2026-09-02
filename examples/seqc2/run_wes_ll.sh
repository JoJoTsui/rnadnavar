#!/usr/bin/env bash
# Run the SEQC2 WES_LL tumor-normal pair (BAM input, DNA-only, consensus).
#
# Usage:
#   bash run_wes_ll.sh             (uses config/runner.yaml defaults)
#   bash run_wes_ll.sh --dry-run   (print the nextflow command, don't run)
#
# Extra args are passed through to scripts/run_pipeline.py.

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

python3 "$HERE/scripts/run_pipeline.py" \
    --config "$HERE/config/runner.yaml" \
    "$@"
