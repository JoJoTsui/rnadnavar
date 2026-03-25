#!/usr/bin/env bash
# Run a single patient.
#
# Usage:
#   PROJECT=PRJNA298376 PATIENT=4060 bash run_single.sh
#   bash run_single.sh   (uses defaults below)

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

PROJECT=${PROJECT:-PRJNA298376}
PATIENT=${PATIENT:-4060}

python3 "$HERE/scripts/run_batch_from_json.py" \
    --config "$HERE/config/runner.yaml" \
    --project "$PROJECT" \
    --patient "$PATIENT"
