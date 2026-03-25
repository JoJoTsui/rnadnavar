#!/usr/bin/env bash
# Run set4 — Gastric / Lung / Pancreatic / Rectal.

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

python3 "$HERE/scripts/run_batch_from_json.py" \
    --config "$HERE/config/runner.yaml" \
    --set 4
