#!/usr/bin/env bash
# Run set3 — Ampullary / Bile Duct / Cholangiocarcinoma / Esophageal / Melanoma.

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

python3 "$HERE/scripts/run_batch_from_json.py" \
    --config "$HERE/config/runner.yaml" \
    --set 3
