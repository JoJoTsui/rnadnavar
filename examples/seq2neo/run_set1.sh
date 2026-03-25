#!/usr/bin/env bash
# Run set1 — Colorectal cancer (strict), all eligible samples.

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

python3 "$HERE/scripts/run_batch_from_json.py" \
    --config "$HERE/config/runner.yaml" \
    --set 1
