#!/usr/bin/env bash
# Dry-run: preview commands for any filter combination without executing.
# Passes all extra arguments through to run_batch_from_json.py.
#
# Usage:
#   bash dry_run.sh                            # all eligible samples
#   bash dry_run.sh --set 1                    # set1 only
#   bash dry_run.sh --set 2 --status standard
#   bash dry_run.sh --project PRJNA298376 --patient 4060
#   bash dry_run.sh --disease "Pancreatic"

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

python3 "$HERE/scripts/run_batch_from_json.py" \
    --config "$HERE/config/runner.yaml" \
    --dry-run \
    "$@"
