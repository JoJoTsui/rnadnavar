#!/usr/bin/env bash
# Default: preparation only. Shared destinations come from the committed config.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
usage() {
    echo 'Usage: bash examples/seq2neo/run_refined_cohort.sh [--prepare | --execute --approve-pilot]'
}
mode=prepare
approved=false
for arg in "$@"; do
    case "$arg" in
        --prepare) mode=prepare ;;
        --execute) mode=execute ;;
        --approve-pilot) approved=true ;;
        --help|-h) usage; exit 0 ;;
        *) usage >&2; exit 2 ;;
    esac
done
if [[ "$mode" = execute && "$approved" != true ]]; then
    echo 'Full generation requires --approve-pilot after reviewing the pilot warning.' >&2
    exit 2
fi
PYTHON="$ROOT/.venv/bin/python"
CONFIG="$ROOT/examples/seq2neo/config/refined_native_v2_cohort.json"
OUT="$("$PYTHON" -c 'import json,sys; print(json.load(open(sys.argv[1]))["output_root"])' "$CONFIG")"
[[ "$OUT" = /* ]] || { echo 'Require absolute output root' >&2; exit 2; }
mkdir -p "$OUT/logs"
LOG="$(mktemp "$OUT/logs/$mode.XXXXXXXX.log")"
args=(--config "$CONFIG")
if [[ "$mode" = execute ]]; then
    args+=(--all --approve-pilot --execute)
fi
printf 'Mode: %s\nCode: %s\nOutputs: %s\nLog: %s\n' "$mode" "$ROOT" "$OUT" "$LOG"
exec "$PYTHON" -u "$ROOT/examples/seq2neo/scripts/run_refined_cohort.py" "${args[@]}" > "$LOG" 2>&1
