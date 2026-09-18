#!/usr/bin/env bash
# Candidate-only preparation by default. Never grants biological training approval.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
usage() {
    echo 'Usage: bash examples/seq2neo/run_three_class_cohort.sh [--prepare | --pilot | --execute --approve-pilot]'
}
mode=prepare
selected=false
approved=false
for arg in "$@"; do
    case "$arg" in
        --prepare|--pilot|--execute)
            [[ "$selected" = false ]] || { usage >&2; exit 2; }
            mode="${arg#--}"; selected=true ;;
        --approve-pilot) approved=true ;;
        --help|-h) usage; exit 0 ;;
        *) usage >&2; exit 2 ;;
    esac
done
if [[ "$mode" = execute && "$approved" != true ]] || [[ "$mode" != execute && "$approved" = true ]]; then
    echo '--approve-pilot is required only with --execute, after reviewing the new-policy pilot.' >&2
    exit 2
fi
PYTHON="$ROOT/.venv/bin/python"
CONFIG="$ROOT/examples/seq2neo/config/separated_three_class_v2_cohort.json"
OUT="$("$PYTHON" -c 'import json,sys; print(json.load(open(sys.argv[1]))["output_root"])' "$CONFIG")"
[[ "$OUT" = /* ]] || { echo 'Require absolute output root' >&2; exit 2; }
mkdir -p "$OUT/logs"
LOG="$(mktemp "$OUT/logs/$mode.XXXXXXXX.log")"
args=(--config "$CONFIG")
if [[ "$mode" = pilot ]]; then args+=(--pilot --execute); fi
if [[ "$mode" = execute ]]; then args+=(--all --approve-pilot --execute); fi
printf 'Candidate-only mode: %s\nCode: %s\nOutputs: %s\nLog: %s\n' "$mode" "$ROOT" "$OUT" "$LOG"
exec "$PYTHON" -u "$ROOT/examples/seq2neo/scripts/run_refined_cohort.py" "${args[@]}" > "$LOG" 2>&1
