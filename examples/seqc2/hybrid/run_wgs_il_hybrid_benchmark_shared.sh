#!/usr/bin/env bash
# Shared-checkout wrapper for the SEQC2 WGS-IL hybrid workflow + benchmark.
#
# Usage:
#   bash run_wgs_il_hybrid_benchmark_shared.sh [outdir] [pair] [compare_dir] [log]
#
# All arguments are optional. The wrapper resolves defaults inside the shared
# checkout and redirects the complete workflow/benchmark log itself.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"

OUTDIR="${1:-$HERE/output/seqc2.wgs.il.hybrid}"
PAIR="${2:-WGS_IL_T_1_vs_WGS_IL_N_1}"
COMPARE_DIR="${3:-$HERE/comparison/$PAIR}"
LOG="${4:-$HERE/wgs_il_benchmark.nohup.log}"

mkdir -p "$OUTDIR" "$COMPARE_DIR" "$(dirname "$LOG")"

cd "$REPO"
exec bash "$HERE/run_wgs_il_hybrid_benchmark.sh" \
  "$OUTDIR" "$PAIR" "$COMPARE_DIR" \
  >"$LOG" 2>&1
