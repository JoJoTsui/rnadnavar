#!/usr/bin/env bash
# Shared-checkout wrapper for the HG008 WGS hybrid workflow + benchmark.
#
# Usage:
#   bash run_hg008_wgs_hybrid_benchmark_shared.sh [outdir] [pair] [compare_dir] [log]
#
# All arguments are optional. The underlying launcher validates the HG008
# truth/reference inputs and requires the final realignment-rescue VEP VCF.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"

OUTDIR="${1:-$HERE/output/hg008.wgs.hybrid}"
PAIR="${2:-HG008_T_1_vs_HG008_N_1}"
COMPARE_DIR="${3:-$HERE/comparison/$PAIR}"
LOG="${4:-$HERE/hg008_wgs_benchmark.nohup.log}"

mkdir -p "$OUTDIR" "$COMPARE_DIR" "$(dirname "$LOG")"

cd "$REPO"
export HG008_PIPELINE_OUTDIR="$OUTDIR"
export HG008_PAIR="$PAIR"
export HG008_COMPARE_DIR="$COMPARE_DIR"
exec bash "$HERE/run_hg008_wgs_hybrid_benchmark.sh" \
  >"$LOG" 2>&1
