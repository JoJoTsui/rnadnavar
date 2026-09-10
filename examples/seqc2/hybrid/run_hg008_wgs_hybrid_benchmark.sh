#!/usr/bin/env bash
# Run the HG008 paired-BAM/RNA-FASTQ hybrid workflow, then benchmark its VCFs.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTDIR="${1:-$HERE/output/hg008.wgs.hybrid}"
PAIR="${2:-HG008_T_1_vs_HG008_N_1}"
COMPARE="${3:-$HERE/comparison/$PAIR}"
bash "$HERE/run_hg008_wgs_hybrid.sh" --outdir "$OUTDIR"
HG008_PIPELINE_OUTDIR="$OUTDIR" HG008_PAIR="$PAIR" HG008_COMPARE_DIR="$COMPARE" \
  bash "$HERE/../scripts/run_hg008_wgs_benchmark.sh"
