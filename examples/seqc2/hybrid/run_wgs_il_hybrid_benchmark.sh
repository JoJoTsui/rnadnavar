#!/usr/bin/env bash
# Run the paired-BAM/RNA-FASTQ hybrid workflow, then benchmark its VCFs.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTDIR="${1:-$HERE/output/seqc2.wgs.il.hybrid}"
PAIR="${2:-WGS_IL_T_1_vs_WGS_IL_N_1}"
COMPARE="${3:-$HERE/comparison/$PAIR}"
bash "$HERE/run_wgs_il_hybrid.sh" --outdir "$OUTDIR"
bash "$HERE/../scripts/run_wgs_benchmark.sh" "$OUTDIR" "$PAIR" "$COMPARE"
