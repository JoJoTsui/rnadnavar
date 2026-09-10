#!/usr/bin/env bash
# Run the HG008 paired-BAM/RNA-FASTQ hybrid workflow, then benchmark its VCFs.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTDIR="${1:-$HERE/output/hg008.wgs.hybrid}"
PAIR="${2:-HG008_T_1_vs_HG008_N_1}"
COMPARE="${3:-$HERE/comparison/$PAIR}"
bash "$HERE/run_hg008_wgs_hybrid.sh" --outdir "$OUTDIR"

# Require a final rescue VCF; caller-only comparisons are not sufficient for
# the hybrid benchmark acceptance surface. Prefer the annotated realignment
# artifact, then fall back to the first filtered rescue VCF.
if [ -z "${RESCUE_VCF:-}" ]; then
  RESCUE_VCF="$(find "$OUTDIR/vcf_realignment/rescue" "$OUTDIR/rescue" -type f \
    -name '*stripped.vep.vcf.gz' 2>/dev/null | sort | head -n 1 || true)"
  if [ -z "$RESCUE_VCF" ]; then
    RESCUE_VCF="$(find "$OUTDIR/rescue" -type f -name '*.filtered.vcf.gz' \
      2>/dev/null | sort | head -n 1 || true)"
  fi
fi
[ -n "$RESCUE_VCF" ] && [ -f "$RESCUE_VCF" ] || { echo "ERROR: final rescue VCF not found under $OUTDIR" >&2; exit 1; }
export RESCUE_VCF
HG008_PIPELINE_OUTDIR="$OUTDIR" HG008_PAIR="$PAIR" HG008_COMPARE_DIR="$COMPARE" \
  bash "$HERE/../scripts/run_hg008_wgs_benchmark.sh"
