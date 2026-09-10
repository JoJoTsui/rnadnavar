#!/usr/bin/env bash
# Run the paired-BAM/RNA-FASTQ hybrid workflow, then benchmark its VCFs.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTDIR="${1:-$HERE/output/seqc2.wgs.il.hybrid}"
PAIR="${2:-WGS_IL_T_1_vs_WGS_IL_N_1}"
COMPARE="${3:-$HERE/comparison/$PAIR}"
bash "$HERE/run_wgs_il_hybrid.sh" --outdir "$OUTDIR"

# Require a final rescue VCF; caller-only comparisons are not sufficient for
# the hybrid benchmark acceptance surface. Prefer the annotated realignment
# artifact, then fall back to the first filtered rescue VCF.
if [ -z "${RESCUE_VCF:-}" ]; then
  RESCUE_VCF="$(find "$OUTDIR/rescue" "$OUTDIR/vcf_realignment/rescue" -type f \
    \( -name '*stripped.vep.vcf.gz' -o -name '*.filtered.vcf.gz' \) \
    2>/dev/null | sort | head -n 1 || true)"
fi
[ -n "$RESCUE_VCF" ] && [ -f "$RESCUE_VCF" ] || { echo "ERROR: final rescue VCF not found under $OUTDIR" >&2; exit 1; }
export RESCUE_VCF
bash "$HERE/../scripts/run_wgs_benchmark.sh" "$OUTDIR" "$PAIR" "$COMPARE"
