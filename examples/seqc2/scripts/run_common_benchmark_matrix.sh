#!/usr/bin/env bash
# Run the reproducible UKB/MedExome benchmark matrix for an existing output.
# This wrapper only consumes caller/rescue VCFs; it never maps or calls variants.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
OUTDIR="${1:?pipeline output directory}"
PAIR="${2:?tumor-normal pair}"
BASE="${3:?comparison output root}"
MED="${MEDEXOME_BED:-$ROOT/data/SeqCap_EZ_MedExome_hg38_empirical_targets.authoritative.bed}"
UKB="${UKB_BED:-/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/bio_db/intervals/ukb.pad50.broad.pad50.union.bed}"

for name in ukb medexome; do
  target="$UKB"; [ "$name" = medexome ] && target="$MED"
  mkdir -p "$BASE/$name"
  TARGET_BED="$target" NORMALIZE_ALL="${NORMALIZE_ALL:-1}" \
    bash "$HERE/run_benchmark.sh" "$OUTDIR" "$PAIR" "$BASE/$name"
done
python3 "$HERE/compare_stage_transitions.py" --help >/dev/null
