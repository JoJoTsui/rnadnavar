#!/usr/bin/env bash
# Run and preserve all four dataset/domain cells. Each output directory is an
# existing caller-ready result; this script never invokes mapping or calling.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
OUT_ROOT="${1:?comparison output root}"
WES_OUT="${2:?WES-LL output directory}"
WES_PAIR="${3:?WES-LL pair}"
WGS_OUT="${4:?WGS-IL output directory}"
WGS_PAIR="${5:?WGS-IL pair}"
for dataset in wes_ll wgs_il; do
  out="$WES_OUT"; pair="$WES_PAIR"; [ "$dataset" = wgs_il ] && out="$WGS_OUT" && pair="$WGS_PAIR"
  for domain in ukb medexome; do
    target="${UKB_BED:-/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/bio_db/intervals/ukb.pad50.broad.pad50.union.bed}"
    [ "$domain" = medexome ] && target="${MEDEXOME_BED:-$ROOT/data/SeqCap_EZ_MedExome_hg38_empirical_targets.authoritative.bed}"
    mkdir -p "$OUT_ROOT/$dataset/$domain"
    TARGET_BED="$target" BENCHMARK_MODE="$dataset" \
      bash "$HERE/run_benchmark.sh" "$out" "$pair" "$OUT_ROOT/$dataset/$domain"
  done
done
printf 'matrix complete: %s/{wes_ll,wgs_il}/{ukb,medexome}\n' "$OUT_ROOT"
