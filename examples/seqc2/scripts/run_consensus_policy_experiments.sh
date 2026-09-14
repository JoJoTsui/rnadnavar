#!/usr/bin/env bash
# Generate and benchmark consensus policy candidates from existing caller VCFs.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTDIR="${1:?pipeline output directory}"
PAIR="${2:?tumor-normal pair}"
DEST="${3:?experiment output directory}"
TARGET_BED="${TARGET_BED:?target BED for this benchmark cell}"
TRUTH="${TRUTH:?merged truth VCF for this benchmark cell}"
HC_BED="${HC_BED:?high-confidence truth BED}"
FASTA="${FASTA:?reference FASTA}"
HAPPY_ENV="${HAPPY_ENV:-happy}"
mkdir -p "$DEST"
TMP="$DEST/.caller_inputs"
mkdir -p "$TMP"
trap 'rm -rf "$TMP"' EXIT

src_dir="$OUTDIR/variant_calling"
link_caller() {
  local caller="$1" pattern="$2" source
  source="$(find "$src_dir/$caller/$PAIR" -maxdepth 1 -type f -name "$pattern" | head -1)"
  [ -n "$source" ] || { echo "missing $caller source under $src_dir" >&2; exit 1; }
  source="$(readlink -f "$source")"
  ln -sf "$source" "$TMP/$PAIR.$caller.vcf.gz"
}
link_caller mutect2 '*.mutect2.filtered.vcf.gz'
link_caller deepsomatic '*.deepsomatic.vcf.gz'
link_caller strelka '*.strelka.variants.vcf.gz'

for policy in ordinary native strict; do
  prefix="$DEST/$policy/$PAIR.consensus"
  mkdir -p "${prefix%/*}"
  args=(--input_dir "$TMP" --expected_callers mutect2,deepsomatic,strelka
        --out_prefix "$prefix" --snv_thr 2 --indel_thr 2 --min_alt_support 3)
  [ "$policy" = native ] && args+=(--native-evidence-snv)
  [ "$policy" = strict ] && args+=(--snv_thr 3 --indel_thr 3)
  python3 "$HERE/../../../bin/run_consensus_vcf.py" "${args[@]}"
  query="$prefix.vcf.gz"
  pass="$DEST/$policy/$PAIR.consensus.somatic.vcf.gz"
  bcftools view -i 'FILTER="Somatic"' "$query" | awk 'BEGIN{OFS="\t"} /^#/{print; next} {$7="PASS"; print}' | bgzip -c > "$pass"
  bcftools index -t "$pass"
  micromamba run -n "$HAPPY_ENV" som.py "$TRUTH" "$pass" -R "$HC_BED" -T "$TARGET_BED" -r "$FASTA" -N -o "$DEST/$policy/benchmark"
done

python3 - "$DEST" "$PAIR" "$TARGET_BED" "$TRUTH" "$HC_BED" "$FASTA" <<'PY'
import json
import sys
from pathlib import Path

dest, pair, target, truth, hc, fasta = sys.argv[1:]
Path(dest, "provenance.json").write_text(json.dumps({
    "pair": pair,
    "policies": ["ordinary", "native", "strict"],
    "target_bed": target,
    "truth": truth,
    "truth_bed": hc,
    "reference": fasta,
    "source_outputs_unchanged": True,
}, indent=2, sort_keys=True) + "\n")
PY
