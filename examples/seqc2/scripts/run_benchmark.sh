#!/usr/bin/env bash
# Benchmark pipeline VCFs against the SEQC2 truth set using som.py (hap.py family).
#
# Compares, for a tumor-normal pair run produced by run_wes_ll.sh:
#   consensus   (FILTER == Somatic only, converted to PASS for som.py)
#   rescue      (optional RESCUE_VCF; Somatic records converted to PASS)
#   mutect2     (PASS records of *.mutect2.filtered.vcf.gz)
#   deepsomatic (PASS records of *.deepsomatic.vcf.gz)
#   strelka     (PASS records of *.strelka.variants.vcf.gz)
# against the merged SEQC2 high-confidence sSNV+sINDEL truth, restricted to
# High-Confidence_Regions_v1.2.bed. Metrics (TP/FP/FN/P/R/F1, split by
# SNV/INDEL) are aggregated into a single CSV by aggregate_benchmark.py.
#
# Usage:
#   bash run_benchmark.sh [PIPELINE_OUTDIR] [PAIR] [COMPARE_DIR]
#
#   PIPELINE_OUTDIR  default: <examples/seqc2>/output/seqc2.wes.ll
#   PAIR             default: WES_LL_T_1_vs_WES_LL_N_1
#   COMPARE_DIR      default: <examples/seqc2>/comparison/<PAIR>
#   RESCUE_VCF       optional cross-modality rescue VCF to benchmark as rescue
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXAMPLE_ROOT="$(cd "$HERE/.." && pwd)"

OUTDIR="${1:-${OUTDIR:-$EXAMPLE_ROOT/output/seqc2.wes.ll}}"
PAIR="${2:-${PAIR:-WES_LL_T_1_vs_WES_LL_N_1}}"
OD="${3:-${COMPARE_DIR:-$EXAMPLE_ROOT/comparison/$PAIR}}"
mkdir -p "$OD"

# Prepared benchmark copies are derived artifacts. Never reuse one merely
# because its index exists: changed production VCFs must invalidate stale copies.
fingerprint_source() {
    local source="$1"; sha256sum "$source" | awk -v p="$source" '{print p"\t"$1}'
}
derived_is_current() {
    local source="$1" derived="$2" stamp="$derived.source.sha256"
    [ -f "$derived" ] && [ -f "$derived.tbi" ] && [ -f "$stamp" ] \
        && [ "$(fingerprint_source "$source")" = "$(cat "$stamp")" ]
}
record_derived_source() {
    local source="$1" derived="$2"; fingerprint_source "$source" > "$derived.source.sha256"
}

# SEQC2 data root + references (override via environment if needed)
SEQ2C="${SEQ2C_ROOT:-/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/data/giab/data/seqc2}"
TRUTH_SNV="${TRUTH_SNV:-$SEQ2C/truth/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz}"
TRUTH_INDEL="${TRUTH_INDEL:-$SEQ2C/truth/high-confidence_sINDEL_in_HC_regions_v1.2.1.vcf.gz}"
HC_BED="${HC_BED:-$SEQ2C/truth/High-Confidence_Regions_v1.2.bed}"
TARGET_BED="${TARGET_BED:-/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/bio_db/intervals/ukb.pad50.broad.pad50.union.bed}"
FA="${FASTA:-/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/bio_db/references/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta}"
HAPPY_ENV="${HAPPY_ENV:-happy}"

# Query VCFs from the pipeline output layout
C_VCF="$OUTDIR/consensus/$PAIR/$PAIR.consensus.vcf.gz"
M2_VCF="$OUTDIR/variant_calling/mutect2/$PAIR/$PAIR.mutect2.filtered.vcf.gz"
DS_VCF="$OUTDIR/variant_calling/deepsomatic/$PAIR/$PAIR.deepsomatic.vcf.gz"
S2_VCF="$OUTDIR/variant_calling/strelka/$PAIR/$PAIR.strelka.variants.vcf.gz"
CLAIR_VCF="${CLAIR_VCF:-}"
RESCUE_VCF="${RESCUE_VCF:-}"

for f in "$C_VCF" "$M2_VCF" "$DS_VCF" "$S2_VCF" "$TRUTH_SNV" "$TRUTH_INDEL" "$HC_BED" "$TARGET_BED" "$FA"; do
    [ -f "$f" ] || { echo "ERROR: missing input: $f" >&2; exit 1; }
done

# 1. Merged, indexed truth VCF (cached)
TRUTH="$OD/high-confidence_sSNV+INDEL_in_HC_regions_v1.2.1.vcf.gz"
if [ ! -f "$TRUTH.tbi" ]; then
    echo ">> Building merged truth VCF: $TRUTH"
    bcftools concat --allow-overlaps --remove-duplicates "$TRUTH_SNV" "$TRUTH_INDEL" \
        | bcftools sort -Oz -o "$TRUTH"
    tabix -p vcf "$TRUTH"
fi

# 2. Consensus query: keep FILTER == Somatic only, then rewrite FILTER to PASS
#    in a benchmark-only copy. The production VCF remains unchanged.
C_SOM="$OD/$PAIR.consensus.somatic.vcf.gz"
if ! derived_is_current "$C_VCF" "$C_SOM"; then
    echo ">> Preparing PASS-only consensus benchmark VCF"
    bcftools view -i 'FILTER="Somatic"' "$C_VCF" | awk 'BEGIN{OFS="\t"} /^#/{print; next} {$7="PASS"; print}' | bgzip -c > "$C_SOM"
    bcftools index -t "$C_SOM"
    record_derived_source "$C_VCF" "$C_SOM"
fi

run_som() {  # <name> <query> <extra som.py args...>
    local name="$1" query="$2"; shift 2
    echo ">> som.py: $name"
    micromamba run -n "$HAPPY_ENV" som.py \
        "$TRUTH" "$query" \
        -R "$HC_BED" \
        -T "$TARGET_BED" \
        -o "$OD/$name" \
        -r "$FA" -N "$@"
}

# 3. Per-query benchmark runs (caller VCFs: PASS records only, som.py default)
run_som consensus "$C_SOM"
run_som mutect2 "$M2_VCF"
run_som deepsomatic "$DS_VCF"
run_som strelka "$S2_VCF"

QUERIES=(consensus mutect2 deepsomatic strelka)
if [ -n "$CLAIR_VCF" ]; then
    [ -f "$CLAIR_VCF" ] || { echo "ERROR: missing Clair input: $CLAIR_VCF" >&2; exit 1; }
    run_som clair "$CLAIR_VCF"
    QUERIES+=(clair)
fi

# Optional cross-modality rescue query. Rescue VCFs use EnsembleVar's biological
# FILTER vocabulary, so benchmark only Somatic records in a PASS-only copy.
if [ -n "$RESCUE_VCF" ]; then
    [ -f "$RESCUE_VCF" ] || { echo "ERROR: missing rescue input: $RESCUE_VCF" >&2; exit 1; }
    RESCUE_SOM="$OD/rescue.somatic.vcf.gz"
    if ! derived_is_current "$RESCUE_VCF" "$RESCUE_SOM"; then
        echo ">> Preparing PASS-only rescue benchmark VCF"
        bcftools view -i 'FILTER="Somatic"' "$RESCUE_VCF" \
            | awk 'BEGIN{OFS="\t"} /^#/{print; next} {$7="PASS"; print}' \
            | bgzip -c > "$RESCUE_SOM"
        bcftools index -t "$RESCUE_SOM"
        record_derived_source "$RESCUE_VCF" "$RESCUE_SOM"
    fi
    run_som rescue "$RESCUE_SOM"
    QUERIES+=(rescue)
fi

# 4. Aggregate metrics into one comparison table
python3 "$HERE/aggregate_benchmark.py" \
    --metrics-dir "$OD" \
    --queries "${QUERIES[@]}" \
    --output "$OD/$PAIR.benchmark_comparison.csv"

echo ">> Done. Table: $OD/$PAIR.benchmark_comparison.csv"
