#!/usr/bin/env bash
# Execute the independent HG008 somatic WGS benchmark.
# Set HG008_* variables to the benchmark truth/region/reference and caller VCFs.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export HG008_ROOT="${HG008_ROOT:-/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/data/giab/data/giab_benchmarks/smvar_v0.3}"
export HG008_TRUTH_VCF="${HG008_TRUTH_VCF:-$HG008_ROOT/HG008-T_somatic_smvar_benchmark_v0.3_somatic_tumornormal.vcf.gz}"
export HG008_REGIONS="${HG008_REGIONS:-$HG008_ROOT/HG008-T_somatic_smvar_benchmark_v0.3_all.bed}"
# Use the same GIAB GRCh38 reference as the HG008 pipeline profile so
# normalization and benchmarking operate in one coordinate/reference system.
export HG008_FASTA="${HG008_FASTA:-/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/bio_db/references/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta}"
: "${HG008_PIPELINE_OUTDIR:?Set HG008_PIPELINE_OUTDIR}"
: "${HG008_PAIR:?Set HG008_PAIR}"
: "${HG008_COMPARE_DIR:?Set HG008_COMPARE_DIR}"
unset CLAIR_VCF
export BENCHMARK_MODE=wgs
# The HG008 benchmark BED defines the comparison domain. Use it for both
# hap.py region and target restrictions unless a narrower domain is explicit.
export TRUTH_VCF="$HG008_TRUTH_VCF" HC_BED="$HG008_REGIONS" FASTA="$HG008_FASTA" TARGET_BED="${TARGET_BED:-$HG008_REGIONS}"
exec bash "$HERE/run_benchmark.sh" "$HG008_PIPELINE_OUTDIR" "$HG008_PAIR" "$HG008_COMPARE_DIR"
