#!/usr/bin/env bash
# Execute the independent HG008 somatic WGS benchmark.
# Set HG008_* variables to the benchmark truth/region/reference and caller VCFs.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export HG008_ROOT="${HG008_ROOT:-/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/data/giab/data/giab_benchmarks/smvar_v0.3}"
export HG008_TRUTH_VCF="${HG008_TRUTH_VCF:-$HG008_ROOT/HG008-T_somatic_smvar_benchmark_v0.3_somatic_tumornormal.vcf.gz}"
export HG008_REGIONS="${HG008_REGIONS:-$HG008_ROOT/HG008-T_somatic_smvar_benchmark_v0.3_all.bed}"
export HG008_FASTA="${HG008_FASTA:-/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/bio_db/references/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta}"
: "${HG008_PIPELINE_OUTDIR:?Set HG008_PIPELINE_OUTDIR}"
: "${HG008_PAIR:?Set HG008_PAIR}"
: "${HG008_COMPARE_DIR:?Set HG008_COMPARE_DIR}"
unset CLAIR_VCF
export BENCHMARK_MODE=wgs
export TRUTH_VCF="$HG008_TRUTH_VCF" HC_BED="$HG008_REGIONS" FASTA="$HG008_FASTA" TARGET_BED="${TARGET_BED:-/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/bio_db/intervals/ukb.pad50.broad.pad50.union.bed}"
exec bash "$HERE/run_benchmark.sh" "$HG008_PIPELINE_OUTDIR" "$HG008_PAIR" "$HG008_COMPARE_DIR"
