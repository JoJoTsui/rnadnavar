#!/usr/bin/env bash
# Execute a SEQC2 WGS benchmark with the same caller comparison as run_benchmark.sh.
# Defaults use the SEQC2 truth/reference and UKB BED for both -R and -T; all can be overridden.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SEQ2C_ROOT="${SEQ2C_ROOT:-/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/data/giab/data/seqc2}"
export TRUTH_SNV="${TRUTH_SNV:-$SEQ2C_ROOT/truth/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz}"
export TRUTH_INDEL="${TRUTH_INDEL:-$SEQ2C_ROOT/truth/high-confidence_sINDEL_in_HC_regions_v1.2.1.vcf.gz}"
export HC_BED="${HC_BED:-$SEQ2C_ROOT/truth/High-Confidence_Regions_v1.2.bed}"
export FASTA="${FASTA:-/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/bio_db/references/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta}"
export TARGET_BED="${TARGET_BED:-/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/bio_db/intervals/ukb.pad50.broad.pad50.union.bed}"
unset CLAIR_VCF
export BENCHMARK_MODE=wgs
exec bash "$HERE/run_benchmark.sh" "${1:?pipeline output directory}" "${2:?tumor-normal pair}" "${3:?comparison output directory}"
