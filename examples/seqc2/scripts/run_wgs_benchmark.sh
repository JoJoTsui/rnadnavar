#!/usr/bin/env bash
# Execute a SEQC2 WGS benchmark with the same caller comparison as run_benchmark.sh.
# Required overrides are truth/reference/output-specific; no WES -T target is used.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
: "${TRUTH_SNV:?Set TRUTH_SNV to the SEQC2 WGS SNV truth VCF}"
: "${TRUTH_INDEL:?Set TRUTH_INDEL to the SEQC2 WGS indel truth VCF}"
: "${HC_BED:?Set HC_BED to the WGS evaluation BED}"
: "${FASTA:?Set FASTA to the matching reference FASTA}"
export BENCHMARK_MODE=wgs
exec bash "$HERE/run_benchmark.sh" "${1:?pipeline output directory}" "${2:?tumor-normal pair}" "${3:?comparison output directory}"
