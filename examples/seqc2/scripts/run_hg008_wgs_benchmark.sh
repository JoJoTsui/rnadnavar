#!/usr/bin/env bash
# Execute the independent HG008 somatic WGS benchmark.
# Set HG008_* variables to the benchmark truth/region/reference and caller VCFs.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
: "${HG008_TRUTH_SNV:?Set HG008_TRUTH_SNV}"
: "${HG008_TRUTH_INDEL:?Set HG008_TRUTH_INDEL}"
: "${HG008_REGIONS:?Set HG008_REGIONS}"
: "${HG008_FASTA:?Set HG008_FASTA}"
: "${HG008_PIPELINE_OUTDIR:?Set HG008_PIPELINE_OUTDIR}"
: "${HG008_PAIR:?Set HG008_PAIR}"
: "${HG008_COMPARE_DIR:?Set HG008_COMPARE_DIR}"
export BENCHMARK_MODE=wgs
export TRUTH_SNV="$HG008_TRUTH_SNV" TRUTH_INDEL="$HG008_TRUTH_INDEL" HC_BED="$HG008_REGIONS" FASTA="$HG008_FASTA"
exec bash "$HERE/run_benchmark.sh" "$HG008_PIPELINE_OUTDIR" "$HG008_PAIR" "$HG008_COMPARE_DIR"
