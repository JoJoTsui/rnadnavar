#!/usr/bin/env bash
# Run variant statistics pipeline on seq2neo output.
#
# Usage:
#   bash run_stats.sh                              # all sets, all wises
#   bash run_stats.sh --set 1                      # set 1 only
#   bash run_stats.sh --max-samples 12             # limit to 12 samples
#   bash run_stats.sh --wise set tier caller       # specific wises only
#   bash run_stats.sh --no-pileup                  # skip BAM pileup (faster)
#   bash run_stats.sh --no-bam                     # skip all BAM processing
#
#   # WES coverage (use BED region total as coverage denominator):
#   bash run_stats.sh --bed /path/to/exome_targets.bed
#
#   # Resume from existing parquet files (skip variant parsing + BAM stats):
#   bash run_stats.sh --resume
#
#   # Combined:
#   bash run_stats.sh --set 1 --wise set tier --no-pileup
#   bash run_stats.sh --bed /path/to/exome.bed --resume
#
# All extra arguments are forwarded to the Python CLI.

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MANIFEST="${HERE}/data/processed/sample_manifest.parquet"
OUTDIR="${HERE}/stats/full"

# Use rnadnavar's virtual environment
VENV="${HERE}/../../.venv"

cd "${HERE}/../.."  # rnadnavar repo root for PYTHONPATH=bin

PYTHONPATH="bin" "${VENV}/bin/python" -m vcf_stats.seq2neo.cli \
    --manifest "$MANIFEST" \
    --output-dir "$OUTDIR" \
    --sample-workers 4 \
    --bam-workers 8 \
    "$@"
