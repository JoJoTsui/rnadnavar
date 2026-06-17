#!/usr/bin/env bash
# Run variant statistics pipeline with filtered BAM pileup (excludes NoConsensus).
#
# This runs the FULL pipeline (not --resume) including Strelka TIR fix,
# but skips BAM pileup for NoConsensus variants (~98M positions) to save
# hours of runtime. NoConsensus variants are still included in all
# statistics and visualizations — only the BAM-level pileup is skipped.
#
# Usage:
#   bash run_stats_filtered.sh                    # full run, filtered pileup
#   bash run_stats_filtered.sh --resume           # resume from existing parquet
#   bash run_stats_filtered.sh --max-samples 12   # limit samples
#   bash run_stats_filtered.sh --no-bam           # skip all BAM processing
#
# All extra arguments are forwarded to the Python CLI.

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MANIFEST="${HERE}/data/processed/sample_manifest.tsv"
OUTDIR="${HERE}/stats/full_filtered"

# Use rnadnavar's virtual environment
VENV="${HERE}/../../.venv"

cd "${HERE}/../.."  # rnadnavar repo root for PYTHONPATH=bin

echo "=== seq2neo stats (filtered pileup) ==="
echo "Manifest: ${MANIFEST}"
echo "Output:   ${OUTDIR}"
echo "Pileup:   filtered (NoConsensus excluded from BAM pileup)"
echo ""

PYTHONPATH="bin" "${VENV}/bin/python" -m vcf_stats.seq2neo.cli \
    --manifest "$MANIFEST" \
    --output-dir "$OUTDIR" \
    --pileup-mode filtered \
    --sample-workers 4 \
    --bam-workers 8 \
    --threads 6 \
    "$@"
