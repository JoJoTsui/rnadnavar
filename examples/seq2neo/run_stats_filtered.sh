#!/usr/bin/env bash
# Run variant statistics pipeline excluding NoConsensus variants.
#
# Uses the unified filter pipeline to exclude NoConsensus variants from
# ALL downstream processing: statistics, visualizations, BAM pileup,
# rescue analytics, and threshold sweeps. (The old --pileup-mode filtered
# flag only affected BAM pileup — it has been replaced.)
#
# Usage:
#   bash run_stats_filtered.sh                    # full run, NoConsensus excluded
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

echo "=== seq2neo stats (exclude NoConsensus) ==="
echo "Manifest: ${MANIFEST}"
echo "Output:   ${OUTDIR}"
echo "Filter:   NoConsensus excluded from ALL output"
echo ""

PYTHONPATH="bin" "${VENV}/bin/python" -m vcf_stats.seq2neo.cli \
    --manifest "$MANIFEST" \
    --output-dir "$OUTDIR" \
    --exclude-filters NoConsensus \
    --sample-workers 4 \
    --bam-workers 8 \
    --threads 6 \
    "$@"
