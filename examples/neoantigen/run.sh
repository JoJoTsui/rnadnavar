#!/usr/bin/env bash
# ==============================================================================
# run.sh — single-sample neoantigen workflow launcher
#
# Usage:
#   bash examples/neoantigen/run.sh [--input CSV] [--outdir DIR] [--dry-run]
#
# Defaults (edit the variables below or override via CLI flags):
#   INPUT   — shared test dataset CSV (accessible to all teammates)
#   OUTDIR  — output/COO8801.neoantigen
#
# Examples:
#   # Run with defaults (test dataset)
#   bash examples/neoantigen/run.sh
#
#   # Run with a custom sample
#   bash examples/neoantigen/run.sh \
#       --input /path/to/my_sample.csv \
#       --outdir /path/to/output
#
#   # Preview the nextflow command without executing
#   bash examples/neoantigen/run.sh --dry-run
# ==============================================================================
set -euo pipefail

# ── Self-locate (script works from any working directory) ─────────────────────
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "${HERE}/../.." && pwd)"

# ── Defaults ──────────────────────────────────────────────────────────────────
MAIN_NF="${REPO}/main.nf"
RDV_CONF="${HERE}/neoantigen.shared.config"
# -- shared rnadnavar test path
RDV_TEST_DIR="/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/work/rnadnavar_test"
RDV_TEST_INPUT_DIR="${RDV_TEST_DIR}/C008801/input"
RDV_TEST_INPUT="${RDV_TEST_INPUT_DIR}/test.rdv.shared.csv"

# Test dataset (small, for debugging and tutorial)
PATIENT_ID="COO8801"
INPUT="${RDV_TEST_INPUT}"
OUTDIR="output/${PATIENT_ID}.neoantigen"

# CONDA environment settings
NXF_CONDA_CACHEDIR="/t9k/mnt/joey/nf_conda_envs"
MICROMAMBA_ENV="nextflow"
HTTPS_PROXY="http://10.233.17.241:3128"
DRY_RUN=false

# ── Argument parsing ──────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --input)    INPUT="$2";   shift ;;
        --outdir)   OUTDIR="$2";  shift ;;
        --main-nf)  MAIN_NF="$2"; shift ;;
        --conf)     RDV_CONF="$2"; shift ;;
        --dry-run)  DRY_RUN=true ;;
        -h|--help)
            sed -n '/^# Usage:/,/^# ==/p' "$0" | grep "^#" | sed 's/^# \{0,2\}//'
            exit 0 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
    shift
done

# ── Print resolved config ─────────────────────────────────────────────────────
echo "========================================"
echo " Neoantigen workflow — single sample run"
echo "========================================"
echo "  main.nf : ${MAIN_NF}"
echo "  config  : ${RDV_CONF}"
echo "  input   : ${INPUT}"
echo "  outdir  : ${OUTDIR}"
echo "  dry-run : ${DRY_RUN}"
echo "========================================"

if [[ "${DRY_RUN}" == true ]]; then
    echo ""
    echo "Command that would run:"
    echo ""
    echo "  HTTPS_PROXY=\"${HTTPS_PROXY}\" \\"
    echo "  NXF_CONDA_CACHEDIR=\"${NXF_CONDA_CACHEDIR}\" \\"
    echo "  NXF_CONDA_USEMAMBA=true \\"
    echo "  micromamba run -n ${MICROMAMBA_ENV} nextflow run \\"
    echo "      ${MAIN_NF} \\"
    echo "      -c ${RDV_CONF} \\"
    echo "      --input  ${INPUT} \\"
    echo "      --outdir ${OUTDIR} \\"
    echo "      -offline -with-conda -resume"
    echo ""
    echo "Dry run complete. Nothing was executed."
    exit 0
fi

# ── Run ───────────────────────────────────────────────────────────────────────
HTTPS_PROXY="${HTTPS_PROXY}" \
NXF_CONDA_CACHEDIR="${NXF_CONDA_CACHEDIR}" \
NXF_CONDA_USEMAMBA=true \
micromamba run -n "${MICROMAMBA_ENV}" nextflow run \
    "${MAIN_NF}" \
    -c "${RDV_CONF}" \
    --input  "${INPUT}" \
    --outdir "${OUTDIR}" \
    -offline -with-conda -resume
