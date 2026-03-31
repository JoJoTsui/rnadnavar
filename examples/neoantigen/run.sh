#!/usr/bin/env bash
# ==============================================================================
# run.sh — single-sample neoantigen workflow launcher
#
# Usage:
#   bash examples/neoantigen/run.sh [--input CSV] [--outdir DIR] [--dry-run]
#
# Shared infrastructure (accessible to all teammates on this cluster):
#   Pipeline repo : /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/rnadnavar/
#   Reference DBs : /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/bio_db/
#   Test dataset  : /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/work/rnadnavar_test/
#   Conda envs    : /t9k/mnt/joey/nf_conda_envs/
#
# Defaults (no edits needed for the test dataset):
#   MAIN_NF — resolved from this script's location inside the shared repo
#   RDV_CONF — neoantigen.shared.config in the same directory as this script
#   INPUT   — shared test dataset CSV (COO8801, small subset for debugging)
#   OUTDIR  — output/COO8801.neoantigen (relative to where you run the script)
#
# Params you may need to change:
#   salmon_index — set in neoantigen.shared.config; must point to a pre-built
#                  Salmon v1.11.x index (see "Building the Salmon index" in README)
#   OUTDIR       — override with --outdir for production runs
#   INPUT        — override with --input for your own sample CSV
#
# Examples:
#   # Run with defaults (test dataset, dry run first)
#   bash examples/neoantigen/run.sh --dry-run
#   bash examples/neoantigen/run.sh
#
#   # Run with a custom sample
#   bash examples/neoantigen/run.sh \
#       --input /path/to/my_sample.csv \
#       --outdir /path/to/output
# ==============================================================================
set -euo pipefail

# ── Self-locate: MAIN_NF and RDV_CONF resolve from the repo, not the caller ──
# This means the script works correctly regardless of which directory you run
# it from, as long as the repo is checked out at the shared path below.
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "${HERE}/../.." && pwd)"

# ── Shared paths (same across all teammates on this cluster) ──────────────────
# Pipeline
MAIN_NF="${REPO}/main.nf"
RDV_CONF="${HERE}/neoantigen.shared.config"
# -- shared rnadnavar test path
RDV_TEST_DIR="/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/work/rnadnavar_test"
RDV_TEST_INPUT_DIR="${RDV_TEST_DIR}/C008801/input"
RDV_TEST_INPUT="${RDV_TEST_INPUT_DIR}/test.rdv.shared.csv"

# Test dataset — small COO8801 subset, used for debugging and CI
PATIENT_ID="COO8801"
INPUT="${RDV_TEST_INPUT}"
OUTDIR="output/${PATIENT_ID}.neoantigen"

# Conda / Nextflow environment
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
