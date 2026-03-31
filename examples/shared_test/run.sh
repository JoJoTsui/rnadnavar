#!/usr/bin/env bash
# ==============================================================================
# run.sh — shared test dataset launcher (data labeling workflow)
#
# This is the canonical single-sample test script for validating a fresh
# rnadnavar installation. It uses the shared COO8801 test dataset and the
# shared pipeline config, so no path edits are needed on this cluster.
#
# Usage:
#   bash examples/shared_test/run.sh [--outdir DIR] [--dry-run]
#
# Shared infrastructure (accessible to all teammates):
#   Pipeline repo : /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/rnadnavar/
#   Reference DBs : /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/bio_db/
#   Test dataset  : /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/work/rnadnavar_test/C008801/
#   Conda envs    : /t9k/mnt/joey/nf_conda_envs/
#
# Expected outputs (installation validation):
#   ${OUTDIR}/variant_calling/mutect2/     — Mutect2 VCFs
#   ${OUTDIR}/variant_calling/strelka/     — Strelka2 VCFs
#   ${OUTDIR}/variant_calling/deepsomatic/ — DeepSomatic VCFs
#   ${OUTDIR}/consensus/                   — consensus VCF + MAF
#   ${OUTDIR}/vcf_realignment/             — realigned RNA VCFs
#   ${OUTDIR}/rescue/                      — cross-modality rescued VCFs
#   ${OUTDIR}/pipeline_info/execution_trace_*.txt
# ==============================================================================
set -euo pipefail

# ── Self-locate: MAIN_NF resolves from the repo, not the caller ──────────────
# Works correctly from any working directory as long as the repo is at the shared path.
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "${HERE}/../.." && pwd)"

# ── Paths ─────────────────────────────────────────────────────────────────────
MAIN_NF="${REPO}/main.nf"

# Shared test dataset — small COO8801 subset, used for debugging and CI
RDV_TEST_DIR="/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/work/rnadnavar_test"
RDV_CONF="${RDV_TEST_DIR}/C008801/input/test.rdv.shared.config"
INPUT="${RDV_TEST_DIR}/C008801/input/test.rdv.shared.csv"
OUTDIR="${RDV_TEST_DIR}/output/COO8801.shared"

NXF_CONDA_CACHEDIR="/t9k/mnt/joey/nf_conda_envs"
MICROMAMBA_ENV="nextflow"
HTTPS_PROXY="http://10.233.17.241:3128"
DRY_RUN=false

# ── Argument parsing ──────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --outdir)  OUTDIR="$2"; shift ;;
        --dry-run) DRY_RUN=true ;;
        -h|--help)
            sed -n '/^# Usage:/,/^# ==/p' "$0" | grep "^#" | sed 's/^# \{0,2\}//'
            exit 0 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
    shift
done

echo "========================================"
echo " rnadnavar — shared test run"
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

HTTPS_PROXY="${HTTPS_PROXY}" \
NXF_CONDA_CACHEDIR="${NXF_CONDA_CACHEDIR}" \
NXF_CONDA_USEMAMBA=true \
micromamba run -n "${MICROMAMBA_ENV}" nextflow run \
    "${MAIN_NF}" \
    -c "${RDV_CONF}" \
    --input  "${INPUT}" \
    --outdir "${OUTDIR}" \
    -offline -with-conda -resume
