#!/usr/bin/env bash
# Run the SEQC2 WES_LL tumor-normal pair (BAM input, DNA-only, consensus).
# Usage: bash run_wes_ll.sh [OUTDIR]
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/.." && pwd)"

INPUT="${INPUT:-$HERE/csv/seqc2_wes_ll.csv}"
OUTDIR="${1:-${OUTDIR:-$HERE/output/seqc2.wes.ll}}"

NXF_CONDA_CACHEDIR="${NXF_CONDA_CACHEDIR:-/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/nf_conda_envs}"
export NXF_CONDA_CACHEDIR NXF_CONDA_USEMAMBA=true

micromamba run -n nextflow nextflow run \
    "$REPO/main.nf" \
    -c "$HERE/seqc2.shared.config" \
    --input  "$INPUT" \
    --outdir "$OUTDIR" \
    -offline -with-conda -resume
