#!/usr/bin/env bash
set -euo pipefail

# Keep this launcher in the hybrid directory. Nextflow's task cache is keyed
# by this launch/work directory and the existing .nextflow database.
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "${HERE}/../../.." && pwd)"

export NXF_CONDA_CACHEDIR="${NXF_CONDA_CACHEDIR:-/t9k/mnt/joey/nf_conda_envs}"

exec micromamba run -n nextflow nextflow run "${REPO}/main.nf" \
  -c "${HERE}/seqc2.hybrid.realignment.config" \
  --input "${HERE}/csv/seqc2_wes_ll_hybrid.csv" \
  --outdir "${HERE}/output/seqc2.wes.ll.hybrid.realign.latest" \
  -offline -resume -with-conda "$@"
