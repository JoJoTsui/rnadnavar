#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
HG008_CONFIG="${HG008_CONFIG:-${HERE}/config_hg008_wgs.yaml}"
HG008_REF="${HG008_FASTA:-/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/bio_db/references/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta}"
if [ ! -f "$HG008_REF" ]; then
  echo "ERROR: HG008 reference FASTA not found: $HG008_REF" >&2
  echo "Set HG008_FASTA to the GIAB GRCh38 reference before starting the pipeline." >&2
  exit 2
fi
HG008_FAI="${HG008_FASTA_FAI:-${HG008_REF}.fai}"
if [ ! -f "$HG008_FAI" ]; then
  echo "ERROR: HG008 reference index not found: $HG008_FAI" >&2
  echo "Set HG008_FASTA_FAI to the matching .fai before starting the pipeline." >&2
  exit 2
fi
exec python3 "${HERE}/../scripts/run_pipeline.py" --config "$HG008_CONFIG" \
  --fasta "$HG008_REF" --fasta-fai "$HG008_FAI" "$@"
