#!/usr/bin/env bash
# Re-run the frozen comparisons into NEW directories. Never overwrites originals.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
DEST="${1:?Usage: bash reproduce.sh /absolute/new-output-root}"
[[ "$DEST" = /* && ! -e "$DEST" ]] || { echo 'Require a new absolute output root' >&2; exit 1; }
cd "$ROOT"
"$ROOT/.venv/bin/python" docs/manuscript/frozen_native_gate_20260916/verify.py
for dataset in seqc2_wes seqc2_wgs hg008; do
    "$ROOT/.venv/bin/python" examples/seqc2/scripts/validate_frozen_hybrid_policy.py \
        --manifest "examples/seqc2/hybrid/frozen_${dataset}_validation.json" \
        --outdir "$DEST/$dataset"
    for stage in consensus first realignment; do
        if [[ "$stage" = consensus ]]; then vcf="$DEST/$dataset/refined.vcf.gz";
        else vcf="$DEST/$dataset/$stage/refined.rescue.vcf.gz"; fi
        "$ROOT/.venv/bin/python" examples/seqc2/scripts/audit_refined_label_contract.py \
            --vcf "$vcf" --report "$DEST/$dataset/$stage.structural_audit_v1.json"
    done
done
