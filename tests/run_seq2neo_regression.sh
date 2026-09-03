#!/usr/bin/env bash
# Bounded DN/DT/RT acceptance test that inherits the production seq2neo config.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUTDIR="${SEQ2NEO_REGRESSION_OUTDIR:-$ROOT/.artifacts/seq2neo-regression}"
WORKDIR="${SEQ2NEO_REGRESSION_WORKDIR:-$ROOT/.tmp/seq2neo-regression-work}"
FIXTURE_DIR="${SEQ2NEO_REGRESSION_FIXTURE_DIR:-$OUTDIR/fixtures}"
MINI_FASTA='/t9k/mnt/hdd/work/Vax/rnadnavar/test-datasets/reference/chr7_hg38/GRCh38.d1.vd1.chr7.mini.fa'
GERMLINE_SOURCE='/t9k/mnt/hdd/work/Vax/rnadnavar/test-datasets/reference/chr7_hg38/af-only-gnomad.hg38.chr7.mini2.vcf.gz'
GERMLINE_VCF="$FIXTURE_DIR/af-only-gnomad.hg38.chr7.ref-consistent.vcf.gz"

export NXF_CONDA_CACHEDIR="${NXF_CONDA_CACHEDIR:-/t9k/mnt/joey/nf_conda_envs}"
export NXF_CONDA_USE_MAMBA="${NXF_CONDA_USE_MAMBA:-true}"

mkdir -p "$OUTDIR" "$WORKDIR" "$FIXTURE_DIR"

# The source mini gnomAD records are already bounded to 39 Mb, but its VCF
# header retained the full chr7 length. Reheader a generated copy to the exact
# mini FASTA dictionary; never weaken GATK sequence validation.
if [[ ! -s "$GERMLINE_VCF" || ! -s "$GERMLINE_VCF.tbi" || "$GERMLINE_SOURCE" -nt "$GERMLINE_VCF" ]]; then
    ln -sfn "$MINI_FASTA" "$FIXTURE_DIR/reference.fa"
    samtools faidx "$FIXTURE_DIR/reference.fa"
    tmp_vcf="$GERMLINE_VCF.$$.tmp.vcf.gz"
    bcftools reheader --fai "$FIXTURE_DIR/reference.fa.fai" --output "$tmp_vcf" "$GERMLINE_SOURCE"
    tabix --force --preset vcf "$tmp_vcf"
    mv "$tmp_vcf" "$GERMLINE_VCF"
    mv "$tmp_vcf.tbi" "$GERMLINE_VCF.tbi"
fi

micromamba run -n nextflow nextflow run "$ROOT/main.nf" \
    -c "$ROOT/tests/config/seq2neo_regression.config" \
    --outdir "$OUTDIR" \
    --germline_resource "$GERMLINE_VCF" \
    --germline_resource_tbi "$GERMLINE_VCF.tbi" \
    -work-dir "$WORKDIR" \
    -offline -resume -with-conda

TRACE="$(find "$OUTDIR/pipeline_info" -maxdepth 1 -name 'execution_trace*.txt' -type f -printf '%T@ %p\n' | sort -nr | head -1 | cut -d' ' -f2-)"
if [[ -z "$TRACE" ]]; then
    echo "ERROR: no execution trace found under $OUTDIR/pipeline_info" >&2
    exit 1
fi

"$ROOT/.venv/bin/python" "$ROOT/tests/seq2neo/validate_regression.py" \
    --trace "$TRACE" \
    --outdir "$OUTDIR"
