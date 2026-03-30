# Neoantigen Workflow — Usage Guide

## Overview

The neoantigen preparation workflow is an opt-in branch of the rnadnavar pipeline that:

1. **Harmonizes FORMAT fields** across DNA variant callers (Mutect2, Strelka2, SAGE, DeepSomatic) to a unified Mutect2-compatible schema, and enriches the consensus VCF with per-caller `AD_BY_CALLER` and `AF_BY_CALLER` INFO fields.
2. **Runs Salmon RNA quantification** (`salmon quant` v1.11.4) in quasi-mapping mode on RNA tumor FASTQs (status=2 samples) in parallel with STAR alignment.
3. **Publishes a neoantigen-ready VCF** for each DNA tumor sample (status=1), ready for use with pVACseq, seq2neo, or similar tools.

All new code paths are gated behind `enable_neoantigen_workflow = true`. Existing `examples/seq2neo/` configurations are unaffected.

---

## Quick start

Before first use, set three values in `config/runner.yaml`:

| Key | What to set |
|-----|-------------|
| `main_nf` | absolute path to `rnadnavar/main.nf` on this machine |
| `rdv_conf` | absolute path to `examples/neoantigen/neoantigen.shared.config` on this machine |
| `seq2neo_root` | where data/output/state will be written |

Then set `salmon_index` in `neoantigen.shared.config` to the path of your pre-built Salmon v1.11.x index (see [Building the Salmon index](#building-the-salmon-index) below).

---

## Parameters

### Neoantigen workflow options

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `enable_neoantigen_workflow` | boolean | `false` | Activates the neoantigen preparation branch. When `false`, the pipeline runs identically to the standard data labeling workflow. |
| `neoantigen_input_source` | string | `'consensus'` | VCF source for the neoantigen output. Accepted values: `consensus` (multi-caller consensus VCF) or `mutect2` (Mutect2-filtered VCF). |
| `salmon_index` | string | `null` | **Required when `enable_neoantigen_workflow = true`.** Absolute path to a pre-built Salmon v1.11.x SSHash index directory. |
| `salmon_libtype` | string | `'A'` | Library type passed to `salmon quant --libType`. `'A'` enables automatic detection. |
| `salmon_gc_bias` | boolean | `false` | When `true`, passes `--gcBias` to `salmon quant` to enable fragment-level GC bias correction. Safe to enable; adds a few minutes per sample. |

### Pre-flight validation

The pipeline validates parameters before launching any workflow steps:

- If `enable_neoantigen_workflow = true` and `salmon_index` is not set → error and exit.
- If `salmon_index` path does not exist or is not a directory → error and exit.
- If `neoantigen_input_source` is not `consensus` or `mutect2` → error and exit.

---

## Required reference files

### Salmon v1.11.x index

A pre-built Salmon v1.11.x SSHash-based index directory is required. This index **cannot** be shared with Salmon v1.10.x or earlier — the index format changed in v1.11.x (SSHash-based k-mer index, replacing the prior pufferfish/ccDBG index). Indices built with v1.10.x must be rebuilt.

See [Building the Salmon index](#building-the-salmon-index) below.

---

## Building the Salmon index

Use the provided script to download Gencode v49 references and build a decoy-aware Salmon v1.11.x index:

```bash
bash /path/to/rnadnavar/scripts/build_salmon_index.sh
```

The script:
- Downloads `gencode.v49.transcripts.fa.gz` and `GRCh38.primary_assembly.genome.fa.gz` from the GENCODE v49 EBI FTP via `aria2c`
- Extracts decoy sequence IDs from the genome FASTA
- Concatenates transcriptome + genome into a gentrome file
- Runs `salmon index` with the `--gencode` flag

### Notes on the `--gencode` flag

The `--gencode` flag is **required** when indexing Gencode transcriptomes. Gencode FASTA headers use a pipe-delimited format (e.g., `ENST00000456328.2|ENSG00000223972.5|...`). Without `--gencode`, Salmon cannot correctly parse these headers during index construction. This pipe-delimited format also means the pipeline's transcript name normalization step (which produces `quant.tsv` from `quant.sf`) is necessary for downstream tools.

### Gencode v49 details

| Field | Value |
|-------|-------|
| Release | v49 |
| Genome assembly | GRCh38.p14 |
| Ensembl release | 115 |
| Freeze date | February 2025 |
| Release date | September 2025 |

### Index version incompatibility

Salmon v1.11.x adopts a new SSHash-based k-mer index format. Indices built with v1.10.3 or earlier are **not compatible** with v1.11.x and must be rebuilt using `salmon index` from v1.11.x before use.

---

## Expected outputs

For each sample, the pipeline publishes to `${outdir}/<sample_id>/`:

### Neoantigen VCF (DNA tumor samples, status=1)

```
<outdir>/<sample_id>/neoantigen/
├── <sample_id>.neoantigen.vcf.gz      # bgzip-compressed neoantigen-ready VCF
└── <sample_id>.neoantigen.vcf.gz.tbi  # tabix index
```

### Salmon quantification (RNA tumor samples, status=2)

```
<outdir>/<sample_id>/salmon/
├── quant.sf                  # transcript-level quantification (Gencode pipe-delimited names)
├── quant.tsv                 # quant.sf with Name column normalized to plain transcript IDs
├── lib_format_counts.json    # inferred library format counts
├── cmd_info.json             # exact salmon quant command used (reproducibility)
└── aux_info/                 # auxiliary output directory
    ├── meta_info.json        # run metadata and mapping statistics
    ├── ambig_info.tsv        # ambiguous mapping information
    └── fld.gz                # fragment length distribution
```

`quant.tsv` is the post-processed version of `quant.sf` with the `Name` column normalized: everything after the first `|` is stripped (e.g., `ENST00000456328.2|ENSG00000223972.5|...` → `ENST00000456328.2`). This is the file consumed by downstream neoantigen tools.

---

## Manual nextflow invocation (single sample)

```bash
REPO=/path/to/rnadnavar
DATA=/path/to/neoantigen_data

NXF_CONDA_CACHEDIR="/path/to/nf_conda_envs" \
NXF_CONDA_USEMAMBA=true \
micromamba run -n nextflow nextflow run \
    $REPO/main.nf \
    -c $REPO/examples/neoantigen/neoantigen.shared.config \
    --input  $DATA/runs/csv/<sample>.csv \
    --outdir $DATA/output/<sample> \
    -offline -with-conda -resume
```
