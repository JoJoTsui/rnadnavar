# Neoantigen Workflow — Usage Guide

## Overview

The neoantigen preparation workflow is an opt-in branch of the rnadnavar pipeline that:

1. **Harmonizes FORMAT fields** across DNA variant callers (Mutect2, Strelka2, SAGE, DeepSomatic) to a unified Mutect2-compatible schema, and enriches the consensus VCF with per-caller `AD_BY_CALLER` and `AF_BY_CALLER` INFO fields.
2. **Runs Salmon RNA quantification** (`salmon quant` v1.11.4) in quasi-mapping mode on RNA tumor FASTQs (status=2 samples) in parallel with STAR alignment.
3. **Publishes a neoantigen-ready VCF** for each DNA tumor sample (status=1), ready for use with pVACseq, seq2neo, or similar tools.

All new code paths are gated behind `enable_neoantigen_workflow = true`. Existing `examples/seq2neo/` configurations are unaffected.

---

## Quick start

```bash
# 1. Set your salmon_index path in neoantigen.shared.config
#    (see "Building the Salmon index" below)

# 2. Preview the command (dry run)
bash examples/neoantigen/run.sh --dry-run

# 3. Run with the test dataset
bash examples/neoantigen/run.sh

# 4. Run with your own sample
bash examples/neoantigen/run.sh \
    --input /path/to/my_sample.csv \
    --outdir /path/to/output
```

The `run.sh` script defaults to the small test dataset at
`/t9k/mnt/hdd/work/Vax/sequencing/aim_exp/rdv_test/C008801/input/test.rdv.shared.csv`,
which is the same dataset used for pipeline debugging and CI.

---

## Files in this directory

```
examples/neoantigen/
├── neoantigen.shared.config   # full pipeline config (references + tools + neoantigen params)
├── run.sh                     # single-sample launcher (debug / tutorial)
└── README.usage.md            # this file
```

`neoantigen.shared.config` contains all standard pipeline params (genome, tools, reference
files) inherited from `seq2neo.shared.config`, plus the neoantigen-specific params on top.
The only value you need to change before first use is `salmon_index`.

---

## VCF source selection: `mutect2` vs `consensus`

### Why `mutect2` is the recommended default

For neoantigen prediction tools (pVACseq, seq2neo), the Mutect2-filtered VCF is the correct input:

| Property | Mutect2 VCF | Consensus VCF |
|----------|-------------|---------------|
| Per-sample FORMAT | ✓ GT/AD/AF/DP/GQ/F1R2/F2R1/SB | ✗ FORMAT column is `.` (empty) |
| Variant count | ~91 (filtered somatic) | ~26,000+ (unfiltered union of all callers) |
| pVACseq compatible | ✓ | ✗ (no FORMAT data) |
| Filter status | Post-filtered (PASS + soft-filtered) | All variants from all callers |

The consensus VCF is a multi-caller aggregation designed for data labeling — all variant evidence lives in INFO fields (`VAF_BY_CALLER`, `DP_BY_CALLER`, etc.) rather than per-sample FORMAT. This is intentional and correct for the pipeline's primary purpose, but incompatible with neoantigen tools that expect per-sample genotype data.

### Caller FORMAT field differences

| Field | Mutect2 | DeepSomatic | Strelka2 |
|-------|---------|-------------|----------|
| `GT` | ✓ | ✓ | ✗ |
| `AD` (ref,alt) | ✓ Number=R | ✓ Number=R | ✗ (uses AU/CU/GU/TU or TAR/TIR) |
| `AF` | ✓ (named `AF`) | ✗ (named `VAF`) | ✗ (must compute from base counts) |
| `DP` | ✓ | ✓ | ✓ |
| `GQ` | ✓ | ✓ | ✗ |
| `F1R2/F2R1/SB` | ✓ | ✗ | ✗ |

This is why `FORMAT_HARMONIZER` exists — it normalizes Strelka2 and DeepSomatic FORMAT fields to Mutect2 conventions before the consensus step. However, the consensus VCF itself still has no per-sample FORMAT data by design.

### When to use `consensus`

Use `neoantigen_input_source = 'consensus'` only if your downstream tool can consume the INFO-field-based format and you specifically want the multi-caller union. The consensus VCF with `--neoantigen` adds `AD_BY_CALLER` and `AF_BY_CALLER` INFO fields that encode per-caller allele depths and frequencies.

### Neoantigen-specific options

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `enable_neoantigen_workflow` | boolean | `false` | Activates the neoantigen preparation branch. |
| `neoantigen_input_source` | string | `'consensus'` | VCF source for neoantigen output: `consensus` or `mutect2`. |
| `salmon_index` | string | `null` | **Required.** Absolute path to a pre-built Salmon v1.11.x SSHash index directory. |
| `salmon_libtype` | string | `'A'` | Library type passed to `salmon quant --libType`. `'A'` = automatic detection. |
| `salmon_gc_bias` | boolean | `false` | Pass `--gcBias` to `salmon quant` for GC bias correction. |

All other params (genome, tools, reference files, etc.) are the same as `seq2neo.shared.config`.

### Pre-flight validation

The pipeline validates before launching any workflow steps:

- `enable_neoantigen_workflow = true` and `salmon_index` not set → error and exit.
- `salmon_index` path does not exist or is not a directory → error and exit.
- `neoantigen_input_source` not `consensus` or `mutect2` → error and exit.

---

## Building the Salmon index

Use the provided script to download Gencode v49 references and build a decoy-aware Salmon v1.11.x index:

```bash
bash /path/to/rnadnavar/scripts/build_salmon_index.sh --outdir /path/to/salmon_indices
```

Then set `salmon_index` in `neoantigen.shared.config`:

```groovy
salmon_index = '/path/to/salmon_indices/salmon_index_gencode_v49_salmon_v1.11'
```

### Notes

- The `--gencode` flag is **required** when indexing Gencode transcriptomes. Gencode FASTA headers use pipe-delimited format (e.g., `ENST00000456328.2|ENSG00000223972.5|...`). The pipeline's `QUANT_TSV_NORMALIZE` step strips these suffixes to produce `quant.tsv` for downstream tools.
- **Salmon v1.11.x index incompatibility**: v1.11.x uses a new SSHash-based k-mer index. Indices built with v1.10.x or earlier are **not compatible** and must be rebuilt.
- Gencode v49: GRCh38.p14, Ensembl 115, released September 2025.

---

## Expected outputs

```
${outdir}/<sample_id>/
├── neoantigen/                              # DNA tumor samples (status=1)
│   ├── <sample_id>.neoantigen.vcf.gz
│   └── <sample_id>.neoantigen.vcf.gz.tbi
│
└── salmon/                                  # RNA tumor samples (status=2)
    ├── quant.sf                             # raw Salmon output (Gencode pipe-delimited names)
    ├── quant.tsv                            # normalized transcript IDs (plain ENST IDs)
    ├── lib_format_counts.json
    ├── cmd_info.json
    └── aux_info/
        ├── meta_info.json
        ├── ambig_info.tsv
        └── fld.gz
```

`quant.tsv` is the file consumed by downstream neoantigen tools (pVACseq, seq2neo).

---

## run.sh options

```
--input  CSV     Input samplesheet CSV (default: test dataset)
--outdir DIR     Output directory (default: output/COO8801.neoantigen)
--main-nf PATH   Path to main.nf (default: auto-detected from repo root)
--conf   PATH    Path to config file (default: neoantigen.shared.config)
--dry-run        Print the nextflow command without executing
-h, --help       Show help
```

---

## Manual nextflow invocation

```bash
REPO=/path/to/rnadnavar

HTTPS_PROXY="http://10.233.17.241:3128" \
NXF_CONDA_CACHEDIR="/t9k/mnt/joey/nf_conda_envs" \
NXF_CONDA_USEMAMBA=true \
micromamba run -n nextflow nextflow run \
    $REPO/main.nf \
    -c $REPO/examples/neoantigen/neoantigen.shared.config \
    --input  /path/to/sample.csv \
    --outdir /path/to/output \
    -offline -with-conda -resume
```
