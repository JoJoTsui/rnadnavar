# Neoantigen Workflow — Usage Guide

## Overview

The neoantigen preparation workflow is an opt-in branch of the rnadnavar pipeline that:

1. **Harmonizes FORMAT fields** across DNA variant callers (Mutect2, Strelka2, SAGE, DeepSomatic) to a unified Mutect2-compatible schema, and enriches the consensus VCF with per-caller `AD_BY_CALLER` and `AF_BY_CALLER` INFO fields.
2. **Runs Salmon RNA quantification** (`salmon quant` v1.11.4) in quasi-mapping mode on RNA tumor FASTQs (status=2 samples) in parallel with STAR alignment.
3. **Publishes a neoantigen-ready VCF** for each DNA tumor sample (status=1), ready for use with pVACseq, seq2neo, or similar tools.

All new code paths are gated behind `enable_neoantigen_workflow = true`. Existing `examples/seq2neo/` configurations are unaffected.

---

## Shared infrastructure

All paths in `neoantigen.shared.config` and `run.sh` point to shared locations on this cluster, accessible to all teammates:

| Resource | Shared path |
|----------|-------------|
| Pipeline repo | `/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/rnadnavar/` |
| Reference databases | `/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/bio_db/` |
| Test dataset (COO8801) | `/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/work/rnadnavar_test/C008801/` |
| Conda environments | `/t9k/mnt/joey/nf_conda_envs/` |
| Salmon index (Gencode v49) | `/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/bio_db/salmon/salmon_index_gencode_v49_salmon_v1.11` |

No path edits are needed to run the test dataset — everything resolves from the shared locations above.

---

## Quick start

```bash
# 1. Dry run — preview the nextflow command without executing
bash examples/neoantigen/run.sh --dry-run

# 2. Run with the shared test dataset (COO8801, small subset)
bash examples/neoantigen/run.sh

# 3. Run with your own sample
bash examples/neoantigen/run.sh \
    --input /path/to/my_sample.csv \
    --outdir /path/to/output
```

`run.sh` resolves `MAIN_NF` and `RDV_CONF` from its own location in the repo (via `$BASH_SOURCE`), so it works correctly from any working directory as long as the repo is at the shared path.

---

## Files in this directory

```
examples/neoantigen/
├── neoantigen.shared.config   # full pipeline config (references + tools + neoantigen params)
├── run.sh                     # single-sample launcher (debug / tutorial)
└── README.usage.md            # this file
```

`neoantigen.shared.config` contains all standard pipeline params (genome, tools, reference files) plus the neoantigen-specific params. All reference paths point to the shared `bio_db`. The only param you may need to change is `salmon_index` if you want to use a different index than the shared one.

---

## Params to change per run

Most params are pre-configured for the shared cluster. The ones you may need to adjust:

| Param | Where | When to change |
|-------|-------|----------------|
| `salmon_index` | `neoantigen.shared.config` | Only if you need a different Gencode version or custom transcriptome. The shared index at `bio_db/salmon/` is ready to use. |
| `neoantigen_input_source` | `neoantigen.shared.config` | Change to `'consensus'` only if your downstream tool requires the multi-caller union VCF (see [VCF source selection](#vcf-source-selection-mutect2-vs-consensus) below). |
| `tools` | `neoantigen.shared.config` | Remove tools you don't need (e.g., remove `realignment` for a faster debug run). |
| `resourceLimits` | `neoantigen.shared.config` | Adjust `cpus`/`memory`/`time` for your compute node. |
| `--input` | `run.sh` CLI flag | Your own sample CSV. Default is the shared COO8801 test dataset. |
| `--outdir` | `run.sh` CLI flag | Output directory. Default is `output/COO8801.neoantigen` (relative). |

---

## VCF source selection: `mutect2` vs `consensus`

### Why `mutect2` is the recommended default

For neoantigen prediction tools (pVACseq, seq2neo), the Mutect2-filtered VCF is the correct input:

| Property | Mutect2 VCF | Consensus VCF |
|----------|-------------|---------------|
| Per-sample FORMAT | ✓ GT/AD/AF/DP/GQ/F1R2/F2R1/SB | ✗ FORMAT column is `.` (empty) |
| Variant count | ~91 (filtered somatic) | ~26,000+ (unfiltered union of all callers) |
| pVACseq/seq2neo compatible | ✓ | ✗ (no FORMAT data) |
| Filter status | Post-filtered (PASS + soft-filtered) | All variants from all callers |

The consensus VCF is a multi-caller aggregation designed for data labeling — all variant evidence lives in INFO fields (`VAF_BY_CALLER`, `DP_BY_CALLER`, etc.) rather than per-sample FORMAT. This is intentional for the pipeline's primary purpose, but incompatible with neoantigen tools that expect per-sample genotype data.

### Caller FORMAT field differences

| Field | Mutect2 | DeepSomatic | Strelka2 |
|-------|---------|-------------|----------|
| `GT` | ✓ | ✓ | ✗ |
| `AD` (ref,alt) | ✓ `Number=R` | ✓ `Number=R` | ✗ (uses `AU`/`CU`/`GU`/`TU` or `TAR`/`TIR`) |
| `AF` | ✓ (named `AF`) | ✗ (named `VAF`) | ✗ (must compute from base counts) |
| `DP` | ✓ | ✓ | ✓ |
| `GQ` | ✓ | ✓ | ✗ |
| `F1R2`/`F2R1`/`SB` | ✓ | ✗ | ✗ |

`FORMAT_HARMONIZER` normalizes Strelka2 and DeepSomatic FORMAT fields to Mutect2 conventions before the consensus step. However, the consensus VCF itself still has no per-sample FORMAT data by design — it cannot be aligned to Mutect2 FORMAT because it represents a union across callers, not a single caller's genotype call.

### When to use `consensus`

Use `neoantigen_input_source = 'consensus'` only if your downstream tool can consume the INFO-field-based format and you specifically want the multi-caller union. The consensus VCF with `--neoantigen` adds `AD_BY_CALLER` and `AF_BY_CALLER` INFO fields encoding per-caller allele depths and frequencies.

---

## Parameters

### Neoantigen-specific options

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `enable_neoantigen_workflow` | boolean | `false` | Activates the neoantigen preparation branch. |
| `neoantigen_input_source` | string | `'mutect2'` | VCF source: `mutect2` (recommended) or `consensus`. |
| `salmon_index` | string | shared index | Absolute path to a pre-built Salmon v1.11.x SSHash index directory. |
| `salmon_libtype` | string | `'A'` | Library type passed to `salmon quant --libType`. `'A'` = automatic detection. |
| `salmon_gc_bias` | boolean | `false` | Pass `--gcBias` to `salmon quant` for GC bias correction (~2 min extra/sample). |

### Pre-flight validation

The pipeline validates before launching any workflow steps:

- `enable_neoantigen_workflow = true` and `salmon_index` not set → error and exit.
- `salmon_index` path does not exist or is not a directory → error and exit.
- `neoantigen_input_source` not `consensus` or `mutect2` → error and exit.

---

## Salmon index

A pre-built Gencode v49 index is available at the shared path:

```
/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/bio_db/salmon/salmon_index_gencode_v49_salmon_v1.11
```

This is already set as the default in `neoantigen.shared.config`. No action needed unless you want a different transcriptome.

To build your own index:

```bash
bash /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/rnadnavar/scripts/build_salmon_index.sh \
    --outdir /path/to/output_dir
```

### Notes

- **Salmon v1.11.x index incompatibility**: v1.11.x uses a new SSHash-based k-mer index. Indices built with v1.10.x or earlier are **not compatible** and must be rebuilt.
- The `--gencode` flag is required when indexing Gencode transcriptomes (pipe-delimited FASTA headers). The pipeline's `QUANT_TSV_NORMALIZE` step strips these suffixes to produce `quant.tsv` for downstream tools.
- Gencode v49: GRCh38.p14, Ensembl 115, released September 2025.

---

## Expected outputs

```
${outdir}/neoantigen/
├── <DNA_tumor_id>/                          # status=1 sample (e.g. COO8801DT_vs_COO8801DN)
│   ├── *.filtered.vcf.gz                   # neoantigen-ready VCF (Mutect2-filtered by default)
│   └── *.filtered.vcf.gz.tbi
│
└── <RNA_tumor_id>/                          # status=2 sample (e.g. COO8801RT-LX)
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
--input  CSV     Input samplesheet CSV (default: shared COO8801 test dataset)
--outdir DIR     Output directory (default: output/COO8801.neoantigen)
--main-nf PATH   Path to main.nf (default: resolved from repo via $BASH_SOURCE)
--conf   PATH    Path to config file (default: neoantigen.shared.config in same dir)
--dry-run        Print the nextflow command without executing
-h, --help       Show help
```

---

## Manual nextflow invocation

```bash
HTTPS_PROXY="http://10.233.17.241:3128" \
NXF_CONDA_CACHEDIR="/t9k/mnt/joey/nf_conda_envs" \
NXF_CONDA_USEMAMBA=true \
micromamba run -n nextflow nextflow run \
    /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/rnadnavar/main.nf \
    -c /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/rnadnavar/examples/neoantigen/neoantigen.shared.config \
    --input  /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/work/rnadnavar_test/C008801/input/test.rdv.shared.csv \
    --outdir output/COO8801.neoantigen \
    -offline -with-conda -resume
```
