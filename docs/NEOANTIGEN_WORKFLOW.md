# Neoantigen Workflow — Developer & User Reference

## Overview

The neoantigen preparation workflow is an opt-in branch of the rnadnavar pipeline, gated behind `enable_neoantigen_workflow = false` by default. When enabled, it adds three coordinated pipeline branches without modifying any existing code paths:

1. **FORMAT_HARMONIZER** — normalizes VCF FORMAT fields across DNA variant callers to Mutect2 conventions; runs between `VCF_NORMALIZE` and `VCF_CONSENSUS`.
2. **Salmon quantification** — runs `salmon quant` v1.11.4 in quasi-mapping mode on RNA tumor FASTQs (status=2) in parallel with STAR alignment.
3. **VCF selection** — selects the consensus or Mutect2-filtered VCF for DNA tumor samples (status=1) and publishes it as the neoantigen-ready VCF.

All existing `examples/seq2neo/` configurations and the main data labeling workflow are unaffected when `enable_neoantigen_workflow` is `false`.

---

## Architecture

### Pipeline integration points

The neoantigen workflow integrates into `workflows/rnadnavar.nf` at two points:

**Point 1 — After `VCF_NORMALIZE`, before `VCF_CONSENSUS`**

When `enable_neoantigen_workflow = true`, the `FORMAT_HARMONIZER` process intercepts the normalized per-caller VCF channel and remaps FORMAT fields to Mutect2 conventions. The harmonized VCFs replace the normalized VCFs as input to `VCF_CONSENSUS`. The `--neoantigen` flag is conditionally appended to `run_consensus_vcf.py` via `task.ext.args` in `conf/modules/consensus/vcf_consensus.config`.

**Point 2 — In parallel with `BAM_ALIGN`**

The `NEOANTIGEN_WORKFLOW` subworkflow is invoked after `BAM_ALIGN` completes (FASTQs are available), running `SALMON_QUANT` concurrently with the STAR-based `BAM_PROCESSING` branch. VCF selection and publishing happen after `BAM_PROCESSING` completes.

### Workflow diagram (text description)

```
BAM_ALIGN
  │
  ├──► BAM_PROCESSING
  │       │
  │       └──► VCF_NORMALIZE
  │               │
  │               ├── [enable_neoantigen_workflow=false] ──► VCF_CONSENSUS
  │               │
  │               └── [enable_neoantigen_workflow=true]
  │                       │
  │                       └──► FORMAT_HARMONIZER ──► VCF_CONSENSUS
  │                                                       │
  │                                                       └──► [neoantigen_input_source]
  │                                                               ├── consensus → select status=1 → publish neoantigen/
  │                                                               └── mutect2  → select status=1 → publish neoantigen/
  │
  └──► [enable_neoantigen_workflow=true]
          │
          └──► SALMON_QUANT (status=2 FASTQs only)
                  │
                  └──► QUANT_TSV_NORMALIZE ──► publish salmon/
```

### Backward compatibility gate

All new code paths are wrapped in `if (params.enable_neoantigen_workflow)` blocks. When `false`:
- `FORMAT_HARMONIZER` is never invoked
- `SALMON_QUANT` is never invoked
- `run_consensus_vcf.py` receives no `--neoantigen` flag
- No `neoantigen/` or `salmon/` output directories are created

---

## Parameters

All new parameters are declared in `nextflow.config` under the `// Neoantigen workflow options` comment block and registered in `nextflow_schema.json` under the `"Neoantigen workflow options"` group.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `enable_neoantigen_workflow` | boolean | `false` | Activates the neoantigen preparation branch. When `false`, the pipeline runs identically to the pre-feature state. |
| `neoantigen_input_source` | string | `'mutect2'` | VCF source for neoantigen output: `consensus` or `mutect2`. See [VCF source selection](#vcf-source-selection-mutect2-vs-consensus) below. |
| `salmon_index` | string | `null` | Absolute path to a pre-built Salmon v1.11.x SSHash index directory. Required when `enable_neoantigen_workflow = true`. |
| `salmon_libtype` | string | `'A'` | Library type passed to `salmon quant --libType`. `'A'` enables automatic detection. |
| `salmon_gc_bias` | boolean | `false` | When `true`, passes `--gcBias` to `salmon quant` for fragment-level GC bias correction. |

### Pre-flight validation

Enforced in `workflows/rnadnavar.nf` before any subworkflow is invoked, and additionally by nf-schema via `nextflow_schema.json`:

```groovy
if (params.enable_neoantigen_workflow) {
    if (!params.salmon_index) {
        error "ERROR: --salmon_index is required when enable_neoantigen_workflow is true."
    }
    if (!file(params.salmon_index).isDirectory()) {
        error "ERROR: --salmon_index path '${params.salmon_index}' does not exist or is not a directory."
    }
    if (!['consensus', 'mutect2'].contains(params.neoantigen_input_source)) {
        error "ERROR: --neoantigen_input_source must be 'consensus' or 'mutect2', got '${params.neoantigen_input_source}'."
    }
}
```

---

## Pipeline branches

### Branch 1: FORMAT_HARMONIZER

**File**: `modules/local/format_harmonizer/main.nf`  
**Script**: `bin/harmonize_vcf_format.py`  
**Container**: `mulled-v2-629aec3ba267b06a1efc3ec454c0f09e134f6ee2` (same as `VCF_CONSENSUS`; includes `cyvcf2`)

**When it runs**: After `VCF_NORMALIZE`, before `VCF_CONSENSUS`, only when `enable_neoantigen_workflow = true`.

**What it does**: Reads each per-caller normalized VCF and remaps FORMAT fields to Mutect2 conventions. The caller is identified from `meta.variantcaller` passed via `--caller`. The harmonized VCF replaces the normalized VCF as input to `VCF_CONSENSUS`.

**Input**: `tuple val(meta), path(vcf), path(tbi)` — one normalized per-caller VCF  
**Output**: `tuple val(meta), path("*.harmonized.vcf.gz"), path("*.harmonized.vcf.gz.tbi")`

#### Caller-specific FORMAT remapping rules

The target schema is Mutect2 FORMAT conventions. The following table shows how each caller's fields are mapped to the canonical `GT`, `AD`, `AF`, `DP`, `GQ` fields:

| Canonical Field | Mutect2 | Strelka2 (SNV) | Strelka2 (Indel) | SAGE | DeepSomatic |
|---|---|---|---|---|---|
| `GT` | passthrough | passthrough | passthrough | passthrough | passthrough |
| `AD` (ref,alt) | passthrough | `{REF}U[0]`, `{ALT}U[0]` (base tier-0 counts) | `TAR[0]`, `TIR[0]` | passthrough | passthrough |
| `AF` (allele freq) | passthrough | alt/(ref+alt) from base counts | TIR[0]/(TAR[0]+TIR[0]) | `VF` → `AF` | `VAF` → `AF` |
| `DP` (total depth) | passthrough | sum of all base tier-0 counts | TAR[0]+TIR[0] | passthrough | passthrough |
| `GQ` (genotype qual) | passthrough | passthrough | passthrough | passthrough | passthrough |

**SNV vs Indel detection for Strelka2**: check whether `AU`/`CU`/`GU`/`TU` FORMAT fields are present (SNV) vs `TAR`/`TIR` (Indel).

**Header updates**: The script rewrites the VCF header to declare `AD` (Number=R, Type=Integer), `AF` (Number=A, Type=Float), `DP` (Number=1, Type=Integer), `GQ` (Number=1, Type=Integer) with Mutect2-compatible descriptions. All non-target FORMAT fields are preserved unchanged.

**Missing field handling**: If a required source field is absent and cannot be computed, the script writes `.` for that field and continues. It does not raise an error.

#### Per-caller FORMAT aggregation in consensus VCF

When `enable_neoantigen_workflow = true`, `run_consensus_vcf.py` is invoked with `--neoantigen`, which adds two additional INFO fields to the consensus VCF:

| Field | Format | Example |
|-------|--------|---------|
| `AD_BY_CALLER` | `caller:ref,alt\|caller:ref,alt\|...` | `mutect2:45,12\|strelka:43,11\|deepsomatic:44,12` |
| `AF_BY_CALLER` | `caller:value\|caller:value\|...` | `mutect2:0.2105\|strelka:0.2037\|deepsomatic:0.2143` |

When `enable_neoantigen_workflow = false`, these fields are never written. The existing `VAF_BY_CALLER` field is unchanged in both modes.

---

### Branch 2: Salmon quantification

**Files**:
- `modules/local/salmon_quant/main.nf`
- `modules/local/quant_tsv_normalize/main.nf`
- `conf/modules/salmon/salmon_quant.config`

**When it runs**: In parallel with `BAM_PROCESSING` (after `BAM_ALIGN`), only when `enable_neoantigen_workflow = true`. Only status=2 (RNA tumor) samples are processed.

#### SALMON_QUANT process

**Container/Conda**:
```groovy
conda "bioconda::salmon=1.11.4 conda-forge::libstdcxx-ng>=14"
container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/salmon:1.11.4--h43eeafb_0' :
    'quay.io/biocontainers/salmon:1.11.4--h43eeafb_0' }"
```

**Resource label**: `process_medium` (12 CPUs, 48 GB, 8 h)

**Command construction** (paired-end):
```bash
salmon quant \
    --index <index_dir> \
    --libType A \
    -1 <R1> -2 <R2> \
    -p ${task.cpus} \
    --validateMappings \
    -o <prefix> \
    ${args}
```

Important: `--libType` must appear **before** `-1`/`-2`/`-r` per Salmon's argument parser requirement.

For single-end reads, `-1 <R1> -2 <R2>` is replaced with `-r <R>`.

For multi-lane samples, all R1 files are passed as a space-separated list to `-1` and all R2 files to `-2` (Salmon natively supports multiple input files per read end).

#### QUANT_TSV_NORMALIZE process

Reads `quant.sf` and produces `quant.tsv` by stripping everything after the first `|` in the `Name` column (e.g., `ENST00000456328.2|ENSG00000223972.5|...` → `ENST00000456328.2`). Values without `|` are left unchanged. All other columns are preserved byte-for-byte.

Fails with a non-zero exit if `quant.sf` is empty or has fewer than 2 lines.

**Container**: Same mulled container as `VCF_CONSENSUS`  
**Resource label**: `process_low`

---

### Branch 3: VCF selection

**File**: `subworkflows/local/neoantigen_workflow/main.nf`

**When it runs**: After `BAM_PROCESSING` completes, only when `enable_neoantigen_workflow = true`.

**What it does**: Based on `params.neoantigen_input_source`, selects either the consensus VCF or the Mutect2-filtered VCF for status=1 (DNA tumor) samples and publishes it as the neoantigen-ready VCF.

- `neoantigen_input_source = 'consensus'` → selects the DNA_Consensus_VCF (output of `VCF_CONSENSUS`) for status=1 samples
- `neoantigen_input_source = 'mutect2'` → selects the Mutect2_Filtered_VCF for status=1 samples

If no status=1 samples are present, a warning is logged via `log.warn` and no neoantigen VCF is published.

---

## Output directory structure

```
${outdir}/<sample_id>/
├── neoantigen/                              # DNA tumor samples (status=1)
│   ├── <sample_id>.neoantigen.vcf.gz        # bgzip-compressed neoantigen-ready VCF
│   └── <sample_id>.neoantigen.vcf.gz.tbi   # tabix index
│
└── salmon/                                  # RNA tumor samples (status=2)
    ├── quant.sf                             # transcript-level quantification (Gencode pipe-delimited names)
    ├── quant.tsv                            # quant.sf with Name column normalized to plain transcript IDs
    ├── lib_format_counts.json               # inferred library format counts
    ├── cmd_info.json                        # exact salmon quant command (reproducibility)
    └── aux_info/                            # auxiliary output directory
        ├── meta_info.json                   # run metadata and mapping statistics
        ├── ambig_info.tsv                   # ambiguous mapping information
        └── fld.gz                           # fragment length distribution
```

---

## Salmon v1.11.4 version constraints

This workflow targets **Salmon v1.11.4** (released 2026-03-11). Key points:

- **Index format change (breaking)**: v1.11.x adopts a new SSHash-based k-mer index, replacing the prior pufferfish/ccDBG index. Indices built with v1.10.3 or earlier are **not compatible** and must be rebuilt with `salmon index` using v1.11.x.
- **`salmon quant` CLI is unchanged**: All flags (`--index`, `--libType`, `-1`/`-2`, `-r`, `-p`, `--validateMappings`, `-o`, etc.) are identical to v1.10.x.
- **Output formats unchanged**: `quant.sf`, `meta_info.json`, `lib_format_counts.json`, and `aux_info/` are compatible with v1.10.x output.
- **`salmon alevin` removed**: The single-cell subcommand has been permanently removed in v1.11.x. This workflow uses bulk RNA-seq only (`salmon quant`).

---

## Integration points with the existing pipeline

### Where FORMAT_HARMONIZER slots in

```
VCF_NORMALIZE
    │
    ├── [enable_neoantigen_workflow=false] ──────────────────────► VCF_CONSENSUS
    │
    └── [enable_neoantigen_workflow=true] ──► FORMAT_HARMONIZER ──► VCF_CONSENSUS
```

The `FORMAT_HARMONIZER` process is inserted into the existing normalization-to-consensus channel in `workflows/rnadnavar.nf`. When `enable_neoantigen_workflow = false`, the channel flows directly from `VCF_NORMALIZE` to `VCF_CONSENSUS` as before.

### Parallel execution with BAM_PROCESSING

`SALMON_QUANT` runs concurrently with the STAR-based `BAM_PROCESSING` branch. Both branches receive their inputs from `BAM_ALIGN` (FASTQs for Salmon, aligned BAMs for BAM_PROCESSING). This parallelism means Salmon quantification does not add to the critical path wall-clock time.

### Versions tracking

`NEOANTIGEN_WORKFLOW.out.versions` is mixed into the pipeline's `versions` channel, so Salmon v1.11.4 is recorded in `nf_core_rnadnavar_software_mqc_versions.yml`.

---

## Module configuration

**`conf/modules/salmon/salmon_quant.config`** — included at the end of the `includeConfig` block in `nextflow.config`:

```groovy
process {
    withName: 'SALMON_QUANT' {
        ext.args = {
            def args = "--validateMappings --libType ${params.salmon_libtype}"
            if (params.salmon_gc_bias) args += " --gcBias"
            args
        }
        publishDir = [
            path: { "${params.outdir}/${meta.id}/salmon" },
            mode: params.publish_dir_mode,
            saveAs: { filename ->
                if (filename == 'versions.yml') return null
                filename
            }
        ]
    }
    withName: 'QUANT_TSV_NORMALIZE' {
        publishDir = [
            path: { "${params.outdir}/${meta.id}/salmon" },
            mode: params.publish_dir_mode,
            saveAs: { filename -> filename == 'versions.yml' ? null : filename }
        ]
    }
}
```

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

The consensus VCF is a multi-caller aggregation designed for data labeling. All variant evidence lives in INFO fields (`VAF_BY_CALLER`, `DP_BY_CALLER`, etc.) rather than per-sample FORMAT — intentional for the pipeline's primary purpose, but incompatible with neoantigen tools that expect per-sample genotype data.

### Caller FORMAT field differences

| Field | Mutect2 | DeepSomatic | Strelka2 |
|-------|---------|-------------|----------|
| `GT` | ✓ | ✓ | ✗ |
| `AD` (ref,alt) | ✓ `Number=R` | ✓ `Number=R` | ✗ (uses `AU`/`CU`/`GU`/`TU` or `TAR`/`TIR`) |
| `AF` | ✓ (named `AF`) | ✗ (named `VAF`) | ✗ (must compute from base counts) |
| `DP` | ✓ | ✓ | ✓ |
| `GQ` | ✓ | ✓ | ✗ |
| `F1R2`/`F2R1`/`SB` | ✓ | ✗ | ✗ |

This is why `FORMAT_HARMONIZER` normalizes Strelka2 and DeepSomatic FORMAT fields to Mutect2 conventions before the consensus step. However, the consensus VCF itself still has no per-sample FORMAT data by design — it cannot be aligned to Mutect2 FORMAT because it represents a union across callers, not a single caller's genotype call.

### When to use `consensus`

Use `neoantigen_input_source = 'consensus'` only if your downstream tool can consume the INFO-field-based format and you specifically want the multi-caller union. The consensus VCF with `--neoantigen` adds `AD_BY_CALLER` and `AF_BY_CALLER` INFO fields encoding per-caller allele depths and frequencies.

---

## Error handling

| Condition | Behavior |
|-----------|----------|
| `enable_neoantigen_workflow=true` and `salmon_index=null` | `error()` before any process; pipeline exits |
| `salmon_index` path does not exist or is not a directory | `error()` before any process; pipeline exits |
| `neoantigen_input_source` not in `['consensus', 'mutect2']` | `error()` before any process; pipeline exits |
| `salmon quant` exits non-zero | Nextflow propagates the error; task marked failed; pipeline halts |
| `QUANT_TSV_NORMALIZE` receives empty/malformed `quant.sf` | Script exits non-zero; task marked failed |
| `FORMAT_HARMONIZER` encounters a field that cannot be computed | Writes `.` for that field; does not fail |
| No status=1 samples when `enable_neoantigen_workflow=true` | `log.warn` emitted; neoantigen VCF channel is empty; no output published |

---

## Cross-references

- **Index build script**: `scripts/build_salmon_index.sh` — downloads Gencode v49 references and builds a decoy-aware Salmon v1.11.x index. See the script header for usage.
- **Example configuration**: `examples/neoantigen/` — ready-to-use config files and usage guide.
- **Example usage guide**: `examples/neoantigen/README.usage.md` — user-facing quick start and parameter reference.
- **FORMAT harmonizer script**: `bin/harmonize_vcf_format.py` — Python script implementing the caller-specific FORMAT remapping logic.
- **Consensus script**: `bin/run_consensus_vcf.py` — updated with `--neoantigen` flag to emit `AD_BY_CALLER`/`AF_BY_CALLER`.
- **Salmon module config**: `conf/modules/salmon/salmon_quant.config` — `ext.args` and `publishDir` for `SALMON_QUANT` and `QUANT_TSV_NORMALIZE`.
