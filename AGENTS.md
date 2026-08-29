# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**nf-core/rnadnavar** is a bioinformatics pipeline for RNA and DNA integrated analysis for somatic mutation detection. It uses Nextflow (DSL2/Groovy) for workflow orchestration and implements a consensus-based approach across multiple variant callers (Mutect2, Strelka2, DeepSomatic). The pipeline is designed for cancer research and supports both single-sample and multi-sample analyses.

> **Note on utilized variant callers:** Only **Mutect2, Strelka2, and DeepSomatic** are currently utilized as variant callers. SAGE and Manta modules exist in the repo (`modules/local/sage/`, `modules/nf-core/manta/`, and corresponding subworkflows/configs), but are presently ignored by the active workflow. Set `--tools` without them for a standard run.

Nextflow >= 24.10.5 is required. The nf-core template version is 3.3.2.

## Common Development Commands

### Running Tests

```bash
# Run full test suite
nf-test test --profile debug,test,docker --verbose

# Run a single test file
nf-test test tests/default.nf.test --profile test,docker

# Run a specific module test
nf-test test modules/local/vcf_consensus/tests/ --profile test,docker

# Run with singularity (CI default)
nf-test test --profile test,singularity

# Run integration tests
./tests/run_final_integration_tests.sh
./tests/run_rna_editing_tests.sh
./tests/run_vcf_realignment_tests.sh
```

The nf-test config (`nf-test.config`) sets `testsDir "."` and profile `test`. CI uses nf-test 0.9.2, Nextflow 24.10.x, and Python 3.13. Tests auto-shard across parallel CI jobs.

### Python Unit Tests (uv-managed venv)

The repo's Python environment is uv-managed: `pyproject.toml` + `uv.lock`, with the venv at `.venv/` (created by uv; Python 3.10). Run the Python suites with the venv interpreter:

```bash
.venv/bin/python -m pytest tests/ -q --ignore=tests/test_vcf_stats --ignore=tests/vcf_stats
```

Note: `.venv` contains extra packages beyond the lockfile (e.g. `maturin`, needed to build `bin/vcf_stats/seq2neo/stats_core.so`). Do NOT run bare `uv sync` — it would uninstall those; use `uv sync --inexact` if you must reconcile. New Python tools should stay stdlib-only where possible (see `bin/label_qc.py`).

### Pipeline Execution (Local Development)

```bash
# Local test run (use . for local development, not nf-core/rnadnavar)
nextflow run . -profile test,docker --input assets/samplesheet.csv --outdir results/

# Run with specific tools
nextflow run . -profile docker \
  --input samplesheet.csv \
  --outdir results/ \
  --tools strelka,mutect2,deepsomatic,vep,consensus

# Resume from checkpoint
nextflow run . -profile docker -resume

# Use a params file
nextflow run . -profile docker -params-file params.yml

# Build reference indexes only
nextflow run . -profile docker --build_only_index
```

### Code Quality

```bash
nf-core lint .        # nf-core linting
nextflow run . --help  # View all parameters
```

## Architecture Overview

### Pipeline Stages

The main workflow (`workflows/rnadnavar.nf`) executes these sequential subworkflows:

1. **PREPARE_REFERENCE_AND_INTERVALS** - Build genome indices and prepare interval files
2. **BAM_ALIGN** - Align reads using BWA/BWA-mem2/dragmap (DNA) or STAR (RNA)
3. **BAM_PROCESSING** - GATK preprocessing (MarkDuplicates, SplitNCigarReads, BQSR) + variant calling
4. **VCF_CONSENSUS** - Combine variant calls using consensus approach
5. **MAF_FILTERING** - Filter MAF files with gnomAD/whitelist/blacklist
6. **MAF_FILTERING_RNA** - RNA-specific filtering
7. **PREPARE_REALIGNMENT** / **REALIGNMENT** - Optional realignment for RNA (rescue workflow)
8. **MULTIQC** - Aggregate QC reports

### Subworkflows

All subworkflows live under `subworkflows/local/`. Each pipeline stage maps to multiple subworkflows (e.g., BAM_PROCESSING uses `bam_gatk_preprocessing`, `bam_variant_calling_pre_post_processing`, etc.). Channel-creation subworkflows (`channel_*`) handle the data routing between stages.

### Module Organization

```
modules/
├── nf-core/           # Standard nf-core modules (BWA, STAR, GATK4, VEP, etc.)
└── local/             # Custom modules for rnadnavar-specific logic
    ├── rna_editing_annotation/    # REDIportal annotation
    ├── cosmic_gnomad_annotation/  # COSMIC/gnomAD annotation
    ├── vcf_consensus/             # Consensus variant calling
    ├── vcf_rescue/                # Cross-modality rescue
    ├── maf_filtering/             # MAF file filtering
    ├── maf_rna_filtering/         # RNA-specific MAF filtering
    ├── vcf_filtering/             # VCF-level filtering
    ├── sage/                      # SAGE caller wrapper (currently unused — see note above)
    └── ...                        # Other utilities (vcf2bed, vcf2maf, vt, etc.)
```

### Python Packages

The `bin/` directory contains two major Python packages and standalone scripts:

**`bin/vcf_utils/`** — Core VCF operations library:
- `io_utils.py` — VCF reading/writing with custom INFO fields
- `aggregation.py` — Variant aggregation from multiple callers/modalities
- `tagging.py` — Caller support and modality metadata tagging
- `filters.py` — Filter normalization across callers
- `variant_classifier.py` / `variant_classifier_unified.py` — Somatic/germline/artifact/reference classification
- `rna_editing_core.py` — RNA editing detection logic
- `cosmic_annotator.py` / `gnomad_annotator.py` — Database annotation
- `annotation_utils.py` / `bcftools_annotator.py` — Annotation helpers
- `evidence_tiering.py` / `variant_statistics.py` — Statistical analysis

**`bin/common/`** — Classification and filtering configuration:
- `vcf_config.py` — VCF processing configuration
- `category_matcher.py` / `database_checker.py` — Variant matching and validation
- `tier_config.py` — Tiering rules and thresholds

**`bin/vcf_stats/`** — Statistics, visualization, and validation:
- `comparison.py` / `caller_comparison.py` — Cross-caller comparisons
- `tiering.py` / `tiering_engine.py` — Tier assignment logic
- `visualizer.py` / `igv_reports_wrapper.py` — Report generation
- `workflow.py` — Orchestrates multi-step stats pipelines
- `bam_validator.py` / `vcf_processor.py` — File validation

**Standalone scripts in `bin/`:**
- `run_consensus_vcf.py` — Within-modality consensus variant calling
- `run_rescue_vcf.py` — Cross-modality (DNA/RNA) variant rescue
- `annotate_rna_editing.py` — REDIportal RNA editing annotation
- `annotate_cosmic_gnomad.py` — COSMIC/gnomAD database annotation
- `filter_vcf.py` / `filter_mutations.py` / `filter_rna_mutations.py` / `filter_rescue_vcf.py` — Various filtering pipelines
- `maf2bed.py` — MAF-to-BED conversion for realignment
- `run_consensus.R` — R-based consensus statistics
- `label_qc.py` — label QC gate for truth VCFs: Tier A (VCF-only rules, pure stdlib) plus optional Tier B BAM verification via samtools mpileup subprocess (`--verify-bam`; defaults in `bin/label_qc_config.json`; design doc: `dev_docs/implementation/label_qc_design.md`; tests: `tests/label_qc/`)

### Configuration System

`nextflow.config` is the root config. Modular configs:
- `conf/base.config` — Executor resource defaults
- `conf/modules.config` — Module-level process configs (publishers, resource limits)
- `conf/modules/` — Per-domain module configs: `alignment`, `variant_calling`, `consensus`, `filtering`, `annotate`, `rescue`, etc.
- `conf/test.config` / `conf/test_full.config` — Test profiles
- `conf/igenomes.config` — Reference genome paths
- `conf/vcf_realignment_optimized.config` — Realignment optimization settings
- `conf/legacy_realignment.config` — Legacy MAF-based realignment config

Custom external configs are loaded from `params.custom_config_base` (default: `/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/configs`).

### Important Parameters

- `--step` — Pipeline starting point: `mapping` (default), `variant_calling`, `preprocessing`
- `--tools` — Comma-separated tools. Utilized variant callers: `mutect2,strelka,deepsomatic`. Also valid: `vep,consensus,filtering,realignment,rna_filtering,norm,vcf2maf,preprocessing,rescue`. SAGE/Manta are recognized by config schema but currently ignored by the active workflow.
- `--rna` / `--dna` — Enable RNA/DNA analysis (both default: true)
- `--aligner` — `bwa-mem` (default), `bwa-mem2`, `dragmap`, or STAR for RNA
- `--defaultvariantcallers` — `sage,strelka,mutect2` in `nextflow.config` (used when `--tools` is null), but only strelka/mutect2 are effectively driven; enable `deepsomatic` via `--tools`
- `--rescue_snv_thr` / `--rescue_indel_thr` — Consensus thresholds (default: 2)
- `--realignment_mode` — `vcf` (new default) or `maf` (legacy)
- `--joint_mutect2` — Patient-wise multi-sample somatic calling
- `--enable_rna_annotation` — REDIportal RNA editing annotation
- `--enable_cosmic_gnomad_annotation` — COSMIC/gnomAD annotation (default: true)
- `--trim_fastq` — Enable FastQC + fastp trimming
- `--save_bam_mapped` / `--save_output_as_bam` — Control intermediate file retention
- `--read_length` — Required for STAR index building (default: 76)
- `--wes` — Use WES models in variant callers (default: false)

### Sample Sheet Format

```csv
sample,lane,fastq_1,fastq_2
CONTROL_REP1,LX,path/to/R1.fastq.gz,path/to/R2.fastq.gz
```

## Key Development Workflows

### Adding a New Variant Caller

1. Add module in `modules/nf-core/` or `modules/local/`
2. Create module config in `conf/modules/variant_calling/`
3. Update the subworkflow in `subworkflows/local/bam_variant_calling_pre_post_processing/`
4. Add to `params.tools` and `params.defaultvariantcallers` in `nextflow.config`
5. Add nf-test in the module's `tests/` directory

### Modifying Consensus Logic

Consensus is implemented in `bin/run_consensus_vcf.py` using the `bin/vcf_utils/` package. Key thresholds:
- `--snv_thr` — Minimum callers for SNV consensus (default: 2)
- `--indel_thr` — Minimum callers for indel consensus (default: 2)

See `docs/consensus_logic_explained.md` and `docs/consensus_vcf_rules.md` for detailed rules.

### Working with the Variant Classification System

The unified classification system (`bin/vcf_utils/variant_classifier_unified.py`, configured via `bin/common/vcf_config.py`) classifies variants into four biological categories (Somatic, Germline, Reference, Artifact). Each caller (Strelka, DeepSomatic, Mutect2) has caller-specific classification logic. Categories are used directly as VCF FILTER values. See `bin/CLASSIFICATION_SYSTEM.md` for full documentation.

### RNA Editing Annotation

Uses `bin/annotate_rna_editing.py` with `bin/vcf_utils/rna_editing_core.py` to distinguish somatic mutations from RNA editing events. Key parameters: `--rediportal_vcf`, `--min_rna_support` (default: 2), `--rna_annotation_log_level` (default: INFO).

### Realignment Modes

Two realignment modes exist:
- **vcf** (default) — New VCF-based realignment with optimized error handling. Configured in `conf/vcf_realignment_optimized.config`.
- **maf** (legacy) — Old MAF-based realignment. Configured in `conf/legacy_realignment.config`.

Control with `--realignment_mode`, `--force_legacy_realignment`, or `--disable_realignment_optimization`. VCF optimization features: `--enable_vcf_realignment_optimization`.

## Testing

Uses **nf-test** framework with `.nf.test` files:
- Module tests in `modules/*/tests/`
- Pipeline tests in `tests/`
- Test profiles in `conf/test/` and `conf/test_full.config`

Test data sources:
- `params.modules_testdata_base_path` — nf-core test data
- `params.pipelines_testdata_base_path` — rnadnavar test data

nf-test plugins used: `nft-bam@0.4.0`, `nft-utils@0.0.5`, `nft-vcf@1.0.7`.

## Reference Documentation

The `docs/` directory contains detailed guides on specific pipeline features:
- `consensus_logic_explained.md` / `consensus_vcf_rules.md` — Consensus variant calling rules
- `rescue_workflow.md` / `rescue_vcf_rules.md` — Cross-modality rescue logic
- `VCF_REALIGNMENT.md` — Realignment workflow details
- `MAF_FILTERING_COMPREHENSIVE_GUIDE.md` — MAF filtering documentation
- `RNA_EDITING_ANNOTATION_GUIDE.md` — RNA editing annotation
- `RNADNAVAR_WORKFLOW.md` — Full pipeline workflow documentation

## Branching

- `main` — Latest stable release
- `dev` — Development branch, merge feature branches here
- Feature branches — Named by feature (e.g., `realignment`, `vcf_refactoring`, `wes`)
