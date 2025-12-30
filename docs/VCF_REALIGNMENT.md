# VCF-Based Realignment Workflow

This document describes the VCF-based realignment workflow in the `rnadnavar` pipeline. This workflow is designed to improve variant calling accuracy by realigning RNA reads in regions where variants are detected in DNA or RNA consensus calls.

## Overview

The VCF-based realignment workflow consists of the following steps:

1.  **Preparation (`PREPARE_REALIGNMENT_VCF`)**:
    *   Extracts candidate regions from RNA consensus VCFs (and optionally DNA consensus VCFs).
    *   Extracts RNA reads mapping to these regions from the original CRAM/BAM files.
    *   Realigns these reads using HISAT2.
2.  **RNA Variant Calling (`RNA_REALIGNMENT_WORKFLOW`)**:
    *   Performs GATK preprocessing (MarkDuplicates, SplitNCigarReads) on realigned BAMs.
    *   Performs variant calling on the realigned RNA BAMs using configured callers.
    *   Normalizes variants (VT decompose + BCFTools norm).
    *   Generates a new RNA consensus VCF.
    *   Filters the consensus VCF.
3.  **Second Rescue (`SECOND_RESCUE_WORKFLOW`)**:
    *   Performs cross-modality rescue between the first-round DNA consensus VCF and the realigned RNA consensus VCF.
    *   Applies RNA editing annotation (REDIportal) if enabled.
    *   Applies COSMIC/gnomAD annotation if enabled.
    *   Filters the rescued variants.

## Key Features

*   **Optimized Metadata Handling**: Uses a lightweight `validateMeta` function to ensure data integrity without the overhead of complex sanitization processes.
*   **Robust Channel Joining**: Safely joins VCF and CRAM channels by patient ID to ensure correct sample pairing.
*   **Integrated Filtering**: Filtering is performed automatically within the rescue workflow, ensuring that the final output is ready for analysis.
*   **Comprehensive Annotation**: Includes COSMIC, gnomAD, and RNA editing annotations (REDIportal) with conditional execution based on parameter settings.
*   **Organized Output Structure**: All realignment-related outputs are placed in dedicated `vcf_realignment/` directory to separate them from main workflow outputs.
*   **Suffix Propagation**: Realigned samples maintain `_realign` suffix throughout the workflow for easy identification.

## Usage

To enable VCF-based realignment, set the following parameters in your configuration or command line:

```bash
--tools realignment,rescue
--realignment_mode vcf
```

### Optional Parameters

**RNA Editing Annotation:**
```bash
--enable_rna_annotation true
--rediportal_vcf /path/to/REDIportal.vcf.gz
--rediportal_tbi /path/to/REDIportal.vcf.gz.tbi
--min_rna_support 2  # Minimum number of callers required for RNA annotation
```

**COSMIC/gnomAD Annotation:**
```bash
--enable_cosmic_gnomad_annotation true
--cosmic_database /path/to/cosmic.vcf.gz
--gnomad_database /path/to/gnomad/directory
--cosmic_gnomad_verbose false  # Set to true for detailed annotation logs
```

## Output

The workflow produces the following outputs organized in dedicated directories:

*   **MAF Results**: `results/maf_results/` (Main workflow only - MAF mode)
    *   Contains `vcf2maf`, `consensus`, and `filtered` MAF files.
*   **Realignment Outputs**: `results/vcf_realignment/`
    *   `preprocessing/`: Contains realigned CRAMs/BAMs from MarkDuplicates and SplitNCigarReads with `_realign` suffix.
    *   `variant_calling/`: Contains raw VCFs from variant callers (mutect2, strelka, deepsomatic, etc.) for realigned samples.
    *   `normalized/`: Contains normalized VCFs (VT decompose + BCFTools norm) organized by variant caller.
    *   `consensus/`: Contains consensus VCFs generated from realigned variant calls.
    *   `filtered/`: Contains filtered consensus VCFs for realigned samples.
*   **Rescued VCFs**: `results/rescue/` (Final output)
    *   First-round rescue: `DNA_TUMOR_vs_DNA_NORMAL_rescued_RNA_TUMOR_vs_DNA_NORMAL.*`
    *   Second-round rescue (with realigned RNA): `DNA_TUMOR_vs_DNA_NORMAL_rescued_RNA_TUMOR_realign_vs_DNA_NORMAL.*`
    *   Annotated rescue VCFs with RNA editing and COSMIC/gnomAD annotations (if enabled).
*   **Reports**: `results/reports/vcf_realignment/`
    *   QC reports, metrics, and statistics for realignment processes.
    *   Includes samtools, mosdepth, bcftools stats, and GATK metrics.

### Output File Naming Convention

All realignment outputs include the `_realign` suffix in sample IDs to distinguish them from main workflow outputs:
- Preprocessing: `RNA_TUMOR_realign.cram`
- Variant calling: `RNA_TUMOR_realign_vs_DNA_NORMAL.vcf.gz`
- Consensus: `RNA_TUMOR_realign_vs_DNA_NORMAL.consensus.vcf.gz`
- Filtered: `RNA_TUMOR_realign_vs_DNA_NORMAL.filtered.vcf.gz`
- Second rescue: `DNA_TUMOR_vs_DNA_NORMAL_rescued_RNA_TUMOR_realign_vs_DNA_NORMAL.rescued.vcf.gz`

## Configuration Architecture

The VCF realignment workflow uses a modular configuration approach:

*   **Base Configuration**: `conf/modules/modules.config` - Defines default publishDir patterns for all processes.
*   **Realignment-Specific Overrides**: `conf/modules/prepare_realignment/vcf_realignment.config` - Overrides publishDir for all processes within `RNA_REALIGNMENT_WORKFLOW` and `SECOND_RESCUE_WORKFLOW` to route outputs to `vcf_realignment/` directory.
*   **Optimization Config**: `conf/modules/prepare_realignment/optimization.config` - Contains only PREPARE_REALIGNMENT_VCF-specific overrides (HISAT2, SAMTOOLS).

The configuration loading order in [nextflow.config](nextflow.config) ensures proper precedence:
1. Base modules.config (default patterns)
2. Specific module configs (variant calling, normalization, etc.)
3. vcf_realignment.config (loaded last to override patterns for realignment processes)

This architecture ensures that:
- Realignment outputs are properly organized without affecting main workflow outputs.
- Wildcard selector patterns (`.*:RNA_REALIGNMENT_WORKFLOW:.*:PROCESS_NAME`) catch all nested processes.
- Configuration is maintainable and easy to extend for new processes.

## Troubleshooting

**Issue: RNA annotation not running despite `enable_rna_annotation=true`**

Check the following:
1. Verify REDIportal VCF files are provided and exist:
   ```bash
   --rediportal_vcf /path/to/REDIportal.vcf.gz
   --rediportal_tbi /path/to/REDIportal.vcf.gz.tbi
   ```
2. Check pipeline logs for parameter debugging messages:
   ```
   === Second Rescue Workflow - RNA Annotation Parameters ===
   enable_rna_annotation: true
   rediportal_vcf channel: [/path/to/file]
   ```
3. Check VCF_RESCUE_POST_PROCESSING logs for conditional evaluation:
   ```
   Conditional evaluation: run_rna_annotation = true
   ```

**Issue: Realignment outputs not in `vcf_realignment/` directory**

1. Verify `vcf_realignment.config` is loaded in [nextflow.config](nextflow.config) after other module configs.
2. Check resolved configuration with:
   ```bash
   nextflow config -profile <your_profile> | grep -A5 "RNA_REALIGNMENT_WORKFLOW"
   ```
3. Ensure process names match the wildcard patterns in the config.

**Issue: `_realign` suffix missing in output files**

1. Verify metadata is propagated correctly in [bam_variant_calling/main.nf](subworkflows/local/bam_variant_calling/main.nf):
   - Should use `tumor[1].id` (which contains suffix) instead of `tumor[1].sample`.
2. Check that GATK preprocessing doesn't reset metadata in [bam_gatk_preprocessing/main.nf](subworkflows/local/bam_gatk_preprocessing/main.nf).

**Issue: Second rescue VCF not generated**

1. Verify `rescue` is in the `--tools` parameter:
   ```bash
   --tools realignment,rescue
   ```
2. Check that first-round DNA consensus VCF exists (status <= 1).
3. Check that realigned RNA consensus VCF was generated successfully.

## Verification

Use the provided verification script to check output structure:

```bash
./tests/verify_output_structure.sh
```

This script checks for:
- Proper directory structure (`vcf_realignment/`, `maf_results/`, `rescue/`)
- Presence of `_realign` suffix in filenames
- Existence of both first-round and second-round rescue outputs
- All expected subdirectories within `vcf_realignment/`
