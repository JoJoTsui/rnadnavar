# VCF-Based Realignment Workflow

This document describes the VCF-based realignment workflow in the `rnadnavar` pipeline. This workflow is designed to improve variant calling accuracy by realigning RNA reads in regions where variants are detected in DNA or RNA consensus calls.

## Overview

The VCF-based realignment workflow consists of the following steps:

1.  **Preparation (`PREPARE_REALIGNMENT_VCF`)**:
    *   Extracts candidate regions from RNA consensus VCFs (and optionally DNA consensus VCFs).
    *   Extracts RNA reads mapping to these regions from the original CRAM/BAM files.
    *   Realigns these reads using HISAT2.
2.  **RNA Variant Calling (`RNA_REALIGNMENT_WORKFLOW`)**:
    *   Performs variant calling on the realigned RNA BAMs.
    *   Generates a new RNA consensus VCF.
3.  **Rescue (`SECOND_RESCUE_WORKFLOW`)**:
    *   Performs cross-modality rescue between the first-round DNA consensus VCF and the realigned RNA consensus VCF.
    *   Annotates and filters the rescued variants.

## Key Features

*   **Optimized Metadata Handling**: Uses a lightweight `validateMeta` function to ensure data integrity without the overhead of complex sanitization processes.
*   **Robust Channel Joining**: Safely joins VCF and CRAM channels by patient ID to ensure correct sample pairing.
*   **Integrated Filtering**: Filtering is performed automatically within the rescue workflow, ensuring that the final output is ready for analysis.
*   **Comprehensive Annotation**: Includes COSMIC, gnomAD, and RNA editing annotations (REDIportal).

## Usage

To enable VCF-based realignment, set the following parameters in your configuration or command line:

```bash
--tools realignment,rescue
--realignment_mode vcf
```

## Output

The workflow produces the following outputs:

*   **Realigned BAMs**: `results/realignment/bam/`
*   **Realigned RNA VCFs**: `results/realignment/vcf/`
*   **Rescued VCFs**: `results/rescue/` (Final output)

## Troubleshooting

If you encounter issues:
*   Check the `multiqc_report.html` for QC metrics.
*   Verify that your input samplesheet contains correct `patient`, `sample`, and `status` fields.
*   Ensure that reference files (FASTA, HISAT2 index, etc.) are correctly specified.
