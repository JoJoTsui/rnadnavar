# RNADNAvar Workflow Documentation

This document provides comprehensive documentation for the `rnadnavar` pipeline, covering the complete workflow architecture from alignment through variant calling, consensus generation, rescue, annotation, and realignment.

## Table of Contents

1. [Overview](#overview)
2. [Pipeline Architecture](#pipeline-architecture)
3. [VCF Consensus Workflow](#vcf-consensus-workflow)
4. [VCF Rescue Workflow](#vcf-rescue-workflow)
5. [RNA Editing Annotation](#rna-editing-annotation)
6. [COSMIC/gnomAD Annotation](#cosmicgnomad-annotation)
7. [VCF Filtering](#vcf-filtering)
8. [VCF Realignment Workflow](#vcf-realignment-workflow)
9. [Configuration Reference](#configuration-reference)
10. [Output Structure](#output-structure)
11. [Logging and Debugging](#logging-and-debugging)

---

## Overview

The **rnadnavar** pipeline is a Nextflow-based bioinformatics workflow for integrating RNA and DNA sequencing data to identify somatic variants with cross-modality validation. The pipeline performs:

- Multi-modality alignment (DNA with BWA-MEM, RNA with HISAT2/STAR)
- Variant calling with multiple callers (Mutect2, Strelka, DeepSomatic, SAGE)
- Consensus generation from multiple callers
- Cross-modality rescue (DNA ↔ RNA variant validation)
- RNA editing annotation using REDIportal database
- COSMIC/gnomAD annotation for somatic/germline classification
- Optional VCF-based realignment for improved variant calling

### Key Features

- **Cross-modality Validation**: Leverages both DNA and RNA sequencing to validate somatic variants
- **Multi-caller Consensus**: Combines calls from multiple variant callers for higher confidence
- **RNA Editing Detection**: Annotates known RNA editing sites to distinguish from true mutations
- **Integrated Annotation**: COSMIC and gnomAD databases for variant classification
- **Clean Logging**: All debug logging gated behind `--debug_verbose` flag

---

## Pipeline Architecture

### High-Level Flow

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                            INPUT SAMPLESHEET                                 │
│                    (DNA normal, DNA tumor, RNA tumor)                        │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                         PREPARE_REFERENCE_AND_INTERVALS                      │
│                    (Build genome indices, intervals)                         │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                              BAM_ALIGN                                       │
│           DNA: BWA-MEM    │    RNA: HISAT2/STAR                             │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                            BAM_PROCESSING                                    │
│  ┌─────────────────┐    ┌─────────────────┐    ┌─────────────────────────┐  │
│  │ GATK Preprocess │ -> │ Variant Calling │ -> │ VCF_CONSENSUS_WORKFLOW  │  │
│  │ (MarkDup, BQSR) │    │ (Mutect2,       │    │ (Union consensus)       │  │
│  │                 │    │  Strelka,       │    └───────────┬─────────────┘  │
│  │                 │    │  DeepSomatic)   │                │               │
│  └─────────────────┘    └─────────────────┘                ▼               │
│                                            ┌─────────────────────────────┐  │
│                                            │ VCF_RESCUE_WORKFLOW         │  │
│                                            │ (DNA ↔ RNA cross-rescue)    │  │
│                                            └───────────┬─────────────────┘  │
│                                                        │                    │
│                                            ┌───────────▼─────────────────┐  │
│                                            │ VCF_RESCUE_POST_PROCESSING  │  │
│                                            │ ┌─────────────────────────┐ │  │
│                                            │ │ COSMIC/gnomAD Annotation│ │  │
│                                            │ │ RNA Editing Annotation  │ │  │
│                                            │ │ VCF Rescue Filtering    │ │  │
│                                            │ └─────────────────────────┘ │  │
│                                            └─────────────────────────────┘  │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                        (if --tools realignment)
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                         PREPARE_REALIGNMENT_VCF                              │
│                    (Extract reads, HISAT2 realign)                           │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                        RNA_REALIGNMENT_WORKFLOW                              │
│       (Variant calling + consensus on realigned RNA)                         │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                        SECOND_RESCUE_WORKFLOW                                │
│       (DNA consensus + realigned RNA consensus rescue)                       │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                              MULTIQC                                         │
│                    (Aggregate QC reports)                                    │
└─────────────────────────────────────────────────────────────────────────────┘
```

### Sample Status Codes

The pipeline uses status codes to distinguish sample types:

| Status | Description | Role |
|--------|-------------|------|
| 0 | DNA Normal | Germline reference for variant calling |
| 1 | DNA Tumor | Somatic mutation detection |
| 2 | RNA Tumor | RNA expression + variant validation |

---

## VCF Consensus Workflow

**Location**: `subworkflows/local/vcf_consensus_workflow/main.nf`

### Purpose

Merge variant calls from multiple callers into a consensus VCF using a union approach. Variants are aggregated and annotated with caller support information.

### Process Flow

1. **Group VCFs by Sample**: Groups VCFs by patient, sample, and status
2. **Execute Consensus**: Runs Python consensus script (`bin/run_consensus_vcf.py`)
3. **Separate by Modality**: Filters DNA (status ≤ 1) and RNA (status = 2) VCFs
4. **Invoke Rescue**: Passes DNA and RNA consensus to rescue workflow

### Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `params.rescue_snv_thr` | 2 | SNV caller count threshold for rescue |
| `params.rescue_indel_thr` | 2 | Indel caller count threshold for rescue |
| `params.defaultvariantcallers` | "sage,strelka,mutect2" | Expected caller list |

### Consensus Algorithm

The consensus script (`bin/run_consensus_vcf.py`):

1. Reads variants from all input VCFs using cyvcf2
2. Aggregates variants with same `chrom:pos:ref:alt` key
3. Tracks which callers support each variant
4. Marks variants meeting caller threshold as "consensus"
5. Adds INFO fields:
   - `NCALLERS` - Number of supporting callers
   - `CALLERS` - List of caller names
   - `VAF_*` - Per-caller VAF values

### Output

```
consensus/
└── DNA_TUMOR_vs_DNA_NORMAL/
    ├── DNA_TUMOR_vs_DNA_NORMAL.consensus.vcf.gz
    └── DNA_TUMOR_vs_DNA_NORMAL.consensus.vcf.gz.tbi
```

---

## VCF Rescue Workflow

**Location**: `subworkflows/local/vcf_rescue_workflow/main.nf`

### Purpose

Cross-modality variant rescue - find variants with support from both DNA and RNA sequencing, increasing confidence in somatic variant calls.

### Rescue Logic

The rescue script (`bin/run_rescue_vcf.py`) performs:

1. **Load Modalities**: Reads DNA and RNA consensus VCFs with modality tags
2. **Key-Based Aggregation**: Uses `(chrom, pos, ref, alt)` as unique key
3. **Cross-Reference**: Identifies variants found in both modalities
4. **Evidence Scoring**: Counts caller support from each modality
5. **Output Tagging**: Adds rescue and cross-modality INFO fields

### Output INFO Fields

| Field | Description |
|-------|-------------|
| `N_DNA_CALLERS_SUPPORT` | Number of DNA callers supporting variant |
| `N_RNA_CALLERS_SUPPORT` | Number of RNA callers supporting variant |
| `RESCUED` | Variant was rescued from other modality |
| `CROSS_MODALITY` | Variant found in both DNA and RNA |
| `MODALITY` | Source modality (DNA, RNA, or BOTH) |

### Output

```
rescue/
└── DNA_TUMOR_vs_DNA_NORMAL_rescued_RNA_TUMOR_vs_DNA_NORMAL/
    ├── *.rescued.vcf.gz
    ├── *.rescued.vcf.gz.tbi
    ├── *.rescue.rna_annotated.vcf.gz      # if RNA annotation enabled
    ├── *.rescue.rna_annotated.vcf.gz.tbi
    ├── *.filtered.vcf.gz
    ├── *.filtered.vcf.gz.tbi
    ├── *.filtered.vcf.stripped.vcf.gz
    └── *.filtered.vcf.stripped.vcf.gz.tbi
```

---

## RNA Editing Annotation

**Location**: `modules/local/rna_editing_annotation/main.nf`

### Purpose

Annotate variants with RNA editing information from the REDIportal database to distinguish true somatic mutations from RNA editing events (primarily A-to-I editing).

### Database

The annotation uses REDIportal, a comprehensive database of RNA editing sites:
- Source: http://srv00.recas.ba.infn.it/atlas/index.html
- Format: VCF with editing site annotations

### Algorithm

The annotation script (`bin/annotate_rna_editing.py`) performs:

1. **REDIportal Preparation**: Converts database to proper VCF format
2. **VCF-to-VCF Annotation**: Uses bcftools for exact coordinate/allele matching
3. **Evidence Tiering**: Classifies variants by RNA editing evidence level
4. **FILTER Updates**: Sets `RNAedit` filter for high-evidence variants

### Evidence Tiers

| Tier | Criteria | Description |
|------|----------|-------------|
| HIGH | REDIportal match + canonical A>G/T>C + RNA consensus | Strong RNA editing evidence |
| MEDIUM | REDIportal match + RNA consensus | Moderate evidence |
| LOW | RNA consensus only | Weak evidence |
| NONE | No RNA editing evidence | Likely true mutation |

### Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `params.enable_rna_annotation` | false | Enable RNA editing annotation |
| `params.rediportal_vcf` | null | Path to REDIportal VCF |
| `params.rediportal_tbi` | null | Path to REDIportal index |
| `params.min_rna_support` | 2 | Min RNA callers for evidence |

### Output INFO Fields

| Field | Description |
|-------|-------------|
| `REDI_EVIDENCE` | Evidence level (HIGH/MEDIUM/LOW/NONE) |
| `REDI_CANONICAL` | Is canonical A>G/T>C transition (YES/NO) |
| `REDI_ACCESSION` | REDIportal accession ID |
| `REDI_DB` | Source database |
| `REDI_TYPE` | Editing type |

### Conditional Execution

The process uses a `when` directive to skip execution if database files are not provided:

```groovy
when:
    rediportal_vcf.name != 'NO_FILE' && rediportal_tbi.name != 'NO_FILE'
```

---

## COSMIC/gnomAD Annotation

**Location**: `modules/local/cosmic_gnomad_annotation/main.nf`

### Purpose

Annotate variants with COSMIC somatic mutation database and gnomAD population frequencies to classify variants as somatic, germline, or artifacts.

### Algorithm

The annotation script (`bin/annotate_cosmic_gnomad.py`) performs:

1. **COSMIC Annotation**: Uses bcftools annotate with COSMIC VCF
2. **gnomAD Annotation**: Scatter-gather by chromosome for parallelism
3. **Variant Classification**: Multi-modal classification with FILTER updates

### Classification Logic

| Classification | Criteria |
|----------------|----------|
| Germline | Population frequency > threshold (default 0.001) |
| Somatic | Cross-modality support + COSMIC recurrence ≥ threshold |
| Reference | Low quality or insufficient evidence |
| Artifact | Failed quality filters |

### Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `params.enable_cosmic_gnomad_annotation` | true | Enable COSMIC/gnomAD annotation |
| `params.cosmic_database` | null | Path to COSMIC VCF |
| `params.gnomad_database` | null | Path to gnomAD directory |
| `params.cosmic_gnomad_germline_freq_threshold` | 0.001 | gnomAD frequency for germline |
| `params.cosmic_gnomad_somatic_consensus_threshold` | 3 | Min callers for somatic |
| `params.cosmic_gnomad_cosmic_recurrence_threshold` | 5 | COSMIC count for rescue |

### Output INFO Fields

| Field | Description |
|-------|-------------|
| `COSMIC_ID` | COSMIC mutation ID |
| `COSMIC_CNT` | COSMIC occurrence count |
| `gnomAD_AF` | gnomAD allele frequency |
| `GNOMAD_RESCUE` | Rescued based on low gnomAD frequency |
| `COSMIC_RESCUE` | Rescued based on COSMIC recurrence |

---

## VCF Filtering

**Location**: `subworkflows/local/vcf_rescue_filter/main.nf`

### Purpose

Apply RaVeX filtering logic to rescue VCFs while preserving classifications and RNA editing annotations.

### Filter Script

The filtering script (`bin/filter_rescue_vcf.py`) applies:

1. **Preservation**: Keeps original FILTER values (Somatic, Germline, Reference, Artifact)
2. **Preservation**: Keeps RNA editing annotations (REDI_* INFO fields)
3. **RaVeX Filters**: Adds RaVeX filter flags to INFO field
4. **Stripped Output**: Generates multiallelic-filtered version

### Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `params.gnomad_freq_threshold` | 0.0001 | gnomAD frequency threshold |
| `params.min_alt_reads` | 2 | Minimum alt read support |
| `params.filter_multiallelic` | false | Filter multiallelic sites |

### Output Files

| File | Description |
|------|-------------|
| `*.filtered.vcf.gz` | Filtered VCF with all INFO preserved |
| `*.filtered.vcf.stripped.vcf.gz` | Stripped VCF with multiallelic sites filtered |

---

## VCF Realignment Workflow

For detailed documentation on the VCF-based realignment workflow, see [VCF_REALIGNMENT.md](VCF_REALIGNMENT.md).

### Overview

The realignment workflow improves variant calling accuracy by:

1. Extracting reads mapping to variant regions from RNA CRAM files
2. Realigning these reads with HISAT2
3. Performing variant calling on realigned BAMs
4. Generating new consensus VCF
5. Running second-round rescue with DNA consensus

### Enabling Realignment

```bash
--tools realignment,rescue
--realignment_mode vcf
```

---

## Configuration Reference

### Essential Parameters

```bash
# Input/Output
--input samplesheet.csv
--outdir results/

# Alignment
--aligner bwa-mem                    # DNA aligner (bwa-mem, bwa-mem2)
--star_index /path/to/star           # RNA aligner index

# Variant Calling
--tools deepsomatic,mutect2,strelka,consensus,rescue
--defaultvariantcallers "deepsomatic,strelka,mutect2"

# Consensus/Rescue
--rescue_snv_thr 2                   # SNV caller threshold
--rescue_indel_thr 2                 # Indel caller threshold
```

### RNA Editing Annotation

```bash
--enable_rna_annotation true
--rediportal_vcf /path/to/REDIportal.vcf.gz
--rediportal_tbi /path/to/REDIportal.vcf.gz.tbi
--min_rna_support 2
```

### COSMIC/gnomAD Annotation

```bash
--enable_cosmic_gnomad_annotation true
--cosmic_database /path/to/cosmic.vcf.gz
--gnomad_database /path/to/gnomad/directory
--cosmic_gnomad_germline_freq_threshold 0.001
--cosmic_gnomad_somatic_consensus_threshold 3
--cosmic_gnomad_cosmic_recurrence_threshold 5
```

### Realignment

```bash
--realignment_mode vcf               # 'vcf' (default) or 'maf'
--force_legacy_realignment false     # Force MAF-based realignment
```

### Debugging

```bash
--debug_verbose true                 # Enable verbose logging to stdout
```

---

## Output Structure

### Standard Output

```
results/
├── preprocessing/
│   ├── markduplicates/
│   ├── recal_table/
│   ├── recalibrated/
│   └── splitncigarreads/
├── variant_calling/
│   ├── deepsomatic/
│   ├── mutect2/
│   └── strelka/
├── normalized/
│   ├── deepsomatic/
│   ├── mutect2/
│   └── strelka/
├── annotation/
│   ├── deepsomatic/
│   ├── mutect2/
│   └── strelka/
├── consensus/
│   ├── DNA_TUMOR_vs_DNA_NORMAL/
│   └── RNA_TUMOR_vs_DNA_NORMAL/
├── filtered/
│   ├── DNA_TUMOR_vs_DNA_NORMAL/
│   └── RNA_TUMOR_vs_DNA_NORMAL/
├── rescue/                          # Consolidated first-round rescue
│   └── DNA_TUMOR_vs_DNA_NORMAL_rescued_RNA_TUMOR_vs_DNA_NORMAL/
│       ├── *.rescued.vcf.gz
│       ├── *.rescue.rna_annotated.vcf.gz
│       └── *.filtered.vcf.gz
└── vcf_realignment/                 # Realignment outputs
    ├── preprocessing/
    ├── variant_calling/
    ├── normalized/
    ├── consensus/
    ├── filtered/
    └── rescue/                      # Second-round rescue
        └── DNA_TUMOR_vs_DNA_NORMAL_rescued_RNA_TUMOR_realign_vs_DNA_NORMAL/
```

---

## Logging and Debugging

### Default Mode (`--debug_verbose false`)

- **stdout**: Only Nextflow execution progress (process completion, task counts)
- **.nextflow.log**: Contains all log messages including debug information

### Debug Mode (`--debug_verbose true`)

Enables detailed workflow state logging:

```bash
nextflow run main.nf --debug_verbose true ...
```

Debug output includes:
- Realignment mode selection
- Channel contents and metadata validation
- CRAM to BAM conversion status
- Memory usage monitoring
- Conditional evaluation results

### Error Messages

Error messages (`log.error`) are always displayed regardless of `--debug_verbose` setting:
- Missing required metadata fields
- File validation failures
- Process execution errors

### Viewing Logs

```bash
# View execution log
cat .nextflow.log

# Follow log in real-time
tail -f .nextflow.log

# Search for specific process
grep "VCF_RESCUE" .nextflow.log
```

---

## Version History

### Current Version (1.0dev)

- **VCF Consensus**: Union-based consensus with caller aggregation
- **Cross-modality Rescue**: DNA ↔ RNA variant validation
- **RNA Editing Annotation**: REDIportal-based A-to-I editing detection
- **COSMIC/gnomAD Annotation**: Somatic/germline classification
- **VCF Realignment**: Improved variant calling through read realignment
- **Unified Logging**: All debug logging gated behind `--debug_verbose`
- **Consolidated Outputs**: First and second-round rescue outputs in dedicated directories
