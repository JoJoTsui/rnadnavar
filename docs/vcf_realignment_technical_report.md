# VCF-Based RNA Realignment: Technical Report

## 1. Executive Summary

### 1.1 Purpose

The VCF-based RNA realignment feature provides a targeted approach to re-aligning RNA-seq reads specifically at regions containing RNA consensus variants. This differs from the traditional MAF-based approach which re-aligns all variants from the filtered MAF file.

### 1.2 Key Features

- **Targeted Realignment**: Only RNA consensus variant regions are re-aligned, reducing computational overhead
- **Two-Stage Rescue**: Cross-modality rescue is performed twice—before and after realignment
- **VCF-Only Workflow**: No MAF dependency in the realignment pipeline, simplifying data flow
- **DNA Normal Pairing**: Properly handles somatic variant calling with DNA normal as paired sample
- **Backward Compatibility**: MAF-based realignment is preserved as an optional mode

### 1.3 Impact on Pipeline Outputs

When VCF-based realignment is enabled (`--realignment_mode vcf`), the following additional outputs are generated:

| Output | Description |
|--------|-------------|
| `realigned_rna_consensus_vcf` | RNA consensus VCF from HISAT2-realigned BAM |
| `realigned_filtered_rna_vcf` | Filtered RNA VCF from realigned BAM |
| `second_rescued_vcf` | Second rescued VCF combining DNA (first round) + realigned RNA, with full annotation |

---

## 2. Architecture Overview

### 2.1 Sample Status Codes

The pipeline uses three status codes in the samplesheet:

| Status | Sample Type | Realignment |
|--------|-------------|-------------|
| 0 | DNA normal | **Included** |
| 1 | DNA tumor | Excluded |
| 2 | RNA tumor | **Included** |

### 2.2 Realignment Sample Composition

Both MAF and VCF realignment modes process the same sample types:

- **Included**: RNA (status=2) + DNA normal (status=0)
- **Excluded**: DNA tumor (status=1)

This design is necessary because somatic variant callers (Mutect2, Strelka, etc.) require tumor-normal pairs:
- Realigned RNA serves as the "tumor" sample
- DNA normal provides the germline baseline

### 2.3 Component Interaction Diagram

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                        FIRST ROUND (Standard Pipeline)                       │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                              │
│   BAM_ALIGN (STAR for RNA, BWA for DNA)                                  │
│          ↓                                                                  │
│   BAM_VARIANT_CALLING_PRE_POST_PROCESSING                                 │
│     ├─ GATK preprocessing (MarkDuplicates, SplitNCigarReads, BQSR)        │
│     ├─ Variant calling (Mutect2, Strelka, SAGE)                          │
│     ├─ VCF normalization                                                   │
│     ├─ VEP annotation                                                     │
│     ├─ VCF_CONSENSUS_WORKFLOW → vcf (DNA+RNA consensus)                 │
│     ├─ VCF_RESCUE_WORKFLOW → vcf_rescue (first rescue)                 │
│     ├─ VCF_RESCUE_POST_PROCESSING → rescued vcf                         │
│     ├─ MAF_CONSENSUS → maf                                                │
│     ├─ VCF_FILTERING → filtered_vcf                                      │
│     └─ MAF_FILTERING → filtered_maf                                      │
│                                                                              │
│   STORE FOR SECOND RESCUE:                                               │
│     ├─ first_round_dna_consensus_vcf (status <= 1)                       │
│     └─ first_round_dna_consensus_vcf (used in second rescue)            │
│                                                                              │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                        REALIGNMENT (VCF Mode Enabled)                      │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                              │
│   RNA consensus VCF (from first round, status == 2)                        │
│          ↓                                                                  │
│   PREPARE_REALIGNMENT_VCF                                                 │
│     ├─ VCF2BED (convert VCF to BED regions)                             │
│     ├─ SAMTOOLS_EXTRACT_READ_IDS (extract reads from regions)            │
│     ├─ SAMTOOLS_CONVERT (CRAM to BAM)                                   │
│     ├─ PICARD_FILTERSAMREADS (filter to extracted reads)                 │
│     ├─ BAM_CONVERT_SAMTOOLS (BAM to FASTQ)                              │
│     └─ FASTQ_ALIGN_HISAT2 (realign with HISAT2)                         │
│          ↓                                                                  │
│   Realigned BAM (RNA + DNA normal)                                        │
│          ↓                                                                  │
│   RNA_REALIGNMENT_WORKFLOW                                                │
│     ├─ BAM_GATK_PREPROCESSING                                            │
│     ├─ BAM_VARIANT_CALLING                                               │
│     ├─ VCF_NORMALIZE                                                      │
│     ├─ VCF_ANNOTATE                                                       │
│     ├─ VCF_CONSENSUS_WORKFLOW → realigned_rna_consensus_vcf             │
│     └─ VCF_FILTERING → realigned_filtered_rna_vcf                        │
│          ↓                                                                  │
│   SECOND_RESCUE_WORKFLOW                                                  │
│     ├─ VCF_RESCUE_WORKFLOW (DNA + realigned RNA)                        │
│     └─ VCF_RESCUE_POST_PROCESSING → second_rescued_vcf                   │
│                                                                              │
└─────────────────────────────────────────────────────────────────────────────┘
```

### 2.4 VCF vs MAF Workflow Comparison

| Feature | MAF Mode | VCF Mode |
|---------|----------|----------|
| Input Source | Filtered MAF (all variants) | RNA consensus VCF only |
| Target Regions | All filtered variants | RNA consensus variants |
| DNA Consensus | Once | Once |
| RNA Consensus | Once | **Twice** (before/after) |
| Rescue Steps | Once | **Twice** (before/after) |
| Output Format | MAF + VCF | **VCF-only** |
| MAF_FILTERING_RNA | Yes | **No** (VCF-only) |
| Computational Load | Higher (all variants) | Lower (targeted) |

---

## 3. Detailed Workflow Description

### 3.1 Phase 1: First Round Processing

The first round is the standard pipeline execution that produces initial variant calls from both DNA and RNA samples.

**Input**: Raw FASTQ files

**Steps**:
1. **Alignment**: STAR for RNA, BWA/BWA-mem2 for DNA
2. **GATK Preprocessing**: MarkDuplicates, SplitNCigarReads, BQSR
3. **Variant Calling**: Mutect2, Strelka, SAGE (based on `--tools`)
4. **VCF Normalization**: Left-align and normalize variants
5. **VEP Annotation**: Add functional annotations
6. **VCF Consensus**: Combine calls from multiple callers
7. **Cross-Modality Rescue**: DNA + RNA rescue
8. **VCF Filtering**: Apply filtering criteria

**Key Outputs**:
- `filtered_vcf`: Filtered VCF files (DNA and RNA)
- `vcf_rescue`: First-round rescued VCF
- `filtered_maf`: Filtered MAF files

**Stored for Second Rescue**:
- `first_round_dna_consensus_vcf`: DNA consensus (status <= 1)

### 3.2 Phase 2: Realignment (VCF Mode)

Realignment is an optional second pass that targets RNA-specific variants with improved alignment.

**Input**: RNA consensus VCF (from first round)

**Steps**:
1. **VCF2BED**: Convert VCF positions to BED format for read extraction
2. **SAMTOOLS_EXTRACT_READ_IDS**: Extract read IDs overlapping target regions
3. **CRAM to BAM Conversion**: Convert CRAM to BAM for filtering
4. **PICARD_FILTERSAMREADS**: Filter BAM to extracted reads only
5. **BAM to FASTQ**: Convert filtered BAM to FASTQ
6. **HISAT2 Alignment**: Re-align reads with HISAT2 splice-aware aligner

**Key Outputs**:
- Realigned BAM files (RNA tumor + DNA normal)

### 3.3 Phase 3: RNA Realignment Workflow

Full variant calling pipeline on realigned BAMs.

**Steps**:
1. **BAM_GATK_PREPROCESSING**: Preprocess realigned BAM
2. **BAM_VARIANT_CALLING**: Call variants on realigned data
3. **VCF_NORMALIZE**: Normalize VCF files
4. **VCF_ANNOTATE**: Annotate with VEP
5. **VCF_CONSENSUS_WORKFLOW**: Generate RNA consensus from realigned data
6. **VCF_FILTERING**: Filter RNA consensus VCF

**Key Outputs**:
- `realigned_rna_consensus_vcf`: RNA consensus from realigned BAM
- `realigned_filtered_rna_vcf`: Filtered RNA VCF

### 3.4 Phase 4: Second Rescue

Combine DNA consensus (first round) with realigned RNA consensus for cross-modality rescue.

**Steps**:
1. **VCF_RESCUE_WORKFLOW**: Aggregate DNA + realigned RNA variants
2. **VCF_RESCUE_POST_PROCESSING**: Apply full annotation and filtering
   - RNA editing annotation (REDIportal)
   - COSMIC database annotation
   - gnomAD population frequencies
   - Rescue-specific filtering

**Key Outputs**:
- `second_rescued_vcf`: Final rescued VCF with full annotation

---

## 4. Configuration Guide

### 4.1 Required Parameters

None required. VCF realignment uses the same parameters as standard pipeline execution.

### 4.2 Optional Parameters

| Parameter | Default | Description |
|----------|---------|-------------|
| `--realignment_mode` | `vcf` | Realignment mode: 'vcf' (new) or 'maf' (legacy) |

### 4.3 Custom Configuration Examples

**Example 1: Basic VCF Realignment**

```bash
nextflow run nf-core/rnadnavar -profile docker \
  --input samplesheet.csv \
  --outdir results/ \
  --tools sage,strelka,mutect2,consensus,rescue,realignment
```

**Example 2: Use MAF Mode (Legacy)**

```bash
nextflow run nf-core/rnadnavar -profile docker \
  --input samplesheet.csv \
  --outdir results/ \
  --tools sage,strelka,mutect2,consensus,rescue,realignment \
  --realignment_mode maf
```

**Example 3: Custom Config File**

Create `custom.config`:
```groovy
params {
    realignment_mode = 'vcf'
    save_align_intermeds = true  // Save intermediate realignment files
    publish_dir_mode = 'copy'
}
```

Run with custom config:
```bash
nextflow run nf-core/rnadnavar -profile docker \
  -c custom.config \
  --input samplesheet.csv \
  --tools sage,strelka,mutect2,consensus,rescue,realignment
```

### 4.4 Parameter Reference

#### Global Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--tools` | string | `null` | Comma-separated tools to enable |
| `--realignment_mode` | string | `vcf` | Realignment workflow mode |

#### Database Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--rediportal_vcf` | path | `null` | REDIportal database VCF for RNA editing annotation |
| `--cosmic_database` | path | `null` | COSMIC database VCF |
| `--gnomad_database` | path | `null` | gnomAD database directory |

#### Rescue Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--rescue_snv_thr` | int | `2` | Minimum callers for SNV consensus |
| `--rescue_indel_thr` | int | `2` | Minimum callers for indel consensus |
| `--min_rna_support` | int | `2` | Minimum RNA support threshold |

---

## 5. Output Files

### 5.1 File Naming Conventions

```
{sample_id}.{extension}.{type}

Examples:
- SAMPLE001.consensus.vcf.gz
- SAMPLE001.consensus.filtered.vcf.gz
- SAMPLE001_rescued.vcf.gz
- SAMPLE001_realign.consensus.vcf.gz
- SAMPLE001_realign_rescued.vcf.gz
```

### 5.2 Output Directory Structure

```
results/
├── preprocessing/
│   └── realignment/
│       ├── vcf2bed/
│       │   └── {sample}/
│       │       └── {sample}.bed
│       └── {sample}/
│           ├── {sample}_filtered.bam
│           └── {sample}_filtered.bai
├── vcf/
│   ├── {sample}.consensus.vcf.gz         # First round consensus
│   ├── {sample}.consensus.filtered.vcf.gz
│   ├── {sample}.rescued.vcf.gz            # First round rescue
│   └── {sample}.rescued.filtered.vcf.gz
├── realignment/
│   ├── {sample}_realign.consensus.vcf.gz         # Realigned RNA consensus
│   ├── {sample}_realign.consensus.filtered.vcf.gz
│   └── {sample}_realign_rescued.vcf.gz             # Second rescue
├── rescue/
│   └── {sample}_realign_rescued.rna_editing.csv
├── maf/
│   └── {sample}.maf
├── pipeline_info/
│   └── nf_core_rnadnavar_software_mqc_versions.yml
└── multiqc_report.html
```

### 5.3 File Format Specifications

**VCF Format**: Standard VCF v4.2 with custom INFO fields

| INFO Field | Description |
|------------|-------------|
| `DNA` | Number of DNA callers supporting variant |
| `RNA` | Number of RNA callers supporting variant |
| `MODALITY` | Source modality (DNA, RNA, both) |
| `RESCUED` | Whether variant was rescued (true/false) |
| `CROSS_MODALITY` | Support from both DNA and RNA (true/false) |
| `RNA_EDITING` | RNA editing annotation (if enabled) |
| `COSMIC_COUNT` | COSMIC database occurrence count |
| `GNOMAD_AF` | gnomAD allele frequency |

---

## 6. Usage Examples

### 6.1 Basic Usage

```bash
# VCF-based realignment (new default)
nextflow run nf-core/rnadnavar -profile docker \
  --input samplesheet.csv \
  --outdir results/ \
  --tools sage,strelka,mutect2,consensus,rescue,realignment

# Explicitly specify VCF mode
nextflow run nf-core/rnadnavar -profile docker \
  --input samplesheet.csv \
  --outdir results/ \
  --tools sage,strelka,mutect2,consensus,rescue,realignment \
  --realignment_mode vcf
```

### 6.2 Advanced Configurations

**With RNA Editing Annotation**:

```bash
nextflow run nf-core/rnadnavar -profile docker \
  --input samplesheet.csv \
  --outdir results/ \
  --tools sage,strelka,mutect2,consensus,rescue,realignment \
  --enable_rna_annotation \
  --rediportal_vcf /path/to/rediportal.vcf.gz \
  --rediportal_tbi /path/to/rediportal.vcf.gz.tbi
```

**With COSMIC/gnomAD Annotation**:

```bash
nextflow run nf-core/rnadnavar -profile docker \
  --input samplesheet.csv \
  --outdir results/ \
  --tools sage,strelka,mutect2,consensus,rescue,realignment \
  --enable_cosmic_gnomad_annotation \
  --cosmic_database /path/to/cosmic.vcf.gz \
  --gnomad_database /path/to/gnomad/
```

### 6.3 Common Use Cases

**Use Case 1: Targeted RNA Variant Re-alignment**

```bash
# For samples where RNA alignment quality is suboptimal
# VCF mode targets only RNA consensus regions for re-alignment
nextflow run nf-core/rnadnavar -profile docker \
  --input samplesheet.csv \
  --tools mutect2,strelka,consensus,rescue,realignment
```

**Use Case 2: Comparison Testing**

```bash
# Run VCF mode
nextflow run nf-core/rnadnavar -profile docker \
  --input samplesheet.csv \
  --tools mutect2,strelka,consensus,rescue,realignment \
  --realignment_mode vcf \
  --outdir results_vcf/

# Run MAF mode for comparison
nextflow run nf-core/rnadnavar -profile docker \
  --input samplesheet.csv \
  --tools mutect2,strelka,consensus,rescue,realignment \
  --realignment_mode maf \
  --outdir results_maf/
```

**Use Case 3: Resume from Checkpoint**

```bash
# Resume if pipeline was interrupted
nextflow run nf-core/rnadnavar -profile docker \
  -resume \
  --input samplesheet.csv \
  --tools sage,strelka,mutect2,consensus,rescue,realignment
```

---

## 7. Troubleshooting

### 7.1 Common Issues

#### Issue: "No reads extracted during realignment"

**Symptoms**: VCF2BED produces a BED file, but no reads are extracted.

**Possible Causes**:
1. RNA consensus VCF is empty
2. BED regions don't overlap with any reads in the CRAM file
3. Sample status codes are incorrect

**Solutions**:
```bash
# Check if VCF has variants
bcftools view -i 'TYPE="snp" results/vcf/*.consensus.vcf.gz | wc -l

# Verify sample status
head -1 results/vcf/*.consensus.vcf.gz | grep "^#CHROM"

# Check realignment mode in log
grep "realignment_mode" work/.nextflow.log
```

#### Issue: "Second rescue not running"

**Symptoms**: `second_rescued_vcf` output is missing.

**Possible Causes**:
1. Rescue tool not in `--tools` list
2. Realignment didn't complete successfully
3. No DNA samples available for pairing

**Solutions**:
```bash
# Verify rescue is enabled
echo $tools | grep rescue

# Check realignment workflow status
grep "RNA_REALIGNMENT" work/.nextflow.log

# Check for DNA normal samples
awk -F',' '$8==0' samplesheet.csv
```

#### Issue: "MAF_FILTERING_RNA error in VCF mode"

**Note**: This is expected behavior. In VCF mode, MAF filtering is skipped since the workflow is VCF-only.

**If you need MAF output**: Use `--realignment_mode maf` instead.

### 7.2 Debug Mode

Enable debug logging for detailed troubleshooting:

```bash
nextflow run nf-core/rnadnavar -profile docker,test \
  --input samplesheet.csv \
  --tools sage,strelka,mutect2,consensus,rescue,realignment \
  --realignment_mode vcf \
  -dump-channels
```

### 7.3 Log File Locations

| Log File | Location |
|----------|----------|
| Nextflow log | `work/.nextflow.log` |
| Trace log | `work/trace_*.txt` |
| Command log | `work/.command.log` |
| Module versions | `results/pipeline_info/nf_core_rnadnavar_software_mqc_versions.yml` |

---

## 8. Performance Considerations

### 8.1 Computational Requirements

| Resource | Recommended | Minimum |
|----------|-------------|----------|
| CPU | 16+ cores | 8 cores |
| Memory | 64 GB | 32 GB |
| Disk | 500 GB | 200 GB |
| Runtime | 6-12 hours | Varies by data size |

### 8.2 Runtime Estimates

| Stage | Typical Runtime |
|-------|----------------|
| First Round (standard) | 4-8 hours |
| Realignment (VCF extraction) | 10-30 minutes |
| HISAT2 Alignment | 1-2 hours |
| RNA Realignment Workflow | 2-3 hours |
| Second Rescue | 30-60 minutes |

**Total with Realignment**: 8-14 hours

### 8.3 Resource Optimization

**Reducing Runtime**:

1. **Use Interval Files**: Target specific genomic regions
2. **Reduce Caller Set**: Use fewer variant callers
3. **Adjust Threads**: Increase `--max_cpus` for alignment steps

**Reducing Disk Usage**:

1. **Disable Intermediate Saving**: Don't set `--save_align_intermeds`
2. **Compress Outputs**: Ensure `publish_dir_mode = 'copy'` not 'symlink'

---

## 9. Comparison with MAF Mode

### 9.1 Feature Comparison Table

| Feature | MAF Mode | VCF Mode | Recommendation |
|---------|----------|----------|----------------|
| Input Source | All filtered variants | RNA consensus only | VCF: Faster |
| Target Region Size | Larger | Smaller (targeted) | VCF: Efficient |
| DNA Consensus | Once | Once | Same |
| RNA Consensus | Once | Twice | VCF: More thorough |
| Rescue Steps | Once | Twice | VCF: More comprehensive |
| Output Format | MAF + VCF | VCF-only | MAF: For MAF workflows |
| MAF_FILTERING_RNA | Yes | No | MAF: For MAF output |
| Computational Cost | Higher | Lower | VCF: More efficient |
| Migration Path | N/A | New feature | VCF: Recommended |

### 9.2 When to Use Each Mode

**Use VCF Mode** (Recommended):
- Default for new analyses
- When computational efficiency is important
- When targeted realignment is sufficient
- For RNA-specific variant discovery
- When VCF output is primary requirement

**Use MAF Mode**:
- For backward compatibility with existing workflows
- When MAF output is required
- When re-aligning all filtered variants is desired
- For comparison with historical analyses

### 9.3 Migration Guide

**From MAF to VCF Mode**:

1. **Change default mode** (already done in code):
   ```groovy
   // In nextflow.config
   realignment_mode = 'vcf'
   ```

2. **Update workflows** that depend on MAF output:
   - Change input from `maf/` to `vcf/`
   - Update file naming conventions
   - Verify VCF parsing handles new fields

3. **Update documentation**:
   - Note VCF mode is now default
   - Update MAF mode to "legacy" status

---

## 10. Appendix

### 10.1 Module Reference

#### VCF2BED

**Location**: `modules/local/vcf2bed/main.nf`

**Purpose**: Convert VCF variant positions to BED format for read extraction

**Input**: VCF file + TBI index

**Output**: BED file (chromosome, start, end)

**Key Command**:
```bash
bcftools view -O v -o - ${vcf} | \
awk 'BEGIN{OFS="\t"} !/^#/ {
    start = $2 - 1
    end = $2 + length($4) - 1
    print $1, start, end
}' > ${prefix}.bed
```

#### PREPARE_REALIGNMENT_VCF

**Location**: `subworkflows/local/prepare_realignment_vcf/main.nf`

**Purpose**: Extract reads and re-align with HISAT2 using VCF input

**Input**: RNA consensus VCF, CRAM files

**Output**: Realigned BAM files

#### RNA_REALIGNMENT_WORKFLOW

**Location**: `subworkflows/local/rna_realignment/main.nf`

**Purpose**: Full variant calling pipeline on realigned BAMs

**Input**: Realigned BAM, reference files

**Output**: RNA consensus VCF, filtered VCF

#### SECOND_RESCUE_WORKFLOW

**Location**: `subworkflows/local/second_rescue/main.nf`

**Purpose**: Cross-modality rescue combining DNA (first round) + realigned RNA

**Input**: DNA consensus VCF, realigned RNA consensus VCF

**Output**: Second rescued VCF with full annotation

### 10.2 Error Codes

| Error | Cause | Solution |
|-------|-------|----------|
| `VCF_FILE_NOT_FOUND` | RNA consensus VCF missing | Check first round completed successfully |
| `NO_EXTRACTED_READS` | BED regions don't overlap reads | Check VCF has variants in regions with coverage |
| `RESCUE_FAILED` | No DNA samples available for pairing | Verify samplesheet includes DNA normal (status=0) |
| `HISAT2_ALIGN_FAILED` | HISAT2 index missing | Ensure `--build_only_index` or provide `--hisat2_index` |

### 10.3 Glossary

| Term | Definition |
|------|------------|
| **Cross-modality rescue** | Combining variants from DNA and RNA data sources to increase confidence |
| **Somatic variant calling** | Identifying variants present in tumor but not in normal tissue |
| **Paired-end sequencing** | Both ends of DNA fragments are sequenced |
| **HISAT2** | Splice-aware aligner for RNA-seq data |
| **VCF normalization** | Left-aligning and decomposing multi-allelic sites |
| **VEP annotation** | Adding functional consequences to variants |
| **REDIportal** | Database of known RNA editing events |
| **COSMIC** | Catalogue of somatic mutations in cancer |
| **gnomAD** | Genome aggregation database of variant frequencies |

---

## Version History

| Version | Date | Changes |
|---------|------|---------|
| 1.0.0 | 2025-12-24 | Initial VCF-based RNA realignment implementation |

---

## Contact & Support

For questions or issues related to VCF-based RNA realignment:

1. Check the [nf-core/rnadnavar documentation](https://nf-co.re/rnadnavar)
2. Review the [GitHub issues](https://github.com/nf-core/rnadnavar/issues)
3. Run the validation script: `./tests/validate_vcf_realignment.sh`

---

*Document Version: 1.0.0*
*Last Updated: December 24, 2025*
