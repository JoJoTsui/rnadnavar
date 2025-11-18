# Cross-Modality Rescue Workflow

## Overview

The rescue workflow is a cross-modality variant aggregation feature that identifies variants supported by evidence from both DNA and RNA sequencing data. It operates after within-modality consensus calling and provides an additional layer of variant validation by leveraging complementary information from different sequencing modalities.

## Purpose

Traditional consensus calling operates within a single modality (DNA or RNA), combining results from multiple variant callers. The rescue workflow extends this by:

1. **Cross-Modality Validation**: Identifying variants detected in both DNA and RNA samples
2. **Evidence Aggregation**: Combining variant information across modalities with proper tagging
3. **Rescue Statistics**: Providing metrics on cross-modality support and rescue effectiveness

## When to Use Rescue

The rescue workflow is beneficial when:

- You have paired DNA and RNA samples from the same patient
- You want to identify variants with cross-modality support
- You need to validate RNA-specific variants with DNA evidence
- You want to rescue variants that may not meet consensus thresholds within a single modality but have support across modalities

## Workflow Architecture

```
┌─────────────────────────────────────────────────────────────┐
│                    Variant Calling                          │
│  DNA: Mutect2, Strelka, DeepSomatic                        │
│  RNA: Mutect2, Strelka, DeepSomatic                        │
└────────────────┬────────────────────────────────────────────┘
                 │
                 ▼
┌─────────────────────────────────────────────────────────────┐
│              Within-Modality Consensus                      │
│  DNA Consensus: Aggregate DNA callers (status ≤ 1)         │
│  RNA Consensus: Aggregate RNA callers (status = 2)         │
└────────────────┬────────────────────────────────────────────┘
                 │
                 ▼
┌─────────────────────────────────────────────────────────────┐
│              Cross-Modality Rescue                          │
│  - Aggregate DNA + RNA consensus variants                   │
│  - Tag with modality information                            │
│  - Mark rescued variants (present in both)                  │
│  - Generate cross-modality statistics                       │
└────────────────┬────────────────────────────────────────────┘
                 │
                 ▼
┌─────────────────────────────────────────────────────────────┐
│                    Filtering                                │
└─────────────────────────────────────────────────────────────┘
```


## Enabling Rescue in the Pipeline

### Configuration

Add `rescue` to the `--tools` parameter:

```bash
nextflow run nf-core/rnadnavar \
    --input samplesheet.csv \
    --outdir results \
    --tools mutect2,strelka,deepsomatic,consensus,rescue \
    --genome GRCh38 \
    -profile docker
```

### Requirements

For rescue to run, you need:

1. **Paired DNA and RNA samples** with the same `patient` ID in the samplesheet
2. **DNA samples** with `status` 0 (normal) or 1 (tumor)
3. **RNA samples** with `status` 2 (tumor RNA)
4. **Consensus enabled** in tools (rescue runs after consensus)

### Samplesheet Example

```csv
patient,status,sample,lane,fastq_1,fastq_2
patient1,0,DNA_NORMAL,L001,/path/to/dna_normal_R1.fastq.gz,/path/to/dna_normal_R2.fastq.gz
patient1,1,DNA_TUMOR,L001,/path/to/dna_tumor_R1.fastq.gz,/path/to/dna_tumor_R2.fastq.gz
patient1,2,RNA_TUMOR,L001,/path/to/rna_tumor_R1.fastq.gz,/path/to/rna_tumor_R2.fastq.gz
```

In this example:
- DNA consensus will aggregate variants from DNA_TUMOR (status=1)
- RNA consensus will aggregate variants from RNA_TUMOR (status=2)
- Rescue will cross DNA and RNA consensus results for patient1


## Rescue Algorithm

### Step 1: Modality Separation

The pipeline automatically separates samples by modality based on the `status` field:

- **DNA modality**: status ≤ 1 (normal DNA = 0, tumor DNA = 1)
- **RNA modality**: status = 2 (tumor RNA)

### Step 2: Within-Modality Consensus

Before rescue, consensus calling is performed within each modality:

**DNA Consensus:**
```
DNA_TUMOR.mutect2.vcf  ┐
DNA_TUMOR.strelka.vcf  ├─→ DNA_TUMOR.consensus.vcf
DNA_TUMOR.deepsomatic.vcf ┘
```

**RNA Consensus:**
```
RNA_TUMOR.mutect2.vcf  ┐
RNA_TUMOR.strelka.vcf  ├─→ RNA_TUMOR.consensus.vcf
RNA_TUMOR.deepsomatic.vcf ┘
```

### Step 3: Cross-Modality Aggregation

The rescue script (`run_rescue_vcf.py`) aggregates variants across modalities:

1. **Read DNA consensus VCF** with modality tag 'DNA'
2. **Read RNA consensus VCF** with modality tag 'RNA'
3. **Read individual caller VCFs** from both modalities (optional, for detailed statistics)
4. **Aggregate all variants** using the same variant key (chrom:pos:ref:alt)
5. **Tag variants** with modality information
6. **Mark rescued variants** (present in both DNA and RNA)

### Step 4: Modality Tagging

Each variant is tagged with:

- **MODALITIES**: Which modalities detected the variant (DNA, RNA, or both)
- **CALLERS_BY_MODALITY**: Callers grouped by modality (e.g., "DNA:mutect2,strelka|RNA:mutect2")
- **DNA_SUPPORT**: Number of DNA callers supporting the variant
- **RNA_SUPPORT**: Number of RNA callers supporting the variant
- **CROSS_MODALITY**: YES if variant has both DNA and RNA support
- **RESCUED**: YES if variant is present in both DNA and RNA consensus

### Step 5: Statistics Generation

The rescue workflow generates comprehensive statistics:

- Total variants in rescue output
- DNA-only variants (detected only in DNA)
- RNA-only variants (detected only in RNA)
- Cross-modality variants (detected in both)
- Number of rescued variants
- Rescue rate (percentage of variants rescued)
- SNV vs indel breakdown


## Output Files

### Rescue VCF

**Location**: `results/variant_calling/vcf_rescue_workflow/{patient}/`

**Filename**: `{DNA_sample}_rescued_{RNA_sample}.rescued.vcf.gz`

**Content**: All variants from DNA and RNA consensus with modality annotations

### Modality Prefix Format

**Important**: In rescue VCFs, all caller names are prefixed with their modality (DNA_ or RNA_) to clearly distinguish the source of evidence. This applies to all caller-related INFO fields.

**Example INFO fields**:
```
##INFO=<ID=N_CALLERS,Number=1,Type=Integer,Description="Total number of variant callers (excludes consensus)">
##INFO=<ID=N_DNA_CALLERS,Number=1,Type=Integer,Description="Number of DNA variant callers">
##INFO=<ID=N_RNA_CALLERS,Number=1,Type=Integer,Description="Number of RNA variant callers">
##INFO=<ID=CALLERS,Number=.,Type=String,Description="List of all callers with modality prefix">
##INFO=<ID=CALLERS_SUPPORT,Number=.,Type=String,Description="Supporting callers with modality prefix">
##INFO=<ID=FILTERS_ORIGINAL,Number=.,Type=String,Description="Original filters with modality prefix (format: MODALITY_caller:filter|...)">
##INFO=<ID=GT_BY_CALLER,Number=.,Type=String,Description="Genotypes with modality prefix (format: MODALITY_caller:GT|...)">
##INFO=<ID=DP_BY_CALLER,Number=.,Type=String,Description="Depth with modality prefix (format: MODALITY_caller:DP|...)">
##INFO=<ID=VAF_BY_CALLER,Number=.,Type=String,Description="VAF with modality prefix (format: MODALITY_caller:VAF|...)">
##INFO=<ID=MODALITIES,Number=.,Type=String,Description="Modalities where variant was detected">
##INFO=<ID=CALLERS_BY_MODALITY,Number=.,Type=String,Description="Callers grouped by modality">
##INFO=<ID=DNA_SUPPORT,Number=1,Type=Integer,Description="Number of DNA callers supporting this variant">
##INFO=<ID=RNA_SUPPORT,Number=1,Type=Integer,Description="Number of RNA callers supporting this variant">
##INFO=<ID=CROSS_MODALITY,Number=1,Type=String,Description="Whether variant has cross-modality support">
##INFO=<ID=RESCUED,Number=1,Type=String,Description="Variant included via cross-modality consensus">
##INFO=<ID=DP_DNA_MEAN,Number=1,Type=Float,Description="Mean depth across DNA callers">
##INFO=<ID=DP_RNA_MEAN,Number=1,Type=Float,Description="Mean depth across RNA callers">
##INFO=<ID=VAF_DNA_MEAN,Number=1,Type=Float,Description="Mean VAF across DNA callers">
##INFO=<ID=VAF_RNA_MEAN,Number=1,Type=Float,Description="Mean VAF across RNA callers">
```

### Example Variant Record

```
chr1  12345  .  A  G  .  PASS  N_CALLERS=3;N_DNA_CALLERS=2;N_RNA_CALLERS=1;CALLERS=DNA_mutect2|DNA_strelka|RNA_deepsomatic;CALLERS_SUPPORT=DNA_mutect2|RNA_deepsomatic;FILTERS_ORIGINAL=DNA_mutect2:PASS|RNA_deepsomatic:LowQuality;GT_BY_CALLER=DNA_mutect2:0/1|RNA_deepsomatic:0/1;DP_BY_CALLER=DNA_mutect2:100|RNA_deepsomatic:50;VAF_BY_CALLER=DNA_mutect2:0.45|RNA_deepsomatic:0.38;MODALITIES=DNA|RNA;CALLERS_BY_MODALITY=DNA:mutect2,strelka|RNA:deepsomatic;DNA_SUPPORT=2;RNA_SUPPORT=1;CROSS_MODALITY=YES;RESCUED=YES;DP_DNA_MEAN=150.5;DP_RNA_MEAN=85.0;VAF_DNA_MEAN=0.25;VAF_RNA_MEAN=0.28
```

This variant:
- Was detected by 3 variant callers total (N_CALLERS=3)
- Has 2 DNA callers (DNA_mutect2, DNA_strelka) and 1 RNA caller (RNA_deepsomatic)
- Was detected in both DNA and RNA modalities
- Has support from 2 DNA callers (mutect2, strelka) and 1 RNA caller (deepsomatic)
- Is marked as rescued (RESCUED=YES)
- Has cross-modality support (CROSS_MODALITY=YES)
- Shows similar VAF in DNA (0.25) and RNA (0.28)
- All caller names are prefixed with their modality (DNA_ or RNA_)


## Interpreting Rescue Results

### Variant Categories

Rescued VCF files contain three categories of variants:

1. **DNA-only variants** (MODALITIES=DNA, CROSS_MODALITY=NO)
   - Detected in DNA consensus but not RNA
   - May be in non-expressed genes or low RNA coverage regions
   - Still valuable for somatic variant analysis

2. **RNA-only variants** (MODALITIES=RNA, CROSS_MODALITY=NO)
   - Detected in RNA consensus but not DNA
   - Could be RNA editing events, artifacts, or low DNA coverage
   - Require careful interpretation

3. **Cross-modality variants** (MODALITIES=DNA|RNA, CROSS_MODALITY=YES, RESCUED=YES)
   - Detected in both DNA and RNA consensus
   - Highest confidence variants
   - Validated across sequencing modalities
   - Likely true somatic mutations in expressed genes

### Filtering Recommendations

For high-confidence variant sets, filter for:

```bash
# Extract only rescued variants (cross-modality support)
bcftools view -i 'INFO/RESCUED="YES"' rescued.vcf.gz > high_confidence.vcf

# Extract variants with minimum support
bcftools view -i 'INFO/DNA_SUPPORT>=2 && INFO/RNA_SUPPORT>=1' rescued.vcf.gz > supported.vcf

# Extract PASS variants with cross-modality support
bcftools view -f PASS -i 'INFO/CROSS_MODALITY="YES"' rescued.vcf.gz > pass_cross_modality.vcf
```

### Statistics Interpretation

Example rescue statistics output:

```
- Statistics:
  - Total variants: 15,234
  - SNVs: 13,500 (rescued: 8,200)
  - Indels: 1,734 (rescued: 950)

- Modality Support:
  - DNA only: 4,500
  - RNA only: 1,584
  - Cross-modality: 9,150

- Rescue Effectiveness:
  - Total rescued: 9,150
  - Rescue rate: 60.1%
```

**Interpretation:**
- 60.1% of variants have cross-modality support (high confidence)
- 4,500 DNA-only variants may be in non-expressed regions
- 1,584 RNA-only variants need careful review
- 9,150 variants validated across both modalities


## Configuration Parameters

### Rescue-Specific Parameters

The rescue workflow can be configured in `conf/modules/consensus/vcf_consensus_workflow.config`:

```groovy
withName: 'VCF_RESCUE' {
    ext.snv_thr = 2        // Minimum callers for SNV consensus
    ext.indel_thr = 2      // Minimum callers for indel consensus
}
```

### Parameter Descriptions

- **snv_thr**: Minimum number of callers required for a SNV to be included in rescue output
  - Default: 2
  - Lower values increase sensitivity but may include more false positives
  - Higher values increase specificity but may miss true variants

- **indel_thr**: Minimum number of callers required for an indel to be included in rescue output
  - Default: 2
  - Indels are typically harder to call consistently
  - Consider using same or lower threshold than SNVs

### Example Custom Configuration

Create a custom config file `rescue_config.config`:

```groovy
process {
    withName: 'VCF_RESCUE' {
        ext.snv_thr = 1        // More sensitive for SNVs
        ext.indel_thr = 2      // Standard for indels
        cpus = 4
        memory = 16.GB
    }
}
```

Run with custom config:

```bash
nextflow run nf-core/rnadnavar \
    --input samplesheet.csv \
    --outdir results \
    --tools mutect2,strelka,deepsomatic,consensus,rescue \
    -c rescue_config.config \
    -profile docker
```


## Use Cases and Examples

### Use Case 1: Validating RNA Variants with DNA Evidence

**Scenario**: You have RNA-seq data with potential somatic variants but want DNA validation.

**Approach**:
1. Run variant calling on both DNA and RNA samples
2. Enable consensus and rescue in tools
3. Filter rescue output for RESCUED=YES variants

**Benefit**: High-confidence variants validated across modalities

### Use Case 2: Identifying Expressed Somatic Mutations

**Scenario**: You want to find somatic mutations that are actually expressed.

**Approach**:
1. Run DNA variant calling to identify all somatic mutations
2. Run RNA variant calling to identify expressed variants
3. Use rescue to find intersection (CROSS_MODALITY=YES)

**Benefit**: Focus on functionally relevant mutations

### Use Case 3: Comprehensive Variant Discovery

**Scenario**: You want maximum sensitivity for variant discovery.

**Approach**:
1. Set low thresholds (snv_thr=1, indel_thr=1)
2. Include all variants from rescue output
3. Use modality tags for downstream filtering

**Benefit**: Capture all variants with any level of support

### Example Workflow

```bash
# Step 1: Run pipeline with rescue
nextflow run nf-core/rnadnavar \
    --input samplesheet.csv \
    --outdir results \
    --tools mutect2,strelka,deepsomatic,consensus,rescue \
    --genome GRCh38 \
    -profile docker

# Step 2: Extract high-confidence rescued variants
bcftools view -i 'INFO/RESCUED="YES" && INFO/DNA_SUPPORT>=2 && INFO/RNA_SUPPORT>=1' \
    results/variant_calling/vcf_rescue_workflow/patient1/DNA_TUMOR_rescued_RNA_TUMOR.rescued.vcf.gz \
    > patient1_high_confidence.vcf

# Step 3: Annotate with VEP (if enabled)
# Annotation is performed automatically if vep is in tools

# Step 4: Convert to MAF for downstream analysis
# Conversion is performed automatically by the pipeline
```


## Technical Implementation

### Nextflow Workflow Structure

The rescue workflow is implemented as a subworkflow:

**File**: `subworkflows/local/vcf_rescue_workflow/main.nf`

**Key Steps**:
1. Cross DNA and RNA consensus VCFs by patient ID
2. Group individual caller VCFs by patient and modality
3. Invoke VCF_RESCUE module with all inputs
4. Output rescued VCF with modality annotations

### Python Script

**File**: `bin/run_rescue_vcf.py`

**Key Functions** (from vcf_utils package):
- `read_variants_from_vcf()`: Read VCF with modality tagging
- `aggregate_variants()`: Combine variants across modalities
- `mark_rescued_variants()`: Identify cross-modality variants
- `tag_variant_with_modality()`: Add modality metadata
- `compute_rescue_statistics()`: Generate statistics
- `write_union_vcf()`: Write output with modality fields

### VCF Utils Package

The rescue workflow uses the `vcf_utils` Python package for core functionality:

```python
from vcf_utils.aggregation import read_variants_from_vcf, aggregate_variants
from vcf_utils.tagging import mark_rescued_variants, tag_variant_with_modality
from vcf_utils.statistics import compute_rescue_statistics, print_statistics
from vcf_utils.io_utils import write_union_vcf
```

See `bin/vcf_utils/README.md` for detailed package documentation.


## Comparison with Consensus Workflow

### Consensus Workflow (Within-Modality)

**Purpose**: Aggregate variants from multiple callers within the same modality

**Input**: Multiple VCF files from different callers (all DNA or all RNA)

**Output**: Single consensus VCF with caller support information

**Example**:
```
DNA_TUMOR.mutect2.vcf  ┐
DNA_TUMOR.strelka.vcf  ├─→ DNA_TUMOR.consensus.vcf
DNA_TUMOR.deepsomatic.vcf ┘
```

**INFO Fields**: N_CALLERS, CALLERS, PASSES_CONSENSUS, DP_MEAN, VAF_MEAN, etc.

### Rescue Workflow (Cross-Modality)

**Purpose**: Aggregate variants across DNA and RNA modalities

**Input**: DNA consensus VCF + RNA consensus VCF (+ optional individual caller VCFs)

**Output**: Single rescued VCF with modality annotations

**Example**:
```
DNA_TUMOR.consensus.vcf ┐
                        ├─→ DNA_TUMOR_rescued_RNA_TUMOR.rescued.vcf
RNA_TUMOR.consensus.vcf ┘
```

**INFO Fields**: All consensus fields + MODALITIES, CALLERS_BY_MODALITY, DNA_SUPPORT, RNA_SUPPORT, CROSS_MODALITY, RESCUED, DP_DNA_MEAN, DP_RNA_MEAN, VAF_DNA_MEAN, VAF_RNA_MEAN

### Key Differences

| Feature | Consensus | Rescue |
|---------|-----------|--------|
| Scope | Within-modality | Cross-modality |
| Input | Multiple callers, same modality | Consensus from each modality |
| Modality tracking | No | Yes |
| Cross-validation | No | Yes |
| Output size | Typically smaller | Typically larger |
| Confidence level | Caller agreement | Modality agreement |


## Best Practices

### 1. Always Run Consensus Before Rescue

Rescue operates on consensus results, not individual caller outputs. Ensure `consensus` is in your tools list before `rescue`:

```bash
--tools mutect2,strelka,deepsomatic,consensus,rescue  # Correct
--tools mutect2,strelka,deepsomatic,rescue            # Incorrect
```

### 2. Use Appropriate Thresholds

- For high-confidence analysis: Use higher thresholds (snv_thr=2, indel_thr=2)
- For discovery analysis: Use lower thresholds (snv_thr=1, indel_thr=1)
- Consider your downstream filtering strategy

### 3. Interpret Modality-Specific Variants Carefully

- **DNA-only variants**: May be in non-expressed genes or low RNA coverage regions
- **RNA-only variants**: Could be RNA editing, artifacts, or low DNA coverage
- **Cross-modality variants**: Highest confidence, validated across modalities

### 4. Leverage Modality-Specific Statistics

Use DP_DNA_MEAN, DP_RNA_MEAN, VAF_DNA_MEAN, VAF_RNA_MEAN to:
- Assess coverage differences between modalities
- Identify potential allelic imbalance
- Filter variants with insufficient coverage

### 5. Combine with Downstream Filtering

Rescue provides comprehensive variant sets. Apply additional filters based on:
- RESCUED status (YES for cross-modality support)
- DNA_SUPPORT and RNA_SUPPORT counts
- FILTER status (PASS vs FAIL)
- Modality-specific depth and VAF thresholds

### 6. Document Your Analysis

When publishing results using rescue workflow:
- Report rescue rate and cross-modality statistics
- Specify thresholds used (snv_thr, indel_thr)
- Describe filtering criteria applied to rescue output
- Mention vcf_utils package version


## Troubleshooting

### Issue: Rescue workflow not running

**Symptoms**: Pipeline completes but no rescue output

**Possible causes**:
1. `rescue` not in `--tools` parameter
2. `consensus` not in `--tools` parameter (rescue requires consensus)
3. Missing paired DNA and RNA samples with same patient ID
4. Incorrect status values in samplesheet

**Solution**:
- Verify tools parameter: `--tools mutect2,strelka,consensus,rescue`
- Check samplesheet has paired DNA (status ≤ 1) and RNA (status = 2) samples
- Ensure patient IDs match exactly between DNA and RNA samples

### Issue: No rescued variants (RESCUED=NO for all variants)

**Symptoms**: Rescue VCF generated but no variants marked as rescued

**Possible causes**:
1. No overlapping variants between DNA and RNA consensus
2. Different variant calling sensitivity between modalities
3. Low RNA coverage or expression

**Solution**:
- Check DNA and RNA consensus VCFs separately for variant counts
- Review variant calling parameters for each modality
- Consider lowering consensus thresholds
- Verify RNA sample quality and coverage

### Issue: Too many variants in rescue output

**Symptoms**: Rescue VCF much larger than expected

**Possible causes**:
1. Low consensus thresholds (snv_thr=1, indel_thr=1)
2. Including all variants from both modalities

**Solution**:
- Increase thresholds in configuration
- Filter rescue output for RESCUED=YES variants only
- Apply additional quality filters

### Issue: Memory errors during rescue

**Symptoms**: Process fails with out-of-memory error

**Possible causes**:
1. Very large VCF files
2. Insufficient memory allocation

**Solution**:
- Increase memory in process configuration
- Use custom config with higher memory limits:
  ```groovy
  process {
      withName: 'VCF_RESCUE' {
          memory = 32.GB
      }
  }
  ```


## Frequently Asked Questions

### Q: Can I run rescue without consensus?

**A**: No, rescue requires consensus to be run first. The rescue workflow operates on consensus VCF files from each modality.

### Q: What if I only have DNA or only RNA samples?

**A**: Rescue requires paired DNA and RNA samples. If you only have one modality, use the consensus workflow instead.

### Q: How does rescue handle multi-sample patients?

**A**: Rescue crosses DNA and RNA samples by patient ID. If a patient has multiple DNA or RNA samples, all combinations will be processed.

### Q: Can I rescue variants from more than two modalities?

**A**: Currently, rescue is designed for DNA-RNA pairs. For additional modalities, you would need to run rescue multiple times or extend the workflow.

### Q: What's the difference between RESCUED and CROSS_MODALITY?

**A**: 
- **RESCUED=YES**: Variant is present in both DNA and RNA consensus VCFs
- **CROSS_MODALITY=YES**: Variant has support from callers in both modalities (may include individual caller support)

Both typically have the same value, but CROSS_MODALITY provides more granular information.

### Q: How do I cite the rescue workflow?

**A**: Cite the rnadnavar pipeline and mention the rescue workflow feature:

```
We used the nf-core/rnadnavar pipeline (version X.X.X) with the cross-modality 
rescue workflow to identify variants with support from both DNA and RNA sequencing 
data. Rescue analysis was performed using the vcf_utils package (version 1.0.0).
```

### Q: Can I use rescue with other variant callers?

**A**: Yes, as long as the callers produce standard VCF output. You may need to adjust the pipeline configuration to include additional callers.

### Q: What's the computational cost of rescue?

**A**: Rescue is relatively lightweight compared to variant calling. It typically requires:
- CPU: 1-4 cores
- Memory: 8-16 GB
- Time: Minutes to tens of minutes depending on variant count


## References and Further Reading

### Related Documentation

- [Consensus Logic Explained](consensus_logic_explained.md): Details on within-modality consensus calling
- [VCF Utils Package](../bin/vcf_utils/README.md): Technical documentation for the underlying Python package
- [Pipeline Usage](usage.md): General pipeline usage and parameters
- [Pipeline Output](output.md): Description of all pipeline outputs

### Key Concepts

- **Modality**: The type of sequencing data (DNA or RNA)
- **Consensus**: Agreement among multiple variant callers within a modality
- **Rescue**: Cross-modality validation of variants
- **Variant Key**: Unique identifier for a variant (chrom:pos:ref:alt)

### External Resources

- [GATK Best Practices for RNA-seq Variant Calling](https://gatk.broadinstitute.org/hc/en-us/articles/360035531192)
- [VCF Format Specification](https://samtools.github.io/hts-specs/VCFv4.3.pdf)
- [Nextflow Documentation](https://www.nextflow.io/docs/latest/index.html)

## Version History

- **1.0.0**: Initial implementation of cross-modality rescue workflow
  - Support for DNA-RNA variant aggregation
  - Modality tagging and tracking
  - Rescue statistics generation
  - Integration with vcf_utils package

## Support

For questions, issues, or feature requests related to the rescue workflow:

1. Check this documentation and the troubleshooting section
2. Search existing [GitHub issues](https://github.com/nf-core/rnadnavar/issues)
3. Join the [nf-core Slack](https://nf-co.re/join/slack) and ask in the #rnadnavar channel
4. Open a new [GitHub issue](https://github.com/nf-core/rnadnavar/issues/new) with:
   - Pipeline version
   - Command used
   - Error message or unexpected behavior
   - Relevant log files
