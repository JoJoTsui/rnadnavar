# VCF Statistics Analysis - Comprehensive Plan & Implementation

## Overview

This notebook provides a complete solution for analyzing VCF files from a multi-tool, multi-modal variant calling pipeline with DNA and RNA sequencing data. The system now supports both **standard** and **realignment** workflows, enabling comprehensive comparison of RNA variant calling with and without realignment.

## Architecture & Design

### 1. **Data Pipeline Structure**

Your pipeline has the following hierarchy:

```
Variant Calling (Per-tool, per-modality)
    ├── DeepSomatic (DNA & RNA)
    ├── Mutect2 (DNA & RNA)
    └── Strelka (DNA & RNA)
           ↓
Normalization (bcftools norm)
           ↓
Annotation (VEP)
           ↓
Consensus (Within modality)
    ├── DNA Consensus (DeepSomatic + Mutect2 + Strelka on DNA)
    └── RNA Consensus (DeepSomatic + Mutect2 + Strelka on RNA)
           ↓
Rescue (Cross-modality)
    └── DNA rescued by RNA
           ↓
Annotation Stages
    ├── COSMIC/gnomAD annotation
    ├── RNA editing annotation
    └── Final filtering
```

### 1.1 **Realignment Workflow Support**

The pipeline now includes an optional **realignment workflow** that applies **only to RNA samples**:

```
Standard RNA Workflow:
    RNA reads → Alignment → Variant calling → Consensus → Rescue

Realignment Workflow (RNA only):
    DNA consensus variants → Extract RNA reads → Realign → Variant calling → Consensus → Rescue
```

**Key Points:**
- **Realignment only applies to RNA-tumor samples** (not DNA)
- DNA samples are analyzed only in the standard workflow
- RNA samples are analyzed in both standard and realignment workflows
- Realignment workflow follows the same stages as standard: normalized → consensus → rescue → annotation stages
- Output structure: `vcf_realignment/` directory parallel to standard workflow directories

**Directory Structure:**
```
pipeline_output/
├── normalized/                    # Standard workflow
│   ├── deepsomatic/
│   │   ├── DNA_TUMOR_vs_DNA_NORMAL/
│   │   └── RNA_TUMOR_vs_DNA_NORMAL/
│   └── mutect2/
│       ├── DNA_TUMOR_vs_DNA_NORMAL/
│       └── RNA_TUMOR_vs_DNA_NORMAL/
├── consensus/
│   ├── DNA_TUMOR_vs_DNA_NORMAL/
│   └── RNA_TUMOR_vs_DNA_NORMAL/
├── rescue/
├── annotation/
│   ├── cosmic_gnomad/
│   └── rna_editing/
├── filtered/
└── vcf_realignment/               # Realignment workflow (RNA only)
    ├── normalized/
    │   ├── deepsomatic/
    │   │   └── RNA_TUMOR_realign_vs_DNA_NORMAL/
    │   └── mutect2/
    │       └── RNA_TUMOR_realign_vs_DNA_NORMAL/
    ├── consensus/
    │   └── RNA_TUMOR_realign_vs_DNA_NORMAL/
    ├── rescue/
    ├── annotation/
    │   ├── cosmic_gnomad/
    │   └── rna_editing/
    └── filtered/
```

### 2. **Implementation Components**

#### A. **WorkflowManager Class** (NEW)
- **Purpose**: Automatically detect and manage standard and realignment workflows
- **Capabilities**:
  - Detect which workflows are present (standard, realignment, or both)
  - Provide unified configuration for both workflow types
  - Enable seamless workflow-aware file discovery
- **Key Methods**:
  - `detect_workflows()`: Identify available workflows
  - `get_workflow_config()`: Get configuration for a specific workflow
  - `get_all_configs()`: Get all workflow configurations

#### B. **VCFFileDiscovery Class** (ENHANCED)
- **Purpose**: Automatically discovers all VCF files in the pipeline output
- **Categories**: variant_calling, normalized, annotated, consensus, rescue, filtered
- **Output**: Organized dictionary of all VCF file paths
- **New Features**:
  - Workflow-aware discovery for both standard and realignment
  - `discover_vcfs_for_workflow()`: Discover files for a specific workflow
  - `discover_all_workflows()`: Discover files for all detected workflows
  - `discover_bam_files_for_workflow()`: Discover alignment files per workflow

#### C. **WorkflowComparator Class** (NEW)
- **Purpose**: Compare statistics between standard and realignment workflows
- **Focus**: RNA-modality comparisons (realignment only applies to RNA)
- **Key Methods**:
  - `compare_rna_variant_counts()`: Compare RNA variant counts between workflows
  - `compare_rna_category_distribution()`: Compare RNA category percentages
  - `compare_rna_annotation_stages()`: Compare annotation stage results
  - `create_integrative_view()`: Combine DNA + RNA standard + RNA realignment
  - `calculate_realignment_impact()`: Quantify realignment effectiveness
  - `export_comparison_report()`: Export comprehensive comparison results

#### D. **VCFStatisticsExtractor Class**
- **Purpose**: Extract comprehensive statistics using cyvcf2
- **Methods**:
  - `extract_basic_stats()`: Variant counts, types, filter status
  - `extract_info_fields()`: INFO field statistics (DP, AF, TLOD, etc.)
  - `extract_format_fields()`: Sample-level FORMAT statistics
- **Key Metrics**:
  - Total variants, SNPs, INDELs, complex variants
  - Quality score distributions
  - Chromosome distribution
  - Pass/filter rates
  - Tool-specific INFO fields

#### D. **VCFStatisticsExtractor Class**
- **Purpose**: Extract comprehensive statistics using cyvcf2
- **Methods**:
  - `extract_basic_stats()`: Variant counts, types, filter status
  - `extract_info_fields()`: INFO field statistics (DP, AF, TLOD, etc.)
  - `extract_format_fields()`: Sample-level FORMAT statistics
- **Key Metrics**:
  - Total variants, SNPs, INDELs, complex variants
  - Quality score distributions
  - Chromosome distribution
  - Pass/filter rates
  - Tool-specific INFO fields
- **Workflow Support**: Works identically for both standard and realignment VCFs

#### E. **BAMValidator Class** (ENHANCED)
- **Purpose**: Validate variants using pysam against alignment files
- **Capabilities**:
  - Read support validation (ref/alt counts)
  - Variant allele frequency (VAF) calculation
  - Quality filtering (min base quality 20)
  - Supports both BAM and CRAM formats
- **Features**:
  - Handles CRAM files with reference genome
  - Pileup-based variant validation
  - Per-sample validation metrics
- **New Features**:
  - **Comprehensive 4-sample validation**: DNA_TUMOR, DNA_NORMAL, RNA_TUMOR, RNA_TUMOR_realign
  - Compare VAF between RNA standard and RNA realignment
  - Identify variants with improved support in realignment
  - Workflow-aware BAM file discovery

#### F. **StatisticsAggregator Class** (ENHANCED)
- **Purpose**: Aggregate and compare statistics across VCFs
- **Methods**:
  - `create_variant_count_summary()`: Aggregate counts
  - `create_quality_summary()`: Quality metrics
  - `compare_tools_by_modality()`: Tool performance
  - `compare_consensus_to_individual()`: Consensus analysis
  - `create_info_field_summary()`: INFO field analysis

#### F. **StatisticsAggregator Class** (ENHANCED)
- **Purpose**: Aggregate and compare statistics across VCFs
- **Methods**:
  - `create_variant_count_summary()`: Aggregate counts
  - `create_quality_summary()`: Quality metrics
  - `compare_tools_by_modality()`: Tool performance
  - `compare_consensus_to_individual()`: Consensus analysis
  - `create_info_field_summary()`: INFO field analysis
- **New Features**:
  - `create_workflow_comparison_summary()`: Compare standard vs realignment
  - `create_rna_stage_comparison_summary()`: RNA-focused stage comparison
  - Workflow-aware aggregation and reporting

#### G. **VCFVisualizer Class** (ENHANCED)
#### G. **VCFVisualizer Class** (ENHANCED)
- **Purpose**: Create interactive visualizations with Plotly
- **Plots**:
  1. Variant counts by tool (grouped bar chart)
  2. Quality score distributions (box plots)
  3. Variant type distribution (pie charts)
  4. Consensus vs individual tools (bar + line chart)
  5. Filter status (stacked bar chart)
- **New Features**:
  - `plot_rna_workflow_comparison()`: Side-by-side RNA standard vs realignment
  - `plot_rna_stage_progression_comparison()`: RNA variant progression through stages
  - `plot_rna_annotation_impact_comparison()`: Heatmap of annotation differences
  - `plot_integrative_view()`: Comprehensive DNA + RNA standard + RNA realignment view

### 3. **Analysis Workflow**

```
1. Configure & Discover
   ├── Set BASE_DIR, discover all VCFs and alignments
   └── Detect available workflows (standard, realignment, or both)

2. Extract Statistics
   ├── Process all VCF files with cyvcf2
   ├── Extract basic stats, INFO, FORMAT fields
   ├── Process both standard and realignment workflows
   └── Store in organized dictionary

3. BAM Validation (Enhanced)
   ├── Map VCFs to alignment files
   ├── Validate variant support across all 4 samples:
   │   ├── DNA_TUMOR (standard)
   │   ├── DNA_NORMAL (standard)
   │   ├── RNA_TUMOR (standard)
   │   └── RNA_TUMOR_realign (realignment)
   ├── Calculate VAF from pileups
   └── Compare RNA standard vs realignment support

4. Aggregate & Compare
   ├── Create summary tables
   ├── Compare tools within modalities
   ├── Compare consensus to individuals
   ├── Analyze rescue effectiveness
   └── NEW: Compare standard vs realignment workflows (RNA-focused)

5. Visualize
   ├── Generate interactive plots
   ├── Quality distributions
   ├── Performance comparisons
   └── NEW: RNA workflow comparison visualizations

6. Export Results
   ├── CSV files for each summary table
   ├── Validation results
   ├── Summary report
   └── NEW: Workflow comparison reports
```

## Key Features

### 1. **Realignment Workflow Support** (NEW)
- **Automatic workflow detection**: System detects standard, realignment, or both workflows
- **RNA-focused comparisons**: Compare RNA variant calling with and without realignment
- **Integrative view**: Comprehensive analysis of DNA + RNA standard + RNA realignment
- **Realignment impact quantification**: Measure effectiveness of realignment on RNA variant calling
- **Backward compatible**: Works seamlessly with datasets that don't include realignment

### 2. **Multi-Tool Comparison**
### 2. **Multi-Tool Comparison**
- Compare DeepSomatic, Mutect2, and Strelka
- Evaluate sensitivity and specificity
- Identify tool-specific biases

### 3. **Multi-Modal Analysis**
- DNA vs RNA variant detection
- Modality-specific characteristics
- Cross-modality validation

### 4. **Consensus Analysis**
- Agreement across tools
- Retention rates
- High-confidence variant identification

### 5. **Rescue Analysis**
- Cross-modality variant recovery
- RNA-specific variants added to DNA consensus
- Quantify rescue effectiveness

### 6. **Quality Control**
- Quality score distributions
- Filter status tracking
- Pass/fail rates

### 7. **BAM-Level Validation** (ENHANCED)
- Independent verification using alignments
- Read-level support validation
- VAF consistency checks
- **NEW: 4-sample comprehensive validation**:
  - DNA_TUMOR, DNA_NORMAL (standard workflow)
  - RNA_TUMOR (standard workflow)
  - RNA_TUMOR_realign (realignment workflow)
- **NEW: RNA realignment VAF comparison**

## Statistics Extracted

### Basic Statistics
- Total variants, SNPs, INDELs, complex variants
- Pass/filter counts and rates
- Chromosome distribution
- Quality score statistics (mean, median, percentiles)

### INFO Field Statistics
Common fields analyzed:
- **DP**: Total depth
- **AF**: Allele frequency
- **TLOD**: Tumor LOD score (Mutect2)
- **NLOD**: Normal LOD score (Mutect2)
- **MBQ**: Median base quality
- **MMQ**: Median mapping quality
- **MPOS**: Median distance from read end
- **GERMQ**: Germline quality score
- **STRANDQ**: Strand bias quality

### FORMAT Field Statistics
Per-sample metrics:
- **DP**: Sample depth
- **AD**: Allelic depth (ref, alt)
- **AF**: Sample allele frequency
- **GQ**: Genotype quality

## Visualizations

1. **Tool Comparison Charts**
   - Side-by-side variant counts
   - SNP vs INDEL ratios
   - Modality-specific performance

2. **Quality Distributions**
   - Box plots by tool and modality
   - Identify quality cutoff thresholds
   - Compare tool confidence levels

3. **Variant Type Analysis**
   - SNP/INDEL proportions
   - Modality-specific patterns
   - Complex variant identification

4. **Consensus Analysis**
   - Individual tool contributions
   - Retention rates visualization
   - Agreement patterns

5. **Rescue Effectiveness**
   - Before/after rescue comparison
   - Variants added by RNA
   - Cross-modality benefits

6. **Realignment Workflow Comparisons** (NEW)
   - **RNA workflow comparison**: Side-by-side RNA standard vs realignment
   - **RNA stage progression**: Variant count changes through processing stages
   - **RNA annotation impact**: Heatmap showing annotation stage differences
   - **Integrative view**: Comprehensive DNA + RNA standard + RNA realignment visualization

## RNA-Focused Comparison Approach

### Why RNA-Focused?

Realignment **only applies to RNA-tumor samples**, not DNA samples. Therefore, meaningful comparisons focus on RNA modality:

- **Primary comparison**: RNA standard vs RNA realignment
- **DNA context**: DNA results provide baseline (no realignment equivalent)
- **Integrative view**: Combines all three for complete picture

### Comparison Workflow

```python
# 1. Detect workflows
workflow_manager = WorkflowManager(BASE_DIR)
workflows = workflow_manager.detect_workflows()
print(f"Detected workflows: {workflows}")

# 2. Discover files for all workflows
discovery = VCFFileDiscovery(BASE_DIR, workflow_manager)
all_workflows = discovery.discover_all_workflows()

# 3. Extract statistics for both workflows
standard_stats = process_all_vcfs(all_workflows['standard'])
realignment_stats = process_all_vcfs(all_workflows['realignment'])

# 4. Create comparator
comparator = WorkflowComparator(standard_stats, realignment_stats)

# 5. Generate RNA-focused comparisons
rna_counts = comparator.compare_rna_variant_counts()
rna_categories = comparator.compare_rna_category_distribution()
rna_stages = comparator.compare_rna_annotation_stages()

# 6. Create integrative view (DNA + RNA standard + RNA realignment)
integrative = comparator.create_integrative_view(
    dna_stats=extract_dna_stats(standard_stats),
    rna_standard_stats=extract_rna_stats(standard_stats),
    rna_realignment_stats=realignment_stats
)

# 7. Calculate realignment impact
impact = comparator.calculate_realignment_impact()
print(f"RNA variants added by realignment: {impact['rna_variants_added']}")
print(f"RNA variants removed by realignment: {impact['rna_variants_removed']}")
print(f"Realignment improvement: {impact['realignment_improvement']:.2f}%")
```

### Interpretation Guidelines

**Positive Realignment Impact:**
- More variants called in realignment → Improved sensitivity
- Higher VAF in realignment → Better alignment quality
- More somatic variants → Better signal-to-noise ratio

**Negative Realignment Impact:**
- Fewer variants in realignment → Potential over-filtering
- Lower VAF in realignment → Alignment issues
- More artifacts → Need parameter tuning

**Key Metrics to Monitor:**
1. **Variant count changes**: Track across all stages
2. **Category distribution shifts**: Somatic vs Germline vs Artifact
3. **VAF improvements**: Compare RNA standard vs realignment
4. **Annotation stage impacts**: cosmic_gnomad, rna_editing, filtered_rescue

### Example Comparison Output

```
RNA Workflow Comparison Summary:
================================
Stage: filtered_rescue
  Somatic:
    Standard RNA: 1,234 variants
    Realignment RNA: 1,456 variants
    Difference: +222 variants (+18.0%)
  
  Germline:
    Standard RNA: 567 variants
    Realignment RNA: 543 variants
    Difference: -24 variants (-4.2%)
  
  RNA_Edit:
    Standard RNA: 89 variants
    Realignment RNA: 102 variants
    Difference: +13 variants (+14.6%)

Realignment Impact:
  - Variants added: 235
  - Variants removed: 13
  - Net improvement: +222 variants
  - Realignment effectiveness: 18.0% increase in somatic calls
```

## Comprehensive BAM Validation with 4 Samples

### Validation Strategy

The enhanced BAM validator now validates variants against **all four samples**:

1. **DNA_TUMOR** (standard workflow BAM)
2. **DNA_NORMAL** (standard workflow BAM)
3. **RNA_TUMOR** (standard workflow BAM)
4. **RNA_TUMOR_realign** (realignment workflow BAM)

### Usage Example

```python
# Discover BAM files for both workflows
standard_bams = discovery.discover_bam_files_for_workflow(
    workflow_manager.get_workflow_config(WorkflowType.STANDARD)
)
realignment_bams = discovery.discover_bam_files_for_workflow(
    workflow_manager.get_workflow_config(WorkflowType.REALIGNMENT)
)

# Combine all BAM files
all_bams = {
    "DNA_TUMOR": standard_bams["DNA_TUMOR"],
    "DNA_NORMAL": standard_bams["DNA_NORMAL"],
    "RNA_TUMOR": standard_bams["RNA_TUMOR"],
    "RNA_TUMOR_realign": realignment_bams["RNA_TUMOR_realign"]
}

# Validate final filtered VCF against all samples
validator = RealignmentBAMValidator(
    filtered_vcf_path=filtered_vcf,
    bam_files=all_bams
)

validation_results = validator.validate_all_samples(max_variants=100)

# Results include VAF for all 4 samples
print(validation_results[['CHROM', 'POS', 'REF', 'ALT', 
                          'DNA_TUMOR_VAF', 'DNA_NORMAL_VAF',
                          'RNA_TUMOR_VAF', 'RNA_TUMOR_realign_VAF']])
```

### Validation Output

```
Validation Results (sample):
============================
CHROM  POS      REF  ALT  DNA_TUMOR_VAF  DNA_NORMAL_VAF  RNA_TUMOR_VAF  RNA_TUMOR_realign_VAF
chr1   12345    A    G    0.45           0.02            0.38           0.52
chr2   67890    C    T    0.52           0.01            0.41           0.48
chr3   11111    G    A    0.38           0.00            0.29           0.44

Key Findings:
- RNA realignment shows improved VAF in 85% of variants
- Average VAF increase: +0.08 (21% improvement)
- Variants with >10% VAF improvement: 42
```

## Output Files

### CSV Files (in `vcf_statistics_output/`)
1. `variant_count_summary.csv` - All variant counts
2. `quality_summary.csv` - Quality statistics
3. `tool_comparison.csv` - Tool performance metrics
4. `consensus_comparison.csv` - Consensus analysis
5. `bam_validation_results.csv` - Validation results
6. **NEW: `rna_workflow_comparison.csv`** - RNA standard vs realignment comparison
7. **NEW: `integrative_view.csv`** - DNA + RNA standard + RNA realignment
8. **NEW: `realignment_impact.csv`** - Realignment effectiveness metrics
9. **NEW: `bam_validation_4samples.csv`** - Comprehensive 4-sample validation

### Report
`summary_report.txt` - Comprehensive text summary (includes realignment analysis)

## Customization Options

### Change Analysis Scope
```python
# Analyze different VCF categories
categories = ['consensus', 'rescue']  # Focus on final outputs

# Change tools analyzed
TOOLS = ['mutect2', 'deepsomatic']  # Skip Strelka

# Different modalities
MODALITIES = ['DNA_TUMOR_vs_DNA_NORMAL']  # DNA only

# Analyze specific workflow
workflow_config = workflow_manager.get_workflow_config(WorkflowType.REALIGNMENT)
realignment_vcfs = discovery.discover_vcfs_for_workflow(workflow_config)
```

### Adjust Validation Parameters
```python
# Validate more variants
validator.validate_variants(vcf, bam_map, max_variants=500)

# Change base quality threshold
# Modify in BAMValidator class: min_base_quality=30

# Comprehensive 4-sample validation
validator = RealignmentBAMValidator(filtered_vcf, all_bams)
results = validator.validate_all_samples(max_variants=200)
```

### Add Custom INFO Fields
```python
# Analyze specific INFO field
custom_stats = aggregator.create_info_field_summary('POPAF')
```

### Custom Visualizations
```python
# Create custom plots
import plotly.graph_objects as go

fig = go.Figure()
# Add your custom visualization

# Or use new comparison visualizations
visualizer.plot_rna_workflow_comparison(rna_standard_stats, rna_realignment_stats)
visualizer.plot_integrative_view(dna_stats, rna_standard_stats, rna_realignment_stats)
```

### Workflow-Specific Analysis
```python
# Analyze only realignment workflow
if WorkflowType.REALIGNMENT in workflow_manager.detect_workflows():
    realignment_config = workflow_manager.get_workflow_config(WorkflowType.REALIGNMENT)
    realignment_vcfs = discovery.discover_vcfs_for_workflow(realignment_config)
    realignment_stats = process_all_vcfs(realignment_vcfs)
    # Analyze realignment-specific results
```

## Technical Requirements

### Dependencies
- `cyvcf2` - Efficient VCF parsing
- `pysam` - BAM/CRAM file handling
- `pandas` - Data manipulation
- `numpy` - Numerical operations
- `plotly` - Interactive visualizations
- `seaborn`, `matplotlib` - Additional plotting

### Reference Genome
For CRAM file validation, set reference path:
```python
REFERENCE_FASTA = "/path/to/Homo_sapiens_assembly38.fasta"
```

## Best Practices

1. **Run in Order**: Execute notebook cells sequentially
2. **Check Outputs**: Verify file discovery before analysis
3. **Memory Management**: Large VCFs may require limiting variants analyzed
4. **Validation Sampling**: Validate representative subset of variants
5. **Export Results**: Always save results before closing notebook

## Advanced Usage

### Parallel Processing
For large datasets, consider parallelizing extraction:
```python
from multiprocessing import Pool

def process_vcf_wrapper(vcf_path):
    extractor = VCFStatisticsExtractor(vcf_path)
    return extractor.extract_all_stats()

with Pool(4) as p:
    results = p.map(process_vcf_wrapper, vcf_list)
```

### Custom Filtering
Add custom filters during extraction:
```python
# Only analyze high-quality variants
if variant.QUAL > 30 and variant.FILTER == "PASS":
    # Process variant
```

### Integration with Pipeline
Export results for downstream analysis:
```python
# Export to pipeline-compatible format
results_df.to_json('vcf_stats.json', orient='records')
```

## Troubleshooting

### Common Issues

1. **CRAM files not opening**
   - Ensure reference genome path is correct
   - Check CRAM index (.crai) exists

2. **Memory errors**
   - Reduce max_variants parameter
   - Process VCFs individually

3. **Missing INFO fields**
   - Not all tools use same INFO fields
   - Check VCF header for available fields

4. **Empty visualizations**
   - Verify data extraction completed
   - Check for empty DataFrames

5. **Realignment workflow not detected** (NEW)
   - Verify `vcf_realignment/` directory exists
   - Check that realignment VCF files are present
   - Ensure directory structure matches expected format

6. **BAM validation fails for realignment** (NEW)
   - Verify realignment BAM files exist in `vcf_realignment/preprocessing/`
   - Check BAM index files (.bai or .crai) are present
   - Ensure sample naming matches expected format

7. **Comparison plots show no data** (NEW)
   - Verify both workflows were processed successfully
   - Check that RNA samples exist in both workflows
   - Ensure statistics extraction completed for both workflows

## Performance Notes

- VCF extraction: ~1-10s per 1000 variants (depends on INFO fields)
- BAM validation: ~5-30s per 100 variants (depends on coverage)
- Visualization: Nearly instant with Plotly
- Export: <1s for typical datasets

## Future Enhancements

Potential additions:
1. Variant overlap analysis (Venn diagrams)
2. Functional annotation statistics
3. Genomic region enrichment
4. Mutation signature analysis
5. Allele frequency comparisons with population databases
6. Automated report generation (PDF/HTML)
7. **Machine learning-based realignment quality prediction**
8. **Automated realignment parameter optimization**
9. **Cross-sample realignment comparison** (multiple patients)

## Citation

If using this notebook for publications, cite the tools:
- cyvcf2: https://github.com/brentp/cyvcf2
- pysam: https://github.com/pysam-developers/pysam
- Plotly: https://plotly.com/python/

## Support

For issues or questions:
1. Check the Usage Guide in the notebook
2. Review this README
3. Inspect individual function docstrings
4. Check tool documentation (cyvcf2, pysam)

---

**Created**: November 2025  
**Updated**: January 2026 (Added realignment workflow support)  
**Version**: 2.0  
**Status**: Production Ready
