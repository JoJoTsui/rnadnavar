# VCF Statistics Analysis - Comprehensive Plan & Implementation

## Overview

This notebook provides a complete solution for analyzing VCF files from a multi-tool, multi-modal variant calling pipeline with DNA and RNA sequencing data.

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
```

### 2. **Implementation Components**

#### A. **VCFFileDiscovery Class**
- **Purpose**: Automatically discovers all VCF files in the pipeline output
- **Categories**: variant_calling, normalized, annotated, consensus, rescue, filtered
- **Output**: Organized dictionary of all VCF file paths

#### B. **VCFStatisticsExtractor Class**
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

#### C. **BAMValidator Class**
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

#### D. **StatisticsAggregator Class**
- **Purpose**: Aggregate and compare statistics across VCFs
- **Methods**:
  - `create_variant_count_summary()`: Aggregate counts
  - `create_quality_summary()`: Quality metrics
  - `compare_tools_by_modality()`: Tool performance
  - `compare_consensus_to_individual()`: Consensus analysis
  - `create_info_field_summary()`: INFO field analysis

#### E. **VCFVisualizer Class**
- **Purpose**: Create interactive visualizations with Plotly
- **Plots**:
  1. Variant counts by tool (grouped bar chart)
  2. Quality score distributions (box plots)
  3. Variant type distribution (pie charts)
  4. Consensus vs individual tools (bar + line chart)
  5. Filter status (stacked bar chart)

### 3. **Analysis Workflow**

```
1. Configure & Discover
   └── Set BASE_DIR, discover all VCFs and alignments

2. Extract Statistics
   ├── Process all VCF files with cyvcf2
   ├── Extract basic stats, INFO, FORMAT fields
   └── Store in organized dictionary

3. BAM Validation
   ├── Map VCFs to alignment files
   ├── Validate variant support
   └── Calculate VAF from pileups

4. Aggregate & Compare
   ├── Create summary tables
   ├── Compare tools within modalities
   ├── Compare consensus to individuals
   └── Analyze rescue effectiveness

5. Visualize
   ├── Generate interactive plots
   ├── Quality distributions
   └── Performance comparisons

6. Export Results
   ├── CSV files for each summary table
   ├── Validation results
   └── Summary report
```

## Key Features

### 1. **Multi-Tool Comparison**
- Compare DeepSomatic, Mutect2, and Strelka
- Evaluate sensitivity and specificity
- Identify tool-specific biases

### 2. **Multi-Modal Analysis**
- DNA vs RNA variant detection
- Modality-specific characteristics
- Cross-modality validation

### 3. **Consensus Analysis**
- Agreement across tools
- Retention rates
- High-confidence variant identification

### 4. **Rescue Analysis**
- Cross-modality variant recovery
- RNA-specific variants added to DNA consensus
- Quantify rescue effectiveness

### 5. **Quality Control**
- Quality score distributions
- Filter status tracking
- Pass/fail rates

### 6. **BAM-Level Validation**
- Independent verification using alignments
- Read-level support validation
- VAF consistency checks

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

## Output Files

### CSV Files (in `vcf_statistics_output/`)
1. `variant_count_summary.csv` - All variant counts
2. `quality_summary.csv` - Quality statistics
3. `tool_comparison.csv` - Tool performance metrics
4. `consensus_comparison.csv` - Consensus analysis
5. `bam_validation_results.csv` - Validation results

### Report
`summary_report.txt` - Comprehensive text summary

## Customization Options

### Change Analysis Scope
```python
# Analyze different VCF categories
categories = ['consensus', 'rescue']  # Focus on final outputs

# Change tools analyzed
TOOLS = ['mutect2', 'deepsomatic']  # Skip Strelka

# Different modalities
MODALITIES = ['DNA_TUMOR_vs_DNA_NORMAL']  # DNA only
```

### Adjust Validation Parameters
```python
# Validate more variants
validator.validate_variants(vcf, bam_map, max_variants=500)

# Change base quality threshold
# Modify in BAMValidator class: min_base_quality=30
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
**Version**: 1.0
**Status**: Production Ready
