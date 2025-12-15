# RNA Editing Annotation Guide

## Overview

The nf-core/rnadnavar pipeline includes comprehensive RNA editing annotation functionality that helps distinguish true somatic mutations from RNA editing events in RNA-seq data. This feature integrates seamlessly with the existing workflow to provide enhanced variant classification and filtering capabilities.

## What is RNA Editing Annotation?

RNA editing is a post-transcriptional modification process where specific nucleotides in RNA molecules are chemically altered, most commonly adenosine-to-inosine (A-to-I) editing. In RNA-seq variant calling, these editing events can appear as false-positive somatic mutations, making it crucial to identify and properly annotate them.

The RNA editing annotation module:
- Annotates rescue VCF files with known RNA editing sites from the REDIportal database
- Provides evidence tiering based on RNA/DNA caller support
- Adds comprehensive annotation fields for downstream analysis
- Integrates seamlessly with existing filtering workflows

## Key Features

### 1. REDIportal Database Integration
- Uses the comprehensive REDIportal database of known RNA editing sites
- Supports automatic conversion from text format to VCF format
- Handles coordinate system matching (hg19/GRCh37 vs hg38/GRCh38)
- Provides automatic tabix indexing

### 2. Evidence Tiering System
The annotation system classifies variants based on RNA/DNA caller support:

- **HIGH**: Strong evidence from multiple RNA callers (≥ min_rna_support threshold)
- **MEDIUM**: Moderate evidence from some RNA callers
- **LOW**: Limited evidence from few RNA callers
- **NONE**: No supporting evidence from RNA callers

### 3. Comprehensive Annotation Fields
Added INFO fields include:
- `REDI_GENE`: Gene name from REDIportal database
- `REDI_REGION`: Genomic region type (e.g., exonic, intronic, UTR)
- `REDI_STRAND`: Strand information for the editing site
- `REDI_EDIT_TYPE`: Type of RNA editing (e.g., A-to-G, C-to-T)
- `REDI_EVIDENCE`: Evidence tier classification
- `REDI_RNA_SUPPORT`: Number of RNA callers supporting the variant

### 4. FILTER Field Updates
- Adds `RNAedit` to the FILTER field for confirmed RNA editing sites
- Preserves existing filter annotations
- Maintains compatibility with downstream filtering tools

## Configuration and Parameters

### Basic Configuration

Enable RNA editing annotation with minimal configuration:

```bash
nextflow run nf-core/rnadnavar \
    --input samplesheet.csv \
    --outdir results \
    --enable_rna_annotation \
    --rediportal_vcf /path/to/REDIportal_hg38_v3.vcf.gz \
    -profile docker
```

### Parameter Reference

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `enable_rna_annotation` | boolean | `false` | Enable RNA editing annotation |
| `rediportal_vcf` | string | `null` | Path to REDIportal VCF database |
| `rediportal_tbi` | string | `null` | Path to tabix index (auto-generated if missing) |
| `min_rna_support` | integer | `2` | Minimum RNA caller support threshold |
| `rna_annotation_memory` | string | `4.GB` | Memory allocation for annotation |
| `rna_annotation_cpus` | integer | `2` | CPU allocation for annotation |
| `rna_annotation_args` | string | `null` | Additional arguments for annotation script |

### Advanced Configuration

For high-throughput or large-scale analyses:

```bash
nextflow run nf-core/rnadnavar \
    --input samplesheet.csv \
    --outdir results \
    --enable_rna_annotation \
    --rediportal_vcf /path/to/REDIportal_hg38_v3.vcf.gz \
    --min_rna_support 3 \
    --rna_annotation_memory 16.GB \
    --rna_annotation_cpus 8 \
    --rna_annotation_log_level DEBUG \
    --rna_annotation_enable_metrics \
    -profile docker
```

## REDIportal Database Setup

### Downloading REDIportal Data

The REDIportal database can be downloaded from the official repository:

```bash
# For human hg38/GRCh38
wget http://srv00.recas.ba.infn.it/webshare/ATLAS/donwload/TABLE1_hg38.txt.gz

# For human hg19/GRCh37
wget http://srv00.recas.ba.infn.it/webshare/ATLAS/donwload/TABLE1_hg19.txt.gz
```

### Converting to VCF Format

Use the provided conversion script to convert the text format to VCF:

```bash
# Convert to VCF format
python scripts/convert_rediportal_to_vcf.py \
    --input TABLE1_hg38.txt.gz \
    --output REDIportal_hg38_v3.vcf.gz \
    --genome hg38

# Create tabix index
tabix -p vcf REDIportal_hg38_v3.vcf.gz
```

### Database Validation

Verify the converted database:

```bash
# Check VCF format
bcftools view -h REDIportal_hg38_v3.vcf.gz | head -20

# Check variant count
bcftools view -H REDIportal_hg38_v3.vcf.gz | wc -l

# Validate tabix index
tabix -l REDIportal_hg38_v3.vcf.gz
```

## Workflow Integration

### Pipeline Flow

The RNA editing annotation integrates seamlessly into the existing workflow:

```
VCF Rescue Generation
         ↓
RNA Editing Annotation (if enabled)
         ↓
VCF Rescue Filtering
         ↓
Final Output
```

### Conditional Processing

- When `enable_rna_annotation = false`: Rescue VCF flows directly to filtering
- When `enable_rna_annotation = true`: Rescue VCF is annotated before filtering
- All existing workflow functionality is preserved
- Backward compatibility is maintained

### Output Files

The annotation process generates:
- **Annotated VCF**: `{sample}.rescue.rna_annotated.vcf.gz`
- **Tabix Index**: `{sample}.rescue.rna_annotated.vcf.gz.tbi`
- **Processing Logs**: Detailed annotation statistics and performance metrics
- **Version Information**: Software versions used in annotation

## Performance Considerations

### Resource Requirements

**Minimum Requirements:**
- Memory: 4 GB
- CPUs: 2
- Disk: 2x input VCF size for temporary files

**Recommended for Large Datasets:**
- Memory: 8-16 GB
- CPUs: 4-8
- Disk: 5x input VCF size for optimal performance

### Performance Optimization

1. **Memory Allocation**: Increase memory for large VCF files
   ```bash
   --rna_annotation_memory 16.GB
   ```

2. **CPU Utilization**: Use multiple CPUs for parallel processing
   ```bash
   --rna_annotation_cpus 8
   ```

3. **Disk I/O**: Use fast storage for temporary files
   ```bash
   --rna_annotation_args "--temp-dir /fast/storage/tmp"
   ```

### Benchmarking Results

Typical performance metrics:
- **Small VCF** (< 1K variants): 1-2 minutes, 2 GB memory
- **Medium VCF** (1K-10K variants): 5-10 minutes, 4 GB memory
- **Large VCF** (> 10K variants): 15-30 minutes, 8-16 GB memory

## Quality Control and Validation

### Annotation Statistics

The pipeline provides comprehensive statistics:
- Total variants processed
- Variants annotated with RNA editing sites
- Evidence tier distribution
- Processing time and resource usage

### Validation Checks

Automatic validation includes:
- Input VCF format validation
- REDIportal database accessibility
- Coordinate system compatibility
- Output file integrity

### Quality Metrics

Key quality indicators:
- **Annotation Rate**: Percentage of variants annotated
- **Evidence Distribution**: Balance across evidence tiers
- **Processing Efficiency**: Variants processed per minute
- **Resource Utilization**: Memory and CPU usage patterns

## Troubleshooting

### Common Issues

1. **Missing REDIportal Database**
   ```
   Error: REDIportal VCF file not found
   Solution: Provide valid --rediportal_vcf path
   ```

2. **Format Validation Errors**
   ```
   Error: Invalid VCF format
   Solution: Use conversion script to generate proper VCF
   ```

3. **Memory Issues**
   ```
   Error: Out of memory during annotation
   Solution: Increase --rna_annotation_memory
   ```

4. **Coordinate Mismatch**
   ```
   Warning: No annotations found
   Solution: Ensure genome builds match (hg19 vs hg38)
   ```

### Debugging Steps

1. **Enable Verbose Logging**
   ```bash
   --rna_annotation_log_level DEBUG
   --rna_annotation_verbose
   ```

2. **Check Input Files**
   ```bash
   bcftools view -h input.vcf.gz
   bcftools view -h rediportal.vcf.gz
   ```

3. **Validate Resources**
   ```bash
   --rna_annotation_enable_metrics
   --rna_annotation_save_reports
   ```

4. **Test with Subset**
   ```bash
   bcftools view -r chr1:1-1000000 input.vcf.gz > test_subset.vcf
   ```

## Best Practices

### Database Management
- Use the latest REDIportal version for your genome build
- Validate database format before running the pipeline
- Keep databases in a centralized, accessible location
- Document database versions and sources

### Parameter Tuning
- Start with default parameters for initial runs
- Adjust `min_rna_support` based on your variant calling strategy
- Scale resources based on input data size
- Monitor performance metrics for optimization

### Quality Assurance
- Review annotation statistics for each run
- Validate output VCF format and content
- Check evidence tier distributions for reasonableness
- Compare results across similar samples

### Integration with Downstream Analysis
- Consider RNA editing annotations in variant filtering
- Use evidence tiers for prioritization
- Document annotation parameters in analysis reports
- Maintain annotation provenance information

## Advanced Usage

### Custom REDIportal Databases
Create custom databases for specific research needs:

```bash
# Filter REDIportal for specific regions
bcftools view -r chr1,chr2,chr3 REDIportal_hg38_v3.vcf.gz > custom_rediportal.vcf
bgzip custom_rediportal.vcf
tabix -p vcf custom_rediportal.vcf.gz
```

### Batch Processing
Process multiple samples efficiently:

```bash
# Use sample sheet with multiple samples
nextflow run nf-core/rnadnavar \
    --input multi_sample_sheet.csv \
    --enable_rna_annotation \
    --rediportal_vcf REDIportal_hg38_v3.vcf.gz \
    --rna_annotation_cpus 4 \
    -profile docker
```

### Integration with Other Tools
Combine with other annotation tools:

```bash
# Run with VEP annotation
nextflow run nf-core/rnadnavar \
    --input samplesheet.csv \
    --enable_rna_annotation \
    --rediportal_vcf REDIportal_hg38_v3.vcf.gz \
    --tools "vep,norm,consensus,filtering" \
    -profile docker
```

## Support and Resources

### Documentation
- [Pipeline Usage Guide](usage.md)
- [Parameter Documentation](https://nf-co.re/rnadnavar/parameters)
- [Output Documentation](output.md)
- [RNA Editing Logging Guide](RNA_EDITING_LOGGING_GUIDE.md)

### Community Support
- [nf-core Slack #rnadnavar channel](https://nfcore.slack.com/channels/rnadnavar)
- [GitHub Issues](https://github.com/nf-core/rnadnavar/issues)
- [nf-core Help](https://nf-co.re/help)

### External Resources
- [REDIportal Database](http://srv00.recas.ba.infn.it/atlas/index.html)
- [RNA Editing Research](https://www.nature.com/articles/s41576-018-0049-z)
- [GATK Best Practices](https://gatk.broadinstitute.org/hc/en-us/articles/360035531192)

## Citation

If you use the RNA editing annotation functionality in your research, please cite:

- The nf-core/rnadnavar pipeline
- The REDIportal database: Picardi E, D'Erchia AM, Lo Giudice C, Pesole G. REDIportal: a comprehensive database of A-to-I RNA editing events in humans. Nucleic Acids Res. 2017.
- Any additional tools used in your specific analysis workflow