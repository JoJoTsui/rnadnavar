# COSMIC/gnomAD Annotation System - Usage Guide

## Overview

The COSMIC/gnomAD annotation system provides high-performance, reliable variant annotation using
population databases (gnomAD) and cancer mutation catalogs (COSMIC). This enhanced system features
automatic performance optimization, comprehensive output management, and robust error handling.

## Key Features

- **High Performance**: Auto-optimized parallel processing with <5 minute target for typical datasets
- **Reliable Classification**: Working FILTER reclassification with validation and statistics
- **Comprehensive Output**: Intermediate files preserved with consistent naming patterns
- **Robust Error Handling**: Immediate exit on failures with actionable diagnostics
- **Resource Optimization**: Automatic worker configuration based on system resources

## Quick Start

### 1. Basic Usage with Auto-Optimization

```bash
# Full annotation with auto-optimized performance
python bin/annotate_cosmic_gnomad.py \
    --input input.vcf.gz \
    --cosmic cosmic_database.vcf.gz \
    --gnomad /path/to/gnomad/directory \
    --output annotated_output.vcf.gz \
    --verbose \
    --stats-output annotation_stats.json

# COSMIC annotation only
python bin/annotate_cosmic_gnomad.py \
    --input input.vcf.gz \
    --cosmic cosmic_database.vcf.gz \
    --output annotated_output.vcf.gz

# gnomAD annotation only with custom thresholds
python bin/annotate_cosmic_gnomad.py \
    --input input.vcf.gz \
    --gnomad /path/to/gnomad/directory \
    --output annotated_output.vcf.gz \
    --germline-freq-threshold 0.005 \
    --somatic-consensus-threshold 2
```

### 2. Real-World Example with Actual Database Paths

```bash
# Example with real database paths (update paths as needed)
python bin/annotate_cosmic_gnomad.py \
    --input /path/to/your/rescued.vcf.gz \
    --cosmic /t9k/mnt/joey/bio_db/COSMIC/ucsc_chr_cosmic_GenomeScreensMutant_Normal_v103_GRCh38.vcf.gz \
    --gnomad /t9k/mnt/joey/bio_db/gnomAD/exomes \
    --output /path/to/output/annotated.vcf.gz \
    --verbose \
    --stats-output /path/to/output/stats.json
```

### 3. Nextflow Pipeline Integration

```bash
# Enable COSMIC/gnomAD annotation in pipeline with auto-optimization
nextflow run main.nf \
    --enable_cosmic_gnomad_annotation true \
    --cosmic_database /path/to/cosmic.vcf.gz \
    --gnomad_database /path/to/gnomad/directory \
    --cosmic_gnomad_verbose true
```

## Configuration

### Database Preparation

#### COSMIC Database
```bash
# Prepare COSMIC database for annotation
python scripts/prepare_cosmic_database.py \
    --input raw_cosmic.vcf \
    --output prepared_cosmic.vcf.gz
```

#### gnomAD Database
The gnomAD database should be organized as chromosome-split VCF files:
```
gnomad_directory/
├── gnomad.exomes.v4.1.sites.chr1.vcf.bgz
├── gnomad.exomes.v4.1.sites.chr2.vcf.bgz
├── ...
└── gnomad.exomes.v4.1.sites.chrX.vcf.bgz
```

### Pipeline Parameters

Key parameters for Nextflow pipeline:

- `enable_cosmic_gnomad_annotation`: Enable/disable annotation (default: true)
- `cosmic_database`: Path to COSMIC VCF database file
- `gnomad_database`: Path to gnomAD database directory
- `cosmic_gnomad_verbose`: Enable verbose logging (default: false)
- `cosmic_gnomad_germline_freq_threshold`: Germline frequency threshold (default: 0.001)
- `cosmic_gnomad_somatic_consensus_threshold`: Somatic consensus threshold (default: 3)
- `cosmic_gnomad_cosmic_recurrence_threshold`: COSMIC recurrence threshold (default: 5)

## Advanced Usage

### Custom Classification Thresholds

```bash
python bin/annotate_cosmic_gnomad.py \
    --input input.vcf.gz \
    --cosmic cosmic.vcf.gz \
    --gnomad /path/to/gnomad \
    --output output.vcf.gz \
    --germline-freq-threshold 0.005 \
    --somatic-consensus-threshold 2 \
    --cosmic-recurrence-threshold 10
```

### Performance Optimization

```bash
# Auto-optimized performance (recommended)
python bin/annotate_cosmic_gnomad.py \
    --input large_input.vcf.gz \
    --cosmic cosmic.vcf.gz \
    --gnomad /path/to/gnomad \
    --output output.vcf.gz \
    --verbose
    # No --workers argument = auto-optimization based on system resources

# Manual worker configuration (if needed)
python bin/annotate_cosmic_gnomad.py \
    --input large_input.vcf.gz \
    --cosmic cosmic.vcf.gz \
    --gnomad /path/to/gnomad \
    --output output.vcf.gz \
    --workers 16 \
    --verbose
```

### Statistics and Monitoring

```bash
# Generate detailed statistics
python bin/annotate_cosmic_gnomad.py \
    --input input.vcf.gz \
    --cosmic cosmic.vcf.gz \
    --gnomad /path/to/gnomad \
    --output output.vcf.gz \
    --verbose \
    --stats-output annotation_stats.json
```

## Output Files

### Comprehensive Output Organization
The system now generates a complete set of output files with consistent naming:

```
output_directory/
├── sample.cosmic.vcf.gz              # COSMIC-annotated intermediate
├── sample.cosmic.vcf.gz.tbi          # COSMIC index
├── sample.gnomad.vcf.gz              # gnomAD-annotated intermediate  
├── sample.gnomad.vcf.gz.tbi          # gnomAD index
├── sample.final.vcf.gz               # Final classified output
├── sample.final.vcf.gz.tbi           # Final index
├── sample.cosmic_gnomad_annotated.vcf.gz  # Main output (copy of final)
├── sample.cosmic_gnomad_annotated.vcf.gz.tbi
├── sample.annotation_stats.json      # Detailed statistics
└── sample.performance_metrics.json   # Performance data (if --verbose)
```

### Annotated VCF Fields
The annotated VCF files contain additional INFO fields and updated FILTER classifications:

**INFO Fields:**
- `GNOMAD_AF`: Alternative allele frequency from gnomAD
- `GNOMAD_faf95_popmax`: Filtering allele frequency 95% (population maximum)
- `COSMIC_CNT`: COSMIC genome screen sample count
- `COSMIC_ID`: COSMIC mutation identifier
- `Somatic_Rescue`: Flag indicating cross-modal rescue classification

**FILTER Classifications:**
- `Germline`: High population frequency variants (>0.001 by default)
- `Somatic`: High-confidence somatic variants (unanimous caller support or rescue)
- `Artifact`: Low-quality or conflicting evidence variants
- Original filters preserved when evidence is insufficient

### Statistics Report (JSON)
When `--stats-output` is used, a comprehensive JSON report is generated:

```json
{
    "summary": {
        "total_processing_time_seconds": 1234.5,
        "fallback_occurred": false,
        "errors_count": 0,
        "warnings_count": 2
    },
    "annotations": {
        "cosmic_annotations_added": 15420,
        "gnomad_annotations_added": 89234,
        "total_annotations_added": 104654
    },
    "classification": {
        "total_variants_classified": 12500,
        "classification_counts": {
            "Germline": 8200,
            "Somatic": 3100,
            "Artifact": 1200
        }
    }
}
```

## Performance Monitoring

### Expected Performance Targets
- **gnomAD annotation**: <5 minutes for typical rescue VCF files (≤50,000 variants)
- **COSMIC annotation**: <2 minutes for most datasets
- **Classification**: <1 minute for most datasets
- **Total pipeline**: <10 minutes for typical workflows

### Performance Optimization Tips
1. **Use auto-optimization**: Omit `--workers` argument for automatic resource detection
2. **Enable verbose logging**: Use `--verbose` to monitor performance bottlenecks
3. **Check system resources**: Ensure adequate CPU cores and memory (2GB per worker recommended)
4. **Optimize I/O**: Use fast storage (SSD) for database files and temporary processing

## Troubleshooting

### Common Issues and Solutions

1. **Performance Issues**
   ```bash
   # Check system resources
   python -c "import psutil; print(f'CPU cores: {psutil.cpu_count()}, Memory: {psutil.virtual_memory().total/1024**3:.1f}GB')"
   
   # Use auto-optimization (recommended)
   python bin/annotate_cosmic_gnomad.py ... # (no --workers argument)
   
   # Manual optimization for high-memory systems
   python bin/annotate_cosmic_gnomad.py ... --workers 16
   ```

2. **Missing Dependencies**
   ```bash
   # Install required tools
   conda install bcftools htslib pysam
   # or
   apt install bcftools tabix python3-pysam
   ```

3. **Database Access Issues**
   ```bash
   # Check COSMIC database
   bcftools view -h /path/to/cosmic.vcf.gz | head -5
   tabix -l /path/to/cosmic.vcf.gz
   
   # Check gnomAD directory
   ls -la /path/to/gnomad/exomes/*.chr*.vcf.bgz | head -5
   ```

4. **Classification Not Working**
   ```bash
   # Check FILTER statistics before and after
   bcftools query -f '%FILTER\n' input.vcf.gz | sort | uniq -c
   bcftools query -f '%FILTER\n' output.vcf.gz | sort | uniq -c
   
   # Enable verbose logging for detailed classification info
   python bin/annotate_cosmic_gnomad.py ... --verbose
   ```

5. **Memory or Disk Space Issues**
   ```bash
   # Check available resources
   df -h .  # Check disk space
   free -h  # Check memory
   
   # Reduce workers if needed
   python bin/annotate_cosmic_gnomad.py ... --workers 2
   ```

### Error Messages and Solutions

**"CRITICAL: No annotation databases provided"**
- Solution: Provide at least one database with `--cosmic` or `--gnomad`

**"CRITICAL: Input VCF file not found"**
- Solution: Check file path and ensure the VCF file exists and is readable

**"CRITICAL: COSMIC database index not found"**
- Solution: Create tabix index: `tabix -p vcf /path/to/cosmic.vcf.gz`

**"CRITICAL: No chromosome-split VCF files found in gnomAD directory"**
- Solution: Ensure gnomAD directory contains files like `gnomad.exomes.v4.1.sites.chr1.vcf.bgz`

**"gnomAD annotation failed"**
- Solution: Check gnomAD directory structure and file accessibility
- Verify chromosome files are properly indexed

**"Variant classification failed"**
- Solution: Ensure pysam library is installed: `pip install pysam`
- Check VCF format compliance

### Validation and Cleanup

```bash
# Validate complete integration
python bin/validate_cosmic_gnomad_integration.py --verbose

# Clean up temporary files and optimize
python bin/cleanup_and_optimize.py --cleanup --optimize
```

## Support and Documentation

- Requirements: `.kiro/specs/cosmic-gnomad-annotation/requirements.md`
- Design: `.kiro/specs/cosmic-gnomad-annotation/design.md`
- Implementation: `.kiro/specs/cosmic-gnomad-annotation/tasks.md`

For additional support, check the pipeline logs and use verbose mode for detailed diagnostics.