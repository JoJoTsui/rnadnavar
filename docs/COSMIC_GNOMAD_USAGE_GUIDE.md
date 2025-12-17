# COSMIC/gnomAD Annotation System - Usage Guide

## Overview

The COSMIC/gnomAD annotation system provides comprehensive variant annotation using
population databases (gnomAD) and cancer mutation catalogs (COSMIC). This guide
covers installation, configuration, and usage.

## Quick Start

### 1. Basic Usage

```bash
# Full annotation with both databases
python bin/annotate_cosmic_gnomad.py \
    --input input.vcf.gz \
    --cosmic cosmic_database.vcf.gz \
    --gnomad /path/to/gnomad/directory \
    --output annotated_output.vcf.gz

# COSMIC annotation only
python bin/annotate_cosmic_gnomad.py \
    --input input.vcf.gz \
    --cosmic cosmic_database.vcf.gz \
    --output annotated_output.vcf.gz

# gnomAD annotation only
python bin/annotate_cosmic_gnomad.py \
    --input input.vcf.gz \
    --gnomad /path/to/gnomad/directory \
    --output annotated_output.vcf.gz
```

### 2. Nextflow Pipeline Integration

```bash
# Enable COSMIC/gnomAD annotation in pipeline
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
# Increase parallel workers for large datasets
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

### Annotated VCF
The main output is an annotated VCF file with additional INFO fields:

- `GNOMAD_FAF95`: Population frequency (filtering allele frequency 95%)
- `GNOMAD_AF`: Alternative allele frequency
- `COSMIC_CNT`: COSMIC genome screen sample count
- `COSMIC_ID`: COSMIC mutation identifier
- `Somatic_Rescue`: Flag indicating cross-modal rescue (when applicable)

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

## Troubleshooting

### Common Issues

1. **Missing Dependencies**
   ```bash
   # Install required tools
   conda install bcftools htslib
   # or
   apt install bcftools tabix
   ```

2. **Permission Errors**
   ```bash
   # Fix file permissions
   python bin/cleanup_and_optimize.py --fix-permissions
   ```

3. **Memory Issues**
   ```bash
   # Reduce parallel workers
   python bin/annotate_cosmic_gnomad.py ... --workers 2
   ```

4. **Database Format Issues**
   ```bash
   # Validate and prepare COSMIC database
   python scripts/prepare_cosmic_database.py --input cosmic.vcf --output cosmic_prepared.vcf.gz
   ```

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