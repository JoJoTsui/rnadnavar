# Variant Classification System

## Overview

This document describes the unified variant classification system implemented across the rnadnavar pipeline. The system classifies variants from different callers (Strelka, DeepSomatic, Mutect2) into four biological categories and applies standardized FILTER values.

## Classification Categories

All variants are classified into one of four biological categories:

1. **Somatic**: True somatic mutations (tumor-specific)
2. **Germline**: Germline variants (present in normal tissue)
3. **Reference**: Reference calls (no variant detected)
4. **Artifact**: Low-quality variants or artifacts

## Unified FILTER Values

Each classification maps to a standardized FILTER value:

| Classification | FILTER Value | Description |
|---------------|--------------|-------------|
| Somatic | PASS | All filters passed - high-confidence somatic variant |
| Germline | GERMLINE | Variant detected in normal sample - likely germline |
| Reference | RefCall | Reference call - no variant detected |
| Artifact | LowQuality | Low quality variant - failed quality filters |

## Caller-Specific Logic

### Strelka Classification

Strelka variants are classified based on three criteria:
- **FILTER field**: Original caller filter status
- **NT field (INFO/NT)**: Indicates reference vs variant status (ref/het/hom)
- **Normal depth**: Depth of coverage in normal sample (≥2 required for rescue)

Classification rules:
1. `PASS + NT=ref` → **Somatic** (tumor variant, normal reference)
2. `PASS + NT≠ref` → **Germline** (variant in both samples)
3. Failed filters + NT=ref + Normal DP≥2 → **Somatic** (rescue low-quality somatic)
4. Failed filters + NT≠ref → **Germline** (germline fails QC)
5. Failed filters + NT=ref + Normal DP<2 → **Artifact** (insufficient normal coverage)
6. `NT=ref (always)` → **Reference** (explicit reference call)

### DeepSomatic Classification

DeepSomatic uses explicit FILTER labels:
- `PASS` → **Somatic**
- `GERMLINE` → **Germline**
- `RefCall` → **Reference**
- Any other filter → **Artifact**

### Mutect2 Classification

Mutect2 classification is primarily based on FILTER:
- `PASS` → **Somatic**
- Filters containing "germline" (case-insensitive) → **Germline**
- Any other filter → **Artifact**

Note: Mutect2 typically doesn't produce explicit reference calls.

## Module Architecture

### `vcf_utils/classification.py`

Core classification module with all classification logic:

**Classification Functions:**
- `classify_strelka_variant(filter_val, nt_val, normal_dp)` - Strelka-specific logic
- `classify_deepsomatic_variant(filter_val)` - DeepSomatic-specific logic
- `classify_mutect2_variant(filter_val, info_dict=None)` - Mutect2-specific logic
- `classify_variant_from_record(variant, caller_name, sample_indices)` - Universal classifier for cyvcf2 objects
- `classify_variant_from_dict(variant_dict, caller_name)` - Classifier for dict representations

**Helper Functions:**
- `get_sample_indices(vcf_obj, caller_name)` - Identifies tumor/normal sample positions
- `normalize_filter_value(classification)` - Maps classification to unified FILTER
- `get_unified_filter_headers()` - Returns standardized FILTER definitions
- `get_classification_info_headers()` - Returns classification INFO field definitions

### `vcf_utils/aggregation.py`

Integration of classification into variant reading:

**Modified Function:**
```python
def read_variants_from_vcf(vcf_path, caller_name, modality=None, 
                          exclude_refcall=False, exclude_germline=False,
                          classify_variants=True)
```

**Changes:**
- Added `classify_variants=True` parameter (enabled by default)
- Calls `get_sample_indices()` for Strelka to identify tumor/normal
- Calls `classify_variant_from_record()` for each variant
- Uses `normalize_filter_value()` for unified FILTER assignment
- Stores 'classification' field in variant data dict

### `vcf_utils/io_utils.py`

VCF header creation with classification metadata:

**Modified Function:**
```python
def create_output_header(template_header, sample_name=None)
```

**Changes:**
- Added INFO fields: `VC`, `VC_CALLERS`, `VC_CONSENSUS`
- Added unified FILTER definitions: PASS, GERMLINE, RefCall, LowQuality
- Retained legacy filters for backwards compatibility

## Usage in Pipeline

### Consensus VCF Generation (`run_consensus_vcf.py`)

Classification is automatically enabled when reading variants:

```python
variants = read_variants_from_vcf(
    vcf_path, 
    caller, 
    modality=None,
    exclude_refcall=args.exclude_refcall,
    exclude_germline=args.exclude_germline
    # classify_variants=True by default
)
```

Each variant in the returned list includes:
- `classification`: Biological category (Somatic/Germline/Reference/Artifact)
- `filter`: Unified FILTER value (PASS/GERMLINE/RefCall/LowQuality)

### Statistics Notebook (`vcf_statistics.ipynb`)

The notebook uses classification for variant statistics:

```python
# Extract basic statistics with classification
stats = extractor.extract_basic_stats()

# Display classification summary
classification_summary = aggregator.create_classification_summary()
print(classification_summary)
```

Output includes counts for each classification category per caller.

## VCF Output Format

### INFO Fields

**VC** (Variant Classification):
- Type: String
- Number: 1
- Description: Biological classification (Somatic/Germline/Reference/Artifact)
- Example: `VC=Somatic`

**VC_CALLERS** (Classification by Callers):
- Type: String
- Number: .
- Description: Classification from each caller that detected the variant
- Example: `VC_CALLERS=strelka:Somatic|mutect2:Somatic`

**VC_CONSENSUS** (Consensus Classification):
- Type: String
- Number: 1
- Description: Consensus biological classification across all callers
- Example: `VC_CONSENSUS=Somatic`

### FILTER Values

**Standardized Filters:**
- `PASS`: High-confidence somatic variant (Classification: Somatic)
- `GERMLINE`: Likely germline variant (Classification: Germline)
- `RefCall`: Reference call (Classification: Reference)
- `LowQuality`: Failed quality filters (Classification: Artifact)

**Legacy Filters (retained for compatibility):**
- `LowDepth`, `StrandBias`, `Artifact`, `NoConsensus`, `LowEvidenceScore`

## Testing

### Module Import Test

```bash
cd /t9k/mnt/hdd/work/Vax/pipeline/rnadnavar/bin
python3 -c "from vcf_utils.classification import *; print('✓ Import successful')"
```

### Classification Function Test

```python
from vcf_utils.classification import (
    classify_strelka_variant,
    classify_deepsomatic_variant,
    classify_mutect2_variant,
    normalize_filter_value
)

# Test Strelka classification
classification = classify_strelka_variant('PASS', 'ref', None)
unified_filter = normalize_filter_value(classification)
print(f"Strelka: {classification} → {unified_filter}")

# Test DeepSomatic classification
classification = classify_deepsomatic_variant('PASS')
unified_filter = normalize_filter_value(classification)
print(f"DeepSomatic: {classification} → {unified_filter}")

# Test Mutect2 classification
classification = classify_mutect2_variant('PASS')
unified_filter = normalize_filter_value(classification)
print(f"Mutect2: {classification} → {unified_filter}")
```

### Integration Test

```python
from vcf_utils.aggregation import read_variants_from_vcf

# Read variants with classification enabled
variants = read_variants_from_vcf(
    'input.vcf',
    'strelka',
    classify_variants=True
)

# Check classification results
for var in variants[:5]:
    print(f"{var['chrom']}:{var['pos']} - {var['classification']} ({var['filter']})")
```

## Implementation Notes

### Strelka Rescue Logic

Strelka has special "rescue" logic for low-quality variants:
- Variants that fail QC but have NT=ref (normal is reference)
- AND normal depth ≥ 2 (sufficient normal coverage)
- Are rescued and classified as **Somatic**

This prevents discarding true somatic variants that failed QC due to technical reasons but have good normal sample evidence.

### Sample Index Detection

For Strelka, the classification requires identifying tumor vs normal samples:
- Searches for sample names containing "tumor" or "normal" (case-insensitive)
- Fallback: assumes first sample is tumor, second is normal
- Sample indices are detected once per VCF file for efficiency

### Backwards Compatibility

The system maintains backwards compatibility:
- Legacy FILTER values are retained in VCF headers
- Classification is enabled by default but can be disabled with `classify_variants=False`
- Existing code that doesn't use classification will continue to work

## Future Enhancements

Potential improvements to the classification system:

1. **Consensus Classification**: Aggregate classifications from multiple callers
2. **VAF-based Classification**: Use allele frequencies to refine germline detection
3. **Population Database Integration**: Use gnomAD or similar for germline confirmation
4. **Machine Learning**: Train classifier on known somatic/germline variants
5. **Quality Score Adjustment**: Boost/penalize variants based on classification confidence

## References

- Strelka2: [Kim et al., Nature Methods 2018](https://www.nature.com/articles/s41592-018-0051-x)
- DeepSomatic: Google Health variant caller
- Mutect2: [GATK Best Practices](https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2)

## Change Log

### 2024 - Initial Implementation
- Created unified classification system
- Implemented caller-specific classification logic
- Integrated classification into variant reading pipeline
- Added unified FILTER scheme
- Added classification INFO fields to VCF output
- Updated statistics notebook to display classification summaries
