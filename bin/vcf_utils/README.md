# VCF Utils Package

A Python package providing reusable functions for VCF variant aggregation, tagging, statistics generation, and I/O operations. Designed to support both within-modality consensus calling and cross-modality variant rescue workflows.

## Overview

The `vcf_utils` package extracts and modularizes core functionality from the rnadnavar pipeline's consensus and rescue workflows. It provides a clean, reusable API for:

- **Variant Aggregation**: Combining variants from multiple callers and modalities
- **Variant Tagging**: Adding caller support and modality metadata
- **Statistics Generation**: Computing consensus and rescue effectiveness metrics
- **VCF I/O**: Reading and writing VCF files with custom INFO fields
- **Filter Normalization**: Standardizing filter names across different callers

## Package Structure

```
vcf_utils/
├── __init__.py           # Package initialization and documentation
├── aggregation.py        # Variant aggregation and genotype extraction
├── tagging.py           # Variant tagging and annotation
├── statistics.py        # Statistics generation and reporting
├── io_utils.py          # VCF I/O operations and header management
├── filters.py           # Filter normalization and categorization
└── README.md            # This file
```

## Installation

The package is part of the rnadnavar pipeline and is automatically available when the pipeline is installed. No separate installation is required.

## Dependencies

- `cyvcf2`: Fast VCF reading
- `pysam`: VCF writing and header management
- Python standard library: `collections`, `statistics`, `pathlib`, `re`


## Quick Start

### Within-Modality Consensus

Aggregate variants from multiple callers within the same modality (DNA or RNA):

```python
from vcf_utils.aggregation import read_variants_from_vcf, aggregate_variants
from vcf_utils.statistics import compute_consensus_statistics, print_statistics
from vcf_utils.io_utils import write_union_vcf
from cyvcf2 import VCF

# Read variants from each caller
mutect2_vars = read_variants_from_vcf('mutect2.vcf.gz', 'mutect2')
strelka_vars = read_variants_from_vcf('strelka.vcf.gz', 'strelka')
deepsomatic_vars = read_variants_from_vcf('deepsomatic.vcf.gz', 'deepsomatic',
                                          exclude_refcall=True, exclude_germline=True)

# Prepare collections for aggregation
collections = [
    ('mutect2', mutect2_vars, None),
    ('strelka', strelka_vars, None),
    ('deepsomatic', deepsomatic_vars, None)
]

# Aggregate with consensus thresholds
aggregated = aggregate_variants(collections, snv_threshold=2, indel_threshold=2)

# Generate and print statistics
stats = compute_consensus_statistics(aggregated, snv_threshold=2, indel_threshold=2)
print_statistics(stats, 'consensus')

# Write output VCF
template = VCF('mutect2.vcf.gz')
write_union_vcf(aggregated, template, 'sample123', 'consensus.vcf.gz',
                'vcf.gz', ['mutect2', 'strelka', 'deepsomatic'])
```


### Cross-Modality Rescue

Aggregate variants across DNA and RNA modalities to identify cross-modality support:

```python
from vcf_utils.aggregation import read_variants_from_vcf, aggregate_variants
from vcf_utils.tagging import mark_rescued_variants, tag_variant_with_modality
from vcf_utils.statistics import compute_rescue_statistics, print_statistics
from vcf_utils.io_utils import write_union_vcf
from cyvcf2 import VCF

# Read DNA consensus VCF
dna_consensus = read_variants_from_vcf('dna_consensus.vcf.gz', 'consensus', modality='DNA')

# Read RNA consensus VCF
rna_consensus = read_variants_from_vcf('rna_consensus.vcf.gz', 'consensus', modality='RNA')

# Aggregate across modalities
collections = [
    ('consensus', dna_consensus, 'DNA'),
    ('consensus', rna_consensus, 'RNA')
]
aggregated = aggregate_variants(collections, snv_threshold=1, indel_threshold=1)

# Mark rescued variants (present in both DNA and RNA)
dna_keys = set(dna_consensus.keys())
rna_keys = set(rna_consensus.keys())
aggregated = mark_rescued_variants(aggregated, dna_keys, rna_keys)

# Tag with modality information
modality_map = {'consensus': 'DNA'}  # Build complete map for all callers
for vkey, data in aggregated.items():
    tag_variant_with_modality(data, modality_map)

# Generate rescue statistics
stats = compute_rescue_statistics(aggregated, dna_keys, rna_keys)
print_statistics(stats, 'rescue')

# Write output with modality information
template = VCF('dna_consensus.vcf.gz')
write_union_vcf(aggregated, template, 'sample123', 'rescued.vcf.gz',
                'vcf.gz', ['consensus'], modality_map)
```


## Module Documentation

### aggregation.py

Core variant aggregation logic for combining variants from multiple callers and modalities.

**Key Functions:**

- `extract_genotype_info(variant, caller)`: Extract GT, DP, AD, VAF, GQ from a variant
  - Special handling for Strelka's non-standard format fields
  - Returns dict with genotype information

- `aggregate_genotypes(genotypes_by_caller, callers_order)`: Aggregate genotype info across callers
  - Computes consensus genotype (majority vote)
  - Calculates mean/min/max for DP and VAF
  - Returns aggregated statistics dict

- `read_variants_from_vcf(vcf_path, caller_name, modality=None, exclude_refcall=False, exclude_germline=False)`: Read variants from VCF file
  - Supports optional modality tagging
  - Can exclude RefCall and GERMLINE variants
  - Returns dict mapping variant_key to variant_data

- `aggregate_variants(variant_collections, snv_threshold=2, indel_threshold=2)`: Aggregate variants from multiple collections
  - Applies consensus thresholds separately for SNVs and indels
  - Aggregates genotype information automatically
  - Returns dict of aggregated variant data

**Example:**
```python
from vcf_utils.aggregation import read_variants_from_vcf, aggregate_variants

# Read variants
vars1 = read_variants_from_vcf('caller1.vcf.gz', 'caller1', modality='DNA')
vars2 = read_variants_from_vcf('caller2.vcf.gz', 'caller2', modality='DNA')

# Aggregate
collections = [('caller1', vars1, 'DNA'), ('caller2', vars2, 'DNA')]
aggregated = aggregate_variants(collections, snv_threshold=2, indel_threshold=2)
```


### tagging.py

Variant tagging and annotation functions for adding metadata to variants.

**Key Functions:**

- `tag_variant_with_callers(variant_data, all_callers)`: Add caller support information
  - Adds n_callers, support_callers, n_support_callers fields
  - Returns updated variant_data dict

- `tag_variant_with_modality(variant_data, modality_map)`: Add modality tracking
  - Groups callers by modality (DNA/RNA)
  - Adds modalities, dna_callers, rna_callers fields
  - Returns updated variant_data dict

- `mark_rescued_variants(variant_data, dna_variants, rna_variants)`: Mark cross-modality variants
  - Sets rescued, dna_support, rna_support, cross_modality flags
  - Returns updated variant_data dict

- `compute_unified_filter(variant_data)`: Compute overall filter status
  - Uses majority vote across callers
  - Returns (status, filter_tags) tuple

**Example:**
```python
from vcf_utils.tagging import mark_rescued_variants, tag_variant_with_modality

# Mark rescued variants
dna_keys = set(dna_variants.keys())
rna_keys = set(rna_variants.keys())
aggregated = mark_rescued_variants(aggregated, dna_keys, rna_keys)

# Tag with modality
modality_map = {'caller1': 'DNA', 'caller2': 'RNA'}
for vkey, data in aggregated.items():
    tag_variant_with_modality(data, modality_map)
```


### statistics.py

Statistics generation for consensus and rescue operations.

**Key Functions:**

- `compute_consensus_statistics(variant_data, snv_threshold, indel_threshold)`: Compute consensus stats
  - Counts total variants, SNVs, indels
  - Counts variants passing consensus
  - Counts single vs multi-caller support
  - Returns statistics dict

- `compute_rescue_statistics(variant_data, dna_variants, rna_variants)`: Compute rescue stats
  - Counts DNA-only, RNA-only, cross-modality variants
  - Counts rescued variants
  - Calculates rescue rate
  - Returns statistics dict

- `print_statistics(stats, operation_type='consensus')`: Print formatted statistics
  - Supports 'consensus' and 'rescue' operation types
  - Prints to stdout in human-readable format

**Example:**
```python
from vcf_utils.statistics import compute_consensus_statistics, print_statistics

# Compute and print consensus statistics
stats = compute_consensus_statistics(aggregated, snv_threshold=2, indel_threshold=2)
print_statistics(stats, 'consensus')

# Output:
# - Statistics:
#   - Total variants: 1,234
#   - SNVs: 1,100 (consensus: 950)
#   - Indels: 134 (consensus: 120)
#   - Single caller: 200
#   - Multiple callers: 1,034
```


### io_utils.py

VCF I/O operations and header management.

**Key Functions:**

- `get_caller_name(filename)`: Extract caller name from VCF filename
  - Handles .variants. and .consensus. patterns
  - Returns caller name string

- `normalize_chromosome(chrom)`: Normalize chromosome names
  - Removes 'chr' prefix for consistent sorting
  - Returns normalized chromosome string

- `variant_key(variant, use_cyvcf2=True)`: Create unique variant identifier
  - Format: "chrom:pos:ref:alt"
  - Works with cyvcf2.Variant or dict
  - Returns variant key string

- `is_snv(ref, alt_list)`: Determine if variant is SNV
  - Checks if all alleles are single nucleotides
  - Returns boolean

- `create_output_header(template_header, sample_name, include_rescue_fields=False)`: Create VCF header
  - Adds all consensus INFO fields
  - Optionally adds rescue-specific fields
  - Returns pysam.VariantHeader

- `write_union_vcf(variant_data, template_header, sample_name, out_file, output_format, all_callers, modality_map=None)`: Write VCF file
  - Writes aggregated variants with all INFO fields
  - Sorts variants by genomic position
  - Optionally includes modality information
  - Returns number of variants written

**Example:**
```python
from vcf_utils.io_utils import create_output_header, write_union_vcf
from cyvcf2 import VCF

# Create header
template = VCF('input.vcf.gz')
header = create_output_header(template, 'sample123', include_rescue_fields=True)

# Write VCF
n_written = write_union_vcf(aggregated, template, 'sample123', 'output.vcf.gz',
                            'vcf.gz', ['caller1', 'caller2'], modality_map)
print(f"Wrote {n_written} variants")
```


### filters.py

Filter normalization and categorization for standardizing filter names across callers.

**Constants:**

- `FILTER_MAP`: Dict mapping caller-specific filter names to normalized names
  - Handles Mutect2, Strelka, DeepSomatic filters
  - Maps to standardized filter names

- `FILTER_CATEGORIES`: Dict grouping normalized filters into categories
  - Categories: quality, depth, bias, germline, artifact, technical, reference

**Key Functions:**

- `normalize_filter(filter_str)`: Normalize filter string
  - Handles PASS, missing, and compound filters
  - Maps to standardized names using FILTER_MAP
  - Returns normalized filter string

- `categorize_filter(normalized_filter)`: Categorize normalized filter
  - Groups filters into major categories
  - Returns category string or 'Other'

**Example:**
```python
from vcf_utils.filters import normalize_filter, categorize_filter

# Normalize Mutect2 filter
filter_str = "germline;base_qual"
normalized = normalize_filter(filter_str)
print(normalized)  # "Germline;LowBaseQuality"

# Categorize
category = categorize_filter(normalized)
print(category)  # "germline;technical"

# Normalize Strelka filter
strelka_filter = "LowDepth;LowEVS"
print(normalize_filter(strelka_filter))  # "LowDepth;LowEvidenceScore"
```


## Design Decisions

### Modular Architecture

The package is organized into five focused modules, each with a single responsibility:

1. **aggregation.py**: Core variant combining logic
2. **tagging.py**: Metadata annotation
3. **statistics.py**: Metrics computation
4. **io_utils.py**: File operations
5. **filters.py**: Filter standardization

This separation enables:
- Easy testing of individual components
- Reuse across different workflows (consensus, rescue, filtering)
- Clear dependencies between modules

### Variant Key Design

Variants are identified using a composite key: `"chrom:pos:ref:alt"`

- Normalized chromosome names (no 'chr' prefix) for consistent sorting
- Includes all alternate alleles for multi-allelic variants
- Works with both cyvcf2 and dict representations

### Genotype Aggregation Strategy

Genotype information is aggregated using:
- **Consensus genotype**: Majority vote across callers
- **Statistics**: Mean, min, max for DP and VAF
- **Per-caller tracking**: Ordered lists maintain caller correspondence

### Filter Normalization

Different callers use different filter naming conventions. The package normalizes these to:
- Standardized names (e.g., 'Germline', 'LowDepth', 'StrandBias')
- Major categories (quality, depth, bias, germline, artifact, technical, reference)
- Unified filter status based on majority vote

### Modality Tracking

For cross-modality rescue, the package tracks:
- Modality tags ('DNA' or 'RNA') for each caller
- Modality-specific statistics (DP_DNA_MEAN, VAF_RNA_MEAN, etc.)
- Cross-modality support flags (rescued, dna_support, rna_support)


## VCF INFO Fields

### Consensus VCF Fields

All consensus VCF files include these INFO fields:

| Field | Type | Description |
|-------|------|-------------|
| N_CALLERS | Integer | Total number of aggregated callers |
| CALLERS | String | Pipe-separated list of all callers |
| N_SUPPORT_CALLERS | Integer | Number of callers that detected this variant |
| CALLERS_SUPPORT | String | Pipe-separated list of supporting callers |
| FILTERS_ORIGINAL | String | Original filter values from each caller |
| FILTERS_NORMALIZED | String | Normalized filter categories |
| FILTERS_CATEGORY | String | Filter categories |
| UNIFIED_FILTER | String | Unified filter status (PASS/FAIL) |
| PASSES_CONSENSUS | String | YES if variant passes consensus threshold |
| RESCUED | String | YES if included via cross-modality consensus |
| QUAL_MEAN/MIN/MAX | Float | Quality score statistics |
| CONSENSUS_GT | String | Consensus genotype across callers |
| GT_BY_CALLER | String | Pipe-separated genotypes from each caller |
| DP_MEAN/MIN/MAX | Float/Integer | Depth statistics |
| DP_BY_CALLER | String | Pipe-separated depth values |
| VAF_MEAN/MIN/MAX | Float | VAF statistics |
| VAF_BY_CALLER | String | Pipe-separated VAF values |

### Rescue VCF Additional Fields

Rescue VCF files include all consensus fields plus:

| Field | Type | Description |
|-------|------|-------------|
| MODALITIES | String | Pipe-separated modalities (DNA\|RNA) |
| CALLERS_BY_MODALITY | String | Callers grouped by modality (DNA:caller1,caller2\|RNA:caller3) |
| DNA_SUPPORT | Integer | Number of DNA callers supporting variant |
| RNA_SUPPORT | Integer | Number of RNA callers supporting variant |
| CROSS_MODALITY | String | YES if variant has both DNA and RNA support |
| DP_DNA_MEAN | Float | Mean depth across DNA callers |
| DP_RNA_MEAN | Float | Mean depth across RNA callers |
| VAF_DNA_MEAN | Float | Mean VAF across DNA callers |
| VAF_RNA_MEAN | Float | Mean VAF across RNA callers |


## Usage in Pipeline Scripts

### run_consensus_vcf.py

The consensus script uses vcf_utils for within-modality consensus:

```python
from vcf_utils.aggregation import read_variants_from_vcf, aggregate_variants
from vcf_utils.statistics import compute_consensus_statistics, print_statistics
from vcf_utils.io_utils import write_union_vcf

# Read variants from each caller
variant_collections = []
for caller, vcf_path in vcf_files.items():
    variants = read_variants_from_vcf(vcf_path, caller, 
                                     exclude_refcall=args.exclude_refcall,
                                     exclude_germline=args.exclude_germline)
    variant_collections.append((caller, variants, None))

# Aggregate
aggregated = aggregate_variants(variant_collections, args.snv_thr, args.indel_thr)

# Statistics
stats = compute_consensus_statistics(aggregated, args.snv_thr, args.indel_thr)
print_statistics(stats, 'consensus')

# Write output
write_union_vcf(aggregated, template_header, sample_name, out_file, 
                args.output_format, all_callers)
```

### run_rescue_vcf.py

The rescue script uses vcf_utils for cross-modality rescue:

```python
from vcf_utils.aggregation import read_variants_from_vcf, aggregate_variants
from vcf_utils.tagging import mark_rescued_variants, tag_variant_with_modality
from vcf_utils.statistics import compute_rescue_statistics, print_statistics
from vcf_utils.io_utils import write_union_vcf

# Read DNA and RNA consensus VCFs
dna_consensus = read_variants_from_vcf(args.dna_consensus, 'consensus', modality='DNA')
rna_consensus = read_variants_from_vcf(args.rna_consensus, 'consensus', modality='RNA')

# Read individual caller VCFs
all_collections = [('consensus', dna_consensus, 'DNA'), 
                   ('consensus', rna_consensus, 'RNA')]
# ... add DNA and RNA caller VCFs ...

# Aggregate across modalities
aggregated = aggregate_variants(all_collections, args.snv_thr, args.indel_thr)

# Mark rescued variants
dna_keys = set(dna_consensus.keys())
rna_keys = set(rna_consensus.keys())
aggregated = mark_rescued_variants(aggregated, dna_keys, rna_keys)

# Tag with modality
modality_map = {caller: mod for caller, _, mod in all_collections}
for vkey, data in aggregated.items():
    tag_variant_with_modality(data, modality_map)

# Statistics
stats = compute_rescue_statistics(aggregated, dna_keys, rna_keys)
print_statistics(stats, 'rescue')

# Write output with modality information
write_union_vcf(aggregated, template_header, sample_name, out_file,
                args.output_format, all_callers, modality_map)
```


## Testing

The package is tested through:

1. **Unit Tests**: Test individual functions in isolation
2. **Integration Tests**: Test complete workflows (consensus, rescue)
3. **Regression Tests**: Verify backward compatibility with original implementation

Example unit test:

```python
def test_aggregate_genotypes():
    genotypes = {
        'caller1': {'GT': '0/1', 'DP': 100, 'VAF': 0.25},
        'caller2': {'GT': '0/1', 'DP': 120, 'VAF': 0.27},
        'caller3': {'GT': '0/1', 'DP': 110, 'VAF': 0.26}
    }
    callers = ['caller1', 'caller2', 'caller3']
    
    agg = aggregate_genotypes(genotypes, callers)
    
    assert agg['consensus_gt'] == '0/1'
    assert agg['dp_mean'] == 110.0
    assert 0.25 <= agg['vaf_mean'] <= 0.27
```

## Performance Considerations

- **Streaming**: VCF reading uses cyvcf2 for fast, memory-efficient parsing
- **Lazy Evaluation**: Variants are processed on-demand where possible
- **Efficient Data Structures**: Uses dicts and sets for O(1) lookups
- **Minimal Copying**: Variant data is modified in-place when safe

## Contributing

When adding new functionality:

1. Add functions to the appropriate module based on responsibility
2. Include comprehensive docstrings with Args, Returns, and Examples
3. Update this README with usage examples
4. Add unit tests for new functions
5. Ensure backward compatibility with existing workflows

## License

This package is part of the rnadnavar pipeline and follows the same license.

## Version History

- **1.0.0**: Initial release with consensus and rescue support
  - Extracted from run_consensus_vcf.py
  - Added modality tracking for rescue workflow
  - Comprehensive documentation and examples
