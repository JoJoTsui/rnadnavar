# Variant Consensus and Cross-Modality Rescue Rules

## Overview
This document defines the unified rules for variant consensus calling within modalities (DNA/RNA) and cross-modality rescue operations used in the rnadnavar pipeline.

## Core Concepts

### Variant Classification
- **SNV**: Single nucleotide variants (point mutations)
- **Indel**: Insertions and deletions
- **Modality**: DNA or RNA sequencing data source
- **Caller**: Individual variant calling algorithm (e.g., Mutect2, Strelka2, DeepSomatic)

### Variant Key Generation
Variants are uniquely identified by:
```
variant_key = (chromosome, position, reference_allele, alternate_allele)
```

## Within-Modality Consensus Rules

### Consensus Thresholds
- **SNV Consensus Threshold**: Minimum number of callers required (default: 2)
- **Indel Consensus Threshold**: Minimum number of callers required (default: 2)

### Consensus Classification Logic
```python
def classify_consensus_variant(variant_data, snv_threshold, indel_threshold):
    caller_count = len(variant_data['callers'])
    variant_type = 'SNV' if is_snv(variant_data) else 'INDEL'
    
    if variant_type == 'SNV':
        is_consensus = caller_count >= snv_threshold
    else:  # INDEL
        is_consensus = caller_count >= indel_threshold
    
    return is_consensus
```

### Variant Filtering Rules
1. **RefCall Exclusion**: Variants marked as "RefCall" (DeepSomatic reference calls) can be excluded
2. **Germline Exclusion**: Variants marked as "GERMLINE" can be excluded
3. **Quality Filtering**: Applied at individual caller level before consensus

## Cross-Modality Rescue Rules

### Modality Assignment
- Each variant caller is assigned to either DNA or RNA modality
- Consensus VCFs inherit the modality of their constituent callers

### Cross-Modality Support Detection
```python
def has_cross_modality_support(variant_data, modality_map):
    dna_callers = [c for c in variant_data['callers'] if modality_map.get(c) == 'DNA']
    rna_callers = [c for c in variant_data['callers'] if modality_map.get(c) == 'RNA']
    
    has_dna_support = len(dna_callers) > 0
    has_rna_support = len(rna_callers) > 0
    
    return has_dna_support and has_rna_support
```

### Rescue Classification Logic
```python
def classify_rescued_variant(variant_key, variant_data, dna_consensus_keys, rna_consensus_keys):
    in_dna_consensus = variant_key in dna_consensus_keys
    in_rna_consensus = variant_key in rna_consensus_keys
    has_cross_modality = has_cross_modality_support(variant_data, modality_map)
    
    # Rescued: Has cross-modality support but failed consensus in at least one modality
    is_rescued = has_cross_modality and not (in_dna_consensus and in_rna_consensus)
    
    return {
        'cross_modality': has_cross_modality,
        'rescued': is_rescued,
        'in_dna_consensus': in_dna_consensus,
        'in_rna_consensus': in_rna_consensus
    }
```

## Variant Categories

### Within-Modality Categories
1. **Consensus Variants**: Meet caller threshold within modality
2. **Sub-consensus Variants**: Below caller threshold within modality

### Cross-Modality Categories
1. **DNA-only Variants**: Supported only by DNA callers
2. **RNA-only Variants**: Supported only by RNA callers
3. **Cross-modality Variants**: Supported by both DNA and RNA callers
4. **Rescued Variants**: Cross-modality variants that failed consensus in ≥1 modality

### Final Classification Matrix
```
| DNA Consensus | RNA Consensus | Cross-Modality | Category |
|---------------|---------------|----------------|----------|
| Yes           | Yes           | Yes            | Dual Consensus |
| Yes           | No            | Yes            | DNA Consensus + Rescued |
| No            | Yes           | Yes            | RNA Consensus + Rescued |
| No            | No            | Yes            | Rescued Only |
| Yes           | No            | No             | DNA-only Consensus |
| No            | Yes           | No             | RNA-only Consensus |
| No            | No            | No             | Sub-consensus |
```

## Aggregation Rules

### Caller Information Aggregation
For each variant, aggregate:
- **Caller List**: All callers supporting the variant
- **Genotype Information**: GT, AD, DP, VAF from each caller
- **Quality Metrics**: QUAL, FILTER status
- **Modality Counts**: Number of DNA vs RNA callers

### Genotype Consensus Rules
```python
def aggregate_genotypes(caller_genotypes):
    # Priority order for genotype selection
    gt_priority = {'1/1': 3, '0/1': 2, '1/0': 2, '0/0': 1, './.': 0}
    
    # Select highest priority genotype
    consensus_gt = max(caller_genotypes, key=lambda gt: gt_priority.get(gt, 0))
    
    # Aggregate depth and allelic depth
    total_dp = sum(dp for dp in depths if dp is not None)
    total_ad = aggregate_allelic_depths(allelic_depths)
    
    return consensus_gt, total_dp, total_ad
```

## Output Tagging Rules

### INFO Field Tags
- **CALLERS**: Comma-separated list of supporting callers
- **CALLER_COUNT**: Total number of supporting callers
- **DNA_CALLERS**: DNA-specific supporting callers
- **RNA_CALLERS**: RNA-specific supporting callers
- **CROSS_MODALITY**: Boolean flag for cross-modality support
- **RESCUED**: Boolean flag for rescued variants

### FILTER Field Rules
- **PASS**: Meets consensus threshold and quality criteria
- **LowSupport**: Below consensus threshold but has cross-modality support
- **SingleModality**: Supported by only one modality

## Mode-Specific Behaviors

### Union Mode (run_consensus_vcf.py)
- Keeps ALL variants from all callers
- Aggregates information across callers within same modality
- Applies consensus thresholds for classification only

### Rescue Mode (run_rescue_vcf.py)
- **Consensus-only Mode**: Merges only DNA and RNA consensus VCFs
- **Full Rescue Mode**: Includes individual callers for potential rescue
- Identifies variants with cross-modality support that failed single-modality consensus

## Quality Control Rules

### Minimum Requirements
- At least one supporting caller per variant
- Valid genomic coordinates
- Non-empty reference and alternate alleles

### Exclusion Criteria
- Variants in excluded regions (if specified)
- Variants below quality thresholds (caller-specific)
- RefCall variants (optional)
- Germline variants (optional)

## Statistical Reporting

### Consensus Statistics
- Total variants per modality
- Consensus vs sub-consensus counts
- SNV vs Indel breakdown
- Caller overlap statistics

### Rescue Statistics
- Total rescued variants
- Rescue rate (rescued / total cross-modality)
- Modality-specific rescue contributions
- Category distribution matrix

## Implementation Notes

### Performance Considerations
- Use variant keys for efficient lookups
- Batch process variants by chromosome
- Memory-efficient data structures for large VCFs

### Extensibility
- Modular design allows adding new modalities
- Configurable thresholds per variant type
- Pluggable filtering and annotation modules