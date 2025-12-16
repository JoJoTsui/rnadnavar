# Variant Consensus and Cross-Modality Rescue Rules

## Overview
This document defines the exact rules and logic for variant consensus calling within modalities (DNA/RNA) and cross-modality rescue operations as implemented in the rnadnavar pipeline scripts.

## Core Concepts

### Variant Classification
- **SNV**: Single nucleotide variants where both REF and all ALT alleles are single bases
- **Indel**: Insertions, deletions, or complex variants (non-SNV)
- **Modality**: DNA or RNA sequencing data source
- **Caller**: Individual variant calling algorithm (mutect2, strelka, deepsomatic)
- **Consensus Caller**: Aggregated results from within-modality consensus (dna_consensus, rna_consensus)

### Variant Key Generation
Variants are uniquely identified by normalized genomic coordinates:
```python
def variant_key(variant):
    chrom = normalize_chromosome(variant.CHROM)  # Remove 'chr' prefix
    pos = variant.POS  # 1-based position
    ref = variant.REF
    alts = ','.join(variant.ALT) if variant.ALT else '.'
    return f"{chrom}:{pos}:{ref}:{alts}"
```

### SNV Detection Logic
```python
def is_snv(ref, alt_list):
    if len(ref) != 1:
        return False
    for alt in alt_list:
        if len(alt) != 1 or alt == '*':
            return False
    return True
```

## Biological Classification System

### Classification Categories
All variants are classified into exactly one biological category:
1. **Somatic**: High-confidence somatic variants specific to tumor
2. **Germline**: Germline variants detected in normal sample  
3. **Reference**: Reference calls (no variant detected)
4. **Artifact**: Low quality variants or technical artifacts

### Caller-Specific Classification Rules

#### Strelka Classification
```python
def classify_strelka_variant(filter_val, nt_val, normal_dp):
    # 1. Somatic: PASS filter (or None/.)
    if filter_val == "PASS" or filter_val is None or filter_val == ".":
        return "Somatic"
    
    # 2. Germline: Filter failed but NT indicates variant in Normal
    if nt_val in ["het", "hom"]:
        if normal_dp >= 2:  # Sufficient normal depth
            return "Germline"
        return "Artifact"
    
    # 3. Reference: Filter failed but NT indicates Normal is Reference
    if nt_val == "ref":
        if normal_dp >= 2:  # Sufficient normal depth
            return "Reference"
        return "Artifact"
    
    # 4. Artifact: Everything else
    return "Artifact"
```

#### DeepSomatic Classification
```python
def classify_deepsomatic_variant(filter_val):
    if filter_val is None or filter_val == "." or filter_val == "PASS":
        return "Somatic"
    
    filter_upper = str(filter_val).upper()
    if "GERMLINE" in filter_upper:
        return "Germline"
    if "REFCALL" in filter_upper or "REF_CALL" in filter_upper:
        return "Reference"
    
    return "Artifact"  # All other filters
```

#### Mutect2 Classification
```python
def classify_mutect2_variant(filter_val, info_dict=None):
    if filter_val is None or filter_val == "." or filter_val == "PASS":
        return "Somatic"
    
    filter_upper = str(filter_val).upper()
    if "GERMLINE" in filter_upper or "NORMAL_ARTIFACT" in filter_upper:
        return "Germline"
    
    return "Artifact"  # All other filters
```

#### Consensus Caller Classification
For consensus VCFs, the FILTER field already contains the biological category:
```python
def classify_consensus_variant(filter_val):
    if filter_val in ["Somatic", "Germline", "Reference", "Artifact"]:
        return filter_val
    elif filter_val == "NoConsensus":
        return "Artifact"  # Shouldn't appear in consensus VCFs being rescued
    elif filter_val is None or filter_val == "." or filter_val == "PASS":
        return "Somatic"  # Legacy PASS value
    else:
        return "Artifact"  # Unknown filter
```

## Within-Modality Consensus Rules

### Consensus Thresholds
- **SNV Consensus Threshold**: Minimum number of callers required (default: 2)
- **Indel Consensus Threshold**: Minimum number of callers required (default: 2)

### Consensus Determination Logic
```python
def determine_consensus(variant_data, snv_threshold, indel_threshold):
    unique_callers = len(set(variant_data['callers']))
    
    if variant_data['is_snv']:
        passes_consensus = unique_callers >= snv_threshold
    else:  # Indel
        passes_consensus = unique_callers >= indel_threshold
    
    return passes_consensus
```

### Unified Filter Computation (Consensus Mode)
For within-modality consensus, use majority vote across individual callers:
```python
def compute_unified_filter_consensus(variant_data):
    # Exclude consensus callers, use only individual callers
    individual_filters = [
        variant_data['filters_normalized'][i] 
        for i, caller in enumerate(variant_data['callers']) 
        if not caller.endswith('_consensus')
    ]
    
    if not individual_filters:
        return 'Artifact'
    
    # Count classifications
    classification_counts = Counter(individual_filters)
    max_count = max(classification_counts.values())
    most_common = [cls for cls, count in classification_counts.items() if count == max_count]
    
    # Break ties using priority: Somatic > Germline > Reference > Artifact
    priority = ['Somatic', 'Germline', 'Reference', 'Artifact']
    for cls in priority:
        if cls in most_common:
            return cls
    
    return most_common[0]  # Fallback
```

### Variant Filtering Rules
1. **RefCall Exclusion**: Skip variants with 'RefCall' in FILTER field (DeepSomatic)
2. **Germline Exclusion**: Skip variants with 'GERMLINE' in FILTER field
3. **Quality Filtering**: Applied at individual caller level before aggregation

## Cross-Modality Rescue Rules (Full Rescue Mode Only)

### Modality Assignment
- Each caller is assigned to DNA or RNA modality via `modality_map`
- Caller names are prefixed with modality in output: `DNA_mutect2`, `RNA_strelka`
- Consensus callers: `DNA_consensus`, `RNA_consensus`

### Cross-Modality Support Detection
```python
def has_cross_modality_support(variant_data, modality_map):
    dna_callers = [c for c in variant_data['callers'] if modality_map.get(c) == 'DNA']
    rna_callers = [c for c in variant_data['callers'] if modality_map.get(c) == 'RNA']
    
    return len(dna_callers) > 0 and len(rna_callers) > 0
```

### Rescue Classification Logic
```python
def mark_rescued_variants(variant_data, dna_consensus_keys, rna_consensus_keys):
    for vkey, data in variant_data.items():
        in_dna_consensus = vkey in dna_consensus_keys
        in_rna_consensus = vkey in rna_consensus_keys
        has_cross_modality = has_cross_modality_support(data, modality_map)
        
        # Cross-modality: supported by both DNA and RNA callers
        data['cross_modality'] = has_cross_modality
        
        # Rescued: has cross-modality support (regardless of consensus status)
        data['rescued'] = has_cross_modality
        
        data['dna_support'] = in_dna_consensus
        data['rna_support'] = in_rna_consensus
    
    return variant_data
```

### Unified Filter Computation (Rescue Mode)
For cross-modality rescue, use deterministic DNA/RNA consensus rule:
```python
def compute_unified_filter_rescue(variant_data):
    # Get consensus labels from DNA and RNA if present
    dna_label = None
    rna_label = None
    
    for i, caller in enumerate(variant_data['callers']):
        if caller.endswith('_consensus'):
            if 'DNA' in caller.upper():
                dna_label = variant_data['filters_normalized'][i]
            elif 'RNA' in caller.upper():
                rna_label = variant_data['filters_normalized'][i]
    
    # Apply rescue consensus rules:
    # 1) If DNA and RNA have same label, use it
    # 2) If one modality has label and other doesn't, use the labeled modality
    # 3) If different labels, take DNA modality's consensus label
    if dna_label and rna_label:
        if dna_label == rna_label:
            return dna_label
        else:
            return dna_label  # DNA takes priority
    elif dna_label and not rna_label:
        return dna_label
    elif rna_label and not dna_label:
        return rna_label
    else:
        return 'Artifact'  # No consensus labels present
```

## Genotype Aggregation Rules

### Genotype Information Extraction
For each caller, extract:
- **GT**: Genotype (0/1, 1/1, etc.)
- **DP**: Total depth
- **AD**: Allelic depths (ref,alt)
- **VAF**: Variant allele frequency
- **GQ**: Genotype quality

### Special Strelka Handling
```python
def extract_strelka_genotype(variant):
    # Use TAR/TIR for indels, AU/CU/GU/TU for SNVs
    # Calculate AD and VAF from tier1 counts
    # Use SGT field for genotype if GT missing
```

### Genotype Consensus Rules
```python
def aggregate_genotypes(genotypes_by_caller, callers_order):
    # Consensus genotype: most common GT across callers
    gt_counts = Counter([info['GT'] for info in genotypes_by_caller.values() if info['GT']])
    consensus_gt = max(gt_counts, key=gt_counts.get) if gt_counts else None
    
    # Statistics: mean, min, max for DP and VAF
    dp_values = [info['DP'] for info in genotypes_by_caller.values() if info['DP'] is not None]
    vaf_values = [info['VAF'] for info in genotypes_by_caller.values() if info['VAF'] is not None]
    
    return {
        'consensus_gt': consensus_gt,
        'dp_mean': mean(dp_values) if dp_values else None,
        'dp_min': min(dp_values) if dp_values else None,
        'dp_max': max(dp_values) if dp_values else None,
        'vaf_mean': mean(vaf_values) if vaf_values else None,
        'vaf_min': min(vaf_values) if vaf_values else None,
        'vaf_max': max(vaf_values) if vaf_values else None,
        'gt_by_caller': [genotypes_by_caller.get(c, {}).get('GT', '.') for c in callers_order],
        'dp_by_caller': [genotypes_by_caller.get(c, {}).get('DP') for c in callers_order],
        'vaf_by_caller': [genotypes_by_caller.get(c, {}).get('VAF') for c in callers_order]
    }
```

## Output INFO Field Tags

### Consensus Mode INFO Fields
- **N_CALLERS**: Total number of variant callers in analysis (excludes consensus)
- **CALLERS**: All variant callers (pipe-separated, excludes consensus)
- **N_SUPPORT_CALLERS**: Number of callers that detected this variant (excludes consensus)
- **CALLERS_SUPPORT**: Callers that detected this variant (pipe-separated, excludes consensus)
- **N_CONSENSUS_SUPPORT**: Number of consensus callers that detected this variant
- **CONSENSUS_SUPPORT**: Consensus callers that detected this variant (pipe-separated)
- **FILTERS_ORIGINAL**: Original filter values by caller (format: caller:filter|...)
- **FILTERS_NORMALIZED**: Normalized biological classifications by caller
- **FILTERS_CATEGORY**: Same as FILTERS_NORMALIZED (biological categories)
- **UNIFIED_FILTER**: Overall biological classification (majority vote)
- **PASSES_CONSENSUS**: Whether variant meets consensus threshold (YES/NO)
- **RESCUED**: Whether variant included via cross-modality consensus (YES/NO)
- **QUAL_MEAN/MIN/MAX**: Quality score statistics
- **CONSENSUS_GT**: Consensus genotype across callers
- **GT_BY_CALLER**: Genotypes by caller (format: caller:GT|...)
- **DP_MEAN/MIN/MAX**: Depth statistics
- **DP_BY_CALLER**: Depth values by caller (format: caller:DP|...)
- **VAF_MEAN/MIN/MAX**: VAF statistics  
- **VAF_BY_CALLER**: VAF values by caller (format: caller:VAF|...)

### Rescue Mode Additional INFO Fields
- **N_DNA_CALLERS**: Total number of DNA callers in analysis
- **N_RNA_CALLERS**: Total number of RNA callers in analysis
- **N_DNA_CALLERS_SUPPORT**: Number of DNA callers that detected this variant
- **N_RNA_CALLERS_SUPPORT**: Number of RNA callers that detected this variant
- **UNIFIED_FILTER_DNA**: DNA-specific biological classification
- **UNIFIED_FILTER_RNA**: RNA-specific biological classification
- **PASSES_CONSENSUS_DNA**: Whether variant passes DNA consensus (YES/NO)
- **PASSES_CONSENSUS_RNA**: Whether variant passes RNA consensus (YES/NO)
- **MODALITIES**: Modalities where variant detected (pipe-separated)
- **CALLERS_BY_MODALITY**: Callers grouped by modality (format: DNA:caller1,caller2|RNA:caller3)
- **DNA_SUPPORT**: Number of DNA callers supporting variant
- **RNA_SUPPORT**: Number of RNA callers supporting variant
- **CROSS_MODALITY**: Whether variant has cross-modality support (YES/NO)
- **DP_DNA_MEAN**: Mean depth across DNA callers
- **DP_RNA_MEAN**: Mean depth across RNA callers
- **VAF_DNA_MEAN**: Mean VAF across DNA callers
- **VAF_RNA_MEAN**: Mean VAF across RNA callers

## Output FILTER Field Rules

### Biological Category FILTER Values
- **Somatic**: High-confidence somatic variant specific to tumor
- **Germline**: Germline variant detected in normal sample
- **Reference**: Reference call - no variant detected
- **Artifact**: Low quality variant or technical artifact
- **NoConsensus**: Does not meet consensus threshold (overrides biological classification)

### FILTER Assignment Logic
1. **Consensus Mode**: Use unified biological classification from majority vote
2. **Rescue Mode**: Use deterministic DNA/RNA consensus rule
3. **NoConsensus Override**: If variant doesn't meet consensus threshold, set FILTER=NoConsensus

## Statistical Reporting

### Consensus Statistics
```python
def compute_consensus_statistics(variant_data, snv_threshold, indel_threshold):
    return {
        'total_variants': len(variant_data),
        'snvs': count_snvs,
        'indels': count_indels,
        'snvs_consensus': count_snvs_passing_consensus,
        'indels_consensus': count_indels_passing_consensus,
        'single_caller': count_single_caller_variants,
        'multi_caller': count_multi_caller_variants
    }
```

### Rescue Statistics
```python
def compute_rescue_statistics(variant_data, dna_variants, rna_variants):
    return {
        'total_variants': len(variant_data),
        'dna_only': count_dna_only_variants,
        'rna_only': count_rna_only_variants,
        'cross_modality': count_cross_modality_variants,
        'rescued': count_rescued_variants,
        'snvs': count_snvs,
        'indels': count_indels,
        'snvs_rescued': count_snvs_rescued,
        'indels_rescued': count_indels_rescued,
        'rescue_rate': (rescued / total) * 100
    }
```

## Implementation Notes

### Caller Name Prefixing
- In rescue mode, caller names are prefixed with modality: `DNA_mutect2`, `RNA_strelka`
- Consensus callers: `DNA_consensus`, `RNA_consensus`
- Individual caller INFO fields exclude consensus callers to avoid double-counting

### Chromosome Sorting
- Chromosomes normalized by removing 'chr' prefix
- Sorting order: 1-22, X=23, Y=24, M/MT=25, then contig order from VCF header

### Error Handling
- Missing genotype fields default to None
- Classification errors default to "Artifact"
- Unknown filters kept but sanitized (spaces → underscores)

### Performance Optimizations
- Variants sorted by genomic position before writing
- Progress printed every 10,000 variants
- Memory-efficient aggregation using defaultdict

## Recent Fixes Applied

### Fixed Issues (Resolved)

1. **✅ RESCUED Flag Logic** - Fixed to use actual `rescued` status from variant data instead of checking for 'consensus' string in caller names
2. **✅ Consensus Caller Detection** - Standardized naming to `DNA_consensus`/`RNA_consensus` and fixed detection logic
3. **✅ Modality Prefix Handling** - Added safe prefixing function to prevent double-prefixing (e.g., `DNA_DNA_mutect2`)
4. **✅ Filter Field Handling** - Replaced semicolon substitution with proper VCF escaping using URL encoding
5. **✅ Input Validation** - Added comprehensive parameter validation for both scripts
6. **✅ Error Handling** - Improved genotype extraction error handling with specific exception types
7. **✅ Rescue Rate Calculation** - Fixed to calculate rate relative to cross-modality variants instead of total variants

## Remaining Issues (Lower Priority)

### 5. Inconsistent VAF Calculation
**Location**: `extract_genotype_info()` function
**Issue**: VAF calculation from AD may not handle edge cases:
```python
if info['VAF'] is None and info['AD'] is not None:
    try:
        ad_values = [int(x) for x in info['AD'].split(',')]
        if len(ad_values) >= 2 and sum(ad_values) > 0:
            info['VAF'] = ad_values[1] / sum(ad_values)
    except Exception:
        pass
```
**Problem**: Assumes AD format is always "ref,alt" but some callers may have different formats or multiple alt alleles.

### 6. Potential Memory Issue with Large VCFs
**Location**: `aggregate_variants()` function
**Issue**: All variants are loaded into memory simultaneously using `defaultdict`.

**Problem**: For very large VCFs (millions of variants), this could cause memory issues.

**Recommendation**: Consider streaming processing for large datasets.

### 7. Chromosome Sorting Edge Case
**Location**: `write_union_vcf()` function
**Issue**: Chromosome sorting logic has potential issues:
```python
try:
    if chrom.isdigit():
        chrom_idx = int(chrom)
    elif chrom == 'X':
        chrom_idx = 23
    # ...
except:
    pass
```
**Problem**: The try-catch is too broad and may hide legitimate errors. Also, doesn't handle non-standard chromosome names well.

## Implementation Notes

### Standardized Naming Convention
- **Consensus callers**: `DNA_consensus`, `RNA_consensus` (standardized)
- **Individual callers**: `DNA_mutect2`, `RNA_strelka` (with safe prefixing)
- **Detection logic**: Uses standardized patterns for reliable identification

### Enhanced Error Handling
- **Specific exceptions**: Catch `KeyError`, `IndexError`, `ValueError`, `TypeError` separately
- **Graceful degradation**: Continue processing with warnings for non-critical errors
- **Informative messages**: Include context about which caller/file caused issues

### VCF Standard Compliance
- **Proper escaping**: Use URL encoding (`%3D`, `%2C`) instead of character replacement
- **Standard delimiters**: Preserve semicolons as compound filter delimiters
- **Valid format**: Maintain VCF specification compliance