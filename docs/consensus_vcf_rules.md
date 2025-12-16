# Consensus VCF Script Rules and Logic (run_consensus_vcf.py)

## Overview
The `run_consensus_vcf.py` script performs within-modality variant consensus calling by aggregating variants from multiple callers within the same modality (DNA or RNA) and determining which variants have sufficient caller support to be considered consensus variants.

## Script Workflow

### 1. Input Processing
```python
# Command line arguments
--input_dir: Directory containing input VCF files
--out_prefix: Prefix for output files
--snv_thr: Number of callers required for SNV consensus (default: 2)
--indel_thr: Number of callers required for indel consensus (default: 2)
--output_format: Output format (vcf, vcf.gz, bcf, default: vcf.gz)
--sample_name: Sample name (optional, uses first sample from VCF if not specified)
--exclude_refcall: Exclude variants marked as RefCall (DeepSomatic reference calls)
--exclude_germline: Exclude variants marked as GERMLINE
```

### 2. VCF File Discovery
```python
# Find VCF files in input directory
for vcf_path in sorted(input_dir.glob("*.vcf*")):
    if vcf_path.suffix in ['.vcf', '.gz'] or vcf_path.name.endswith('.vcf.gz'):
        caller = get_caller_name(str(vcf_path))  # Extract caller name from filename
        vcf_files[caller] = str(vcf_path)
```

### 3. Caller Name Extraction Logic
```python
def get_caller_name(filename):
    name = Path(filename).name
    # Handle .variants. pattern (e.g., sample.strelka.variants.vcf.gz)
    if '.variants.' in name:
        parts = name.split('.')
        for i, part in enumerate(parts):
            if part == 'variants' and i > 0:
                return parts[i-1]  # Return caller before 'variants'
    # Handle .consensus. pattern
    if '.consensus.' in name:
        return 'consensus'
    # Fallback: take second component
    parts = name.split('.')
    return parts[1] if len(parts) > 1 else parts[0]
```

### 4. Variant Reading and Classification
```python
# For each VCF file
for caller, vcf_path in vcf_files.items():
    variants = read_variants_from_vcf(
        vcf_path, 
        caller, 
        modality=None,  # No modality assignment in consensus mode
        exclude_refcall=args.exclude_refcall,
        exclude_germline=args.exclude_germline
    )
```

#### Variant Classification Process
Each variant is classified into biological categories:
- **Strelka**: Uses FILTER + INFO/NT + normal depth
- **DeepSomatic**: Uses explicit FILTER labels (PASS→Somatic, GERMLINE→Germline, RefCall→Reference)
- **Mutect2**: Primarily PASS→Somatic, other filters→Artifact
- **Consensus**: FILTER field already contains biological category

### 5. Variant Aggregation Logic
```python
def aggregate_variants(variant_collections, snv_threshold, indel_threshold):
    # Group variants by genomic position (chrom:pos:ref:alt)
    aggregated = defaultdict(lambda: {
        'CHROM': None, 'POS': None, 'REF': None, 'ALT': None,
        'is_snv': None, 'callers': [], 'filters_original': [],
        'filters_normalized': [], 'filters_category': [],
        'qualities': [], 'genotypes': {}, 'ids': []
    })
    
    # Process each variant collection
    for caller_name, variant_dict, modality in variant_collections:
        for vkey, variant_data in variant_dict.items():
            data = aggregated[vkey]
            # Store variant info and caller-specific data
            data['callers'].append(caller_name)
            data['filters_original'].append(variant_data['filter_original'])
            data['filters_normalized'].append(variant_data['filter_normalized'])
            # ... store other fields
    
    # Apply consensus thresholds
    for vkey, data in aggregated.items():
        n_support = len(set(data['callers']))  # Unique callers
        if data['is_snv']:
            data['passes_consensus'] = n_support >= snv_threshold
        else:  # Indel
            data['passes_consensus'] = n_support >= indel_threshold
        
        # Aggregate genotype information
        data['gt_aggregated'] = aggregate_genotypes(data['genotypes'], data['callers'])
```

### 6. Genotype Aggregation Rules
```python
def aggregate_genotypes(genotypes_by_caller, callers_order):
    # Consensus genotype: most common GT across callers
    gt_counts = defaultdict(int)
    for caller in callers_order:
        info = genotypes_by_caller.get(caller, {})
        if info and 'GT' in info and info['GT'] is not None:
            gt_counts[info['GT']] += 1
    
    consensus_gt = max(gt_counts, key=gt_counts.get) if gt_counts else None
    
    # Calculate statistics for DP and VAF
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

### 7. Unified Filter Computation (Consensus Mode) - UPDATED
```python
def compute_unified_classification_consensus(variant_data, snv_threshold, indel_threshold):
    # 1. Check consensus threshold first
    if not variant_data.get('passes_consensus', False):
        return 'NoConsensus'
    
    # 2. Get individual caller classifications (exclude consensus callers)
    individual_filters = [
        variant_data['filters_normalized'][i] 
        for i, caller in enumerate(variant_data['callers']) 
        if not caller.endswith('_consensus')
    ]
    
    if not individual_filters:
        return 'Artifact'
    
    # 3. NEW: Check for inconsistent classifications if ≥2 callers
    if len(individual_filters) >= 2:
        unique_classifications = set(individual_filters)
        # If multiple different classifications, mark as Artifact
        if len(unique_classifications) > 1:
            return 'Artifact'
    
    # 4. Use majority vote with priority: Somatic > Germline > Reference > Artifact
    classification_counts = Counter(individual_filters)
    max_count = max(classification_counts.values())
    most_common = [cls for cls, count in classification_counts.items() if count == max_count]
    
    priority = ['Somatic', 'Germline', 'Reference', 'Artifact']
    for cls in priority:
        if cls in most_common:
            return cls
    
    return most_common[0]  # Fallback
```

#### Key Changes:
1. **NoConsensus Priority**: Variants not meeting consensus threshold are immediately classified as NoConsensus
2. **Inconsistency Detection**: Variants with ≥2 callers having different classifications are marked as Artifact
3. **Enhanced Artifact Category**: Now includes both low-quality variants AND inconsistent classifications

### 8. Statistics Computation
```python
def compute_consensus_statistics(variant_data, snv_threshold, indel_threshold):
    stats = {
        'total_variants': len(variant_data),
        'snvs': 0, 'indels': 0,
        'snvs_consensus': 0, 'indels_consensus': 0,
        'single_caller': 0, 'multi_caller': 0
    }
    
    for vkey, data in variant_data.items():
        if data['is_snv']:
            stats['snvs'] += 1
            if data.get('passes_consensus', False):
                stats['snvs_consensus'] += 1
        else:
            stats['indels'] += 1
            if data.get('passes_consensus', False):
                stats['indels_consensus'] += 1
        
        unique_callers = len(set(data['callers']))
        if unique_callers == 1:
            stats['single_caller'] += 1
        else:
            stats['multi_caller'] += 1
    
    return stats
```

### 9. Output VCF Generation
```python
def write_union_vcf(variant_data, template_header, sample_name, out_file, 
                    output_format, all_callers, modality_map=None):
    # Create header with consensus INFO fields (no rescue fields)
    output_header = create_output_header(template_header, sample_name, 
                                       include_rescue_fields=False)
    
    # Sort variants by genomic position
    sorted_variants = sorted(variant_data.items(), key=sort_key)
    
    for vkey, data in sorted_variants:
        record = vcf_out.new_record(...)
        
        # Add consensus-specific INFO fields
        record.info['N_CALLERS'] = len(all_callers)
        record.info['CALLERS'] = '|'.join(all_callers)
        record.info['N_SUPPORT_CALLERS'] = len(set(data['callers']))
        record.info['CALLERS_SUPPORT'] = '|'.join(data['callers'])
        record.info['PASSES_CONSENSUS'] = 'YES' if data['passes_consensus'] else 'NO'
        record.info['RESCUED'] = 'NO'  # Always NO in consensus mode
        
        # Set FILTER field to unified biological classification
        record.filter.clear()
        record.filter.add(unified_classification)
        
        # Override with NoConsensus if variant doesn't meet threshold
        if not data['passes_consensus']:
            record.filter.clear()
            record.filter.add('NoConsensus')
```

## Key Consensus Mode Characteristics

### Union Behavior
- **Keeps ALL variants** from all callers regardless of consensus status
- Variants below threshold are marked with `PASSES_CONSENSUS=NO` and `FILTER=NoConsensus`
- No variants are excluded based on consensus thresholds

### Single Modality Processing
- No modality assignment or cross-modality analysis
- All callers treated as same modality
- No modality-specific INFO fields generated

### Biological Classification Priority - UPDATED
1. **Somatic** (highest priority in tie-breaking)
2. **Germline**
3. **Reference** 
4. **Artifact** (includes inconsistent classifications)
5. **NoConsensus** (lowest priority - threshold not met)

### INFO Fields Generated (Consensus Mode)
- Caller tracking: `N_CALLERS`, `CALLERS`, `N_SUPPORT_CALLERS`, `CALLERS_SUPPORT`
- Filter information: `FILTERS_ORIGINAL`, `FILTERS_NORMALIZED`, `UNIFIED_FILTER`
- Consensus status: `PASSES_CONSENSUS` (always `RESCUED=NO`)
- Quality stats: `QUAL_MEAN/MIN/MAX`
- Genotype stats: `CONSENSUS_GT`, `GT_BY_CALLER`, `DP_MEAN/MIN/MAX`, `VAF_MEAN/MIN/MAX`

### FILTER Values Used - UPDATED
- **Somatic**: High-confidence somatic variants
- **Germline**: Germline variants
- **Reference**: Reference calls
- **Artifact**: Low-quality variants OR variants with inconsistent classifications across callers
- **NoConsensus**: Variants below consensus threshold or not fitting other categories