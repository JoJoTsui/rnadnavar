# Rescue VCF Script Rules and Logic (run_rescue_vcf.py)

## Overview
The `run_rescue_vcf.py` script performs cross-modality variant aggregation between DNA and RNA consensus VCF files, identifying variants with support from both modalities and rescuing variants that may have failed consensus in one modality but have cross-modality support.

## Script Workflow

### 1. Input Processing
```python
# Command line arguments
--dna_consensus: Path to DNA consensus VCF file (required)
--rna_consensus: Path to RNA consensus VCF file (required)
--dna_vcf: Path to individual DNA caller VCF file (can be specified multiple times)
--rna_vcf: Path to individual RNA caller VCF file (can be specified multiple times)
--out_prefix: Output file prefix (required)
--output_format: Output format (vcf, vcf.gz, bcf, default: vcf.gz)
--snv_thr: Minimum number of callers for SNV consensus (default: 2)
--indel_thr: Minimum number of callers for indel consensus (default: 2)
--consensus_only: Only merge DNA and RNA consensus VCFs (no individual callers)
```

### 2. Input Validation
```python
# Validate required consensus files exist
dna_consensus_path = Path(args.dna_consensus)
rna_consensus_path = Path(args.rna_consensus)

if not dna_consensus_path.exists():
    print(f"Error: DNA consensus VCF not found: {dna_consensus_path}", file=sys.stderr)
    sys.exit(1)
if not rna_consensus_path.exists():
    print(f"Error: RNA consensus VCF not found: {rna_consensus_path}", file=sys.stderr)
    sys.exit(1)
```

### 3. Template Header and Sample Name Extraction
```python
# Get template header and sample name from DNA consensus VCF
template_vcf = VCF(str(dna_consensus_path))
template_header = template_vcf
sample_name = template_vcf.samples[0] if template_vcf.samples else 'SAMPLE'
```

### 4. Consensus VCF Reading with Modality Assignment
```python
# Read DNA consensus with modality='DNA'
dna_consensus = read_variants_from_vcf(
    str(dna_consensus_path),
    'dna_consensus',
    modality='DNA'
)

# Read RNA consensus with modality='RNA'
rna_consensus = read_variants_from_vcf(
    str(rna_consensus_path),
    'rna_consensus',
    modality='RNA'
)
```

### 5. Individual Caller VCF Processing (Full Rescue Mode Only)
```python
def read_caller_vcfs(vcf_paths, modality_name):
    caller_variants = {}
    if vcf_paths and not args.consensus_only:  # Skip if consensus_only mode
        for vcf_path in vcf_paths:
            caller_name = get_caller_name(Path(vcf_path).name)
            variants = read_variants_from_vcf(vcf_path, caller_name, modality=modality_name)
            caller_variants[caller_name] = variants
    return caller_variants

# Read DNA and RNA individual callers
dna_caller_variants = read_caller_vcfs(args.dna_vcf, 'DNA')
rna_caller_variants = read_caller_vcfs(args.rna_vcf, 'RNA')
```

### 6. Variant Collection Assembly with Modality Prefixing
```python
# Combine consensus and individual caller variant collections
# Use modality prefix as KEY to avoid collisions in modality_map
all_collections = [
    ('DNA_consensus', dna_consensus, 'DNA'),
    ('RNA_consensus', rna_consensus, 'RNA')
]

# Add individual DNA callers WITH prefix as key
for caller_name, variants in dna_caller_variants.items():
    all_collections.append((f'DNA_{caller_name}', variants, 'DNA'))

# Add individual RNA callers WITH prefix as key
for caller_name, variants in rna_caller_variants.items():
    all_collections.append((f'RNA_{caller_name}', variants, 'RNA'))
```

### 7. Variant Aggregation (Same as Consensus Mode)
```python
variant_data = aggregate_variants(
    all_collections,
    args.snv_thr,
    args.indel_thr
)
# Uses same aggregation logic as consensus mode
# Groups variants by genomic position and applies thresholds
```

### 8. Modality Map Creation
```python
# Create modality_map from all input sources
modality_map = {}
for caller_name, _, modality in all_collections:
    modality_map[caller_name] = modality

# Example modality_map:
# {
#     'DNA_consensus': 'DNA',
#     'RNA_consensus': 'RNA', 
#     'DNA_mutect2': 'DNA',
#     'DNA_strelka': 'DNA',
#     'RNA_deepsomatic': 'RNA'
# }
```

### 9. Modality Tagging
```python
def tag_variant_with_modality(variant_data, modality_map):
    for vkey, data in variant_data.items():
        modalities = set()
        dna_callers = []
        rna_callers = []
        caller_modality_map = {}
        
        for caller in data.get('callers', []):
            modality = modality_map.get(caller, 'UNKNOWN')
            modalities.add(modality)
            caller_modality_map[caller] = modality
            
            if modality == 'DNA':
                dna_callers.append(caller)
            elif modality == 'RNA':
                rna_callers.append(caller)
        
        data['modalities'] = modalities
        data['dna_callers'] = dna_callers
        data['rna_callers'] = rna_callers
        data['caller_modality_map'] = caller_modality_map
```

### 10. Rescue Variant Marking
```python
def mark_rescued_variants(variant_data, dna_consensus_keys, rna_consensus_keys):
    for vkey, data in variant_data.items():
        in_dna_consensus = vkey in dna_consensus_keys
        in_rna_consensus = vkey in rna_consensus_keys
        
        # Check cross-modality support
        has_dna_support = len([c for c in data['callers'] if modality_map.get(c) == 'DNA']) > 0
        has_rna_support = len([c for c in data['callers'] if modality_map.get(c) == 'RNA']) > 0
        has_cross_modality = has_dna_support and has_rna_support
        
        # Set support flags
        data['dna_support'] = in_dna_consensus
        data['rna_support'] = in_rna_consensus
        data['cross_modality'] = has_cross_modality
        
        # Mark as rescued if has cross-modality support
        data['rescued'] = has_cross_modality
    
    return variant_data
```

### 11. Rescue Statistics Computation
```python
def compute_rescue_statistics(variant_data, dna_variants, rna_variants):
    stats = {
        'total_variants': len(variant_data),
        'dna_only': 0, 'rna_only': 0, 'cross_modality': 0,
        'rescued': 0, 'snvs': 0, 'indels': 0,
        'snvs_rescued': 0, 'indels_rescued': 0
    }
    
    for vkey, data in variant_data.items():
        # Count modality support
        dna_support = data.get('dna_support', False)
        rna_support = data.get('rna_support', False)
        
        if dna_support and rna_support:
            stats['cross_modality'] += 1
        elif dna_support:
            stats['dna_only'] += 1
        elif rna_support:
            stats['rna_only'] += 1
        
        # Count rescued variants
        if data.get('rescued', False):
            stats['rescued'] += 1
            if data['is_snv']:
                stats['snvs_rescued'] += 1
            else:
                stats['indels_rescued'] += 1
        
        # Count variant types
        if data['is_snv']:
            stats['snvs'] += 1
        else:
            stats['indels'] += 1
    
    # Calculate rescue rate
    if stats['total_variants'] > 0:
        stats['rescue_rate'] = (stats['rescued'] / stats['total_variants']) * 100
    else:
        stats['rescue_rate'] = 0.0
    
    return stats
```

### 12. Unified Filter Computation (Rescue Mode) - UPDATED
```python
def compute_unified_classification_rescue(variant_data, modality_map):
    # Get consensus labels from DNA and RNA if present
    dna_label = None
    rna_label = None
    
    for i, caller in enumerate(variant_data['callers']):
        if caller.endswith('_consensus'):
            if 'DNA' in caller.upper():
                dna_label = variant_data['filters_normalized'][i]
            elif 'RNA' in caller.upper():
                rna_label = variant_data['filters_normalized'][i]
    
    # Apply NEW rescue consensus rules:
    if dna_label and rna_label:
        # Both DNA and RNA consensus available
        if dna_label == rna_label:
            # Same classification - use it
            return dna_label
        elif dna_label == 'Artifact' and rna_label == 'Artifact':
            # Both are Artifact
            return 'Artifact'
        elif dna_label != 'Artifact' and rna_label != 'Artifact':
            # Different non-Artifact classifications - inconsistent
            return 'Artifact'
        else:
            # One is Artifact, other is not - take DNA priority
            return dna_label
    elif dna_label and not rna_label:
        # Only DNA consensus available
        return dna_label
    elif rna_label and not dna_label:
        # Only RNA consensus available
        return rna_label
    else:
        # No consensus labels present
        return 'NoConsensus'
```

#### Key Changes:
1. **Enhanced Artifact Logic**: Cross-modality inconsistency (e.g., DNA=Somatic, RNA=Germline) → Artifact
2. **NoConsensus Category**: No consensus labels available → NoConsensus (instead of Artifact)
3. **Consistent Artifact Handling**: Both DNA and RNA having Artifact → Artifact

### 13. Output VCF Generation with Rescue Fields
```python
def write_union_vcf(variant_data, template_header, sample_name, out_file, 
                    output_format, all_callers, modality_map):
    # Create header with rescue fields enabled
    output_header = create_output_header(template_header, sample_name, 
                                       include_rescue_fields=True)
    
    for vkey, data in sorted_variants:
        record = vcf_out.new_record(...)
        
        # Add all consensus INFO fields PLUS rescue-specific fields
        
        # Modality-specific caller counts
        total_dna_callers = sum(1 for c in all_callers if modality_map.get(c) == 'DNA')
        total_rna_callers = sum(1 for c in all_callers if modality_map.get(c) == 'RNA')
        record.info['N_DNA_CALLERS'] = total_dna_callers
        record.info['N_RNA_CALLERS'] = total_rna_callers
        
        # Modality-specific support counts
        dna_callers_support = sum(1 for c in actual_callers_in_variant if modality_map.get(c) == 'DNA')
        rna_callers_support = sum(1 for c in actual_callers_in_variant if modality_map.get(c) == 'RNA')
        record.info['N_DNA_CALLERS_SUPPORT'] = dna_callers_support
        record.info['N_RNA_CALLERS_SUPPORT'] = rna_callers_support
        
        # Modality information
        modalities = set()
        dna_callers = []
        rna_callers = []
        
        for caller in data['callers']:
            if not is_consensus_caller(caller):  # Skip consensus callers
                modality = modality_map.get(caller, 'UNKNOWN')
                modalities.add(modality)
                if modality == 'DNA':
                    dna_callers.append(caller)
                elif modality == 'RNA':
                    rna_callers.append(caller)
        
        record.info['MODALITIES'] = '|'.join(sorted(modalities)) if modalities else '.'
        record.info['DNA_SUPPORT'] = len(set(dna_callers))
        record.info['RNA_SUPPORT'] = len(set(rna_callers))
        record.info['CROSS_MODALITY'] = 'YES' if len(modalities) > 1 else 'NO'
        
        # Modality-specific statistics
        # Calculate DNA and RNA mean DP/VAF separately
        
        # Set RESCUED flag based on actual rescued status
        record.info['RESCUED'] = 'YES' if data.get('rescued', False) else 'NO'
```

## Key Rescue Mode Characteristics

### Cross-Modality Analysis
- **Processes both DNA and RNA** consensus VCFs with modality assignment
- **Identifies cross-modality variants** supported by both DNA and RNA callers
- **Rescues variants** that have cross-modality support regardless of individual consensus status

### Modality Prefixing System
- Caller names prefixed with modality: `DNA_mutect2`, `RNA_strelka`
- Consensus callers: `DNA_consensus`, `RNA_consensus`
- Prevents naming collisions in aggregation

### Rescue Logic
```python
# A variant is "rescued" if it has cross-modality support
rescued = has_dna_support and has_rna_support
```

### Two Operating Modes

#### Consensus-Only Mode (`--consensus_only`)
- Only merges DNA and RNA consensus VCFs
- No individual caller VCFs processed
- Final counts cannot exceed DNA + RNA consensus sums
- Ensures conservative rescue approach

#### Full Rescue Mode (Default)
- Includes individual caller VCFs from both modalities
- May rescue variants that failed consensus in one modality
- Can identify additional cross-modality variants

### Unified Filter Priority (Rescue Mode) - UPDATED
1. If DNA and RNA consensus agree → use agreed classification
2. If both DNA and RNA are 'Artifact' → Artifact
3. If DNA and RNA have different non-Artifact classifications → **Artifact** (inconsistent)
4. If only one modality has consensus → use that classification  
5. If one is Artifact and other is not → **DNA takes priority**
6. If no consensus labels → **NoConsensus** (instead of Artifact)

### Additional INFO Fields (Rescue Mode Only)
- **Modality tracking**: `MODALITIES`, `CALLERS_BY_MODALITY`
- **Support counts**: `DNA_SUPPORT`, `RNA_SUPPORT`, `CROSS_MODALITY`
- **Consensus flags**: `PASSES_CONSENSUS_DNA`, `PASSES_CONSENSUS_RNA`
- **Unified filters**: `UNIFIED_FILTER_DNA`, `UNIFIED_FILTER_RNA`
- **Statistics**: `DP_DNA_MEAN`, `DP_RNA_MEAN`, `VAF_DNA_MEAN`, `VAF_RNA_MEAN`

### Rescue Statistics Reported
- **Total variants**: All aggregated variants
- **Modality distribution**: DNA-only, RNA-only, cross-modality
- **Rescue effectiveness**: Total rescued, rescue rate percentage
- **Variant type breakdown**: SNVs vs indels rescued

### Key Differences from Consensus Mode
1. **Modality awareness**: Tracks DNA vs RNA caller support
2. **Cross-modality rescue**: Identifies variants with both DNA and RNA support
3. **Deterministic classification**: Uses DNA/RNA consensus hierarchy instead of majority vote
4. **Extended INFO fields**: Includes modality-specific metadata
5. **Rescue statistics**: Reports cross-modality effectiveness metrics