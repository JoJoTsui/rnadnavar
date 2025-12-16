# Rescue VCF Script Issues and Bug Analysis

## Critical Issues

### 1. Consensus Caller Detection Bug
**Location**: `write_union_vcf()` function in `io_utils.py`
**Severity**: Critical
**Issue**: 
```python
# Current implementation (INCORRECT)
has_dna_consensus = any(c == 'dna_consensus' for c in data['callers'])
has_rna_consensus = any(c == 'rna_consensus' for c in data['callers'])
```
**Problem**: 
- Looks for `dna_consensus`/`rna_consensus` but actual caller names are `DNA_consensus`/`RNA_consensus`
- This causes consensus flags to always be False
- Results in incorrect `PASSES_CONSENSUS_DNA`/`PASSES_CONSENSUS_RNA` values

**Correct Implementation**:
```python
has_dna_consensus = any('DNA' in c.upper() and c.endswith('_consensus') for c in data['callers'])
has_rna_consensus = any('RNA' in c.upper() and c.endswith('_consensus') for c in data['callers'])
```

### 2. Modality Prefix Double-Prefixing Risk
**Location**: `main()` function in `run_rescue_vcf.py`
**Severity**: High
**Issue**: 
```python
# Risk of double prefixing
all_collections.append((f'DNA_{caller_name}', variants, 'DNA'))
```
**Problem**: 
- If caller is already named with modality prefix (e.g., `DNA_mutect2`), creates `DNA_DNA_mutect2`
- Could happen if VCF files are pre-processed with modality prefixes
- Breaks modality mapping logic

**Recommendation**:
```python
def safe_prefix_caller(caller_name, modality):
    if caller_name.startswith(f'{modality}_'):
        return caller_name  # Already prefixed
    return f'{modality}_{caller_name}'
```

### 3. Inconsistent Consensus Caller Naming
**Location**: Multiple functions
**Severity**: High
**Issue**: The `is_consensus_caller()` function checks multiple patterns:
```python
def is_consensus_caller(caller):
    return caller in ['dna_consensus', 'rna_consensus', 'DNA_consensus', 'RNA_consensus'] or caller.endswith('_consensus')
```
**Problem**: 
- Rescue script uses `DNA_consensus`/`RNA_consensus` but some checks look for `dna_consensus`/`rna_consensus`
- Inconsistent naming throughout codebase
- Causes logic errors in consensus detection

**Recommendation**: Standardize on one naming convention throughout

### 4. Rescue Flag Logic Inconsistency
**Location**: `mark_rescued_variants()` vs `write_union_vcf()`
**Severity**: High
**Issue**: Two different definitions of "rescued":
```python
# In mark_rescued_variants()
data['rescued'] = has_cross_modality  # Any cross-modality support

# In write_union_vcf() (INCORRECT)
record.info['RESCUED'] = 'YES' if 'consensus' in set(data['callers']) else 'NO'
```
**Problem**: 
- `mark_rescued_variants()` correctly identifies cross-modality variants
- `write_union_vcf()` uses incorrect logic checking for 'consensus' string
- Results in wrong RESCUED flag in output VCF

**Correct Implementation**:
```python
record.info['RESCUED'] = 'YES' if data.get('rescued', False) else 'NO'
```

## Significant Issues

### 5. Modality Map Key Collision Handling
**Location**: `main()` function
**Severity**: Medium
**Issue**: Potential key collisions in modality_map
```python
# If same caller appears in both DNA and RNA collections
modality_map[caller_name] = modality  # Later assignment overwrites earlier
```
**Problem**: 
- If same caller name appears in both DNA and RNA VCF lists, later assignment overwrites
- Could cause incorrect modality assignment
- Silent failure mode

**Recommendation**: Validate for duplicates and use prefixed names consistently

### 6. Template Header Source Assumption
**Location**: `main()` function
**Issue**: 
```python
# Always uses DNA consensus as template
template_vcf = VCF(str(dna_consensus_path))
```
**Problem**: 
- Assumes DNA consensus VCF has complete/correct header
- If DNA consensus is malformed, entire rescue fails
- No fallback to RNA consensus header

**Recommendation**: Add header validation and fallback logic

### 7. Sample Name Consistency
**Location**: `main()` function
**Issue**: No validation that sample names match across DNA and RNA VCFs
**Problem**: 
- May process VCFs from different samples
- Could lead to incorrect cross-modality analysis
- Silent failure mode

**Recommendation**: 
```python
# Validate sample consistency
dna_samples = set(VCF(dna_consensus_path).samples)
rna_samples = set(VCF(rna_consensus_path).samples)
if dna_samples != rna_samples:
    print(f"Warning: Sample mismatch - DNA: {dna_samples}, RNA: {rna_samples}")
```

### 8. Unified Filter DNA Priority Logic
**Location**: `write_union_vcf()` function
**Issue**: DNA priority rule may not be appropriate in all cases
```python
# If different labels, take DNA modality's consensus label
if dna_label and rna_label:
    if dna_label == rna_label:
        unified_classification = dna_label
    else:
        unified_classification = dna_label  # DNA takes priority
```
**Problem**: 
- Always prioritizes DNA over RNA may not be biologically appropriate
- No consideration of confidence scores or quality metrics
- Hardcoded priority without user control

**Recommendation**: Make priority configurable or use quality-based tie-breaking

## Moderate Issues

### 9. Individual Caller Exclusion Logic
**Location**: `write_union_vcf()` function
**Issue**: Complex logic to exclude consensus callers from various INFO fields
```python
# Multiple places with similar logic
if not is_consensus_caller(caller):  # Skip consensus callers
```
**Problem**: 
- Repeated logic scattered throughout function
- Easy to miss exclusions in new INFO fields
- Maintenance burden

**Recommendation**: Centralize caller filtering logic

### 10. Modality-Specific Statistics Calculation
**Location**: `write_union_vcf()` function
**Issue**: Complex nested loops for DNA/RNA statistics
```python
# Calculate DNA statistics
dna_dp_values = []
for i, caller in enumerate(data['callers']):
    if modality_map.get(caller) == 'DNA':
        if i < len(agg['dp_by_caller']) and agg['dp_by_caller'][i] is not None:
            dna_dp_values.append(agg['dp_by_caller'][i])
```
**Problem**: 
- Inefficient nested loops
- Index-based access prone to errors
- Difficult to maintain and extend

**Recommendation**: Pre-compute modality-specific data structures

### 11. Error Handling in VCF Reading
**Location**: `main()` function
**Issue**: Limited error handling for VCF file reading
```python
# No try-catch around VCF reading
dna_consensus = read_variants_from_vcf(str(dna_consensus_path), 'dna_consensus', modality='DNA')
```
**Problem**: 
- Malformed VCF files cause script crash
- No graceful degradation
- Poor error messages for users

**Recommendation**: Add comprehensive error handling with informative messages

### 12. Memory Usage with Large Individual Caller VCFs
**Location**: `main()` function
**Issue**: All individual caller VCFs loaded into memory simultaneously
**Problem**: 
- In full rescue mode with many individual callers, memory usage can be excessive
- No streaming or chunked processing
- May cause out-of-memory errors

**Recommendation**: Implement streaming processing or memory usage monitoring

## Minor Issues

### 13. Progress Reporting Inconsistency
**Location**: Various functions
**Issue**: Inconsistent progress reporting across different phases
**Problem**: 
- Some phases report progress, others don't
- Different reporting formats
- Difficult to track overall progress

### 14. Hardcoded String Constants
**Location**: Multiple functions
**Issue**: Hardcoded modality strings ('DNA', 'RNA') throughout code
**Problem**: 
- Difficult to extend to additional modalities
- Typo-prone
- No central configuration

**Recommendation**: Define constants at module level

### 15. Consensus-Only Mode Validation
**Location**: `main()` function
**Issue**: Limited validation in consensus-only mode
```python
if args.consensus_only and args.dna_vcf:
    print(f"Skipping {len(args.dna_vcf)} DNA caller VCFs (consensus_only mode)")
```
**Problem**: 
- Only prints warning, doesn't validate mode consistency
- User may expect individual callers to be processed
- Potential confusion about mode behavior

## Logic Issues

### 16. Cross-Modality Definition Inconsistency
**Location**: `mark_rescued_variants()` vs actual usage
**Issue**: Different interpretations of "cross-modality"
```python
# In mark_rescued_variants()
data['cross_modality'] = has_dna_support and has_rna_support

# But "support" is defined as consensus presence, not caller presence
```
**Problem**: 
- Cross-modality based on consensus presence, not actual caller support
- May miss variants with individual caller cross-modality support
- Inconsistent with rescue goal

### 17. Rescue Rate Calculation
**Location**: `compute_rescue_statistics()` function
**Issue**: 
```python
stats['rescue_rate'] = (stats['rescued'] / stats['total_variants']) * 100
```
**Problem**: 
- Rescue rate should be relative to cross-modality variants, not total variants
- Current calculation underestimates rescue effectiveness
- Misleading metric for evaluation

**Correct Calculation**:
```python
if stats['cross_modality'] > 0:
    stats['rescue_rate'] = (stats['rescued'] / stats['cross_modality']) * 100
else:
    stats['rescue_rate'] = 0.0
```

## Recommendations for Improvement

### Critical Priority
1. **Fix consensus caller detection** - Breaks core functionality
2. **Fix RESCUED flag logic** - Incorrect output metadata
3. **Standardize caller naming** - Prevents logic errors
4. **Fix modality prefix handling** - Prevents double-prefixing

### High Priority
5. **Add input validation** - Sample consistency, file format checks
6. **Improve error handling** - Graceful failure with informative messages
7. **Fix rescue rate calculation** - Correct statistical reporting

### Medium Priority
8. **Optimize memory usage** - Handle large datasets better
9. **Centralize caller filtering** - Reduce code duplication
10. **Add comprehensive logging** - Better debugging and monitoring

### Low Priority
11. **Performance optimization** - Reduce nested loops and redundant calculations
12. **Code documentation** - Add more detailed comments
13. **Configuration flexibility** - Make hardcoded values configurable