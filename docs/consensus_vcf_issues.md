# Consensus VCF Script Issues and Bug Analysis

## Critical Issues

### 1. Incorrect RESCUED Flag Logic
**Location**: `write_union_vcf()` function in `io_utils.py`
**Severity**: High
**Issue**: 
```python
# Current implementation (INCORRECT)
record.info['RESCUED'] = 'YES' if 'consensus' in set(data['callers']) else 'NO'
```
**Problem**: 
- Checks if string 'consensus' appears anywhere in caller names
- In consensus mode, this should always be 'NO' since no rescue operation occurs
- Logic is fundamentally flawed for consensus workflow

**Correct Implementation**:
```python
# For consensus mode, RESCUED should always be 'NO'
record.info['RESCUED'] = 'NO'
```

### 2. Genotype Extraction Error Handling
**Location**: `extract_genotype_info()` function in `aggregation.py`
**Severity**: Medium
**Issue**: Overly broad exception handling masks specific errors
```python
try:
    # Complex Strelka parsing logic
    # ... many operations that could fail
except Exception as e:
    print(f"Warning: Error extracting genotype info from {caller}: {e}")
```
**Problems**:
- Catches all exceptions, including programming errors
- May hide legitimate bugs in parsing logic
- Continues processing with potentially invalid data

**Recommendation**:
```python
try:
    # Specific operations
except (KeyError, IndexError, ValueError, TypeError) as e:
    print(f"Warning: Data format error in {caller}: {e}")
except Exception as e:
    print(f"Error: Unexpected error in {caller}: {e}")
    raise  # Re-raise unexpected errors
```

### 3. Filter Field Semicolon Replacement
**Location**: `write_union_vcf()` function
**Severity**: Medium
**Issue**: 
```python
# Replace semicolons with commas in filter values
filter_val = str(data['filters_original'][i]).replace(';', ',')
```
**Problems**:
- Semicolons are standard VCF delimiters for compound filters
- Replacement may break downstream VCF parsing tools
- Creates non-standard VCF format

**Recommendation**: Use proper VCF escaping or keep original format

### 4. Chromosome Sorting Edge Cases
**Location**: `write_union_vcf()` function
**Severity**: Low
**Issue**: 
```python
try:
    if chrom.isdigit():
        chrom_idx = int(chrom)
    elif chrom == 'X':
        chrom_idx = 23
    # ...
except:
    pass  # Too broad exception handling
```
**Problems**:
- Overly broad exception handling
- Doesn't handle non-standard chromosome names well
- May silently fail on legitimate chromosome names

### 5. Missing Input Validation
**Location**: `main()` function
**Severity**: Medium
**Issues**:
- No validation that consensus thresholds are reasonable (> 0, <= total callers)
- No check for minimum number of VCF files
- No validation of VCF file format/readability before processing
- No sample name consistency checks across VCFs

**Recommendations**:
```python
# Add validation
if args.snv_thr <= 0 or args.indel_thr <= 0:
    print("Error: Consensus thresholds must be > 0")
    sys.exit(1)

if len(vcf_files) < args.snv_thr and len(vcf_files) < args.indel_thr:
    print(f"Warning: Only {len(vcf_files)} VCF files found, but thresholds require {max(args.snv_thr, args.indel_thr)}")
```

## Moderate Issues

### 6. Memory Usage for Large VCFs
**Location**: `aggregate_variants()` function
**Severity**: Medium
**Issue**: All variants loaded into memory simultaneously
```python
aggregated = defaultdict(lambda: {...})  # Stores all variants in memory
```
**Problem**: For very large VCFs (millions of variants), could cause memory issues
**Recommendation**: Consider streaming processing or chunked processing for large datasets

### 7. Inconsistent Caller Name Handling
**Location**: `get_caller_name()` function
**Severity**: Low
**Issue**: Multiple fallback strategies may produce inconsistent results
```python
# Fallback: take second component
parts = name.split('.')
return parts[1] if len(parts) > 1 else parts[0]
```
**Problem**: May not work reliably for all filename patterns
**Recommendation**: Define strict naming conventions and validate against them

### 8. Genotype Consensus Tie-Breaking
**Location**: `aggregate_genotypes()` function
**Severity**: Low
**Issue**: No explicit tie-breaking for genotype consensus
```python
consensus_gt = max(gt_counts, key=gt_counts.get) if gt_counts else None
```
**Problem**: If multiple genotypes have same count, result is arbitrary
**Recommendation**: Implement priority-based tie-breaking (e.g., 1/1 > 0/1 > 0/0)

### 9. Quality Score Handling
**Location**: Multiple functions
**Severity**: Low
**Issue**: No validation of quality score ranges or formats
**Problem**: Invalid quality scores may cause downstream issues
**Recommendation**: Add quality score validation and normalization

## Minor Issues

### 10. Progress Reporting Granularity
**Location**: `write_union_vcf()` function
**Severity**: Low
**Issue**: Progress reported every 10,000 variants may be too frequent for small datasets
**Recommendation**: Adaptive progress reporting based on dataset size

### 11. Hardcoded Constants
**Location**: Multiple functions
**Severity**: Low
**Issue**: Magic numbers and hardcoded values scattered throughout code
**Examples**:
- Strelka normal depth threshold: `>= 2`
- Progress reporting interval: `10000`
**Recommendation**: Define constants at module level

### 12. Error Message Clarity
**Location**: Multiple functions
**Severity**: Low
**Issue**: Some error messages lack context about which file or variant caused the issue
**Recommendation**: Include more context in error messages for debugging

## Potential Logic Issues

### 13. SNV Detection Logic
**Location**: `is_snv()` function
**Severity**: Low
**Issue**: May not handle all edge cases correctly
```python
def is_snv(ref, alt_list):
    if len(ref) != 1:
        return False
    for alt in alt_list:
        if len(alt) != 1 or alt == '*':  # '*' indicates deletion
            return False
    return True
```
**Potential Issue**: The '*' check may not be comprehensive for all VCF formats

### 14. Variant Key Collision Risk
**Location**: `variant_key()` function
**Severity**: Low
**Issue**: Potential for key collisions with complex variants
```python
return f"{chrom}:{pos}:{ref}:{alts}"
```
**Problem**: May not uniquely identify all variant types (e.g., complex structural variants)

## Recommendations for Improvement

### High Priority
1. **Fix RESCUED flag logic** - Critical for correct output
2. **Improve error handling** - Use specific exception types
3. **Add input validation** - Prevent runtime errors
4. **Fix filter field handling** - Maintain VCF standard compliance

### Medium Priority
5. **Memory optimization** - Handle large datasets better
6. **Enhance logging** - Better debugging and monitoring
7. **Standardize caller naming** - Reduce inconsistencies

### Low Priority
8. **Code documentation** - Add more inline comments
9. **Unit testing** - Add comprehensive test coverage
10. **Performance profiling** - Identify bottlenecks for optimization