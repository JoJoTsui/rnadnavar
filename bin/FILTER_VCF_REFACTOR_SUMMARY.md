# filter_vcf.py Refactoring Summary

## Overview
Updated `filter_vcf.py` to use shared functions from the `vcf_utils` module, improving code consistency and maintainability across the consensus/rescue workflow.

## Changes Made

### 1. Added vcf_utils Imports
- Imported `normalize_filter` from `vcf_utils.filters`
- Imported `normalize_chromosome` and `variant_key` from `vcf_utils.io_utils`
- Added proper path setup to ensure vcf_utils is accessible

### 2. Updated apply_filters() Function

#### Whitelist Processing
- **Before**: Used simple string concatenation for variant IDs
- **After**: Uses `normalize_chromosome()` for consistent chromosome naming
- **Benefit**: Ensures whitelist matching works correctly regardless of chromosome naming convention (chr1 vs 1)

#### Blacklist Processing
- **Before**: Used raw chromosome names from blacklist file
- **After**: Uses `normalize_chromosome()` for consistent chromosome naming
- **Benefit**: Ensures blacklist matching works correctly across different naming conventions

#### Variant Key Generation
- **Before**: Used manual string concatenation: `f"{chrom}:{pos}:{ref}:{alt}"`
- **After**: Uses shared `variant_key()` function with `use_cyvcf2=True`
- **Benefit**: Consistent variant identification across all scripts in the workflow

#### Filter Normalization
- **Before**: Direct comparison of variant.FILTER with allowed filters
- **After**: Uses `normalize_filter()` to standardize filter values before comparison
- **Benefit**: Handles different filter naming conventions from various callers (e.g., GERMLINE vs germline)

## Maintained Functionality
All existing filtering logic remains unchanged:
- Minimum alt reads filtering
- gnomAD frequency filtering
- Blacklist region filtering
- Noncoding region filtering
- IG/pseudogene filtering
- Homopolymer context filtering
- Variant caller filter handling

## Testing
Created and ran validation tests to verify:
- ✓ vcf_utils imports work correctly
- ✓ normalize_filter() handles various filter formats
- ✓ normalize_chromosome() normalizes chromosome names
- ✓ variant_key() generates consistent variant identifiers

## Benefits
1. **Consistency**: Uses same normalization logic as consensus and rescue workflows
2. **Maintainability**: Filter mapping updates in vcf_utils.filters automatically apply here
3. **Reliability**: Reduces risk of inconsistent variant matching due to naming differences
4. **Code Reuse**: Eliminates duplicate logic across multiple scripts

## Requirements Satisfied
- Requirement 4.4: "THE System SHALL provide shared utility functions for filter normalization"
- Imports and reuses shared functions from vcf_utils
- Maintains existing functionality while improving code quality
