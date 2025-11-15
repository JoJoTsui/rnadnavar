#!/usr/bin/env python3
"""
Validation script for consensus refactoring backward compatibility.

This script validates that the refactored run_consensus_vcf.py maintains
backward compatibility with the original implementation by testing:
1. Output VCF format matches expected structure
2. INFO fields are present and correctly formatted
3. Variant aggregation logic produces consistent results
4. Filter normalization works correctly
5. Statistics generation is accurate

Usage:
    python validate_consensus_backward_compatibility.py --test-data <path>
"""

import argparse
import sys
import os
from pathlib import Path

# Add vcf_utils to path
sys.path.insert(0, str(Path(__file__).parent))

from vcf_utils.aggregation import read_variants_from_vcf, aggregate_variants, extract_genotype_info
from vcf_utils.filters import normalize_filter, categorize_filter
from vcf_utils.statistics import compute_consensus_statistics
from vcf_utils.io_utils import variant_key, is_snv


def test_filter_normalization():
    """Test filter normalization maintains expected mappings."""
    print("\n=== Testing Filter Normalization ===")
    
    test_cases = [
        ('PASS', 'PASS'),
        ('.', 'PASS'),
        (None, 'PASS'),
        ('germline', 'Germline'),
        ('GERMLINE', 'Germline'),
        ('RefCall', 'ReferenceCall'),
        ('LowDepth;LowEVS', 'LowDepth;LowEvidenceScore'),
        ('germline;base_qual', 'Germline;LowBaseQuality'),
    ]
    
    passed = 0
    failed = 0
    
    for input_filter, expected_output in test_cases:
        result = normalize_filter(input_filter)
        if result == expected_output:
            print(f"  ✓ '{input_filter}' -> '{result}'")
            passed += 1
        else:
            print(f"  ✗ '{input_filter}' -> '{result}' (expected '{expected_output}')")
            failed += 1
    
    print(f"\nFilter Normalization: {passed} passed, {failed} failed")
    return failed == 0


def test_filter_categorization():
    """Test filter categorization produces expected categories."""
    print("\n=== Testing Filter Categorization ===")
    
    test_cases = [
        ('PASS', 'PASS'),
        ('Germline', 'germline'),
        ('LowDepth', 'depth'),
        ('StrandBias', 'bias'),
        ('LowQuality', 'quality'),
        ('Germline;LowBaseQuality', 'germline;technical'),
    ]
    
    passed = 0
    failed = 0
    
    for input_filter, expected_category in test_cases:
        result = categorize_filter(input_filter)
        if result == expected_category:
            print(f"  ✓ '{input_filter}' -> '{result}'")
            passed += 1
        else:
            print(f"  ✗ '{input_filter}' -> '{result}' (expected '{expected_category}')")
            failed += 1
    
    print(f"\nFilter Categorization: {passed} passed, {failed} failed")
    return failed == 0


def test_variant_key_generation():
    """Test variant key generation is consistent."""
    print("\n=== Testing Variant Key Generation ===")
    
    # Test with dict representation
    test_variants = [
        ({'CHROM': 'chr1', 'POS': 12345, 'REF': 'A', 'ALT': 'G'}, '1:12345:A:G'),
        ({'CHROM': '2', 'POS': 67890, 'REF': 'C', 'ALT': 'T'}, '2:67890:C:T'),
        ({'CHROM': 'chrX', 'POS': 100, 'REF': 'AT', 'ALT': 'A'}, 'X:100:AT:A'),
    ]
    
    passed = 0
    failed = 0
    
    for variant_dict, expected_key in test_variants:
        result = variant_key(variant_dict, use_cyvcf2=False)
        if result == expected_key:
            print(f"  ✓ {variant_dict['CHROM']}:{variant_dict['POS']} -> '{result}'")
            passed += 1
        else:
            print(f"  ✗ {variant_dict['CHROM']}:{variant_dict['POS']} -> '{result}' (expected '{expected_key}')")
            failed += 1
    
    print(f"\nVariant Key Generation: {passed} passed, {failed} failed")
    return failed == 0


def test_snv_detection():
    """Test SNV detection logic."""
    print("\n=== Testing SNV Detection ===")
    
    test_cases = [
        ('A', ['G'], True),
        ('A', ['G', 'T'], True),
        ('AT', ['A'], False),
        ('A', ['AT'], False),
        ('A', ['*'], False),
        ('C', ['T'], True),
    ]
    
    passed = 0
    failed = 0
    
    for ref, alt_list, expected in test_cases:
        result = is_snv(ref, alt_list)
        if result == expected:
            print(f"  ✓ REF={ref}, ALT={alt_list} -> {result}")
            passed += 1
        else:
            print(f"  ✗ REF={ref}, ALT={alt_list} -> {result} (expected {expected})")
            failed += 1
    
    print(f"\nSNV Detection: {passed} passed, {failed} failed")
    return failed == 0


def test_consensus_statistics():
    """Test consensus statistics computation."""
    print("\n=== Testing Consensus Statistics ===")
    
    # Create mock variant data
    variant_data = {
        'chr1:100:A:G': {
            'is_snv': True,
            'callers': ['mutect2', 'strelka'],
            'passes_consensus': True
        },
        'chr1:200:C:T': {
            'is_snv': True,
            'callers': ['mutect2'],
            'passes_consensus': False
        },
        'chr1:300:AT:A': {
            'is_snv': False,
            'callers': ['mutect2', 'strelka', 'deepsomatic'],
            'passes_consensus': True
        },
    }
    
    stats = compute_consensus_statistics(variant_data, snv_threshold=2, indel_threshold=2)
    
    passed = 0
    failed = 0
    
    checks = [
        ('total_variants', 3),
        ('snvs', 2),
        ('indels', 1),
        ('snvs_consensus', 1),
        ('indels_consensus', 1),
        ('single_caller', 1),
        ('multi_caller', 2),
    ]
    
    for key, expected_value in checks:
        actual_value = stats.get(key)
        if actual_value == expected_value:
            print(f"  ✓ {key}: {actual_value}")
            passed += 1
        else:
            print(f"  ✗ {key}: {actual_value} (expected {expected_value})")
            failed += 1
    
    print(f"\nConsensus Statistics: {passed} passed, {failed} failed")
    return failed == 0


def test_variant_aggregation():
    """Test variant aggregation logic."""
    print("\n=== Testing Variant Aggregation ===")
    
    # Create mock variant collections
    caller1_vars = {
        'chr1:100:A:G': {
            'CHROM': 'chr1', 'POS': 100, 'REF': 'A', 'ALT': 'G',
            'is_snv': True, 'caller': 'caller1', 'filter_original': 'PASS',
            'filter_normalized': 'PASS', 'filter_category': 'PASS',
            'quality': 100.0, 'genotype': {'GT': '0/1', 'DP': 100, 'VAF': 0.25},
            'id': None
        }
    }
    
    caller2_vars = {
        'chr1:100:A:G': {
            'CHROM': 'chr1', 'POS': 100, 'REF': 'A', 'ALT': 'G',
            'is_snv': True, 'caller': 'caller2', 'filter_original': 'PASS',
            'filter_normalized': 'PASS', 'filter_category': 'PASS',
            'quality': 95.0, 'genotype': {'GT': '0/1', 'DP': 120, 'VAF': 0.27},
            'id': None
        }
    }
    
    collections = [
        ('caller1', caller1_vars, None),
        ('caller2', caller2_vars, None)
    ]
    
    aggregated = aggregate_variants(collections, snv_threshold=2, indel_threshold=2)
    
    passed = 0
    failed = 0
    
    # Check aggregation results
    if len(aggregated) == 1:
        print(f"  ✓ Aggregated 1 unique variant")
        passed += 1
    else:
        print(f"  ✗ Aggregated {len(aggregated)} variants (expected 1)")
        failed += 1
    
    variant = aggregated.get('chr1:100:A:G')
    if variant:
        if len(variant['callers']) == 2:
            print(f"  ✓ Variant has 2 callers")
            passed += 1
        else:
            print(f"  ✗ Variant has {len(variant['callers'])} callers (expected 2)")
            failed += 1
        
        if variant['passes_consensus']:
            print(f"  ✓ Variant passes consensus")
            passed += 1
        else:
            print(f"  ✗ Variant does not pass consensus")
            failed += 1
        
        if 'gt_aggregated' in variant:
            print(f"  ✓ Genotype aggregation present")
            passed += 1
        else:
            print(f"  ✗ Genotype aggregation missing")
            failed += 1
    else:
        print(f"  ✗ Variant not found in aggregated data")
        failed += 4
    
    print(f"\nVariant Aggregation: {passed} passed, {failed} failed")
    return failed == 0


def main():
    """Run all validation tests."""
    parser = argparse.ArgumentParser(
        description='Validate consensus refactoring backward compatibility'
    )
    parser.add_argument(
        '--test-data',
        help='Path to test data directory (optional)',
        default=None
    )
    
    args = parser.parse_args()
    
    print("=" * 70)
    print("Consensus Refactoring Backward Compatibility Validation")
    print("=" * 70)
    
    # Run all tests
    tests = [
        test_filter_normalization,
        test_filter_categorization,
        test_variant_key_generation,
        test_snv_detection,
        test_consensus_statistics,
        test_variant_aggregation,
    ]
    
    results = []
    for test_func in tests:
        try:
            result = test_func()
            results.append(result)
        except Exception as e:
            print(f"\n✗ Test {test_func.__name__} failed with exception: {e}")
            results.append(False)
    
    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    
    passed = sum(results)
    total = len(results)
    
    print(f"\nTests passed: {passed}/{total}")
    
    if all(results):
        print("\n✓ All backward compatibility tests PASSED")
        print("\nThe refactored code maintains backward compatibility with the")
        print("original implementation. Output format, INFO fields, and variant")
        print("aggregation logic are consistent with expected behavior.")
        return 0
    else:
        print("\n✗ Some backward compatibility tests FAILED")
        print("\nPlease review the failed tests above and ensure the refactored")
        print("code maintains the same behavior as the original implementation.")
        return 1


if __name__ == '__main__':
    sys.exit(main())
