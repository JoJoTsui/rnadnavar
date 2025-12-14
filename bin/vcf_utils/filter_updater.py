#!/usr/bin/env python3
"""
FILTER Column Update Logic for RNA Editing Annotation

This module implements the FILTER column update logic that processes annotated variants
and updates FILTER values based on evidence tiering results. It loops through all
annotated variants for evidence evaluation, updates FILTER to "RNAedit" for variants
meeting high evidence criteria, and preserves original FILTER values for variants
not meeting criteria.

Requirements Satisfied: 5.1, 5.2, 5.3, 5.4

Author: RNA Editing Enhancement Pipeline
Date: 2025-12-14
"""

import logging
from typing import Dict, Any, Optional, List, Tuple
from pathlib import Path

logger = logging.getLogger(__name__)

class FilterUpdater:
    """
    FILTER column updater for RNA editing variants.
    
    This class implements the FILTER column update logic that:
    1. Loops through all annotated variants for evidence evaluation
    2. Updates FILTER to "RNAedit" for variants meeting high evidence criteria
    3. Preserves original FILTER values for variants not meeting criteria
    4. Ensures complete processing of all variants in input VCF
    """
    
    def __init__(self):
        """Initialize FILTER updater."""
        self.stats = {
            'total_variants_processed': 0,
            'filter_updates_to_rnaedit': 0,
            'filter_values_preserved': 0,
            'original_filter_distribution': {},
            'final_filter_distribution': {}
        }
        
        logger.info("FILTER updater initialized")
    
    def should_update_filter_to_rnaedit(self, evidence_tier: str, 
                                       has_rediportal_match: bool,
                                       rna_support: int, min_rna_support: int) -> bool:
        """
        Determine if FILTER should be updated to "RNAedit" based on evidence criteria.
        
        FILTER Update Rules (Requirements 5.1, 5.2):
        - Update to "RNAedit" if:
          - N_RNA_CALLERS_SUPPORT >= min_rna_support (RNA consensus)
          - Has exact REDIportal match
          - Evidence tier is HIGH, MEDIUM, or LOW
        
        Args:
            evidence_tier: Evidence tier (HIGH, MEDIUM, LOW, NONE)
            has_rediportal_match: Whether variant has exact REDIportal match
            rna_support: Number of RNA callers supporting the variant
            min_rna_support: Minimum RNA support threshold
            
        Returns:
            True if FILTER should be updated to "RNAedit", False to preserve original
        """
        # Check RNA consensus requirement
        if rna_support < min_rna_support:
            return False
        
        # Check REDIportal match requirement
        if not has_rediportal_match:
            return False
        
        # Check evidence tier requirement
        if evidence_tier not in ['HIGH', 'MEDIUM', 'LOW']:
            return False
        
        return True
    
    def update_variant_filter(self, original_filter: str, evidence_tier: str,
                            has_rediportal_match: bool, rna_support: int,
                            min_rna_support: int) -> Tuple[str, bool]:
        """
        Update variant FILTER based on evidence criteria.
        
        Args:
            original_filter: Original FILTER value from variant
            evidence_tier: Evidence tier (HIGH, MEDIUM, LOW, NONE)
            has_rediportal_match: Whether variant has exact REDIportal match
            rna_support: Number of RNA callers supporting the variant
            min_rna_support: Minimum RNA support threshold
            
        Returns:
            Tuple of (new_filter_value, was_updated)
        """
        self.stats['total_variants_processed'] += 1
        
        # Track original filter distribution
        if original_filter not in self.stats['original_filter_distribution']:
            self.stats['original_filter_distribution'][original_filter] = 0
        self.stats['original_filter_distribution'][original_filter] += 1
        
        # Determine if FILTER should be updated
        should_update = self.should_update_filter_to_rnaedit(
            evidence_tier, has_rediportal_match, rna_support, min_rna_support
        )
        
        if should_update:
            new_filter = "RNAedit"
            self.stats['filter_updates_to_rnaedit'] += 1
            was_updated = True
        else:
            new_filter = original_filter
            self.stats['filter_values_preserved'] += 1
            was_updated = False
        
        # Track final filter distribution
        if new_filter not in self.stats['final_filter_distribution']:
            self.stats['final_filter_distribution'][new_filter] = 0
        self.stats['final_filter_distribution'][new_filter] += 1
        
        # Log individual decisions for debugging
        self.log_filter_decision(original_filter, new_filter, evidence_tier, 
                               has_rediportal_match, rna_support, was_updated)
        
        # Log progress at regular intervals
        self.log_processing_progress(self.stats['total_variants_processed'])
        
        return new_filter, was_updated
    
    def process_all_variants(self, variants_data: List[Dict[str, Any]], 
                           min_rna_support: int) -> List[Dict[str, Any]]:
        """
        Process all variants for FILTER updates (Requirement 5.4).
        
        This method ensures complete processing of all variants in input VCF
        by looping through all annotated variants for evidence evaluation.
        
        Args:
            variants_data: List of variant dictionaries with evidence data
            min_rna_support: Minimum RNA support threshold
            
        Returns:
            List of variant dictionaries with updated FILTER information
        """
        logger.info(f"Processing {len(variants_data)} variants for FILTER updates...")
        
        updated_variants = []
        
        for variant_data in variants_data:
            # Extract required fields
            original_filter = variant_data.get('original_filter', 'PASS')
            evidence_tier = variant_data.get('evidence_tier', 'NONE')
            has_rediportal_match = variant_data.get('has_rediportal_match', False)
            rna_support = variant_data.get('rna_support', 0)
            
            # Update FILTER
            new_filter, was_updated = self.update_variant_filter(
                original_filter, evidence_tier, has_rediportal_match, 
                rna_support, min_rna_support
            )
            
            # Create updated variant data
            updated_variant = variant_data.copy()
            updated_variant['new_filter'] = new_filter
            updated_variant['filter_was_updated'] = was_updated
            
            updated_variants.append(updated_variant)
        
        logger.info(f"✓ FILTER processing completed: {self.stats['filter_updates_to_rnaedit']} updates, "
                   f"{self.stats['filter_values_preserved']} preserved")
        
        return updated_variants
    
    def get_statistics(self) -> Dict[str, Any]:
        """
        Get FILTER update statistics.
        
        Returns:
            Dictionary with FILTER update statistics
        """
        stats = self.stats.copy()
        
        # Calculate percentages
        total = stats['total_variants_processed']
        if total > 0:
            stats['filter_update_rate'] = (stats['filter_updates_to_rnaedit'] / total) * 100
            stats['filter_preservation_rate'] = (stats['filter_values_preserved'] / total) * 100
        else:
            stats['filter_update_rate'] = 0.0
            stats['filter_preservation_rate'] = 0.0
        
        return stats
    
    def log_statistics(self) -> None:
        """Log comprehensive FILTER update statistics."""
        stats = self.get_statistics()
        
        logger.info("=== FILTER Update Statistics ===")
        logger.info(f"Total variants processed: {stats['total_variants_processed']:,}")
        logger.info(f"FILTER updates to RNAedit: {stats['filter_updates_to_rnaedit']:,} ({stats['filter_update_rate']:.1f}%)")
        logger.info(f"FILTER values preserved: {stats['filter_values_preserved']:,} ({stats['filter_preservation_rate']:.1f}%)")
        
        logger.info("Original FILTER distribution:")
        for filter_val, count in sorted(stats['original_filter_distribution'].items(), key=lambda x: x[1], reverse=True):
            percentage = (count / stats['total_variants_processed'] * 100) if stats['total_variants_processed'] > 0 else 0
            logger.info(f"  {filter_val}: {count:,} ({percentage:.1f}%)")
        
        logger.info("Final FILTER distribution:")
        for filter_val, count in sorted(stats['final_filter_distribution'].items(), key=lambda x: x[1], reverse=True):
            percentage = (count / stats['total_variants_processed'] * 100) if stats['total_variants_processed'] > 0 else 0
            logger.info(f"  {filter_val}: {count:,} ({percentage:.1f}%)")
        
        # Log FILTER transition analysis
        if stats['filter_updates_to_rnaedit'] > 0:
            logger.info("FILTER transitions to RNAedit:")
            transitions = {}
            
            # Calculate transitions (this is a simplified version - in practice we'd track this)
            for original_filter in stats['original_filter_distribution']:
                if original_filter != 'RNAedit':  # Only count transitions TO RNAedit
                    # Estimate transitions based on the difference
                    original_count = stats['original_filter_distribution'][original_filter]
                    final_count = stats['final_filter_distribution'].get(original_filter, 0)
                    transition_count = max(0, original_count - final_count)
                    if transition_count > 0:
                        transitions[original_filter] = transition_count
            
            for original_filter, count in sorted(transitions.items(), key=lambda x: x[1], reverse=True):
                percentage = (count / stats['filter_updates_to_rnaedit'] * 100) if stats['filter_updates_to_rnaedit'] > 0 else 0
                logger.info(f"  {original_filter} → RNAedit: {count:,} ({percentage:.1f}%)")
        
        # Log efficiency metrics
        if stats['total_variants_processed'] > 0:
            logger.debug("=== FILTER Update Decision Details ===")
            logger.debug("Update criteria: RNA consensus + REDIportal match + evidence tier HIGH/MEDIUM/LOW")
            logger.debug("Preserve criteria: No RNA consensus OR no REDIportal match OR evidence tier NONE")
    
    def log_filter_decision(self, original_filter: str, new_filter: str, evidence_tier: str, 
                           has_rediportal_match: bool, rna_support: int, was_updated: bool):
        """Log individual FILTER update decisions for debugging."""
        if logger.isEnabledFor(logging.DEBUG) and was_updated:
            logger.debug(f"FILTER update: {original_filter} → {new_filter}")
            logger.debug(f"  Evidence tier: {evidence_tier}, RNA support: {rna_support}")
            logger.debug(f"  REDIportal match: {has_rediportal_match}")
    
    def log_processing_progress(self, processed_count: int, interval: int = 10000):
        """Log FILTER processing progress at regular intervals."""
        if processed_count % interval == 0:
            logger.info(f"FILTER update progress: {processed_count:,} variants processed")
            
            # Log current update rate
            current_stats = self.get_statistics()
            if current_stats['total_variants_processed'] > 0:
                logger.info(f"  Current update rate: {current_stats['filter_update_rate']:.1f}%")


class FilterUpdateValidator:
    """
    Validator for FILTER update results.
    
    This class provides validation methods to ensure FILTER update
    results are consistent and follow the expected rules.
    """
    
    @staticmethod
    def validate_filter_update_logic(original_filter: str, new_filter: str,
                                   evidence_tier: str, has_rediportal_match: bool,
                                   rna_support: int, min_rna_support: int) -> bool:
        """
        Validate that FILTER update follows the expected logic.
        
        Args:
            original_filter: Original FILTER value
            new_filter: New FILTER value after update
            evidence_tier: Evidence tier
            has_rediportal_match: Whether variant has REDIportal match
            rna_support: RNA caller support
            min_rna_support: Minimum RNA support threshold
            
        Returns:
            True if logic is consistent, False otherwise
        """
        # Check if update criteria are met
        should_update = (
            rna_support >= min_rna_support and
            has_rediportal_match and
            evidence_tier in ['HIGH', 'MEDIUM', 'LOW']
        )
        
        if should_update:
            return new_filter == "RNAedit"
        else:
            return new_filter == original_filter
    
    @staticmethod
    def validate_complete_processing(input_variant_count: int, 
                                   processed_variant_count: int) -> bool:
        """
        Validate that all variants were processed (Requirement 5.4).
        
        Args:
            input_variant_count: Number of input variants
            processed_variant_count: Number of processed variants
            
        Returns:
            True if all variants were processed, False otherwise
        """
        return input_variant_count == processed_variant_count
    
    @staticmethod
    def validate_filter_values(filter_value: str) -> bool:
        """
        Validate that FILTER value is valid.
        
        Args:
            filter_value: FILTER value to validate
            
        Returns:
            True if valid, False otherwise
        """
        # Common valid FILTER values
        valid_filters = {
            'PASS', 'RNAedit', 'FAIL', 'LowQual', 'weak_evidence', 
            'clustered_events', 'multiallelic', 'normal_artifact',
            'orientation', 'slippage', 'strand_bias', 'base_qual',
            'map_qual', 'fragment', 'contamination'
        }
        
        # Allow any non-empty string as FILTER value (VCF spec allows custom filters)
        return isinstance(filter_value, str) and len(filter_value.strip()) > 0


def create_filter_updater() -> FilterUpdater:
    """
    Factory function to create FILTER updater.
    
    Returns:
        Configured FilterUpdater instance
    """
    return FilterUpdater()


if __name__ == '__main__':
    # Simple testing functionality
    import json
    
    # Test FILTER updater
    updater = create_filter_updater()
    
    # Test cases
    test_cases = [
        {
            'name': 'HIGH evidence - should update to RNAedit',
            'original_filter': 'PASS',
            'evidence_tier': 'HIGH',
            'has_rediportal_match': True,
            'rna_support': 3,
            'min_rna_support': 2,
            'expected_filter': 'RNAedit',
            'expected_updated': True
        },
        {
            'name': 'MEDIUM evidence - should update to RNAedit',
            'original_filter': 'PASS',
            'evidence_tier': 'MEDIUM',
            'has_rediportal_match': True,
            'rna_support': 2,
            'min_rna_support': 2,
            'expected_filter': 'RNAedit',
            'expected_updated': True
        },
        {
            'name': 'NONE evidence - should preserve original',
            'original_filter': 'PASS',
            'evidence_tier': 'NONE',
            'has_rediportal_match': False,
            'rna_support': 1,
            'min_rna_support': 2,
            'expected_filter': 'PASS',
            'expected_updated': False
        },
        {
            'name': 'No REDIportal match - should preserve original',
            'original_filter': 'LowQual',
            'evidence_tier': 'HIGH',
            'has_rediportal_match': False,
            'rna_support': 3,
            'min_rna_support': 2,
            'expected_filter': 'LowQual',
            'expected_updated': False
        },
        {
            'name': 'Insufficient RNA support - should preserve original',
            'original_filter': 'PASS',
            'evidence_tier': 'HIGH',
            'has_rediportal_match': True,
            'rna_support': 1,
            'min_rna_support': 2,
            'expected_filter': 'PASS',
            'expected_updated': False
        }
    ]
    
    print("Testing FILTER Updater:")
    print("=" * 50)
    
    for test_case in test_cases:
        new_filter, was_updated = updater.update_variant_filter(
            test_case['original_filter'],
            test_case['evidence_tier'],
            test_case['has_rediportal_match'],
            test_case['rna_support'],
            test_case['min_rna_support']
        )
        
        filter_correct = new_filter == test_case['expected_filter']
        update_correct = was_updated == test_case['expected_updated']
        status = "✓" if (filter_correct and update_correct) else "✗"
        
        print(f"{status} {test_case['name']}")
        print(f"   Original: {test_case['original_filter']} → New: {new_filter} (expected: {test_case['expected_filter']})")
        print(f"   Updated: {was_updated} (expected: {test_case['expected_updated']})")
        print(f"   RNA support: {test_case['rna_support']}, REDIportal match: {test_case['has_rediportal_match']}")
        print(f"   Evidence tier: {test_case['evidence_tier']}")
        print()
    
    # Print statistics
    print("FILTER Update Statistics:")
    updater.log_statistics()
    
    # Test validator
    print("\nTesting FILTER Update Validator:")
    print("=" * 50)
    
    validator = FilterUpdateValidator()
    
    # Test validation functions
    validation_tests = [
        ('Valid update logic (should update)', 
         validator.validate_filter_update_logic('PASS', 'RNAedit', 'HIGH', True, 3, 2), True),
        ('Valid preserve logic (no REDIportal match)', 
         validator.validate_filter_update_logic('PASS', 'PASS', 'HIGH', False, 3, 2), True),
        ('Invalid update logic (should preserve but updated)', 
         validator.validate_filter_update_logic('PASS', 'RNAedit', 'NONE', False, 1, 2), False),
        ('Valid complete processing', 
         validator.validate_complete_processing(100, 100), True),
        ('Invalid complete processing', 
         validator.validate_complete_processing(100, 95), False),
        ('Valid FILTER value', 
         validator.validate_filter_values('RNAedit'), True),
        ('Invalid FILTER value', 
         validator.validate_filter_values(''), False)
    ]
    
    for test_name, result, expected in validation_tests:
        status = "✓" if result == expected else "✗"
        print(f"{status} {test_name}: {result} (expected: {expected})")