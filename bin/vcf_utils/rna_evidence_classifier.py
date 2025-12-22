#!/usr/bin/env python3
"""
RNA Evidence Classifier for RNA Editing Annotation

** NOT THE SAME AS HYBRID TIERING (CxDy SYSTEM) **
This module is specifically for RNA editing evidence classification (HIGH/MEDIUM/LOW/NONE).
It is separate from the hybrid tiering engine (C1-C7 caller tiers + D0-D1 database evidence).

This module implements the evidence tiering system that classifies RNA editing
confidence based on RNA/DNA caller support and updates FILTER columns accordingly.
It extracts N_RNA_CALLERS_SUPPORT and N_DNA_CALLERS_SUPPORT from variant INFO fields,
implements RNA consensus detection, RNA-only variant detection, and creates evidence
tier assignment logic (HIGH/MEDIUM/LOW/NONE).

For hybrid tiering (CxDy), see: bin/common/tiering_engine.py

Requirements Satisfied: 4.1, 4.2, 4.3, 4.4

Author: RNA Editing Enhancement Pipeline
Date: 2025-12-13
"""

import logging
from typing import Dict, Any, Optional, List, Tuple
from pathlib import Path

logger = logging.getLogger(__name__)

class EvidenceTieringProcessor:
    """
    Evidence tiering processor for RNA editing variants.
    
    This class implements the evidence tiering system that:
    1. Extracts N_RNA_CALLERS_SUPPORT and N_DNA_CALLERS_SUPPORT from variant INFO fields
    2. Implements RNA consensus detection (N_RNA_CALLERS_SUPPORT >= 2)
    3. Implements RNA-only variant detection (N_DNA_CALLERS_SUPPORT = 0)
    4. Creates evidence tier assignment logic (HIGH/MEDIUM/LOW/NONE)
    """
    
    def __init__(self, min_rna_support: int = 2):
        """
        Initialize evidence tiering processor.
        
        Args:
            min_rna_support: Minimum RNA caller support threshold for consensus
        """
        self.min_rna_support = min_rna_support
        self.stats = {
            'total_variants_processed': 0,
            'rna_consensus_variants': 0,
            'rna_only_variants': 0,
            'evidence_levels': {'HIGH': 0, 'MEDIUM': 0, 'LOW': 0, 'NONE': 0},
            'filter_updates': 0
        }
        
        logger.info(f"Evidence tiering processor initialized with min_rna_support={min_rna_support}")
    
    def extract_caller_support(self, variant_info: Dict[str, Any]) -> Tuple[int, int]:
        """
        Extract N_RNA_CALLERS_SUPPORT and N_DNA_CALLERS_SUPPORT from variant INFO fields.
        
        Args:
            variant_info: Dictionary containing variant INFO fields
            
        Returns:
            Tuple of (rna_support, dna_support)
        """
        try:
            rna_support = int(variant_info.get('N_RNA_CALLERS_SUPPORT', 0))
            dna_support = int(variant_info.get('N_DNA_CALLERS_SUPPORT', 0))
            
            # Validate values are non-negative
            if rna_support < 0:
                logger.warning(f"Invalid negative RNA support value: {rna_support}, using 0")
                rna_support = 0
            if dna_support < 0:
                logger.warning(f"Invalid negative DNA support value: {dna_support}, using 0")
                dna_support = 0
                
            return rna_support, dna_support
            
        except (ValueError, TypeError) as e:
            logger.warning(f"Error extracting caller support values: {e}")
            return 0, 0
    
    def detect_rna_consensus(self, rna_support: int) -> bool:
        """
        Implement RNA consensus detection (N_RNA_CALLERS_SUPPORT >= min_rna_support).
        
        Args:
            rna_support: Number of RNA callers supporting the variant
            
        Returns:
            True if RNA consensus is detected, False otherwise
        """
        return rna_support >= self.min_rna_support
    
    def detect_rna_only_variant(self, dna_support: int) -> bool:
        """
        Implement RNA-only variant detection (N_DNA_CALLERS_SUPPORT = 0).
        
        Args:
            dna_support: Number of DNA callers supporting the variant
            
        Returns:
            True if variant is RNA-only, False otherwise
        """
        return dna_support == 0
    
    def assign_evidence_tier(self, rna_support: int, dna_support: int, 
                           has_rediportal_match: bool) -> str:
        """
        Create evidence tier assignment logic (HIGH/MEDIUM/LOW/NONE).
        
        Evidence Tier Assignment Rules:
        - HIGH: RNA consensus (>=2) + RNA-only (DNA=0) + REDIportal match
        - MEDIUM: RNA consensus (>=2) + has DNA support + REDIportal match  
        - LOW: RNA consensus (>=2) + REDIportal match (any DNA support level)
        - NONE: No RNA consensus OR no REDIportal match
        
        Args:
            rna_support: Number of RNA callers supporting the variant
            dna_support: Number of DNA callers supporting the variant
            has_rediportal_match: Whether variant has exact REDIportal match
            
        Returns:
            Evidence tier string (HIGH, MEDIUM, LOW, NONE)
        """
        # Prerequisites for any evidence tier
        if not has_rediportal_match:
            return 'NONE'
        
        if not self.detect_rna_consensus(rna_support):
            return 'NONE'
        
        # RNA consensus is present and REDIportal match exists
        rna_only = self.detect_rna_only_variant(dna_support)
        
        if rna_only:
            return 'HIGH'  # RNA-only with consensus and REDIportal match
        else:
            return 'MEDIUM'  # Has DNA support but still has RNA consensus and REDIportal match
    
    def should_update_filter(self, evidence_tier: str) -> bool:
        """
        Determine if FILTER column should be updated based on evidence tier.
        
        Args:
            evidence_tier: Evidence tier (HIGH, MEDIUM, LOW, NONE)
            
        Returns:
            True if FILTER should be updated to "RNAedit", False to preserve original
        """
        return evidence_tier in ['HIGH', 'MEDIUM', 'LOW']
    
    def process_variant(self, variant_info: Dict[str, Any], 
                       has_rediportal_match: bool) -> Dict[str, Any]:
        """
        Process a single variant for evidence tiering.
        
        Args:
            variant_info: Dictionary containing variant INFO fields
            has_rediportal_match: Whether variant has exact REDIportal match
            
        Returns:
            Dictionary with evidence tiering results
        """
        self.stats['total_variants_processed'] += 1
        
        # Extract caller support
        rna_support, dna_support = self.extract_caller_support(variant_info)
        
        # Detect RNA consensus
        rna_consensus = self.detect_rna_consensus(rna_support)
        if rna_consensus:
            self.stats['rna_consensus_variants'] += 1
        
        # Detect RNA-only variants
        rna_only = self.detect_rna_only_variant(dna_support)
        if rna_only:
            self.stats['rna_only_variants'] += 1
        
        # Assign evidence tier
        evidence_tier = self.assign_evidence_tier(rna_support, dna_support, has_rediportal_match)
        self.stats['evidence_levels'][evidence_tier] += 1
        
        # Determine if FILTER should be updated
        update_filter = self.should_update_filter(evidence_tier)
        if update_filter:
            self.stats['filter_updates'] += 1
        
        result = {
            'rna_support': rna_support,
            'dna_support': dna_support,
            'rna_consensus': rna_consensus,
            'rna_only': rna_only,
            'evidence_tier': evidence_tier,
            'update_filter': update_filter,
            'new_filter': 'RNAedit' if update_filter else None
        }
        
        # Log individual classification decisions for debugging
        self.log_classification_decision(variant_info, result)
        
        # Log progress at regular intervals
        self.log_processing_progress(self.stats['total_variants_processed'])
        
        return result
    
    def get_statistics(self) -> Dict[str, Any]:
        """
        Get processing statistics.
        
        Returns:
            Dictionary with processing statistics
        """
        stats = self.stats.copy()
        
        # Calculate percentages
        total = stats['total_variants_processed']
        if total > 0:
            stats['rna_consensus_rate'] = (stats['rna_consensus_variants'] / total) * 100
            stats['rna_only_rate'] = (stats['rna_only_variants'] / total) * 100
            stats['filter_update_rate'] = (stats['filter_updates'] / total) * 100
        else:
            stats['rna_consensus_rate'] = 0.0
            stats['rna_only_rate'] = 0.0
            stats['filter_update_rate'] = 0.0
        
        return stats
    
    def log_statistics(self) -> None:
        """Log comprehensive processing statistics."""
        stats = self.get_statistics()
        
        logger.info("=== Evidence Tiering Statistics ===")
        logger.info(f"Total variants processed: {stats['total_variants_processed']:,}")
        logger.info(f"RNA consensus variants (>={self.min_rna_support} callers): {stats['rna_consensus_variants']:,} ({stats['rna_consensus_rate']:.1f}%)")
        logger.info(f"RNA-only variants (DNA=0): {stats['rna_only_variants']:,} ({stats['rna_only_rate']:.1f}%)")
        logger.info(f"FILTER updates to RNAedit: {stats['filter_updates']:,} ({stats['filter_update_rate']:.1f}%)")
        
        logger.info("Evidence level distribution:")
        for level, count in stats['evidence_levels'].items():
            percentage = (count / stats['total_variants_processed'] * 100) if stats['total_variants_processed'] > 0 else 0
            logger.info(f"  {level}: {count:,} ({percentage:.1f}%)")
        
        # Log detailed classification decisions for debugging
        if stats['total_variants_processed'] > 0:
            logger.debug("=== Evidence Tiering Decision Details ===")
            logger.debug(f"Classification criteria (min_rna_support={self.min_rna_support}):")
            logger.debug("  HIGH: RNA consensus + RNA-only + REDIportal match")
            logger.debug("  MEDIUM: RNA consensus + has DNA support + REDIportal match")
            logger.debug("  LOW: RNA consensus + REDIportal match (fallback)")
            logger.debug("  NONE: No RNA consensus OR no REDIportal match")
            
            # Calculate and log decision pathway statistics
            rna_consensus_with_rediportal = stats['filter_updates']  # These are the ones that got updated
            rna_consensus_without_rediportal = stats['rna_consensus_variants'] - rna_consensus_with_rediportal
            
            if rna_consensus_without_rediportal > 0:
                logger.info(f"RNA consensus variants without REDIportal match: {rna_consensus_without_rediportal:,}")
            
            # Log efficiency metrics
            if stats['rna_consensus_variants'] > 0:
                rediportal_match_rate = (rna_consensus_with_rediportal / stats['rna_consensus_variants']) * 100
                logger.info(f"REDIportal match rate among RNA consensus variants: {rediportal_match_rate:.1f}%")
    
    def log_classification_decision(self, variant_info: Dict[str, Any], result: Dict[str, Any]):
        """Log individual classification decisions for debugging."""
        if logger.isEnabledFor(logging.DEBUG):
            chrom = variant_info.get('CHROM', 'unknown')
            pos = variant_info.get('POS', 'unknown')
            
            logger.debug(f"Classification decision for {chrom}:{pos}")
            logger.debug(f"  RNA support: {result['rna_support']}, DNA support: {result['dna_support']}")
            logger.debug(f"  RNA consensus: {result['rna_consensus']}, RNA-only: {result['rna_only']}")
            logger.debug(f"  Evidence tier: {result['evidence_tier']}")
            logger.debug(f"  FILTER update: {result['update_filter']}")
    
    def log_processing_progress(self, processed_count: int, interval: int = 10000):
        """Log processing progress at regular intervals."""
        if processed_count % interval == 0:
            logger.info(f"Evidence tiering progress: {processed_count:,} variants processed")
            
            # Log current statistics
            current_stats = self.get_statistics()
            if current_stats['total_variants_processed'] > 0:
                logger.info(f"  Current RNA consensus rate: {current_stats['rna_consensus_rate']:.1f}%")
                logger.info(f"  Current FILTER update rate: {current_stats['filter_update_rate']:.1f}%")


class EvidenceTieringValidator:
    """
    Validator for evidence tiering results.
    
    This class provides validation methods to ensure evidence tiering
    results are consistent and follow the expected rules.
    """
    
    @staticmethod
    def validate_caller_support_values(rna_support: int, dna_support: int) -> bool:
        """
        Validate caller support values are non-negative integers.
        
        Args:
            rna_support: RNA caller support value
            dna_support: DNA caller support value
            
        Returns:
            True if valid, False otherwise
        """
        try:
            return (isinstance(rna_support, int) and rna_support >= 0 and
                   isinstance(dna_support, int) and dna_support >= 0)
        except Exception:
            return False
    
    @staticmethod
    def validate_evidence_tier(evidence_tier: str) -> bool:
        """
        Validate evidence tier is one of the expected values.
        
        Args:
            evidence_tier: Evidence tier string
            
        Returns:
            True if valid, False otherwise
        """
        return evidence_tier in ['HIGH', 'MEDIUM', 'LOW', 'NONE']
    
    @staticmethod
    def validate_tiering_logic(rna_support: int, dna_support: int, 
                             has_rediportal_match: bool, evidence_tier: str,
                             min_rna_support: int) -> bool:
        """
        Validate that evidence tier follows the expected logic.
        
        Args:
            rna_support: RNA caller support value
            dna_support: DNA caller support value
            has_rediportal_match: Whether variant has REDIportal match
            evidence_tier: Assigned evidence tier
            min_rna_support: Minimum RNA support threshold
            
        Returns:
            True if logic is consistent, False otherwise
        """
        # NONE tier validation
        if not has_rediportal_match or rna_support < min_rna_support:
            return evidence_tier == 'NONE'
        
        # HIGH tier validation (RNA consensus + RNA-only + REDIportal match)
        if rna_support >= min_rna_support and dna_support == 0 and has_rediportal_match:
            return evidence_tier == 'HIGH'
        
        # MEDIUM tier validation (RNA consensus + has DNA + REDIportal match)
        if rna_support >= min_rna_support and dna_support > 0 and has_rediportal_match:
            return evidence_tier == 'MEDIUM'
        
        return False
    
    @staticmethod
    def validate_filter_update_logic(evidence_tier: str, update_filter: bool) -> bool:
        """
        Validate that filter update decision follows the expected logic.
        
        Args:
            evidence_tier: Evidence tier
            update_filter: Whether filter should be updated
            
        Returns:
            True if logic is consistent, False otherwise
        """
        if evidence_tier in ['HIGH', 'MEDIUM', 'LOW']:
            return update_filter is True
        elif evidence_tier == 'NONE':
            return update_filter is False
        else:
            return False


def create_evidence_tiering_processor(min_rna_support: int = 2) -> EvidenceTieringProcessor:
    """
    Factory function to create evidence tiering processor.
    
    Args:
        min_rna_support: Minimum RNA caller support threshold
        
    Returns:
        Configured EvidenceTieringProcessor instance
    """
    return EvidenceTieringProcessor(min_rna_support=min_rna_support)


if __name__ == '__main__':
    # Simple testing functionality
    import json
    
    # Test evidence tiering processor
    processor = create_evidence_tiering_processor(min_rna_support=2)
    
    # Test cases
    test_cases = [
        {
            'name': 'HIGH evidence (RNA-only with consensus)',
            'variant_info': {'N_RNA_CALLERS_SUPPORT': 3, 'N_DNA_CALLERS_SUPPORT': 0},
            'has_rediportal_match': True,
            'expected_tier': 'HIGH'
        },
        {
            'name': 'MEDIUM evidence (RNA consensus with DNA support)',
            'variant_info': {'N_RNA_CALLERS_SUPPORT': 2, 'N_DNA_CALLERS_SUPPORT': 1},
            'has_rediportal_match': True,
            'expected_tier': 'MEDIUM'
        },
        {
            'name': 'NONE evidence (no RNA consensus)',
            'variant_info': {'N_RNA_CALLERS_SUPPORT': 1, 'N_DNA_CALLERS_SUPPORT': 0},
            'has_rediportal_match': True,
            'expected_tier': 'NONE'
        },
        {
            'name': 'NONE evidence (no REDIportal match)',
            'variant_info': {'N_RNA_CALLERS_SUPPORT': 3, 'N_DNA_CALLERS_SUPPORT': 0},
            'has_rediportal_match': False,
            'expected_tier': 'NONE'
        }
    ]
    
    print("Testing Evidence Tiering Processor:")
    print("=" * 50)
    
    for test_case in test_cases:
        result = processor.process_variant(
            test_case['variant_info'], 
            test_case['has_rediportal_match']
        )
        
        status = "✓" if result['evidence_tier'] == test_case['expected_tier'] else "✗"
        print(f"{status} {test_case['name']}")
        print(f"   RNA support: {result['rna_support']}, DNA support: {result['dna_support']}")
        print(f"   RNA consensus: {result['rna_consensus']}, RNA-only: {result['rna_only']}")
        print(f"   Evidence tier: {result['evidence_tier']} (expected: {test_case['expected_tier']})")
        print(f"   Update filter: {result['update_filter']}")
        print()
    
    # Print statistics
    print("Processing Statistics:")
    processor.log_statistics()
    
    # Test validator
    print("\nTesting Evidence Tiering Validator:")
    print("=" * 50)
    
    validator = EvidenceTieringValidator()
    
    # Test validation functions
    validation_tests = [
        ('Valid caller support', validator.validate_caller_support_values(3, 1), True),
        ('Invalid negative RNA support', validator.validate_caller_support_values(-1, 1), False),
        ('Valid evidence tier', validator.validate_evidence_tier('HIGH'), True),
        ('Invalid evidence tier', validator.validate_evidence_tier('INVALID'), False),
        ('Valid HIGH tier logic', validator.validate_tiering_logic(3, 0, True, 'HIGH', 2), True),
        ('Invalid HIGH tier logic', validator.validate_tiering_logic(1, 0, True, 'HIGH', 2), False),
        ('Valid filter update logic', validator.validate_filter_update_logic('HIGH', True), True),
        ('Invalid filter update logic', validator.validate_filter_update_logic('NONE', True), False)
    ]
    
    for test_name, result, expected in validation_tests:
        status = "✓" if result == expected else "✗"
        print(f"{status} {test_name}: {result} (expected: {expected})")