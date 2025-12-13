#!/usr/bin/env python3
"""
RNA Editing Core Functions

This module provides core functionality for RNA editing evidence classification
and biological categorization. It implements the four-tier evidence matrix
and FILTER update logic for RNA editing annotation.

SAFETY FEATURES:
- Comprehensive input validation
- Defensive programming practices
- Extensive error handling and logging
- No file system operations (pure functions)

Author: RNA Editing Enhancement Pipeline
Date: 2025-12-13
"""

import logging
from typing import Dict, Any, Union, Optional

logger = logging.getLogger(__name__)

# Evidence levels in order of strength
EVIDENCE_LEVELS = ['VERY_HIGH', 'HIGH', 'MEDIUM', 'LOW', 'NONE']

# Canonical RNA editing transitions
CANONICAL_TRANSITIONS = {
    ('A', 'G'): True,  # A-to-G editing
    ('T', 'C'): True,  # T-to-C editing (reverse complement of A-to-G)
    ('G', 'A'): False, # Not canonical
    ('C', 'T'): False, # Not canonical
}

def is_canonical_editing_transition(ref: str, alt: str) -> bool:
    """
    Check if a variant represents a canonical RNA editing transition.
    
    Canonical transitions are A>G and T>C (which are reverse complements).
    
    Args:
        ref: Reference allele
        alt: Alternative allele(s), comma-separated
        
    Returns:
        True if canonical transition, False otherwise
    """
    if not ref or not alt:
        return False
    
    # Handle multiple alternative alleles
    alt_alleles = alt.split(',')
    
    # Check if any alternative allele forms a canonical transition
    for alt_allele in alt_alleles:
        alt_allele = alt_allele.strip()
        if not alt_allele:
            continue
            
        # Only consider single nucleotide variants
        if len(ref) == 1 and len(alt_allele) == 1:
            transition = (ref.upper(), alt_allele.upper())
            if CANONICAL_TRANSITIONS.get(transition, False):
                return True
    
    return False

def meets_rna_support_threshold(variant_data: Dict[str, Any], min_support: int) -> bool:
    """
    Check if variant meets RNA caller support threshold.
    
    Args:
        variant_data: Dictionary containing variant information
        min_support: Minimum required RNA caller support
        
    Returns:
        True if meets threshold, False otherwise
    """
    rna_support = variant_data.get('N_RNA_CALLERS_SUPPORT', 0)
    
    try:
        rna_support = int(rna_support)
        return rna_support >= min_support
    except (ValueError, TypeError):
        logger.warning(f"Invalid RNA support value: {rna_support}")
        return False

def classify_rna_editing_evidence(variant_data: Dict[str, Any], 
                                exact_rediportal_match: bool, 
                                min_rna_support: int) -> str:
    """
    Classify RNA editing evidence level using the four-tier matrix.
    
    Evidence Classification Matrix:
    - VERY_HIGH: Canonical + RNA-only + High RNA VAF + Low DNA VAF + REDIportal + RNA consensus
    - HIGH: Canonical + RNA-only + High RNA VAF + REDIportal + RNA consensus
    - MEDIUM: Canonical + REDIportal + RNA consensus (has DNA presence)
    - LOW: Non-canonical + REDIportal + RNA consensus
    - NONE: No exact REDIportal match OR insufficient RNA support
    
    Args:
        variant_data: Dictionary containing variant information
        exact_rediportal_match: Whether variant has exact REDIportal coordinate+allele match
        min_rna_support: Minimum RNA caller support threshold
        
    Returns:
        Evidence level string (VERY_HIGH, HIGH, MEDIUM, LOW, NONE)
    """
    try:
        # Prerequisites for any evidence classification
        if not exact_rediportal_match:
            return 'NONE'
        
        if not meets_rna_support_threshold(variant_data, min_rna_support):
            return 'NONE'
        
        # Extract variant characteristics
        ref = variant_data.get('REF', '')
        alt = variant_data.get('ALT', '')
        rna_support = variant_data.get('N_RNA_CALLERS_SUPPORT', 0)
        dna_support = variant_data.get('N_DNA_CALLERS_SUPPORT', 0)
        rna_vaf = variant_data.get('VAF_RNA_MEAN', 0.0)
        dna_vaf = variant_data.get('VAF_DNA_MEAN', 0.0)
        
        # Convert to appropriate types
        try:
            rna_support = int(rna_support)
            dna_support = int(dna_support)
            rna_vaf = float(rna_vaf)
            dna_vaf = float(dna_vaf)
        except (ValueError, TypeError) as e:
            logger.warning(f"Invalid numeric values in variant data: {e}")
            return 'LOW'  # Default to LOW if we can't parse values properly
        
        # Check canonical transition
        canonical = is_canonical_editing_transition(ref, alt)
        
        # Determine RNA-only support (no DNA callers support this variant)
        rna_only = (dna_support == 0)
        
        # Define VAF thresholds
        high_rna_vaf = (rna_vaf >= 0.1)  # 10% or higher
        low_dna_vaf = (dna_vaf <= 0.05)  # 5% or lower
        
        # Apply four-tier classification matrix
        if canonical and rna_only and high_rna_vaf and low_dna_vaf:
            return 'VERY_HIGH'
        elif canonical and rna_only and high_rna_vaf:
            return 'HIGH'
        elif canonical:
            return 'MEDIUM'
        else:
            return 'LOW'
            
    except Exception as e:
        logger.error(f"Error in evidence classification: {e}")
        return 'NONE'

def classify_rna_editing_biological_category(evidence_level: str, 
                                           original_classification: str) -> str:
    """
    Determine biological classification (FILTER field) based on evidence level.
    
    Rules:
    - VERY_HIGH, HIGH, MEDIUM, LOW evidence -> RNAedit
    - NONE evidence -> preserve original classification
    
    Args:
        evidence_level: RNA editing evidence level
        original_classification: Original FILTER classification
        
    Returns:
        Biological classification for FILTER field
    """
    if evidence_level in ['VERY_HIGH', 'HIGH', 'MEDIUM', 'LOW']:
        return 'RNAedit'
    else:
        return original_classification

def initialize_category_transition_stats() -> Dict[str, Any]:
    """
    Initialize statistics tracking for category transitions.
    
    Returns:
        Dictionary structure for tracking statistics
    """
    return {
        'category_transitions': {},
        'evidence_distribution': {
            'reclassified': {'VERY_HIGH': 0, 'HIGH': 0, 'MEDIUM': 0, 'LOW': 0},
            'preserved': {'NONE': 0}
        },
        'rediportal_matches': {
            'exact_matches': 0,
            'total_variants': 0
        }
    }

def update_category_transition_stats(stats: Dict[str, Any], 
                                   original_category: str,
                                   final_category: str,
                                   evidence_level: str,
                                   exact_rediportal_match: bool) -> Dict[str, Any]:
    """
    Update category transition statistics.
    
    Args:
        stats: Statistics dictionary to update
        original_category: Original FILTER category
        final_category: Final FILTER category
        evidence_level: RNA editing evidence level
        exact_rediportal_match: Whether variant had exact REDIportal match
        
    Returns:
        Updated statistics dictionary
    """
    # Track total variants
    stats['rediportal_matches']['total_variants'] += 1
    
    # Track REDIportal matches
    if exact_rediportal_match:
        stats['rediportal_matches']['exact_matches'] += 1
    
    # Initialize category if not seen before
    if original_category not in stats['category_transitions']:
        stats['category_transitions'][original_category] = {
            'total': 0,
            'to_rnaedit': 0,
            'percentage': 0.0
        }
    
    # Update category statistics
    stats['category_transitions'][original_category]['total'] += 1
    
    if final_category == 'RNAedit':
        stats['category_transitions'][original_category]['to_rnaedit'] += 1
        stats['evidence_distribution']['reclassified'][evidence_level] += 1
    else:
        stats['evidence_distribution']['preserved'][evidence_level] += 1
    
    return stats

def finalize_category_transition_stats(stats: Dict[str, Any]) -> Dict[str, Any]:
    """
    Finalize and calculate summary statistics.
    
    Args:
        stats: Statistics dictionary to finalize
        
    Returns:
        Finalized statistics with calculated percentages and summaries
    """
    # Calculate percentages for each category
    for category, category_stats in stats['category_transitions'].items():
        if category_stats['total'] > 0:
            category_stats['percentage'] = (category_stats['to_rnaedit'] / category_stats['total']) * 100
    
    # Calculate summary statistics
    total_variants = stats['rediportal_matches']['total_variants']
    exact_matches = stats['rediportal_matches']['exact_matches']
    total_reclassified = sum(stats['evidence_distribution']['reclassified'].values())
    
    stats['summary'] = {
        'total_variants_processed': total_variants,
        'total_reclassified_to_rnaedit': total_reclassified,
        'overall_reclassification_rate': (total_reclassified / total_variants * 100) if total_variants > 0 else 0.0,
        'exact_match_rate': (exact_matches / total_variants * 100) if total_variants > 0 else 0.0
    }
    
    return stats

def validate_variant_data(variant_data: Dict[str, Any]) -> bool:
    """
    Validate variant data structure for evidence classification.
    
    Args:
        variant_data: Dictionary containing variant information
        
    Returns:
        True if valid, False otherwise
    """
    required_fields = ['CHROM', 'POS', 'REF', 'ALT']
    
    for field in required_fields:
        if field not in variant_data:
            logger.warning(f"Missing required field: {field}")
            return False
    
    # Validate numeric fields if present
    numeric_fields = ['N_RNA_CALLERS_SUPPORT', 'N_DNA_CALLERS_SUPPORT', 'VAF_RNA_MEAN', 'VAF_DNA_MEAN']
    
    for field in numeric_fields:
        if field in variant_data:
            try:
                if field.startswith('VAF'):
                    float(variant_data[field])
                else:
                    int(variant_data[field])
            except (ValueError, TypeError):
                logger.warning(f"Invalid numeric value for {field}: {variant_data[field]}")
                return False
    
    return True

def get_evidence_level_description(evidence_level: str) -> str:
    """
    Get human-readable description of evidence level.
    
    Args:
        evidence_level: Evidence level string
        
    Returns:
        Description of the evidence level
    """
    descriptions = {
        'VERY_HIGH': 'Canonical transition, RNA-only support, high RNA VAF (≥10%), low DNA VAF (≤5%), exact REDIportal match, RNA consensus',
        'HIGH': 'Canonical transition, RNA-only support, high RNA VAF (≥10%), exact REDIportal match, RNA consensus',
        'MEDIUM': 'Canonical transition, exact REDIportal match, RNA consensus (has DNA presence)',
        'LOW': 'Non-canonical transition, exact REDIportal match, RNA consensus',
        'NONE': 'No exact REDIportal match OR insufficient RNA support (<2 callers)'
    }
    
    return descriptions.get(evidence_level, 'Unknown evidence level')

def summarize_evidence_classification(variant_data: Dict[str, Any], 
                                    evidence_level: str,
                                    exact_rediportal_match: bool) -> Dict[str, Any]:
    """
    Generate summary of evidence classification for a variant.
    
    Args:
        variant_data: Dictionary containing variant information
        evidence_level: Classified evidence level
        exact_rediportal_match: Whether variant had exact REDIportal match
        
    Returns:
        Summary dictionary with classification details
    """
    ref = variant_data.get('REF', '')
    alt = variant_data.get('ALT', '')
    canonical = is_canonical_editing_transition(ref, alt)
    
    return {
        'variant': f"{variant_data.get('CHROM', '')}:{variant_data.get('POS', '')}{ref}>{alt}",
        'evidence_level': evidence_level,
        'evidence_description': get_evidence_level_description(evidence_level),
        'canonical_transition': canonical,
        'exact_rediportal_match': exact_rediportal_match,
        'rna_support': variant_data.get('N_RNA_CALLERS_SUPPORT', 0),
        'dna_support': variant_data.get('N_DNA_CALLERS_SUPPORT', 0),
        'rna_vaf': variant_data.get('VAF_RNA_MEAN', 0.0),
        'dna_vaf': variant_data.get('VAF_DNA_MEAN', 0.0)
    }

if __name__ == '__main__':
    # Simple testing functionality
    import json
    
    # Test canonical transition detection
    test_cases = [
        ('A', 'G', True),
        ('T', 'C', True),
        ('G', 'A', False),
        ('C', 'T', False),
        ('A', 'T', False),
        ('G', 'C', False)
    ]
    
    print("Testing canonical transition detection:")
    for ref, alt, expected in test_cases:
        result = is_canonical_editing_transition(ref, alt)
        status = "✓" if result == expected else "✗"
        print(f"  {status} {ref}>{alt}: {result} (expected: {expected})")
    
    # Test evidence classification
    print("\nTesting evidence classification:")
    
    test_variant = {
        'CHROM': 'chr1',
        'POS': 12345,
        'REF': 'A',
        'ALT': 'G',
        'N_RNA_CALLERS_SUPPORT': 3,
        'N_DNA_CALLERS_SUPPORT': 0,
        'VAF_RNA_MEAN': 0.15,
        'VAF_DNA_MEAN': 0.02
    }
    
    evidence = classify_rna_editing_evidence(test_variant, True, 2)
    print(f"  Evidence level: {evidence}")
    
    summary = summarize_evidence_classification(test_variant, evidence, True)
    print(f"  Summary: {json.dumps(summary, indent=2)}")