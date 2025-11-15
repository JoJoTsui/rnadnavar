"""
Variant tagging functions for consensus and rescue workflows.

This module provides reusable functions for tagging variants with caller
information, modality metadata, and rescue indicators. These functions add
metadata to variant data dictionaries to support consensus calling and
cross-modality rescue operations.

Functions:
    tag_variant_with_callers: Add caller support information to variant
    tag_variant_with_modality: Add modality tracking to variant
    mark_rescued_variants: Mark variants with cross-modality support
    compute_unified_filter: Compute overall filter status from multiple callers

Example:
    >>> from vcf_utils.tagging import (tag_variant_with_callers, 
    ...                                 tag_variant_with_modality,
    ...                                 mark_rescued_variants)
    >>> 
    >>> # Tag variant with caller information
    >>> variant_data = {'callers': ['mutect2', 'strelka']}
    >>> all_callers = ['mutect2', 'strelka', 'deepsomatic']
    >>> tagged = tag_variant_with_callers(variant_data, all_callers)
    >>> print(f"Support: {tagged['n_support_callers']}/{tagged['n_callers']}")
    Support: 2/3
    >>> 
    >>> # Tag with modality
    >>> modality_map = {'mutect2': 'DNA', 'strelka': 'DNA'}
    >>> tagged = tag_variant_with_modality(variant_data, modality_map)
    >>> print(f"Modalities: {tagged['modalities']}")
    Modalities: {'DNA'}
    >>> 
    >>> # Mark rescued variants
    >>> dna_variants = {'chr1:100:A:G', 'chr1:200:C:T'}
    >>> rna_variants = {'chr1:100:A:G', 'chr1:300:G:A'}
    >>> aggregated = {'chr1:100:A:G': {}, 'chr1:200:C:T': {}, 'chr1:300:G:A': {}}
    >>> marked = mark_rescued_variants(aggregated, dna_variants, rna_variants)
    >>> print(f"Rescued: {marked['chr1:100:A:G']['rescued']}")
    Rescued: True
"""
from collections import Counter


def tag_variant_with_callers(variant_data, all_callers):
    """
    Tag variant with caller support information.
    
    This function adds caller-related metadata to variant data, including
    the list of all callers in the analysis and which callers support
    this specific variant.
    
    Args:
        variant_data: Variant data dict containing 'callers' key with list of supporting callers
        all_callers: List of all callers in the analysis
    
    Returns:
        dict: Updated variant data with caller tags:
            - 'all_callers': List of all callers in analysis
            - 'n_callers': Total number of callers
            - 'support_callers': Set of callers that detected this variant
            - 'n_support_callers': Number of supporting callers
    """
    # Store all callers
    variant_data['all_callers'] = all_callers
    variant_data['n_callers'] = len(all_callers)
    
    # Store supporting callers
    if 'support_callers' not in variant_data:
        variant_data['support_callers'] = set(variant_data.get('callers', []))
    variant_data['n_support_callers'] = len(variant_data['support_callers'])
    
    return variant_data


def tag_variant_with_modality(variant_data, modality_map):
    """
    Tag variant with modality information from modality_map.
    
    This function adds modality-related metadata to variant data, grouping
    callers by their modality (DNA or RNA) and computing modality-specific
    support counts.
    
    Args:
        variant_data: Variant data dict containing 'callers' key
        modality_map: Dict mapping caller names to modality ('DNA' or 'RNA')
    
    Returns:
        dict: Updated variant data with modality tags:
            - 'modalities': Set of modalities where variant was detected
            - 'dna_callers': List of DNA callers that detected this variant
            - 'rna_callers': List of RNA callers that detected this variant
            - 'caller_modality_map': Dict mapping each caller to its modality
    """
    # Initialize modality tracking
    modalities = set()
    dna_callers = []
    rna_callers = []
    caller_modality_map = {}
    
    # Group callers by modality
    for caller in variant_data.get('callers', []):
        modality = modality_map.get(caller, 'UNKNOWN')
        modalities.add(modality)
        caller_modality_map[caller] = modality
        
        if modality == 'DNA':
            dna_callers.append(caller)
        elif modality == 'RNA':
            rna_callers.append(caller)
    
    # Store modality information
    variant_data['modalities'] = modalities
    variant_data['dna_callers'] = dna_callers
    variant_data['rna_callers'] = rna_callers
    variant_data['caller_modality_map'] = caller_modality_map
    
    return variant_data


def mark_rescued_variants(variant_data, dna_variants, rna_variants):
    """
    Mark variants present in both DNA and RNA variant sets.
    
    This function identifies variants that have cross-modality support
    (present in both DNA and RNA consensus results) and marks them with
    rescue flags.
    
    Args:
        variant_data: Aggregated variant data dict (keyed by variant_key)
        dna_variants: Set of DNA variant keys
        rna_variants: Set of RNA variant keys
    
    Returns:
        dict: Updated variant data with rescue flags:
            - 'rescued': True if variant present in both DNA and RNA
            - 'dna_support': True if variant present in DNA
            - 'rna_support': True if variant present in RNA
            - 'cross_modality': True if variant has both DNA and RNA support
    """
    for vkey, data in variant_data.items():
        # Check if variant is in DNA and/or RNA sets
        in_dna = vkey in dna_variants
        in_rna = vkey in rna_variants
        
        # Set support flags
        data['dna_support'] = in_dna
        data['rna_support'] = in_rna
        data['cross_modality'] = in_dna and in_rna
        
        # Mark as rescued if present in both modalities
        data['rescued'] = in_dna and in_rna
    
    return variant_data


def compute_unified_filter(variant_data):
    """
    Compute unified filter status from multiple callers.
    
    This function determines the overall filter status for a variant based
    on the filter values from all supporting callers. It uses a majority
    vote approach: if at least half of the callers report PASS, the unified
    filter is PASS; otherwise, it's FAIL with the most common failure category.
    
    Args:
        variant_data: Variant data dict containing:
            - 'filters_normalized': List of normalized filter strings
            - 'filters_category': List of filter categories
            - 'callers': List of supporting callers
    
    Returns:
        tuple: (unified_filter_status, filter_tags)
            - unified_filter_status: 'PASS' or 'FAIL'
            - filter_tags: List of filter tag strings to add to FILTER column
    """
    filters_normalized = variant_data.get('filters_normalized', [])
    filters_category = variant_data.get('filters_category', [])
    n_support_callers = len(set(variant_data.get('callers', [])))
    
    # Count PASS filters
    pass_count = sum(1 for f in filters_normalized if f == 'PASS')
    
    # Majority vote: if at least half pass, unified filter is PASS
    if pass_count >= n_support_callers / 2:
        return 'PASS', []
    
    # Otherwise, find most common failure category
    non_pass_cats = [c for c in filters_category if c != 'PASS']
    filter_tags = []
    
    if non_pass_cats:
        most_common_cat = Counter(non_pass_cats).most_common(1)[0][0]
        
        # Map category to filter tag
        if 'quality' in most_common_cat.lower():
            filter_tags.append('LowQuality')
        elif 'depth' in most_common_cat.lower():
            filter_tags.append('LowDepth')
        elif 'bias' in most_common_cat.lower():
            filter_tags.append('StrandBias')
        elif 'germline' in most_common_cat.lower():
            filter_tags.append('Germline')
        elif 'artifact' in most_common_cat.lower():
            filter_tags.append('Artifact')
        else:
            filter_tags.append('LowQuality')  # Default fallback
    
    return 'FAIL', filter_tags
