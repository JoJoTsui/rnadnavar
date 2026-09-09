"""
Statistics generation functions for consensus and rescue workflows.

This module provides reusable functions for generating statistics
about variant support, caller coverage, and cross-modality rescue.
"""


def compute_consensus_statistics(variant_data, snv_threshold, indel_threshold):
    """
    Compute statistics for consensus operation.
    
    This function calculates comprehensive statistics about variant support
    across callers, including total counts, variant types, consensus status,
    and caller distribution.
    
    Args:
        variant_data (dict): Aggregated variant data dictionary where keys are
            variant identifiers and values contain variant information including
            'is_snv', 'callers', and 'passes_consensus' fields
        snv_threshold (int): Minimum number of callers required for SNV consensus
        indel_threshold (int): Minimum number of callers required for indel consensus
    
    Returns:
        dict: Statistics summary containing:
            - total_variants: Total number of variants
            - snvs: Number of SNVs
            - indels: Number of indels
            - snvs_consensus: Number of SNVs passing consensus threshold
            - indels_consensus: Number of indels passing consensus threshold
            - single_caller: Number of variants supported by single caller
            - multi_caller: Number of variants supported by multiple callers
    
    Example:
        >>> stats = compute_consensus_statistics(variant_data, snv_thr=2, indel_thr=2)
        >>> print(f"Total variants: {stats['total_variants']}")
    """
    stats = {
        'total_variants': len(variant_data),
        'snvs': 0,
        'indels': 0,
        'snvs_consensus': 0,
        'indels_consensus': 0,
        'single_caller': 0,
        'multi_caller': 0,
    }
    
    for vkey, data in variant_data.items():
        # Count variant types and consensus status
        if data['is_snv']:
            stats['snvs'] += 1
            if data.get('passes_consensus', False):
                stats['snvs_consensus'] += 1
        else:
            stats['indels'] += 1
            if data.get('passes_consensus', False):
                stats['indels_consensus'] += 1
        
        # Count caller distribution
        unique_callers = len(set(data['callers']))
        if unique_callers == 1:
            stats['single_caller'] += 1
        else:
            stats['multi_caller'] += 1
    
    return stats


def compute_rescue_statistics(variant_data, dna_variants, rna_variants):
    """
    Compute statistics for rescue operation.
    
    This function calculates statistics about cross-modality variant support,
    including DNA-only, RNA-only, and cross-modality variants, as well as
    descriptive outcome counts. These are not truth-based accuracy metrics.
    
    Args:
        variant_data (dict): Aggregated variant data dictionary where keys are
            variant identifiers and values contain variant information including
            'dna_support', 'rna_support', 'rescued', and 'is_snv' fields
        dna_variants (dict): DNA consensus records with filter_normalized labels
        rna_variants (dict): RNA consensus records with filter_normalized labels
    
    Returns:
        dict: Rescue statistics summary containing:
            - total_variants: Total number of variants
            - dna_only: Number of variants with DNA support only
            - rna_only: Number of variants with RNA support only
            - cross_modality: Number of variants with both DNA and RNA support
            - rescued: Number of variants rescued by cross-modality support
            - snvs: Number of SNVs
            - indels: Number of indels
            - snvs_rescued: Number of SNVs rescued
            - indels_rescued: Number of indels rescued
            - rescued_union_fraction: rescued / total union records (0..1)
            - dna_somatic_baseline: DNA input Somatic records
            - output_somatic: Final output Somatic records
            - somatic_retained_from_dna: Final Somatic and DNA input Somatic
            - somatic_new_vs_dna: Final Somatic but not DNA input Somatic
            - somatic_lost_from_dna: DNA input Somatic but not final Somatic
    
    Example:
        >>> stats = compute_rescue_statistics(variant_data, dna_records, rna_records)
        >>> print(stats['rescued_union_fraction'])
    """
    stats = {
        'total_variants': len(variant_data),
        'dna_only': 0,
        'rna_only': 0,
        'cross_modality': 0,
        'rescued': 0,
        'snvs': 0,
        'indels': 0,
        'snvs_rescued': 0,
        'indels_rescued': 0,
    }
    if not isinstance(dna_variants, dict) or not isinstance(rna_variants, dict):
        raise TypeError("Rescue outcome statistics require labeled consensus record dictionaries")
    dna_somatic = {key for key, record in dna_variants.items()
                   if record.get('filter_normalized') == 'Somatic'}
    missing_final = [key for key, data in variant_data.items() if 'final_classification' not in data]
    if missing_final:
        raise ValueError("Rescue outcome statistics require finalized output classifications")
    output_somatic = {key for key, data in variant_data.items()
                      if data['final_classification'] == 'Somatic'}
    stats.update({
        'dna_somatic_baseline': len(dna_somatic),
        'output_somatic': len(output_somatic),
        'somatic_retained_from_dna': len(output_somatic & dna_somatic),
        'somatic_new_vs_dna': len(output_somatic - dna_somatic),
        'somatic_lost_from_dna': len(dna_somatic - output_somatic),
    })
    
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
        if data.get('rescued', False) and vkey in output_somatic:
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
    
    stats['rescued_union_fraction'] = (
        stats['rescued'] / stats['total_variants'] if stats['total_variants'] else 0.0
    )
    
    return stats


def print_statistics(stats, operation_type='consensus'):
    """
    Print formatted statistics to stdout.
    
    This function formats and prints statistics in a human-readable format,
    supporting both consensus and rescue operation types.
    
    Args:
        stats (dict): Statistics dictionary from compute_consensus_statistics()
            or compute_rescue_statistics()
        operation_type (str): Type of operation - 'consensus' or 'rescue'
            Default: 'consensus'
    
    Example:
        >>> stats = compute_consensus_statistics(variant_data, 2, 2)
        >>> print_statistics(stats, 'consensus')
        
        >>> rescue_stats = compute_rescue_statistics(variant_data, dna_records, rna_records)
        >>> print_statistics(rescue_stats, 'rescue')
    """
    print("\n- Statistics:")
    print(f"  - Total variants: {stats['total_variants']:,}")
    
    if operation_type == 'consensus':
        # Print consensus-specific statistics
        print(f"  - SNVs: {stats['snvs']:,} (consensus: {stats['snvs_consensus']:,})")
        print(f"  - Indels: {stats['indels']:,} (consensus: {stats['indels_consensus']:,})")
        print(f"  - Single caller: {stats['single_caller']:,}")
        print(f"  - Multiple callers: {stats['multi_caller']:,}")
        
    elif operation_type == 'rescue':
        # Print rescue-specific statistics
        print(f"  - SNVs: {stats['snvs']:,} (rescued: {stats['snvs_rescued']:,})")
        print(f"  - Indels: {stats['indels']:,} (rescued: {stats['indels_rescued']:,})")
        print("\n- Modality Support:")
        print(f"  - DNA only: {stats['dna_only']:,}")
        print(f"  - RNA only: {stats['rna_only']:,}")
        print(f"  - Cross-modality: {stats['cross_modality']:,}")
        print("\n- Rescue Outcomes (descriptive counts; not truth-based accuracy):")
        print(f"  - Total rescued: {stats['rescued']:,}")
        print(f"  - Rescued / union records: {stats['rescued']:,}/{stats['total_variants']:,} "
              f"({stats['rescued_union_fraction']:.2%})")
        print(f"  - DNA input Somatic baseline: {stats['dna_somatic_baseline']:,}")
        print(f"  - Final output Somatic: {stats['output_somatic']:,}")
        print(f"  - Somatic retained from DNA: {stats['somatic_retained_from_dna']:,}")
        print(f"  - Somatic new versus DNA: {stats['somatic_new_vs_dna']:,}")
        print(f"  - Somatic lost from DNA: {stats['somatic_lost_from_dna']:,}")
    
    else:
        raise ValueError(f"Unknown operation_type: {operation_type}. Must be 'consensus' or 'rescue'.")
