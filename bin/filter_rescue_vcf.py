#!/usr/bin/env python3
"""
Rescue VCF filtering script with modality-aware variant classification.

This script processes rescued VCF files and determines variant types based on
modality-specific caller information and filter categories.
"""

import argparse
import sys
from pathlib import Path
from cyvcf2 import VCF, Writer


def argparser():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description='Filter rescued VCF with modality-aware variant classification',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        '--input',
        type=str,
        required=True,
        help='Input rescued VCF file'
    )
    parser.add_argument(
        '--output',
        type=str,
        required=True,
        help='Output filtered VCF file'
    )
    parser.add_argument(
        '--dna_threshold',
        type=int,
        default=2,
        help='Minimum number of DNA callers for consensus'
    )
    parser.add_argument(
        '--rna_threshold',
        type=int,
        default=2,
        help='Minimum number of RNA callers for consensus'
    )
    
    return parser.parse_args()


def determine_modality_variant_type(filters_normalized, filters_category, n_callers):
    """
    Determine variant type for a single modality based on caller filters.
    
    Args:
        filters_normalized: List of normalized filter values for this modality
        filters_category: List of filter categories for this modality
        n_callers: Number of callers that detected this variant in this modality
    
    Returns:
        str: Variant type - 'SOMATIC', 'GERMLINE', 'REFERENCE', or 'ARTIFACT'
    """
    if n_callers == 0:
        return None
    
    # Count filter types
    germline_count = sum(1 for f in filters_normalized if 'Germline' in f or 'germline' in f)
    reference_count = sum(1 for f in filters_normalized if 'ReferenceCall' in f or 'reference' in f)
    pass_count = sum(1 for f in filters_normalized if f == 'PASS')
    
    # Determine variant type by majority vote
    if germline_count >= n_callers / 2:
        return 'GERMLINE'
    elif reference_count >= n_callers / 2:
        return 'REFERENCE'
    elif pass_count >= n_callers / 2:
        return 'SOMATIC'
    else:
        return 'ARTIFACT'


def determine_final_variant_type(dna_type, rna_type, dna_consensus, rna_consensus, 
                                  dna_n_callers, rna_n_callers, dna_threshold, rna_threshold):
    """
    Determine final variant type across modalities.
    
    Rules:
    1. If DNA and RNA match, use the common type
    2. If DNA and RNA don't match:
       - If DNA is consensus (>= threshold), use DNA type
       - Else if RNA is consensus (>= threshold), use RNA type
       - Else use DNA type
    
    Args:
        dna_type: DNA variant type
        rna_type: RNA variant type
        dna_consensus: Whether DNA passes consensus
        rna_consensus: Whether RNA passes consensus
        dna_n_callers: Number of DNA callers
        rna_n_callers: Number of RNA callers
        dna_threshold: DNA consensus threshold
        rna_threshold: RNA consensus threshold
    
    Returns:
        str: Final variant type
    """
    # If only one modality has data, use that
    if dna_type and not rna_type:
        return dna_type
    if rna_type and not dna_type:
        return rna_type
    
    # If both modalities agree, use common type
    if dna_type == rna_type:
        return dna_type
    
    # If they disagree, use consensus logic
    if dna_consensus and dna_n_callers >= dna_threshold:
        return dna_type
    elif rna_consensus and rna_n_callers >= rna_threshold:
        return rna_type
    else:
        # Default to DNA type
        return dna_type if dna_type else rna_type


def parse_modality_filters(filters_str, modality):
    """
    Parse filter string and extract filters for specific modality.
    
    Args:
        filters_str: Filter string in format "DNA_caller:filter|RNA_caller:filter|..."
        modality: 'DNA' or 'RNA'
    
    Returns:
        list: List of filter values for the specified modality
    """
    if not filters_str or filters_str == '.':
        return []
    
    filters = []
    for item in filters_str.split('|'):
        if ':' in item:
            caller_part, filter_val = item.split(':', 1)
            if caller_part.startswith(f'{modality}_'):
                filters.append(filter_val)
    
    return filters


def main():
    """Main filtering workflow."""
    args = argparser()
    
    print("=" * 80)
    print("Rescue VCF Filtering with Modality-Aware Classification")
    print("=" * 80)
    print(f"\nInput: {args.input}")
    print(f"Output: {args.output}")
    print(f"DNA threshold: {args.dna_threshold}")
    print(f"RNA threshold: {args.rna_threshold}")
    
    # Open input VCF
    vcf_in = VCF(args.input)
    
    # Add new INFO fields for modality-specific variant types
    vcf_in.add_info_to_header({
        'ID': 'DNA_VARIANT_TYPE',
        'Number': '1',
        'Type': 'String',
        'Description': 'Variant type determined from DNA callers (SOMATIC, GERMLINE, REFERENCE, ARTIFACT)'
    })
    vcf_in.add_info_to_header({
        'ID': 'RNA_VARIANT_TYPE',
        'Number': '1',
        'Type': 'String',
        'Description': 'Variant type determined from RNA callers (SOMATIC, GERMLINE, REFERENCE, ARTIFACT)'
    })
    vcf_in.add_info_to_header({
        'ID': 'FINAL_VARIANT_TYPE',
        'Number': '1',
        'Type': 'String',
        'Description': 'Final variant type across modalities (SOMATIC, GERMLINE, REFERENCE, ARTIFACT)'
    })
    
    # Add FILTER entries if not present
    if 'SOMATIC' not in vcf_in.header_iter():
        vcf_in.add_filter_to_header({'ID': 'SOMATIC', 'Description': 'Somatic variant'})
    if 'GERMLINE' not in vcf_in.header_iter():
        vcf_in.add_filter_to_header({'ID': 'GERMLINE', 'Description': 'Germline variant'})
    if 'REFERENCE' not in vcf_in.header_iter():
        vcf_in.add_filter_to_header({'ID': 'REFERENCE', 'Description': 'Reference call'})
    if 'ARTIFACT' not in vcf_in.header_iter():
        vcf_in.add_filter_to_header({'ID': 'ARTIFACT', 'Description': 'Likely artifact'})
    
    # Open output VCF
    vcf_out = Writer(args.output, vcf_in)
    
    # Process variants
    processed_count = 0
    for variant in vcf_in:
        # Extract modality-specific information
        filters_normalized = variant.INFO.get('FILTERS_NORMALIZED', '.')
        filters_category = variant.INFO.get('FILTERS_CATEGORY', '.')
        
        dna_n_callers = variant.INFO.get('N_DNA_CALLERS', 0)
        rna_n_callers = variant.INFO.get('N_RNA_CALLERS', 0)
        
        passes_consensus_dna = variant.INFO.get('PASSES_CONSENSUS_DNA', 'NO') == 'YES'
        passes_consensus_rna = variant.INFO.get('PASSES_CONSENSUS_RNA', 'NO') == 'YES'
        
        # Parse modality-specific filters
        dna_filters_normalized = parse_modality_filters(filters_normalized, 'DNA')
        rna_filters_normalized = parse_modality_filters(filters_normalized, 'RNA')
        
        dna_filters_category = parse_modality_filters(filters_category, 'DNA')
        rna_filters_category = parse_modality_filters(filters_category, 'RNA')
        
        # Determine modality-specific variant types
        dna_type = determine_modality_variant_type(
            dna_filters_normalized, dna_filters_category, dna_n_callers
        )
        rna_type = determine_modality_variant_type(
            rna_filters_normalized, rna_filters_category, rna_n_callers
        )
        
        # Determine final variant type
        final_type = determine_final_variant_type(
            dna_type, rna_type,
            passes_consensus_dna, passes_consensus_rna,
            dna_n_callers, rna_n_callers,
            args.dna_threshold, args.rna_threshold
        )
        
        # Update INFO fields
        if dna_type:
            variant.INFO['DNA_VARIANT_TYPE'] = dna_type
        if rna_type:
            variant.INFO['RNA_VARIANT_TYPE'] = rna_type
        if final_type:
            variant.INFO['FINAL_VARIANT_TYPE'] = final_type
        
        # Update FILTER column based on final type
        if final_type:
            if final_type == 'SOMATIC':
                # SOMATIC variants pass - set to empty list (PASS)
                variant.FILTER = []
            else:
                # Add appropriate filter
                variant.FILTER = [final_type]
        
        # Write variant
        vcf_out.write_record(variant)
        processed_count += 1
        
        if processed_count % 10000 == 0:
            print(f"  - Processed {processed_count:,} variants...")
    
    vcf_out.close()
    vcf_in.close()
    
    print(f"\n{'=' * 80}")
    print(f"Filtering completed successfully!")
    print(f"{'=' * 80}")
    print(f"  - Total variants processed: {processed_count:,}")
    print(f"  - Output file: {args.output}")
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
