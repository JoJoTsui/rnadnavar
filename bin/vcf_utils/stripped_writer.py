"""
Unified VCF writer for FORMAT-stripped outputs.

This module provides a raw text VCF writer that completely omits FORMAT/sample
columns, bypassing pysam's requirement to write these columns when samples exist.
"""

import pysam


def write_vcf_stripped(records, input_vcf_path, output_path, use_cyvcf2=False):
    """
    Write VCF WITHOUT FORMAT/sample columns using raw text writing.
    
    This bypasses pysam's VariantFile writer which always emits FORMAT/sample
    columns if samples exist in the header.
    
    Args:
        records: List of (variant, filter_status, filter_list) tuples
                 variant can be cyvcf2.Variant or pysam.VariantRecord
        input_vcf_path: Path to input VCF (for header template)
        output_path: Path to output VCF
        use_cyvcf2: If True, variants are cyvcf2.Variant objects
    
    Returns:
        int: Number of records written
    """
    # Open input to get header
    input_vcf = pysam.VariantFile(input_vcf_path)
    
    # Build header WITHOUT samples
    header_lines = []
    
    # Copy all header lines except sample-related
    for line in str(input_vcf.header).strip().split('\n'):
        # Skip the #CHROM header line - we'll write our own
        if line.startswith('#CHROM'):
            continue
        header_lines.append(line)
    
    # Add biological classification FILTER definitions
    classification_filters = {
        'Somatic': 'Somatic variant',
        'Germline': 'Germline variant',
        'Reference': 'Reference/wildtype',
        'Artifact': 'Artifact/technical error'
    }
    
    for filt, description in classification_filters.items():
        filter_line = f'##FILTER=<ID={filt},Description="{description}">'
        if filter_line not in header_lines:
            header_lines.append(filter_line)
    
    # Add RaVeX filter INFO fields
    if '##INFO=<ID=RaVeX_FILTER' not in '\n'.join(header_lines):
        header_lines.append('##INFO=<ID=RaVeX_FILTER,Number=.,Type=String,Description="RaVeX filter reasons: semicolon-separated list of filter flags">')
    
    filter_flags = {
        'min_alt_reads': 'Variant filtered due to insufficient alternate reads',
        'gnomad': 'Variant filtered due to high gnomAD allele frequency',
        'blacklist': 'Variant in blacklisted region',
        'noncoding': 'Variant in noncoding region',
        'ig_pseudo': 'Variant in immunoglobulin or pseudogene',
        'homopolymer': 'Variant in homopolymer region',
        'vc_filter': 'Variant failed variant caller filters',
        'not_consensus': 'Variant not in consensus',
        'multiallelic': 'Multiallelic site with multiple ALT alleles'
    }
    
    for flag, description in filter_flags.items():
        info_line = f'##INFO=<ID={flag},Number=0,Type=Flag,Description="{description}">'
        if info_line not in header_lines:
            header_lines.append(info_line)
    
    # Add column header WITHOUT FORMAT/sample columns
    header_lines.append('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO')
    
    # Build variant lines
    variant_lines = []
    written_count = 0
    
    for variant, filter_status, filter_list in records:
        # Extract variant fields based on type
        if use_cyvcf2:
            # cyvcf2.Variant
            chrom_str = variant.CHROM
            pos_str = str(variant.POS)
            id_str = variant.ID if variant.ID else '.'
            ref_str = variant.REF
            alt_str = ','.join(variant.ALT) if variant.ALT else '.'
            qual_str = str(variant.QUAL) if variant.QUAL is not None else '.'
            
            # Get original FILTER
            if variant.FILTER:
                filter_str = variant.FILTER
            else:
                filter_str = 'PASS'
            
            # Copy existing INFO fields from cyvcf2
            # cyvcf2 doesn't iterate well over INFO, so we parse the raw INFO string
            info_parts = []
            try:
                # Get raw INFO string from variant
                info_str_raw = str(variant).strip().split('\t')[7]
                if info_str_raw and info_str_raw != '.':
                    # Split and add each INFO field
                    for info_field in info_str_raw.split(';'):
                        if info_field:
                            info_parts.append(info_field)
            except Exception:
                pass
        else:
            # pysam.VariantRecord
            chrom_str = variant.contig
            pos_str = str(variant.pos)
            id_str = variant.id if variant.id else '.'
            ref_str = variant.ref
            alt_str = ','.join(variant.alts) if variant.alts else '.'
            qual_str = str(variant.qual) if variant.qual is not None else '.'
            
            # Get original FILTER
            if hasattr(variant, 'filter') and variant.filter and len(variant.filter) > 0:
                filter_str = ';'.join(variant.filter.keys())
            else:
                filter_str = 'PASS'
            
            # Copy existing INFO fields
            info_parts = []
            for key in variant.info:
                try:
                    value = variant.info[key]
                    if value is True:
                        info_parts.append(key)
                    elif value is not None and value != '':
                        value_str = str(value).replace(';', '%3B').replace('=', '%3D')
                        info_parts.append(f"{key}={value_str}")
                except Exception:
                    pass
        
        # Add RaVeX filter INFO
        if filter_status == "PASS" or not filter_list:
            info_parts.append('RaVeX_FILTER=PASS')
        else:
            # Use semicolon within the value, not as separator between INFO fields
            ravex_filter_value = ",".join(filter_list)  # Use comma to separate multiple filter reasons
            info_parts.append(f'RaVeX_FILTER={ravex_filter_value}')
            # Add individual filter flags as separate INFO fields
            for flag in set(filter_list):  # Use set to avoid duplicates
                info_parts.append(flag)
        
        info_str = ';'.join(info_parts) if info_parts else '.'
        
        # Write variant line (8 columns only - no FORMAT/sample)
        variant_line = f"{chrom_str}\t{pos_str}\t{id_str}\t{ref_str}\t{alt_str}\t{qual_str}\t{filter_str}\t{info_str}"
        variant_lines.append(variant_line)
        written_count += 1
    
    # Write output (BGZF compressed for tabix compatibility)
    with pysam.BGZFile(output_path, 'w') as f:
        # Write header
        for line in header_lines:
            f.write(line.encode('utf-8') + b'\n')
        # Write variants
        if variant_lines:
            for line in variant_lines:
                f.write(line.encode('utf-8') + b'\n')
    
    input_vcf.close()
    
    return written_count
