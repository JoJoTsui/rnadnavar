#!/usr/bin/env python3
"""
Rescue VCF Filtering Script

Applies RaVeX filtering logic to rescue VCF files while preserving the original
FILTER field (Somatic, Germline, Reference, Artifact) from the rescue/consensus stage.

This script:
1. PRESERVES original FILTER field values (never overwrites them)
2. Adds RaVeX filter flags to INFO field
3. Generates dual outputs: standard + stripped (multiallelic-filtered, FORMAT-stripped)
4. Compatible with both rescue and consensus VCFs
"""

import argparse
import pysam
import sys
from pathlib import Path

# Add vcf_utils to path
script_dir = Path(__file__).parent
sys.path.insert(0, str(script_dir))

from vcf_utils.unified_filters import (
    apply_ravex_filters,
    add_filter_info_fields,
    add_classification_filters_to_header,
    is_multiallelic
)
from vcf_utils.io_utils import normalize_chromosome


def argparser():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description='Filter rescue/consensus VCF with RaVeX filtering logic',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        '-i', '--input',
        type=str,
        required=True,
        help='Input VCF file (rescue or consensus)'
    )
    parser.add_argument(
        '-o', '--output',
        type=str,
        required=True,
        help='Output filtered VCF file'
    )
    parser.add_argument(
        '--output_stripped',
        type=str,
        help='Output stripped VCF file (multiallelic-filtered, FORMAT-stripped). '
             'If not provided, will be <output>.stripped.vcf.gz'
    )
    parser.add_argument(
        '-g', '--gnomad_thr',
        type=float,
        default=0.0001,
        help='GnomAD threshold for variants'
    )
    parser.add_argument(
        '--whitelist',
        type=str,
        help='BED file with variants to keep (CHROM POS REF ALT)'
    )
    parser.add_argument(
        '--blacklist',
        type=str,
        help='BED file with regions to remove (CHROM START END)'
    )
    parser.add_argument(
        '--ref',
        type=str,
        help='FASTA reference file to extract context'
    )
    parser.add_argument(
        '--min_alt_reads',
        type=int,
        default=2,
        help='Minimum alt reads'
    )
    parser.add_argument(
        '--filter_multiallelic',
        action='store_true',
        help='Filter out multiallelic sites in stripped output'
    )
    
    return parser.parse_args()


def load_whitelist(whitelist_path):
    """Load whitelist variants from BED file."""
    whitelist_vars = set()
    if not whitelist_path:
        return whitelist_vars
    
    with open(whitelist_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                chrom = normalize_chromosome(parts[0])
                var_id = f"{chrom}:{parts[1]}:{parts[2]}:{parts[3]}"
                whitelist_vars.add(var_id)
    
    return whitelist_vars


def load_blacklist(blacklist_path):
    """Load blacklist regions from BED file."""
    blacklist_regions = []
    if not blacklist_path:
        return blacklist_regions
    
    with open(blacklist_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                chrom = normalize_chromosome(parts[0])
                blacklist_regions.append((chrom, int(parts[1]), int(parts[2])))
    
    return blacklist_regions


def write_filtered_vcf(input_vcf_path, output_path, args, genome, whitelist_vars, 
                       blacklist_regions, strip_format=False, filter_multiallelic=False):
    """
    Write filtered VCF with RaVeX filter tags in INFO field.
    
    CRITICAL: This function PRESERVES the original FILTER field and only adds
    INFO/RaVeX_FILTER tags to indicate which filters were applied.
    
    Args:
        input_vcf_path: Path to input VCF
        output_path: Path to output VCF
        args: Command-line arguments
        genome: pysam.FastaFile for reference (optional)
        whitelist_vars: Set of whitelisted variant IDs
        blacklist_regions: List of blacklist regions
        strip_format: If True, remove FORMAT column
        filter_multiallelic: If True, filter multiallelic sites
    """
    # Open input VCF
    input_vcf = pysam.VariantFile(input_vcf_path)
    
    # Create new header
    new_header = input_vcf.header.copy()
    
    # Strip FORMAT fields if requested
    if strip_format:
        # Remove all FORMAT definitions. Keep sample names in header; we will
        # not write any FORMAT/sample data in records, which results in no
        # FORMAT column being emitted in the output.
        for fmt_key in list(new_header.formats.keys()):
            new_header.formats.remove_header(fmt_key)
    
    # Add biological classification FILTER definitions
    new_header = add_classification_filters_to_header(new_header)
    
    # Add RaVeX filter INFO fields
    new_header = add_filter_info_fields(new_header)
    
    # Create output VCF
    output_vcf = pysam.VariantFile(output_path, 'w', header=new_header)
    
    # Get chromosome order for sorting
    chrom_order = {contig: idx for idx, contig in enumerate(new_header.contigs)}
    
    # Process variants
    processed_count = 0
    filtered_count = 0
    multiallelic_count = 0
    
    for record in input_vcf:
        # Check whitelist first
        chrom = normalize_chromosome(record.chrom)
        vkey = f"{chrom}:{record.pos}:{record.ref}:{record.alts[0] if record.alts else ''}"
        
        if whitelist_vars and vkey in whitelist_vars:
            # Whitelisted - always pass
            filters = []
        else:
            # Apply RaVeX filters
            filters = apply_ravex_filters(
                record, args, genome=genome, 
                blacklist_regions=blacklist_regions,
                use_cyvcf2=False
            )
        
        # Check multiallelic filtering for stripped output
        if filter_multiallelic and is_multiallelic(record.ref, record.alts):
            if strip_format:
                # Skip this variant in stripped output
                multiallelic_count += 1
                continue
            else:
                # Add to filter list but don't skip
                if "multiallelic" not in filters:
                    filters.append("multiallelic")
        
        # Create new record
        new_record = output_vcf.new_record(
            contig=record.contig,
            start=record.start,
            stop=record.stop,
            alleles=record.alleles,
            id=record.id,
            qual=record.qual
        )
        
        # Copy INFO fields
        for key in record.info:
            if key not in new_header.info:
                continue
            try:
                value = record.info[key]
                if value is not None and value != '':
                    new_record.info[key] = value
            except Exception:
                pass
        
        # Copy FORMAT and sample data (if not stripped)
        if not strip_format and record.format:
            for fmt_key in record.format.keys():
                if fmt_key not in new_header.formats:
                    continue
                try:
                    for sample_idx, sample in enumerate(record.samples):
                        value = record.samples[sample_idx][fmt_key]
                        if value is not None:
                            new_record.samples[sample_idx][fmt_key] = value
                except Exception:
                    pass
        
        # CRITICAL: Preserve original FILTER field
        # Copy all FILTER values from input record
        if hasattr(record, 'filter') and record.filter:
            for filt in record.filter:
                # Skip PASS - handled automatically
                if filt == 'PASS':
                    continue
                # Only add if defined in header
                if filt in new_header.filters:
                    new_record.filter.add(filt)
        
        # Set RaVeX filter INFO tags (NOT FILTER field)
        if not filters:
            new_record.info['RaVeX_FILTER'] = "PASS"
        else:
            new_record.info['RaVeX_FILTER'] = ";".join(filters)
            # Set individual filter flags
            for flag in filters:
                if flag in new_header.info:
                    new_record.info[flag] = True
            filtered_count += 1
        
        # Write variant
        output_vcf.write(new_record)
        processed_count += 1
        
        if processed_count % 10000 == 0:
            print(f"  Processed {processed_count:,} variants...")
    
    # Close files
    output_vcf.close()
    input_vcf.close()
    
    return processed_count, filtered_count, multiallelic_count


def main():
    """Main filtering workflow."""
    args = argparser()
    
    # Auto-generate stripped output path if not provided
    if not args.output_stripped:
        output_path = Path(args.output)
        args.output_stripped = str(output_path.parent / f"{output_path.stem}.stripped.vcf.gz")
    
    print("=" * 80)
    print("Rescue/Consensus VCF Filtering with RaVeX Rules")
    print("=" * 80)
    print(f"\nInput:             {args.input}")
    print(f"Output (standard): {args.output}")
    print(f"Output (stripped): {args.output_stripped}")
    print(f"gnomAD threshold:  {args.gnomad_thr}")
    print(f"Min alt reads:     {args.min_alt_reads}")
    print(f"Filter multiallelic: {args.filter_multiallelic}")
    
    # Open genome if provided
    genome = None
    if args.ref:
        genome = pysam.FastaFile(args.ref)
        print(f"Reference:         {args.ref}")
    
    # Load whitelist and blacklist
    whitelist_vars = load_whitelist(args.whitelist)
    blacklist_regions = load_blacklist(args.blacklist)
    
    if whitelist_vars:
        print(f"Whitelist:         {len(whitelist_vars)} variants")
    if blacklist_regions:
        print(f"Blacklist:         {len(blacklist_regions)} regions")
    
    print("\n" + "=" * 80)
    print("Processing Standard Output (all variants, FORMAT preserved)...")
    print("=" * 80)
    
    # Generate standard output
    proc_std, filt_std, _ = write_filtered_vcf(
        args.input, args.output, args, genome, 
        whitelist_vars, blacklist_regions,
        strip_format=False, filter_multiallelic=False
    )
    
    print(f"\n✓ Standard output complete:")
    print(f"  Processed:  {proc_std:,} variants")
    print(f"  Filtered:   {filt_std:,} variants")
    print(f"  Output:     {args.output}")
    
    print("\n" + "=" * 80)
    print("Processing Stripped Output (multiallelic-filtered, FORMAT-stripped)...")
    print("=" * 80)
    
    # Generate stripped output
    proc_strip, filt_strip, multi_skip = write_filtered_vcf(
        args.input, args.output_stripped, args, genome,
        whitelist_vars, blacklist_regions,
        strip_format=True, filter_multiallelic=args.filter_multiallelic
    )
    
    print(f"\n✓ Stripped output complete:")
    print(f"  Processed:       {proc_strip:,} variants")
    print(f"  Filtered:        {filt_strip:,} variants")
    if args.filter_multiallelic:
        print(f"  Multiallelic skipped: {multi_skip:,} variants")
    print(f"  Output:          {args.output_stripped}")
    
    # Index outputs
    print("\n" + "=" * 80)
    print("Indexing outputs...")
    print("=" * 80)
    
    import subprocess
    for vcf_path in [args.output, args.output_stripped]:
        try:
            subprocess.run(['tabix', '-p', 'vcf', vcf_path], 
                          check=True, capture_output=True)
            print(f"  ✓ Indexed: {vcf_path}")
        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            print(f"  ⚠ Warning: Failed to index {vcf_path}")
    
    if genome:
        genome.close()
    
    print("\n" + "=" * 80)
    print("Filtering completed successfully!")
    print("=" * 80)
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
