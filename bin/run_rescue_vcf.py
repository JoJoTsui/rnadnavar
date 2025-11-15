#!/usr/bin/env python3
"""
Cross-modality VCF rescue script.

This script performs cross-modality variant aggregation between DNA and RNA
consensus VCF files, identifying variants with support from both modalities.
"""

import argparse
import sys
from pathlib import Path


def argparser():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description='Cross-modality VCF rescue: aggregate variants from DNA and RNA consensus',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Input files
    parser.add_argument(
        '--dna_consensus',
        type=str,
        required=True,
        help='Path to DNA consensus VCF file'
    )
    parser.add_argument(
        '--rna_consensus',
        type=str,
        required=True,
        help='Path to RNA consensus VCF file'
    )
    parser.add_argument(
        '--dna_vcfs',
        type=str,
        required=True,
        help='Directory containing individual DNA caller VCF files'
    )
    parser.add_argument(
        '--rna_vcfs',
        type=str,
        required=True,
        help='Directory containing individual RNA caller VCF files'
    )
    
    # Output options
    parser.add_argument(
        '--out_prefix',
        type=str,
        required=True,
        help='Output file prefix (e.g., sample.rescued)'
    )
    parser.add_argument(
        '--output_format',
        type=str,
        default='vcf.gz',
        choices=['vcf', 'vcf.gz', 'bcf'],
        help='Output format'
    )
    
    # Consensus thresholds
    parser.add_argument(
        '--snv_thr',
        type=int,
        default=2,
        help='Minimum number of callers for SNV consensus'
    )
    parser.add_argument(
        '--indel_thr',
        type=int,
        default=2,
        help='Minimum number of callers for indel consensus'
    )
    
    return parser.parse_args()


def find_vcf_files(directory):
    """Find all VCF files in a directory."""
    vcf_files = {}
    directory = Path(directory)
    
    # Look for .vcf, .vcf.gz, and .bcf files
    for pattern in ['*.vcf', '*.vcf.gz', '*.bcf']:
        for vcf_path in directory.glob(pattern):
            # Extract caller name from filename
            from vcf_utils.io_utils import get_caller_name
            caller = get_caller_name(vcf_path.name)
            vcf_files[caller] = str(vcf_path)
    
    return vcf_files


def main():
    """Main rescue workflow."""
    args = argparser()
    
    print("=" * 80)
    print("Cross-Modality VCF Rescue")
    print("=" * 80)
    
    # Validate input files
    dna_consensus_path = Path(args.dna_consensus)
    rna_consensus_path = Path(args.rna_consensus)
    dna_vcfs_dir = Path(args.dna_vcfs)
    rna_vcfs_dir = Path(args.rna_vcfs)
    
    if not dna_consensus_path.exists():
        print(f"Error: DNA consensus VCF not found: {dna_consensus_path}", file=sys.stderr)
        sys.exit(1)
    if not rna_consensus_path.exists():
        print(f"Error: RNA consensus VCF not found: {rna_consensus_path}", file=sys.stderr)
        sys.exit(1)
    if not dna_vcfs_dir.is_dir():
        print(f"Error: DNA VCFs directory not found: {dna_vcfs_dir}", file=sys.stderr)
        sys.exit(1)
    if not rna_vcfs_dir.is_dir():
        print(f"Error: RNA VCFs directory not found: {rna_vcfs_dir}", file=sys.stderr)
        sys.exit(1)
    
    print(f"\nInput files:")
    print(f"  - DNA consensus: {dna_consensus_path}")
    print(f"  - RNA consensus: {rna_consensus_path}")
    print(f"  - DNA VCFs directory: {dna_vcfs_dir}")
    print(f"  - RNA VCFs directory: {rna_vcfs_dir}")
    print(f"\nOutput:")
    print(f"  - Prefix: {args.out_prefix}")
    print(f"  - Format: {args.output_format}")
    print(f"\nThresholds:")
    print(f"  - SNV: {args.snv_thr}")
    print(f"  - Indel: {args.indel_thr}")
    
    # Import vcf_utils functions
    from vcf_utils.aggregation import read_variants_from_vcf, aggregate_variants
    from cyvcf2 import VCF
    
    # Step 1: Read DNA consensus VCF with modality='DNA'
    print("\n" + "=" * 80)
    print("Reading DNA consensus VCF...")
    print("=" * 80)
    dna_consensus = read_variants_from_vcf(
        str(dna_consensus_path),
        'dna_consensus',
        modality='DNA'
    )
    print(f"  - Loaded {len(dna_consensus):,} variants from DNA consensus")
    
    # Step 2: Read RNA consensus VCF with modality='RNA'
    print("\n" + "=" * 80)
    print("Reading RNA consensus VCF...")
    print("=" * 80)
    rna_consensus = read_variants_from_vcf(
        str(rna_consensus_path),
        'rna_consensus',
        modality='RNA'
    )
    print(f"  - Loaded {len(rna_consensus):,} variants from RNA consensus")
    
    # Step 3: Read individual DNA caller VCFs with modality='DNA'
    print("\n" + "=" * 80)
    print("Reading individual DNA caller VCFs...")
    print("=" * 80)
    dna_vcf_files = find_vcf_files(dna_vcfs_dir)
    print(f"  - Found {len(dna_vcf_files)} DNA caller VCF files")
    
    dna_callers = []
    for caller, vcf_path in dna_vcf_files.items():
        print(f"  - Reading {caller} from {Path(vcf_path).name}")
        variants = read_variants_from_vcf(vcf_path, caller, modality='DNA')
        print(f"    - Loaded {len(variants):,} variants")
        dna_callers.append((caller, variants, 'DNA'))
    
    # Step 4: Read individual RNA caller VCFs with modality='RNA'
    print("\n" + "=" * 80)
    print("Reading individual RNA caller VCFs...")
    print("=" * 80)
    rna_vcf_files = find_vcf_files(rna_vcfs_dir)
    print(f"  - Found {len(rna_vcf_files)} RNA caller VCF files")
    
    rna_callers = []
    for caller, vcf_path in rna_vcf_files.items():
        print(f"  - Reading {caller} from {Path(vcf_path).name}")
        variants = read_variants_from_vcf(vcf_path, caller, modality='RNA')
        print(f"    - Loaded {len(variants):,} variants")
        rna_callers.append((caller, variants, 'RNA'))
    
    # Step 5: Aggregate all variants using aggregate_variants()
    print("\n" + "=" * 80)
    print("Aggregating variants across modalities...")
    print("=" * 80)
    
    # Combine all variant collections
    all_collections = [
        ('dna_consensus', dna_consensus, 'DNA'),
        ('rna_consensus', rna_consensus, 'RNA')
    ] + dna_callers + rna_callers
    
    print(f"  - Total variant sources: {len(all_collections)}")
    print(f"    - DNA sources: {1 + len(dna_callers)} (1 consensus + {len(dna_callers)} callers)")
    print(f"    - RNA sources: {1 + len(rna_callers)} (1 consensus + {len(rna_callers)} callers)")
    
    variant_data = aggregate_variants(
        all_collections,
        args.snv_thr,
        args.indel_thr
    )
    
    print(f"  - Aggregated {len(variant_data):,} unique variants")
    
    # Store template header and sample name for later use
    template_vcf = VCF(str(dna_consensus_path))
    template_header = template_vcf
    sample_name = template_vcf.samples[0] if template_vcf.samples else 'SAMPLE'
    
    # Step 6: Create modality_map from all input sources
    print("\n" + "=" * 80)
    print("Creating modality map...")
    print("=" * 80)
    
    from vcf_utils.tagging import tag_variant_with_modality, mark_rescued_variants
    
    modality_map = {}
    for caller_name, _, modality in all_collections:
        modality_map[caller_name] = modality
    
    print(f"  - Created modality map for {len(modality_map)} callers/sources")
    print(f"    - DNA: {sum(1 for m in modality_map.values() if m == 'DNA')}")
    print(f"    - RNA: {sum(1 for m in modality_map.values() if m == 'RNA')}")
    
    # Step 7: Tag variants with modality using tag_variant_with_modality()
    print("\n" + "=" * 80)
    print("Tagging variants with modality information...")
    print("=" * 80)
    
    for vkey, data in variant_data.items():
        tag_variant_with_modality(data, modality_map)
    
    print(f"  - Tagged {len(variant_data):,} variants with modality information")
    
    # Step 8: Mark rescued variants using mark_rescued_variants()
    print("\n" + "=" * 80)
    print("Marking rescued variants...")
    print("=" * 80)
    
    # Get variant keys from DNA and RNA consensus
    dna_variant_keys = set(dna_consensus.keys())
    rna_variant_keys = set(rna_consensus.keys())
    
    print(f"  - DNA consensus variants: {len(dna_variant_keys):,}")
    print(f"  - RNA consensus variants: {len(rna_variant_keys):,}")
    
    variant_data = mark_rescued_variants(variant_data, dna_variant_keys, rna_variant_keys)
    
    # Count rescued variants
    rescued_count = sum(1 for data in variant_data.values() if data.get('rescued', False))
    cross_modality_count = sum(1 for data in variant_data.values() if data.get('cross_modality', False))
    
    print(f"  - Rescued variants (cross-modality support): {rescued_count:,}")
    print(f"  - Cross-modality variants: {cross_modality_count:,}")
    
    # Step 9: Compute rescue statistics using compute_rescue_statistics()
    print("\n" + "=" * 80)
    print("Computing rescue statistics...")
    print("=" * 80)
    
    from vcf_utils.statistics import compute_rescue_statistics, print_statistics
    
    stats = compute_rescue_statistics(variant_data, dna_variant_keys, rna_variant_keys)
    
    # Step 10: Print statistics using print_statistics()
    print_statistics(stats, operation_type='rescue')
    
    # Step 11: Write output VCF with modality_map using write_vcf()
    print("\n" + "=" * 80)
    print("Writing rescued VCF...")
    print("=" * 80)
    
    from vcf_utils.io_utils import write_union_vcf
    
    # Determine output file path
    if args.output_format == 'vcf':
        out_file = f"{args.out_prefix}.vcf"
    elif args.output_format == 'vcf.gz':
        out_file = f"{args.out_prefix}.vcf.gz"
    elif args.output_format == 'bcf':
        out_file = f"{args.out_prefix}.bcf"
    else:
        out_file = f"{args.out_prefix}.vcf.gz"
    
    # Get all caller names for the output
    all_caller_names = [caller_name for caller_name, _, _ in all_collections]
    
    # Write VCF with modality_map to enable modality-specific INFO fields
    written_count = write_union_vcf(
        variant_data,
        template_header,
        sample_name,
        out_file,
        args.output_format,
        all_caller_names,
        modality_map=modality_map
    )
    
    print("\n" + "=" * 80)
    print("Rescue workflow completed successfully!")
    print("=" * 80)
    print(f"  - Output file: {out_file}")
    print(f"  - Total variants written: {written_count:,}")
    print(f"  - Rescued variants: {rescued_count:,}")
    print(f"  - Rescue rate: {stats['rescue_rate']:.1f}%")
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
