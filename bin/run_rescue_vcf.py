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
        '--dna_vcf',
        type=str,
        action='append',
        default=[],
        help='Path to individual DNA caller VCF file (can be specified multiple times)'
    )
    parser.add_argument(
        '--rna_vcf',
        type=str,
        action='append',
        default=[],
        help='Path to individual RNA caller VCF file (can be specified multiple times)'
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
    
    # Rescue mode
    parser.add_argument(
        '--consensus_only',
        action='store_true',
        help='Only merge DNA and RNA consensus VCFs (no individual callers). '
             'This ensures final counts do not exceed DNA + RNA consensus counts.'
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
    
    if not dna_consensus_path.exists():
        print(f"Error: DNA consensus VCF not found: {dna_consensus_path}", file=sys.stderr)
        sys.exit(1)
    if not rna_consensus_path.exists():
        print(f"Error: RNA consensus VCF not found: {rna_consensus_path}", file=sys.stderr)
        sys.exit(1)
    
    # Input validation
    if args.snv_thr <= 0 or args.indel_thr <= 0:
        print("ERROR: Consensus thresholds must be > 0", file=sys.stderr)
        sys.exit(1)
    
    # Validate sample consistency between DNA and RNA VCFs
    try:
        from cyvcf2 import VCF
        dna_vcf_temp = VCF(str(dna_consensus_path))
        rna_vcf_temp = VCF(str(rna_consensus_path))
        
        dna_samples = set(dna_vcf_temp.samples) if dna_vcf_temp.samples else set()
        rna_samples = set(rna_vcf_temp.samples) if rna_vcf_temp.samples else set()
        
        if dna_samples and rna_samples and dna_samples != rna_samples:
            print(f"WARNING: Sample mismatch between DNA and RNA VCFs")
            print(f"  DNA samples: {dna_samples}")
            print(f"  RNA samples: {rna_samples}")
    except Exception as e:
        print(f"WARNING: Could not validate sample consistency: {e}")
    
    print("\nInput files:")
    print(f"  - DNA consensus: {dna_consensus_path}")
    print(f"  - RNA consensus: {rna_consensus_path}")
    print(f"  - DNA caller VCFs: {len(args.dna_vcf)}")
    for vcf in args.dna_vcf:
        print(f"    - {vcf}")
    print(f"  - RNA caller VCFs: {len(args.rna_vcf)}")
    for vcf in args.rna_vcf:
        print(f"    - {vcf}")
    print("\nOutput:")
    print(f"  - Prefix: {args.out_prefix}")
    print(f"  - Format: {args.output_format}")
    print("\nMode:")
    if args.consensus_only:
        print("  - Consensus-only merge (DNA + RNA consensus VCFs only)")
        print("  - Final counts will NOT exceed DNA + RNA consensus sums")
    else:
        print("  - Full cross-modality rescue (consensus + individual callers)")
        print("  - May rescue variants that failed consensus in one modality")
    print("\nThresholds:")
    print(f"  - SNV: {args.snv_thr}")
    print(f"  - Indel: {args.indel_thr}")
    
    # Import vcf_utils functions
    from vcf_utils.aggregation import read_variants_from_vcf, aggregate_variants
    from vcf_utils.io_utils import get_caller_name
    from cyvcf2 import VCF
    
    # Get template header and sample name early (before reading all variants)
    template_vcf = VCF(str(dna_consensus_path))
    template_header = template_vcf
    sample_name = template_vcf.samples[0] if template_vcf.samples else 'SAMPLE'
    
    # Step 1: Read DNA consensus VCF with modality='DNA'
    print("\n" + "=" * 80)
    print("Reading DNA consensus VCF...")
    print("=" * 80)
    dna_consensus = read_variants_from_vcf(
        str(dna_consensus_path),
        'DNA_consensus',
        modality='DNA'
    )
    print(f"  - Loaded {len(dna_consensus):,} variants from DNA consensus")
    
    # Step 2: Read RNA consensus VCF with modality='RNA'
    print("\n" + "=" * 80)
    print("Reading RNA consensus VCF...")
    print("=" * 80)
    rna_consensus = read_variants_from_vcf(
        str(rna_consensus_path),
        'RNA_consensus',
        modality='RNA'
    )
    print(f"  - Loaded {len(rna_consensus):,} variants from RNA consensus")
    
    # Helper function to read caller VCFs
    def read_caller_vcfs(vcf_paths, modality_name):
        """Read individual caller VCFs for a modality."""
        caller_variants = {}
        if vcf_paths:
            print(f"\n" + "=" * 80)
            print(f"Reading individual {modality_name} caller VCFs...")
            print("=" * 80)
            for vcf_path in vcf_paths:
                caller_name = get_caller_name(Path(vcf_path).name)
                print(f"  - Reading {caller_name} from {vcf_path}")
                variants = read_variants_from_vcf(vcf_path, caller_name, modality=modality_name)
                caller_variants[caller_name] = variants
                print(f"    - Loaded {len(variants):,} variants")
        return caller_variants
    
    # Step 3: Read individual DNA caller VCFs (skip if consensus_only mode)
    dna_caller_variants = {} if args.consensus_only else read_caller_vcfs(args.dna_vcf, 'DNA')
    if args.consensus_only and args.dna_vcf:
        print(f"\n  - Skipping {len(args.dna_vcf)} DNA caller VCFs (consensus_only mode)")
    
    # Step 4: Read individual RNA caller VCFs (skip if consensus_only mode)
    rna_caller_variants = {} if args.consensus_only else read_caller_vcfs(args.rna_vcf, 'RNA')
    if args.consensus_only and args.rna_vcf:
        print(f"  - Skipping {len(args.rna_vcf)} RNA caller VCFs (consensus_only mode)")
    
    # Step 5: Aggregate all variants (consensus + individual callers)
    print("\n" + "=" * 80)
    print("Aggregating all variants across modalities...")
    print("=" * 80)
    
    # Combine consensus and individual caller variant collections
    # Use modality prefix as KEY to avoid collisions in modality_map
    # but store the original caller name for later use
    all_collections = [
        ('DNA_consensus', dna_consensus, 'DNA'),
        ('RNA_consensus', rna_consensus, 'RNA')
    ]
    
    # Helper function to safely add modality prefix (prevents double-prefixing)
    def safe_prefix_caller(caller_name, modality):
        """Safely add modality prefix to caller name, preventing double-prefixing."""
        if caller_name.startswith(f'{modality}_'):
            return caller_name  # Already prefixed
        return f'{modality}_{caller_name}'
    
    # Add individual DNA callers WITH safe prefix as key to avoid modality_map collisions
    for caller_name, variants in dna_caller_variants.items():
        prefixed_name = safe_prefix_caller(caller_name, 'DNA')
        all_collections.append((prefixed_name, variants, 'DNA'))
    
    # Add individual RNA callers WITH safe prefix as key to avoid modality_map collisions
    for caller_name, variants in rna_caller_variants.items():
        prefixed_name = safe_prefix_caller(caller_name, 'RNA')
        all_collections.append((prefixed_name, variants, 'RNA'))
    
    print(f"  - Total variant sources: {len(all_collections)}")
    print(f"    - DNA sources: {1 + len(dna_caller_variants)} (1 consensus + {len(dna_caller_variants)} callers)")
    print(f"    - RNA sources: {1 + len(rna_caller_variants)} (1 consensus + {len(rna_caller_variants)} callers)")
    
    variant_data = aggregate_variants(
        all_collections,
        args.snv_thr,
        args.indel_thr
    )
    
    print(f"  - Aggregated {len(variant_data):,} unique variants")
    
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
    
    # Print modality map for debugging
    print("\n  Modality map:")
    for caller, modality in sorted(modality_map.items()):
        print(f"    - {caller}: {modality}")
    
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
    out_file = f"{args.out_prefix}.{args.output_format}"
    
    # Collect all caller names (excluding consensus sources)
    all_caller_names = [name for name, _, _ in all_collections if not name.endswith('_consensus')]
    
    print(f"\n  - Individual callers to track: {len(all_caller_names)}")
    for caller in sorted(all_caller_names):
        print(f"    - {caller}")
    
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
