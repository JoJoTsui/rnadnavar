#!/usr/bin/env python3
"""
VCF Union/Consensus Script
Union mode: keeps ALL variants from all callers with aggregated information
Uses cyvcf2 for fast reading and pysam for writing
"""
import argparse
import sys
from pathlib import Path

# Import I/O utilities
from vcf_utils.io_utils import get_caller_name, write_union_vcf
# Import aggregation utilities
from vcf_utils.aggregation import read_variants_from_vcf, aggregate_variants
# Import statistics utilities
from vcf_utils.statistics import compute_consensus_statistics, print_statistics


def argparser():
    parser = argparse.ArgumentParser(
        description="Union all variants from multiple VCF files with caller aggregation"
    )
    parser.add_argument("--input_dir", required=True, help="Directory containing input VCF files")
    parser.add_argument("--out_prefix", required=True, help="Prefix for output files")
    parser.add_argument("--snv_thr", type=int, default=2, help="Number of callers required for SNV consensus (default: 2)")
    parser.add_argument("--indel_thr", type=int, default=2, help="Number of callers required for indel consensus (default: 2)")
    parser.add_argument("--output_format", choices=['vcf', 'vcf.gz', 'bcf'], default='vcf.gz', 
                       help="Output format (default: vcf.gz)")
    parser.add_argument("--sample_name", help="Sample name (if not specified, use first sample from VCF)")
    parser.add_argument("--exclude_refcall", action='store_true', 
                       help="Exclude variants marked as RefCall (DeepSomatic reference calls)")
    parser.add_argument("--exclude_germline", action='store_true', 
                       help="Exclude variants marked as GERMLINE")
    return parser.parse_args()


# Functions moved to vcf_utils modules:
# - get_caller_name() -> vcf_utils.io_utils
# - normalize_chromosome() -> vcf_utils.io_utils
# - variant_key() -> vcf_utils.io_utils
# - is_snv() -> vcf_utils.io_utils
# - extract_genotype_info() -> vcf_utils.aggregation
# - aggregate_genotypes() -> vcf_utils.aggregation
# - read_variants_from_vcf() -> vcf_utils.aggregation
# - aggregate_variants() -> vcf_utils.aggregation
# - compute_consensus_statistics() -> vcf_utils.statistics
# - create_output_header() -> vcf_utils.io_utils
# - write_union_vcf() -> vcf_utils.io_utils


def main():
    args = argparser()
    
    # Find VCF files
    input_dir = Path(args.input_dir)
    vcf_files = {}
    
    print(f"Searching for VCF files in: {input_dir}")
    for vcf_path in sorted(input_dir.glob("*.vcf*")):
        if vcf_path.suffix in ['.vcf', '.gz'] or vcf_path.name.endswith('.vcf.gz'):
            caller = get_caller_name(str(vcf_path))
            vcf_files[caller] = str(vcf_path)
            print(f"  - Found: {caller} -> {vcf_path.name}")
    
    if not vcf_files:
        print(f"ERROR: No VCF files found in {input_dir}")
        sys.exit(1)
    
    # Input validation
    if args.snv_thr <= 0 or args.indel_thr <= 0:
        print("ERROR: Consensus thresholds must be > 0")
        sys.exit(1)
    
    if len(vcf_files) < max(args.snv_thr, args.indel_thr):
        print(f"WARNING: Only {len(vcf_files)} VCF files found, but thresholds require {max(args.snv_thr, args.indel_thr)}")
        print("Some variants may not meet consensus thresholds")
    
    print(f"\nProcessing {len(vcf_files)} VCF files with callers: {', '.join(vcf_files.keys())}")
    print(f"SNV consensus threshold: {args.snv_thr}")
    print(f"Indel consensus threshold: {args.indel_thr}")
    if args.exclude_refcall:
        print("Excluding RefCall variants")
    if args.exclude_germline:
        print("Excluding GERMLINE variants")
    print()
    
    # Read variants from each VCF file using vcf_utils
    print("- Reading variants from VCF files")
    variant_collections = []
    template_header = None
    sample_name = None
    all_callers = list(vcf_files.keys())
    
    from cyvcf2 import VCF
    
    for caller, vcf_path in vcf_files.items():
        print(f"  - Reading {caller}: {vcf_path}")
        
        # Get template header and sample name from first VCF
        if template_header is None:
            vcf = VCF(vcf_path)
            template_header = vcf
            if vcf.samples:
                sample_name = vcf.samples[0]
        
        # Read variants using vcf_utils
        variants = read_variants_from_vcf(
            vcf_path, 
            caller, 
            modality=None,
            exclude_refcall=args.exclude_refcall,
            exclude_germline=args.exclude_germline
        )
        
        print(f"    - Read {len(variants):,} variants from {caller}")
        variant_collections.append((caller, variants, None))
    
    # Override sample name if provided
    if args.sample_name:
        sample_name = args.sample_name
    
    # Aggregate variants using vcf_utils
    print("\n- Aggregating variants across callers")
    variant_data = aggregate_variants(
        variant_collections,
        snv_threshold=args.snv_thr,
        indel_threshold=args.indel_thr
    )
    
    print(f"- Total unique variants: {len(variant_data):,}")
    
    # Compute statistics using vcf_utils
    stats = compute_consensus_statistics(variant_data, args.snv_thr, args.indel_thr)
    print_statistics(stats, operation_type='consensus')
    
    # Write output VCF using vcf_utils
    out_file = f"{args.out_prefix}.{args.output_format}"
    write_union_vcf(
        variant_data, 
        template_header, 
        sample_name, 
        out_file, 
        args.output_format, 
        all_callers,
        modality_map=None
    )
    
    print(f"\n{'='*60}")
    print("DONE! Union VCF created successfully.")
    print(f"Output: {out_file}")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
