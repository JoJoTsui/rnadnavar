#!/usr/bin/env python3
"""
VCF Union/Consensus Script
Union mode: keeps ALL variants from all callers with aggregated information
Uses cyvcf2 for fast reading and pysam for writing

Updated to use unified classification configuration and chromosome filtering.
By default, only canonical chromosomes (1-22, X, Y, M) are included.
Use --include-non-canonical to include all chromosomes.
"""

import argparse
import sys
from pathlib import Path

# Import aggregation utilities
from vcf_utils.aggregation import aggregate_variants, read_variants_from_vcf

# Import chromosome utilities
from vcf_utils.chromosome_utils import get_canonical_chromosome_list

# Import unified configuration
from vcf_utils.classification_config import DEFAULT_THRESHOLDS, validate_thresholds

# Import I/O utilities
from vcf_utils.io_utils import get_caller_name, write_union_vcf

# Import statistics utilities
from vcf_utils.variant_statistics import compute_consensus_statistics, print_statistics


def argparser():
    parser = argparse.ArgumentParser(
        description="Union all variants from multiple VCF files with caller aggregation",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--input_dir", required=True, help="Directory containing input VCF files"
    )
    parser.add_argument(
        "--expected_callers",
        required=True,
        help="Comma-separated caller panel required for this sample/modality; "
        "missing, duplicate, and unexpected caller VCFs are fatal",
    )
    parser.add_argument("--out_prefix", required=True, help="Prefix for output files")
    parser.add_argument(
        "--snv_thr",
        type=int,
        default=DEFAULT_THRESHOLDS["consensus_snv_threshold"],
        help="Number of callers required for SNV consensus",
    )
    parser.add_argument(
        "--indel_thr",
        type=int,
        default=DEFAULT_THRESHOLDS["consensus_indel_threshold"],
        help="Number of callers required for indel consensus",
    )
    parser.add_argument(
        "--min_alt_support",
        type=int,
        default=DEFAULT_THRESHOLDS["consensus_min_alt_support"],
        help="Minimum tumor alt reads for a caller's record to count toward "
        "consensus support (0 disables the floor)",
    )
    parser.add_argument(
        "--preserve-baseline-callers",
        default="",
        help="Optional comma-separated caller names whose Somatic support is retained; opt-in only",
    )
    parser.add_argument(
        "--native-evidence-snv",
        action="store_true",
        help="Opt-in native-evidence policy for SNVs; indels remain threshold-based",
    )
    parser.add_argument(
        "--output_format",
        choices=["vcf", "vcf.gz", "bcf"],
        default="vcf.gz",
        help="Output format",
    )
    parser.add_argument(
        "--sample_name",
        help="Sample name (if not specified, use first sample from VCF)",
    )
    parser.add_argument(
        "--exclude_refcall",
        action="store_true",
        help="Exclude variants marked as RefCall (DeepSomatic reference calls)",
    )
    parser.add_argument(
        "--exclude_germline",
        action="store_true",
        help="Exclude variants marked as GERMLINE",
    )
    parser.add_argument(
        "--include-non-canonical",
        action="store_true",
        dest="include_non_canonical",
        help="Include non-canonical chromosomes (by default, only chr1-22, X, Y, M are included)",
    )
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
    try:
        validate_thresholds(vars(args))
    except ValueError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(2)

    # Find VCF files
    input_dir = Path(args.input_dir)
    vcf_files = {}
    expected_callers = [
        caller.strip().lower()
        for caller in args.expected_callers.split(",")
        if caller.strip()
    ]
    preserve_baseline_callers = {
        caller.strip().lower()
        for caller in args.preserve_baseline_callers.split(",")
        if caller.strip()
    }
    unknown_baseline = preserve_baseline_callers - set(expected_callers)
    if unknown_baseline:
        print(
            f"ERROR: --preserve-baseline-callers contains callers outside expected panel: {sorted(unknown_baseline)}",
            file=sys.stderr,
        )
        sys.exit(2)
    if not expected_callers:
        print("ERROR: --expected_callers must contain at least one caller", file=sys.stderr)
        sys.exit(2)
    if len(expected_callers) != len(set(expected_callers)):
        print(
            f"ERROR: duplicate callers in expected panel: {args.expected_callers}",
            file=sys.stderr,
        )
        sys.exit(2)

    print(f"Searching for VCF files in: {input_dir}")
    for vcf_path in sorted(input_dir.glob("*.vcf*")):
        if vcf_path.suffix in [".vcf", ".gz"] or vcf_path.name.endswith(".vcf.gz"):
            caller = get_caller_name(str(vcf_path)).lower()
            if caller in vcf_files:
                print(
                    f"ERROR: duplicate VCFs discovered for caller {caller}: "
                    f"{Path(vcf_files[caller]).name}, {vcf_path.name}",
                    file=sys.stderr,
                )
                sys.exit(2)
            vcf_files[caller] = str(vcf_path)
            print(f"  - Found: {caller} -> {vcf_path.name}")

    if not vcf_files:
        print(f"ERROR: No VCF files found in {input_dir}")
        sys.exit(1)

    actual_callers = set(vcf_files)
    expected_set = set(expected_callers)
    missing = sorted(expected_set - actual_callers)
    unexpected = sorted(actual_callers - expected_set)
    if missing or unexpected:
        print(
            "ERROR: caller panel is incomplete or inconsistent; "
            f"expected={expected_callers}, actual={sorted(actual_callers)}, "
            f"missing={missing}, unexpected={unexpected}",
            file=sys.stderr,
        )
        sys.exit(2)

    # Input validation
    if args.snv_thr <= 0 or args.indel_thr <= 0:
        print("ERROR: Consensus thresholds must be > 0")
        sys.exit(1)

    if len(expected_callers) < max(args.snv_thr, args.indel_thr):
        print(
            f"ERROR: expected caller panel has {len(expected_callers)} callers, "
            f"but thresholds require {max(args.snv_thr, args.indel_thr)}",
            file=sys.stderr,
        )
        sys.exit(2)

    print(
        f"\nProcessing complete caller panel: {', '.join(expected_callers)}"
    )
    print(f"SNV consensus threshold: {args.snv_thr}")
    print(f"Indel consensus threshold: {args.indel_thr}")
    print(f"Baseline preservation callers: {sorted(preserve_baseline_callers) or 'none (legacy default)'}")
    print(f"Native-evidence SNV policy: {'enabled' if args.native_evidence_snv else 'disabled'}")
    if args.exclude_refcall:
        print("Excluding RefCall variants")
    if args.exclude_germline:
        print("Excluding GERMLINE variants")

    # Chromosome filtering setting
    if args.include_non_canonical:
        print("Including ALL chromosomes (non-canonical included)")
    else:
        canonical_list = get_canonical_chromosome_list()
        print(
            f"Filtering to canonical chromosomes only: {', '.join(canonical_list[:5])}...X, Y, M"
        )
    print()

    # Per-chromosome streaming: read each caller's records for one chromosome,
    # aggregate, write the sorted chunk, then free it before the next
    # chromosome. Whole-genome materialization peaked at 19-43 GB RSS on the
    # cohort (OOM kills in the 78 GB pod cgroup); chunking bounds peak memory
    # at roughly chr1's share. Output is identical to whole-genome processing:
    # write_union_vcf sorts each chunk by (chrom, pos) and chunks are written
    # in union_contig_order(), which mirrors that sort.
    from cyvcf2 import VCF

    from vcf_utils.io_utils import open_union_vcf, union_contig_order

    all_callers = expected_callers

    # Template header and sample name from the first VCF.
    # Note (audit M8a, ticket 07): sample_name is no longer written to the
    # output — the output header carries no sample column because records
    # carry no FORMAT/sample data. It is still parsed here for API
    # compatibility with open_union_vcf/write_union_vcf and --sample_name.
    template_header = VCF(vcf_files[expected_callers[0]])
    sample_name = (
        template_header.samples[0] if template_header.samples else None
    )

    # Override sample name if provided
    if args.sample_name:
        sample_name = args.sample_name

    out_file = f"{args.out_prefix}.{args.output_format}"
    vcf_out = open_union_vcf(
        template_header,
        sample_name,
        out_file,
        args.output_format,
        modality_map=None,
        include_non_canonical=args.include_non_canonical,
    )
    chroms = union_contig_order(
        template_header, include_non_canonical=args.include_non_canonical
    )

    total_stats = {
        "total_variants": 0,
        "snvs": 0,
        "indels": 0,
        "snvs_consensus": 0,
        "indels_consensus": 0,
        "single_caller": 0,
        "multi_caller": 0,
    }
    total_written = 0

    print("- Streaming per chromosome over " + ", ".join(chroms))
    import gc

    for chrom in chroms:
        variant_collections = []
        n_records = 0
        for caller, vcf_path in vcf_files.items():
            variants = read_variants_from_vcf(
                vcf_path,
                caller,
                modality=None,
                exclude_refcall=args.exclude_refcall,
                exclude_germline=args.exclude_germline,
                include_non_canonical=args.include_non_canonical,
                chrom=chrom,
            )
            n_records += len(variants)
            variant_collections.append((caller, variants, None))

        if n_records == 0:
            continue

        variant_data = aggregate_variants(
            variant_collections,
            snv_threshold=args.snv_thr,
            indel_threshold=args.indel_thr,
            min_alt_support=args.min_alt_support,
            preserve_baseline_callers=preserve_baseline_callers,
            native_evidence_snv=args.native_evidence_snv,
        )

        chunk_stats = compute_consensus_statistics(
            variant_data, args.snv_thr, args.indel_thr
        )
        for key in total_stats:
            total_stats[key] += chunk_stats[key]

        total_written += write_union_vcf(
            variant_data,
            template_header,
            sample_name,
            out_file,
            args.output_format,
            all_callers,
            modality_map=None,
            snv_threshold=args.snv_thr,
            indel_threshold=args.indel_thr,
            include_non_canonical=args.include_non_canonical,
            vcf_out=vcf_out,
        )
        print(
            f"  - {chrom}: {len(variant_data):,} unique variants "
            f"(running total {total_written:,})"
        )

        del variant_collections, variant_data
        gc.collect()

    vcf_out.close()
    print(f"- Successfully wrote {total_written:,} variants to {out_file}")

    print_statistics(total_stats, operation_type="consensus")

    print(f"\n{'=' * 60}")
    print("DONE! Union VCF created successfully.")
    print(f"Output: {out_file}")
    if not args.include_non_canonical:
        print("Note: Only canonical chromosomes (1-22, X, Y, M) included")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    main()
