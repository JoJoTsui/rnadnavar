#!/usr/bin/env python3
"""
Cross-modality VCF rescue script.

This script performs cross-modality variant aggregation between DNA and RNA
consensus VCF files, identifying variants with support from both modalities.

Updated to use unified classification configuration for consistent thresholds.
"""

import argparse
import sys
from pathlib import Path

# Import unified configuration
from vcf_utils.classification_config import DEFAULT_THRESHOLDS, validate_thresholds


def argparser():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Cross-modality VCF rescue: aggregate variants from DNA and RNA consensus",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Input files
    parser.add_argument(
        "--dna_consensus",
        type=str,
        required=True,
        help="Path to DNA consensus VCF file",
    )
    parser.add_argument(
        "--rna_consensus",
        type=str,
        required=True,
        help="Path to RNA consensus VCF file",
    )
    parser.add_argument(
        "--dna_vcf",
        type=str,
        action="append",
        default=[],
        help="Path to individual DNA caller VCF file (can be specified multiple times)",
    )
    parser.add_argument(
        "--rna_vcf",
        type=str,
        action="append",
        default=[],
        help="Path to individual RNA caller VCF file (can be specified multiple times)",
    )

    # Output options
    parser.add_argument(
        "--out_prefix",
        type=str,
        required=True,
        help="Output file prefix (e.g., sample.rescued)",
    )
    parser.add_argument(
        "--output_format",
        type=str,
        default="vcf.gz",
        choices=["vcf", "vcf.gz", "bcf"],
        help="Output format",
    )

    # Consensus thresholds - use unified defaults
    parser.add_argument(
        "--snv_thr",
        type=int,
        default=DEFAULT_THRESHOLDS["consensus_snv_threshold"],
        help="Minimum number of callers for SNV consensus",
    )
    parser.add_argument(
        "--indel_thr",
        type=int,
        default=DEFAULT_THRESHOLDS["consensus_indel_threshold"],
        help="Minimum number of callers for indel consensus",
    )
    parser.add_argument(
        "--min_alt_support",
        type=int,
        default=DEFAULT_THRESHOLDS["consensus_min_alt_support"],
        help="Minimum tumor alt reads for a caller's record to count toward "
        "consensus support (0 disables the floor)",
    )

    # Rescue contract (audit M1/M2)
    parser.add_argument(
        "--disable_rescue_promotion",
        action="store_true",
        help="Disable cross-modality promotion: individual DNA+RNA callers "
        "agreeing on Somatic no longer rescue a site that failed "
        "within-modality consensus (legacy NoConsensus outcome)",
    )
    parser.add_argument(
        "--rescue_min_dna_callers",
        type=int,
        default=DEFAULT_THRESHOLDS["rescue_promotion_min_dna_callers"],
        help="Minimum DNA callers agreeing on Somatic for cross-modality promotion",
    )
    parser.add_argument(
        "--rescue_min_rna_callers",
        type=int,
        default=DEFAULT_THRESHOLDS["rescue_promotion_min_rna_callers"],
        help="Minimum RNA callers agreeing on Somatic for cross-modality promotion",
    )
    parser.add_argument(
        "--rescue_veto",
        type=str,
        choices=["dna", "rna", "none"],
        default=DEFAULT_THRESHOLDS["rescue_veto_direction"],
        help="Which modality's Artifact consensus label vetoes the other "
        "modality's non-Artifact evidence ('none' restores the legacy "
        "RNA-first override)",
    )

    # Rescue mode
    parser.add_argument(
        "--consensus_only",
        action="store_true",
        help="Only merge DNA and RNA consensus VCFs (no individual callers). "
        "This ensures final counts do not exceed DNA + RNA consensus counts.",
    )
    parser.add_argument(
        "--verification-json",
        type=str,
        default=None,
        help="Optional DNA tumor/normal verification manifest keyed by chrom:pos:ref:alt",
    )
    parser.add_argument(
        "--alignment-round",
        choices=["first", "realignment"],
        default="first",
        help="Alignment round represented by this rescue output",
    )

    # Chromosome filtering
    parser.add_argument(
        "--include-non-canonical",
        action="store_true",
        default=False,
        dest="include_non_canonical",
        help="Include non-canonical chromosomes (e.g., chrUn_*, *_random, etc.). "
        "By default, only canonical chromosomes (1-22, X, Y, M/MT) are included.",
    )

    return parser.parse_args()


def find_vcf_files(directory):
    """Find all VCF files in a directory."""
    vcf_files = {}
    directory = Path(directory)

    # Look for .vcf, .vcf.gz, and .bcf files
    for pattern in ["*.vcf", "*.vcf.gz", "*.bcf"]:
        for vcf_path in directory.glob(pattern):
            # Extract caller name from filename
            from vcf_utils.io_utils import get_caller_name

            caller = get_caller_name(vcf_path.name)
            vcf_files[caller] = str(vcf_path)

    return vcf_files


def main():
    """Main rescue workflow."""
    args = argparser()
    try:
        validate_thresholds(vars(args))
    except ValueError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(2)

    print("=" * 80)
    print("Cross-Modality VCF Rescue")
    print("=" * 80)

    # Validate input files
    dna_consensus_path = Path(args.dna_consensus)
    rna_consensus_path = Path(args.rna_consensus)

    if not dna_consensus_path.exists():
        print(
            f"Error: DNA consensus VCF not found: {dna_consensus_path}", file=sys.stderr
        )
        sys.exit(1)
    if not rna_consensus_path.exists():
        print(
            f"Error: RNA consensus VCF not found: {rna_consensus_path}", file=sys.stderr
        )
        sys.exit(1)

    verification = {}
    if args.verification_json:
        verification_path = Path(args.verification_json)
        if not verification_path.exists():
            print(f"Error: verification manifest not found: {verification_path}", file=sys.stderr)
            sys.exit(1)
        import json
        with verification_path.open() as fh:
            verification_doc = json.load(fh)
        allowed_statuses = {"confirmed", "rejected", "inconclusive"}
        for row in verification_doc.get("results", []):
            key = ":".join(str(row.get(k, "")) for k in ("chrom", "pos", "ref", "alt"))
            status = row.get("status", "inconclusive")
            if status not in allowed_statuses:
                print(f"Error: invalid verification status {status!r} for {key}", file=sys.stderr)
                sys.exit(2)
            if key in verification:
                print(f"Error: duplicate verification key: {key}", file=sys.stderr)
                sys.exit(2)
            verification[key] = status
        print(f"Loaded DNA verification outcomes: {len(verification)}")

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
            print("WARNING: Sample mismatch between DNA and RNA VCFs")
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
    print("\nRescue contract:")
    if args.disable_rescue_promotion:
        print("  - Cross-modality promotion: DISABLED (legacy NoConsensus)")
    else:
        print(
            f"  - Cross-modality promotion: ON (min DNA callers: "
            f"{args.rescue_min_dna_callers}, min RNA callers: "
            f"{args.rescue_min_rna_callers})"
        )
    print(f"  - Artifact veto direction: {args.rescue_veto}")
    print("\nChromosome Filtering:")
    if args.include_non_canonical:
        print("  - Including ALL chromosomes (non-canonical enabled)")
    else:
        print("  - Canonical chromosomes only (1-22, X, Y, M/MT)")

    # Import vcf_utils functions
    from cyvcf2 import VCF
    from vcf_utils.aggregation import aggregate_variants, read_variants_from_vcf
    from vcf_utils.io_utils import (
        get_caller_name,
        open_union_vcf,
        union_contig_order,
        write_union_vcf,
    )

    # Get template header and sample name early (before reading all variants)
    template_vcf = VCF(str(dna_consensus_path))
    template_header = template_vcf
    sample_name = template_vcf.samples[0] if template_vcf.samples else "SAMPLE"

    # Caller/modality bookkeeping is data-independent — build it once upfront.
    # Helper function to safely add modality prefix (prevents double-prefixing)
    def safe_prefix_caller(caller_name, modality):
        """Safely add modality prefix to caller name, preventing double-prefixing."""
        if caller_name.startswith(f"{modality}_"):
            return caller_name  # Already prefixed
        return f"{modality}_{caller_name}"

    # (prefixed_name, path) pairs; empty in consensus_only mode, matching the
    # legacy behavior of skipping individual-caller reads entirely
    dna_sources = []
    rna_sources = []
    if not args.consensus_only:
        for modality, paths, sources in (
            ("DNA", args.dna_vcf, dna_sources),
            ("RNA", args.rna_vcf, rna_sources),
        ):
            seen = {}
            for vcf_path in paths:
                if not Path(vcf_path).is_file():
                    print(f"Error: {modality} caller VCF not found: {vcf_path}", file=sys.stderr)
                    sys.exit(2)
                caller_name = safe_prefix_caller(get_caller_name(Path(vcf_path).name), modality)
                identity = caller_name.casefold()
                if identity in seen:
                    print(
                        f"Error: duplicate {modality} caller identity {caller_name}: "
                        f"{seen[identity]} and {vcf_path}", file=sys.stderr,
                    )
                    sys.exit(2)
                seen[identity] = vcf_path
                sources.append((caller_name, vcf_path))
    else:
        if args.dna_vcf:
            print(
                f"\n  - Skipping {len(args.dna_vcf)} DNA caller VCFs "
                "(consensus_only mode)"
            )
        if args.rna_vcf:
            print(
                f"\n  - Skipping {len(args.rna_vcf)} RNA caller VCFs "
                "(consensus_only mode)"
            )

    modality_map = {"DNA_consensus": "DNA", "RNA_consensus": "RNA"}
    for prefixed_name, _ in dna_sources:
        modality_map[prefixed_name] = "DNA"
    for prefixed_name, _ in rna_sources:
        modality_map[prefixed_name] = "RNA"

    # Collect all caller names (excluding consensus sources)
    all_caller_names = [name for name, _ in dna_sources + rna_sources]

    print(f"\n  - Total variant sources: {2 + len(dna_sources) + len(rna_sources)}")
    print(f"    - DNA sources: {1 + len(dna_sources)} (1 consensus + {len(dna_sources)} callers)")
    print(f"    - RNA sources: {1 + len(rna_sources)} (1 consensus + {len(rna_sources)} callers)")
    print(f"\n  - Individual callers to track: {len(all_caller_names)}")
    for caller in sorted(all_caller_names):
        print(f"    - {caller}")
    print("\n  Modality map:")
    for caller, modality in sorted(modality_map.items()):
        print(f"    - {caller}: {modality}")

    # Rescue contract config (audit M1/M2) forwarded to the rescue classifier
    rescue_config = {
        "rescue_promotion_enabled": not args.disable_rescue_promotion,
        "rescue_promotion_min_dna_callers": args.rescue_min_dna_callers,
        "rescue_promotion_min_rna_callers": args.rescue_min_rna_callers,
        "rescue_veto_direction": args.rescue_veto,
    }

    from vcf_utils.tagging import mark_rescued_variants, tag_variant_with_modality
    from vcf_utils.variant_statistics import compute_rescue_statistics, print_statistics

    # Per-chromosome streaming: read/aggregate/tag/mark/write one chromosome
    # at a time, then free it. Whole-genome materialization peaked at 20-40 GB
    # RSS on the cohort (OOM kills in the 78 GB pod cgroup). Output is
    # identical to whole-genome processing: write_union_vcf sorts each chunk
    # by (chrom, pos) and chunks are written in union_contig_order(), which
    # mirrors that sort; aggregation, tagging, and rescue marking are all
    # per-site operations.
    out_file = f"{args.out_prefix}.{args.output_format}"
    vcf_out = open_union_vcf(
        template_header,
        sample_name,
        out_file,
        args.output_format,
        modality_map=modality_map,
        include_non_canonical=args.include_non_canonical,
    )
    chroms = union_contig_order(
        template_header, include_non_canonical=args.include_non_canonical
    )

    total_stats = {
        "total_variants": 0,
        "dna_only": 0,
        "rna_only": 0,
        "cross_modality": 0,
        "rescued": 0,
        "snvs": 0,
        "indels": 0,
        "snvs_rescued": 0,
        "indels_rescued": 0,
        "dna_somatic_baseline": 0,
        "output_somatic": 0,
        "somatic_retained_from_dna": 0,
        "somatic_new_vs_dna": 0,
        "somatic_lost_from_dna": 0,
    }
    dna_variant_keys = set()
    rna_variant_keys = set()
    total_written = 0

    print("\n- Streaming per chromosome over " + ", ".join(chroms))
    import gc

    for chrom in chroms:
        dna_consensus = read_variants_from_vcf(
            str(dna_consensus_path),
            "DNA_consensus",
            modality="DNA",
            alignment_round="first",
            include_non_canonical=args.include_non_canonical,
            chrom=chrom,
        )
        rna_consensus = read_variants_from_vcf(
            str(rna_consensus_path),
            "RNA_consensus",
            modality="RNA",
            alignment_round=args.alignment_round,
            include_non_canonical=args.include_non_canonical,
            chrom=chrom,
        )

        all_collections = [
            ("DNA_consensus", dna_consensus, "DNA"),
            ("RNA_consensus", rna_consensus, "RNA"),
        ]
        n_records = len(dna_consensus) + len(rna_consensus)

        for prefixed_name, vcf_path in dna_sources:
            caller_name = get_caller_name(Path(vcf_path).name)
            variants = read_variants_from_vcf(
                vcf_path,
                caller_name,
                modality="DNA",
                alignment_round="first",
                include_non_canonical=args.include_non_canonical,
                chrom=chrom,
            )
            n_records += len(variants)
            all_collections.append((prefixed_name, variants, "DNA"))

        for prefixed_name, vcf_path in rna_sources:
            caller_name = get_caller_name(Path(vcf_path).name)
            variants = read_variants_from_vcf(
                vcf_path,
                caller_name,
                modality="RNA",
                alignment_round=args.alignment_round,
                include_non_canonical=args.include_non_canonical,
                chrom=chrom,
            )
            n_records += len(variants)
            all_collections.append((prefixed_name, variants, "RNA"))

        if n_records == 0:
            continue

        variant_data = aggregate_variants(
            all_collections,
            args.snv_thr,
            args.indel_thr,
            min_alt_support=args.min_alt_support,
        )

        for data in variant_data.values():
            data["alignment_round"] = args.alignment_round
            data["evidence_correlation"] = "correlated_reassessment" if args.alignment_round == "realignment" else "independent_input_round"

        if args.verification_json is not None:
            for data in variant_data.values():
                key = ":".join(str(data.get(k, "")) for k in ("CHROM", "POS", "REF", "ALT"))
                data["dna_verification_status"] = verification.get(key, "inconclusive")

        for data in variant_data.values():
            tag_variant_with_modality(data, modality_map)

        dna_variant_keys |= set(dna_consensus.keys())
        rna_variant_keys |= set(rna_consensus.keys())

        # Pass the record dicts (not just keys): rescue flags are computed
        # from records that PASSED as Somatic, never from mere presence in the
        # union file (audit M3)
        variant_data = mark_rescued_variants(variant_data, dna_consensus, rna_consensus)

        total_written += write_union_vcf(
            variant_data,
            template_header,
            sample_name,
            out_file,
            args.output_format,
            all_caller_names,
            modality_map=modality_map,
            snv_threshold=args.snv_thr,
            indel_threshold=args.indel_thr,
            include_non_canonical=args.include_non_canonical,
            rescue_config=rescue_config,
            vcf_out=vcf_out,
        )
        # The writer finalizes classification and promotion flags.
        chunk_stats = compute_rescue_statistics(
            variant_data, dna_consensus, rna_consensus
        )
        for key in total_stats:
            total_stats[key] += chunk_stats[key]
        print(
            f"  - {chrom}: {len(variant_data):,} unique variants "
            f"(running total {total_written:,})"
        )

        del all_collections, variant_data
        gc.collect()

    vcf_out.close()
    print(f"- Successfully wrote {total_written:,} variants to {out_file}")

    print(f"\n  - DNA consensus variants: {len(dna_variant_keys):,}")
    print(f"  - RNA consensus variants: {len(rna_variant_keys):,}")
    print(f"  - Rescued variants (cross-modality support): {total_stats['rescued']:,}")
    print(f"  - Cross-modality variants: {total_stats['cross_modality']:,}")

    # Recompute the fraction from summed counts, never average chromosome rates.
    total_stats["rescued_union_fraction"] = (
        total_stats["rescued"] / total_stats["total_variants"]
        if total_stats["total_variants"] else 0.0
    )

    print_statistics(total_stats, operation_type="rescue")

    print("\n" + "=" * 80)
    print("Rescue workflow completed successfully!")
    print("=" * 80)
    print(f"  - Output file: {out_file}")
    print(f"  - Total variants written: {total_written:,}")
    print(f"  - Rescued variants: {total_stats['rescued']:,}")
    print(f"  - Rescued / union records: {total_stats['rescued']:,}/{total_stats['total_variants']:,} "
          f"({total_stats['rescued_union_fraction']:.2%})")

    return 0


if __name__ == "__main__":
    sys.exit(main())
