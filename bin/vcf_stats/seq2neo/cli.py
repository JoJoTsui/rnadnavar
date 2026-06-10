#!/usr/bin/env python3
"""
CLI entry point for seq2neo variant statistics.

Usage:
    .venv/bin/python -m vcf_stats.seq2neo.cli \\
        --manifest data/processed/sample_manifest.parquet \\
        --output-dir stats/variant/ \\
        --threads 8
"""

import argparse
import os
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import polars as pl

from .bam_stats import compute_all_bam_stats
from .bam_validation import validate_bam_all_samples as validate_bam
from .caller_parser import join_caller_columns, parse_all_callers
from .manifest_loader import filter_complete, load_manifest
from .rescue_parser import parse_rescue_vcf as _py_parse_rescue
from .rust_vcf import parse_rescue_vcf as _rust_parse_rescue
from .rescue_validator import validate_all_samples, validation_summary
from .rust_bam import pileup_variants
from .tiering_stats import compute_tiers_for_dataframe, tier_summary as compute_tier_summary
from .statistics import (
    compute_all_per_variant,
    dataset_summary,
    disease_summary,
    flag_filter_breakdown,
    gt_concordance,
    sample_summary,
    sample_tier_summary,
    set_summary,
)
from .visualizer import (
    generate_dashboard,
    plot_bam_metrics_bars,
    plot_bam_coverage_violin,
    plot_caller_overlap,
    plot_cosmic_gnomad_annotation,
    plot_cross_modality,
    plot_dna_vs_rna_dp,
    plot_dna_vs_rna_vaf,
    plot_dp_boxplot_per_tier,
    plot_gt_concordance,
    plot_gt_concordance_per_tier,
    plot_per_sample_distribution,
    plot_per_sample_tier_distribution,
    plot_per_tier_cross_sample_vaf,
    plot_ref_alt_dp_scatter,
    plot_tiered_caller_overlap,
    plot_caller_agreement_matrix,
    plot_chromosome_density,
    plot_dna_vs_rna_per_caller,
    plot_filter_distribution,
    plot_per_tier_vaf_boxplot,
    plot_tier_quality_distribution,
    plot_tiered_variant_types,
    plot_ti_tv_ratio,
    plot_redi_evidence,
    plot_vaf_distribution,
    plot_vaf_boxplot_per_tier,
    plot_validation_heatmap,
    plot_variant_type_distribution,
    plot_vc_distribution,
)


def process_single_sample(row: dict, max_workers: int = 1, use_rust: bool = True) -> dict:
    """Process one sample: parse rescue VCF + all caller VCFs + compute stats.

    Returns a dict with 'sample_id', 'df', and 'stats'.
    """
    sample_id = row["sample_id"]
    rescue_path = row["rescue_vcf_path"]
    base_dir = row["base_output_dir"]
    dir_name = row["dir_name"]
    vcf_prefix = row["vcf_prefix"]

    parse_fn = _rust_parse_rescue if use_rust else _py_parse_rescue
    print(f"  [{sample_id}] Parsing rescue VCF ({'rust' if use_rust else 'python'})...")
    rescue_df = parse_fn(rescue_path)
    if rescue_df.is_empty():
        print(f"  [{sample_id}] WARNING: No variants in rescue VCF")
        return {"sample_id": sample_id, "df": None, "stats": None}

    # Build target positions from rescue VCF as (CHROM, POS, REF, ALT) 4-tuples.
    # Using all 4 columns ensures correct matching at multiallelic sites
    # when joining caller FORMAT data. Normalized caller VCFs guarantee
    # consistent REF/ALT representation.
    chroms = rescue_df["CHROM"].to_list()
    poss = rescue_df["POS"].to_list()
    refs = rescue_df["REF"].to_list()
    alts = rescue_df["ALT"].to_list()
    target_positions = set(zip(chroms, poss, refs, alts))

    print(f"  [{sample_id}] Found {len(target_positions)} positions, parsing 6 caller VCFs (max_workers={max_workers})...")
    caller_data = parse_all_callers(base_dir, dir_name, vcf_prefix, target_positions, max_workers=max_workers)

    # Join caller columns onto rescue dataframe
    df = join_caller_columns(rescue_df, caller_data)

    # Add sample metadata
    df = df.with_columns([
        pl.lit(sample_id).alias("sample_id"),
        pl.lit(row["set_number"]).cast(pl.Int64).alias("set_number"),
        pl.lit(row["disease"]).alias("disease"),
        pl.lit(row["disease_normalized"]).alias("disease_normalized"),
    ])

    # Compute per-variant statistics (VAF, means)
    df = compute_all_per_variant(df)

    # Compute CxDy tiers via tiering engine
    df = compute_tiers_for_dataframe(df)

    # Compute sample-level summary
    stats = sample_summary(df, sample_id)
    stats["set_number"] = row["set_number"]
    stats["disease"] = row["disease"]

    print(f"  [{sample_id}] Done: {len(df)} variants")
    return {"sample_id": sample_id, "df": df, "stats": stats}


def main():
    parser = argparse.ArgumentParser(description="Seq2neo variant statistics")
    parser.add_argument("--manifest", required=True, help="Path to sample manifest CSV/Parquet")
    parser.add_argument("--output-dir", required=True, help="Output directory for statistics")
    parser.add_argument("--threads", type=int, default=6,
                        help="Threads for within-sample caller parsing (1-6, default: 6).")
    parser.add_argument("--sample-workers", type=int, default=1,
                        help="Parallel sample processing threads (default: 1). "
                             "Uses ThreadPoolExecutor (safe with htslib). Set to 4-8 for 32-core machines.")
    parser.add_argument("--max-samples", type=int, default=None, help="Limit number of samples")
    parser.add_argument("--set", type=int, default=None, help="Process only this set (1-4)")
    parser.add_argument("--sample-ids", nargs="*", default=None, help="Process specific sample IDs")
    parser.add_argument("--no-validate", action="store_true", help="Skip rescue VCF validation")
    parser.add_argument("--tolerance", type=float, default=0.01, help="Validation tolerance")
    parser.add_argument("--pileup-mode", choices=["all", "filtered"], default="all",
                        help="Variant-wise BAM pileup mode: all variants (default) or exclude NoConsensus")
    parser.add_argument("--no-bam", action="store_true", help="Skip all BAM processing")
    parser.add_argument("--no-pileup", action="store_true", help="Skip variant-wise BAM pileup (whole-genome only)")
    parser.add_argument("--bam-workers", type=int, default=8,
                        help="Threads for parallel BAM stats processing (default: 8). "
                             "Uses ThreadPoolExecutor (safe with htslib).")
    parser.add_argument("--parser", choices=["rust", "python"], default="rust",
                        help="VCF parser: rust (default) or python (cyvcf2 fallback)")
    args = parser.parse_args()

    # Load and filter manifest
    manifest = load_manifest(args.manifest)
    manifest = filter_complete(manifest)

    if args.set:
        manifest = manifest.filter(pl.col("set_number") == args.set)
    if args.sample_ids:
        manifest = manifest.filter(pl.col("sample_id").is_in(args.sample_ids))
    if args.max_samples:
        manifest = manifest.head(args.max_samples)

    if manifest.is_empty():
        print("No samples to process.")
        sys.exit(0)

    use_rust = args.parser == "rust"
    print(f"Processing {len(manifest)} samples (parser={args.parser}, caller_threads={args.threads}, sample_workers={args.sample_workers}, bam_workers={args.bam_workers})")
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    rows = manifest.to_dicts()
    samples_data = {}
    all_stats = []
    bam_stats_df = pl.DataFrame()  # populated below if --no-bam not set

    if args.sample_workers > 1:
        # Parallel sample processing via ThreadPoolExecutor (threads, safe with htslib)
        with ThreadPoolExecutor(max_workers=args.sample_workers) as executor:
            futures = {}
            for i, row in enumerate(rows):
                future = executor.submit(process_single_sample, row, args.threads, use_rust)
                futures[future] = (i, row["sample_id"])

            for future in as_completed(futures):
                i, sid = futures[future]
                print(f"[{i+1}/{len(manifest)}] {sid} - processing...")
                try:
                    result = future.result()
                    if result["df"] is not None:
                        samples_data[result["sample_id"]] = result["df"]
                    if result["stats"] is not None:
                        all_stats.append(result["stats"])
                    print(f"[{i+1}/{len(manifest)}] {sid} - Done")
                except Exception as e:
                    import traceback
                    print(f"[{i+1}/{len(manifest)}] {sid} - ERROR: {e}")
                    traceback.print_exc()
    else:
        # Sequential processing
        for i, row in enumerate(rows):
            sid = row["sample_id"]
            print(f"\n[{i+1}/{len(manifest)}] {sid}")
            try:
                result = process_single_sample(row, max_workers=args.threads, use_rust=use_rust)
                if result["df"] is not None:
                    samples_data[result["sample_id"]] = result["df"]
                if result["stats"] is not None:
                    all_stats.append(result["stats"])
            except Exception as e:
                import traceback
                print(f"  [{sid}] ERROR: {e}")
                traceback.print_exc()

    if not samples_data:
        print("No data processed.")
        sys.exit(1)

    # Combine all per-variant DataFrames
    combined_df = pl.concat(list(samples_data.values()), how="diagonal_relaxed")
    variant_details_path = output_dir / "variant_details.parquet"
    combined_df.write_parquet(str(variant_details_path))
    print(f"Variant details: {variant_details_path} ({len(combined_df)} variants)")

    # BAM statistics (per-sample per-modality)
    if not args.no_bam:
        print("Computing per-sample BAM statistics...")
        bam_stats_df = compute_all_bam_stats(rows, max_workers=args.bam_workers)
        if not bam_stats_df.is_empty():
            bam_stats_df.write_csv(str(output_dir / "bam_stats.csv"))
            print(f"BAM stats: {output_dir / 'bam_stats.csv'}")
    else:
        print("Skipping BAM statistics (--no-bam)")

    # Sample summary
    if all_stats:
        sample_stats_df = pl.DataFrame(all_stats)
        sample_stats_df.write_csv(str(output_dir / "sample_summary.csv"))
        print(f"Sample summary: {output_dir / 'sample_summary.csv'}")

        # Set summary
        set_summary_df = set_summary(sample_stats_df)
        if not set_summary_df.is_empty():
            set_summary_df.write_csv(str(output_dir / "set_summary.csv"))

        # Disease summary
        disease_summary_df = disease_summary(combined_df)
        if not disease_summary_df.is_empty():
            disease_summary_df.write_csv(str(output_dir / "disease_summary.csv"))

        # Tier summary
        tier_summary_df = compute_tier_summary(combined_df)
        if not tier_summary_df.is_empty():
            tier_summary_df.write_csv(str(output_dir / "tier_summary.csv"))
            print(f"Tier summary: {output_dir / 'tier_summary.csv'}")

        # Dataset summary (whole-dataset aggregates)
        ds_summary = dataset_summary(combined_df)
        if ds_summary:
            pl.DataFrame([ds_summary]).write_csv(str(output_dir / "dataset_summary.csv"))
            print(f"Dataset summary: {output_dir / 'dataset_summary.csv'}")

        # Sample-tier summary (Level 4)
        sample_tier_df = sample_tier_summary(combined_df)
        if not sample_tier_df.is_empty():
            sample_tier_df.write_csv(str(output_dir / "sample_tier_summary.csv"))
            print(f"Sample-tier summary: {output_dir / 'sample_tier_summary.csv'}")

    # Caller overlap
    from .statistics import caller_overlap_distribution
    overlap_df = caller_overlap_distribution(combined_df)
    if not overlap_df.is_empty():
        overlap_df.write_csv(str(output_dir / "caller_overlap.csv"))

    # Filter distribution
    from .statistics import filter_distribution, variant_type_distribution, vc_distribution
    filter_df = filter_distribution(combined_df)
    if not filter_df.is_empty():
        filter_df.write_csv(str(output_dir / "filter_distribution.csv"))

    # Variant type distribution
    vt_df = variant_type_distribution(combined_df)
    if not vt_df.is_empty():
        vt_df.write_csv(str(output_dir / "variant_type_distribution.csv"))

    # GT concordance
    concordance_data = gt_concordance(combined_df)
    if concordance_data:
        pl.DataFrame({"agreement_level": list(concordance_data.keys()), "count": list(concordance_data.values())}).write_csv(
            str(output_dir / "gt_concordance.csv")
        )

    # Validation
    if not args.no_validate:
        print("Running rescue VCF validation...")
        report = validate_all_samples(samples_data, args.tolerance)
        if not report.is_empty():
            report.write_csv(str(output_dir / "rescue_validation_report.csv"))
            summary = validation_summary(report)
            if not summary.is_empty():
                summary.write_csv(str(output_dir / "rescue_validation_summary.csv"))
                print(f"  Validation report: {output_dir / 'rescue_validation_report.csv'}")

        # BAM validation
        print("Running BAM validation...")
        bam_report = validate_bam(samples_data)
        if not bam_report.is_empty():
            bam_report.write_csv(str(output_dir / "bam_validation.csv"))
            print(f"  BAM validation: {output_dir / 'bam_validation.csv'}")
    else:
        report = None

    # Visualizations
    print("Generating visualizations...")
    figs = []

    # Add set_number to combined_df if available from manifest
    if "set_number" not in combined_df.columns:
        set_map = {r["sample_id"]: r["set_number"] for r in rows}
        combined_df = combined_df.with_columns(
            pl.col("sample_id").replace_strict(set_map, default=None).alias("set_number")
        )

    figs.append(plot_cosmic_gnomad_annotation(combined_df, str(output_dir)))
    figs.append(plot_gt_concordance(combined_df, str(output_dir)))
    figs.append(plot_vc_distribution(combined_df, str(output_dir)))
    figs.append(plot_caller_overlap(combined_df, str(output_dir)))
    figs.append(plot_vaf_distribution(combined_df, str(output_dir)))
    figs.append(plot_vaf_boxplot_per_tier(combined_df, str(output_dir)))
    figs.append(plot_dp_boxplot_per_tier(combined_df, str(output_dir)))
    figs.append(plot_dna_vs_rna_vaf(combined_df, str(output_dir)))
    figs.append(plot_dna_vs_rna_dp(combined_df, str(output_dir)))
    figs.append(plot_ref_alt_dp_scatter(combined_df, str(output_dir)))
    figs.append(plot_variant_type_distribution(combined_df, str(output_dir)))
    figs.append(plot_ti_tv_ratio(combined_df, str(output_dir)))
    figs.append(plot_cross_modality(combined_df, str(output_dir)))
    figs.append(plot_gt_concordance_per_tier(combined_df, str(output_dir)))
    figs.append(plot_tiered_caller_overlap(combined_df, str(output_dir)))
    figs.append(plot_tiered_variant_types(combined_df, str(output_dir)))
    figs.append(plot_filter_distribution(combined_df, str(output_dir)))
    figs.append(plot_per_tier_vaf_boxplot(combined_df, str(output_dir)))
    figs.append(plot_caller_agreement_matrix(combined_df, str(output_dir)))
    figs.append(plot_chromosome_density(combined_df, str(output_dir)))
    figs.append(plot_dna_vs_rna_per_caller(combined_df, str(output_dir)))
    figs.append(plot_tier_quality_distribution(combined_df, str(output_dir)))
    figs.append(plot_redi_evidence(combined_df, str(output_dir)))

    # BAM charts
    if not args.no_bam and not bam_stats_df.is_empty():
        figs.append(plot_bam_metrics_bars(bam_stats_df, str(output_dir)))
        figs.append(plot_bam_coverage_violin(combined_df, str(output_dir)))

    if all_stats:
        sample_stats_df = pl.DataFrame(all_stats)
        figs.append(plot_per_sample_distribution(sample_stats_df, str(output_dir)))

        # Per-sample per-tier charts
        sample_tier_df = sample_tier_summary(combined_df)
        if not sample_tier_df.is_empty():
            figs.append(plot_per_sample_tier_distribution(sample_tier_df, str(output_dir)))
            figs.append(plot_per_tier_cross_sample_vaf(sample_tier_df, str(output_dir)))

    if not args.no_validate and report is not None:
        figs.append(plot_validation_heatmap(report, str(output_dir)))

    generate_dashboard(figs, str(output_dir))

    print(f"\nAll outputs written to: {output_dir}")
    print("Done.")


if __name__ == "__main__":
    main()
