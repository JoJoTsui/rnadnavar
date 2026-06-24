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
import ctypes
import gc
import os
import sys
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import polars as pl

# ── Allocator utilities ────────────────────────────────────────────────────

def _malloc_trim() -> None:
    """Release free glibc heap pages back to the kernel via malloc_trim(0).

    On non-glibc platforms (musl, macOS), the call is silently skipped.
    Takes <1ms — safe to call frequently.
    """
    try:
        libc = ctypes.CDLL("libc.so.6")
        libc.malloc_trim(0)
    except Exception:
        pass  # Non-glibc platform — skip

# ── Large-sample throttle ──────────────────────────────────────────────────
# thread mode: threading.Semaphore
# spawn mode: multiprocessing.Manager().Semaphore (passed via args)
# Both share the same acquire/release interface.

_LARGE_THRESHOLD = 2_000_000
_THREAD_LARGE_SEM = threading.Semaphore(1)  # used in --process-mode thread

# Columns needed during caller join + tiering. The remaining ~63 rescue INFO
# fields are output-only — they hstack back before parquet write.
_SLIM_COLS = [
    "CHROM", "POS", "REF", "ALT", "FILTER",
    "FILTERS_NORMALIZED", "GNOMAD_AF", "COSMIC_CNT", "REDI_EVIDENCE",
    "N_DNA_CALLERS_SUPPORT", "N_RNA_CALLERS_SUPPORT",
]


def _mem(msg: str) -> None:
    """Log current RSS in GB. Uses /proc/self/status — no external deps."""
    try:
        with open("/proc/self/status") as f:
            for line in f:
                if line.startswith("VmRSS:"):
                    kb = int(line.split()[1])
                    print(f"  [MEM {kb / 1024**2:.1f}GB] {msg}", flush=True)
                    return
    except Exception:
        print(f"  [MEM ?GB] {msg}", flush=True)

from .bam_stats import compute_all_bam_stats, read_and_merge_bed
from .bam_validation import validate_bam_vs_caller as validate_bam_one
from .caller_parser import _parse_one_caller, CALLERS_STRELKA, CALLERS_WITH_GT, parse_pre_norm_multiallelic
from .manifest_loader import CALLER_CONFIGS, filter_complete, load_manifest, get_pre_norm_vcf_path
from .rust_vcf import parse_rescue_vcf as _rust_parse_rescue
from .rescue_validator import validate_all_samples, validate_sample, validation_summary
from .rust_bam import pileup_variants, pileup_variants_multi
from .tiering_stats import compute_tiers_for_dataframe, tier_summary as compute_tier_summary
from .statistics import (
    compute_all_per_variant,
    compute_caller_wise_summary,
    compute_filter_effectiveness_matrix,
    compute_filter_vaf_dp_cross_tab,
    compute_fp_cross_tab,
    compute_low_vaf_rna_support,
    compute_rescue_breakdown,
    compute_rescue_by_caller_support,
    compute_rescue_by_filter,
    compute_rescue_by_tier,
    compute_rescue_cross_tab,
    compute_rescue_vaf_dp,
    compute_somatic_modality,
    compute_dp_threshold_sweep,
    compute_vaf_threshold_sweep,
    compute_wise_summary,
    dataset_summary,
    disease_summary,
    flag_filter_breakdown,
    gt_concordance,
    sample_rescue_summary,
    sample_summary,
    sample_tier_summary,
    set_summary,
    write_tsv,
)
from .visualizer import (
    _sort_chromosomes,
    generate_dashboard,
    plot_allele_balance_scatter,
    plot_bam_coverage_distribution,
    plot_bam_coverage_violin,
    plot_bam_metrics_bars,
    plot_bam_metrics_sample_wise,
    plot_caller_overlap,
    plot_bam_dp_distribution,
    plot_caller_tier_heatmap,
    plot_dp_distribution,
    plot_category_conflict_summary,
    plot_filter_vaf_dp_heatmap,
    plot_fp_cross_tab_heatmap,
    plot_low_vaf_rna_support,
    plot_mean_vaf_per_group,
    plot_mean_dp_per_group,
    plot_multiallelic_classification,
    plot_n_support_callers_dist,
    plot_sample_overview_scatter,
    plot_somatic_modality_pie,
    plot_somatic_modality_bars,
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
    plot_per_tier_dp_boxplot,
    plot_per_tier_vaf_boxplot,
    plot_tier_quality_distribution,
    plot_tiered_variant_types,
    plot_ti_tv_ratio,
    plot_redi_evidence,
    plot_vaf_distribution,
    plot_vaf_boxplot_per_tier,
    plot_vaf_sum_histogram,
    plot_dp_threshold_sweep,
    plot_vaf_threshold_sweep,
    plot_caller_concordance_vs_vaf,
    plot_filter_effectiveness_heatmap,
    plot_database_enrichment_by_tier,
    plot_rescue_breakdown,
    plot_rescue_by_filter,
    plot_rescue_by_tier,
    plot_rescue_caller_support,
    plot_rescue_cross_tab_heatmap,
    plot_rescue_dp_boxplot,
    plot_rescue_rate_trend,
    plot_rescue_sample_distribution,
    plot_rescue_vaf_boxplot,
    plot_validation_heatmap,
    plot_variant_type_distribution,
    plot_vc_distribution,
)


_STRELKA_CALLERS = ["DNA_strelka", "RNA_strelka"]
_DNA_CALLERS = ["DNA_deepsomatic", "DNA_mutect2", "DNA_strelka"]
_RNA_CALLERS = ["RNA_deepsomatic", "RNA_mutect2", "RNA_strelka"]


def build_unified_filter(args, combined_df_columns: set) -> pl.Expr | None:
    """Build a unified polars filter expression from CLI flags.

    Returns a polars expression that evaluates to True when a variant
    should be KEPT, or None when no filter is configured (pass all).

    The RNA leg of --min-vaf uses N_RNA_CALLERS_SUPPORT >= 2 AND
    RNA_DP_mean >= 10 (NOT RNA_VAF_mean), because RNA VAF is
    confounded by allele-specific expression.

    Args:
        args: Parsed argparse Namespace.
        combined_df_columns: Set of column names in combined_df for
            existence checks (columns are optional).

    Returns:
        polars expression combining all active filters, or None.
    """
    if args.no_filter:
        return None

    # Stack conditions with a FALSE baseline: no variant is kept unless
    # at least one keep clause fires.
    keep_expr = pl.lit(False)

    # ── --include-filters ──────────────────────────────────────────────────
    if args.include_filters:
        include_list = list(args.include_filters)
        keep_expr = keep_expr | pl.col("FILTER").is_in(include_list)
    else:
        keep_expr = pl.lit(True)  # include all FILTER values by default

    # ── --exclude-filters ──────────────────────────────────────────────────
    if args.exclude_filters:
        exclude_list = list(args.exclude_filters)
        keep_expr = keep_expr & ~pl.col("FILTER").is_in(exclude_list)

    # ── --min-dna-callers ─────────────────────────────────────────────────
    if args.min_dna_callers > 0:
        keep_expr = keep_expr & (pl.col("N_DNA_CALLERS_SUPPORT") >= args.min_dna_callers)

    # ── --min-rna-callers ─────────────────────────────────────────────────
    if args.min_rna_callers > 0:
        keep_expr = keep_expr & (pl.col("N_RNA_CALLERS_SUPPORT") >= args.min_rna_callers)

    # ── --min-vaf ─────────────────────────────────────────────────────────
    # DNA VAF >= threshold  OR  (>=2 RNA callers support AND RNA_DP_mean >= 10)
    if args.min_vaf and args.min_vaf > 0:
        dna_pass = pl.lit(False)
        if "DNA_VAF_mean" in combined_df_columns:
            dna_pass = dna_pass | (pl.col("DNA_VAF_mean") >= args.min_vaf)
        rna_pass = pl.lit(False)
        if "N_RNA_CALLERS_SUPPORT" in combined_df_columns and "RNA_DP_mean" in combined_df_columns:
            rna_pass = rna_pass | (
                (pl.col("N_RNA_CALLERS_SUPPORT") >= 2)
                & (pl.col("RNA_DP_mean") >= 10)
            )
        keep_expr = keep_expr & (dna_pass | rna_pass)

    # ── --min-dp ──────────────────────────────────────────────────────────
    if args.min_dp and args.min_dp > 0:
        dp_pass = pl.lit(False)
        if "DNA_DP_mean" in combined_df_columns:
            dp_pass = dp_pass | (pl.col("DNA_DP_mean") >= args.min_dp)
        if "RNA_DP_mean" in combined_df_columns:
            dp_pass = dp_pass | (pl.col("RNA_DP_mean") >= args.min_dp)
        keep_expr = keep_expr & dp_pass

    # ── --max-gnomad-af ───────────────────────────────────────────────────
    if args.max_gnomad_af is not None:
        keep_expr = keep_expr & (
            pl.col("GNOMAD_AF").is_null()
            | (pl.col("GNOMAD_AF") <= args.max_gnomad_af)
        )

    # ── --min-evidence-tier ───────────────────────────────────────────────
    if args.min_evidence_tier is not None:
        evidence_order = {
            "cross_modality": 0,
            "dna_confident": 1,
            "rna_rescued": 2,
            "low_confidence": 3,
        }
        min_rank = evidence_order[args.min_evidence_tier]
        # Build an OR chain: keep if rank <= min_rank
        tier_pass = pl.lit(False)
        for tier_name, rank in evidence_order.items():
            if rank <= min_rank:
                tier_pass = tier_pass | (pl.col("modality_evidence_caller") == tier_name)
        keep_expr = keep_expr & tier_pass

    # ── --exclude-multiallelic-conflict ───────────────────────────────────
    if args.exclude_multiallelic_conflict:
        if "flag_category_conflict" in combined_df_columns:
            keep_expr = keep_expr & ~pl.col("flag_category_conflict")

    return keep_expr


def _repair_strelka_columns(lf: pl.LazyFrame) -> pl.LazyFrame:
    """Detect and repair the Strelka TAR/TIR inversion in pre-existing parquet.

    Old parquet files have AD_ALT == TAR (ref counts) instead of TIR (alt counts).
    This function detects the inversion by sampling and applies lazy transforms
    to fix AD_ALT, AD_REF, VAF, and all derived mean columns.
    """
    cols = lf.collect_schema().names()
    # Check if raw TAR/TIR columns exist (needed for repair)
    tar_col = "DNA_strelka_TAR"
    tir_col = "DNA_strelka_TIR"
    ad_alt_col = "DNA_strelka_AD_ALT"
    if tar_col not in cols or tir_col not in cols or ad_alt_col not in cols:
        print("  Strelka repair: raw TAR/TIR columns not in parquet — skipping")
        return lf

    # Detect inversion: sample rows where both TAR and AD_ALT are non-null
    sample = (
        lf.filter(pl.col(tar_col).is_not_null() & pl.col(ad_alt_col).is_not_null())
        .select([tar_col, tir_col, ad_alt_col])
        .head(200)
        .collect()
    )
    if sample.is_empty():
        print("  Strelka repair: no non-null TAR/AD_ALT rows — skipping")
        return lf

    n_match_tar = (sample[ad_alt_col] == sample[tar_col]).sum()
    n_match_tir = (sample[ad_alt_col] == sample[tir_col]).sum()

    if n_match_tir > n_match_tar:
        print(f"  Strelka repair: AD_ALT already == TIR ({n_match_tir}/{sample.height}) — no repair needed")
        return lf

    print(f"  Strelka repair: DETECTED inversion — AD_ALT == TAR for {n_match_tar}/{sample.height} rows")
    print("  Repairing: AD_ALT ← TIR, AD_REF ← TAR, recomputing VAF and means...")

    # Apply repair for each Strelka caller
    repair_exprs = []
    for caller in _STRELKA_CALLERS:
        tar = f"{caller}_TAR"
        tir = f"{caller}_TIR"
        dp = f"{caller}_DP"
        ad_alt = f"{caller}_AD_ALT"
        ad_ref = f"{caller}_AD_REF"
        vaf = f"{caller}_VAF"

        if tir in cols:
            repair_exprs.append(pl.col(tir).alias(ad_alt))
        if tar in cols:
            repair_exprs.append(pl.col(tar).alias(ad_ref))
        if tir in cols and dp in cols:
            repair_exprs.append(
                pl.when(pl.col(dp).is_not_null() & (pl.col(dp) > 0))
                .then((pl.col(tir).cast(pl.Float64) / pl.col(dp).cast(pl.Float64)).clip(0.0, 1.0))
                .otherwise(None)
                .alias(vaf)
            )

    if repair_exprs:
        lf = lf.with_columns(repair_exprs)

    # Recompute mean columns across all callers (not just Strelka)
    mean_exprs = []
    for prefix, callers in [("DNA", _DNA_CALLERS), ("RNA", _RNA_CALLERS)]:
        vaf_cols = [f"{c}_VAF" for c in callers if f"{c}_VAF" in cols]
        alt_cols = [f"{c}_AD_ALT" for c in callers if f"{c}_AD_ALT" in cols]
        ref_cols = [f"{c}_AD_REF" for c in callers if f"{c}_AD_REF" in cols]
        if vaf_cols:
            mean_exprs.append(pl.mean_horizontal(vaf_cols).alias(f"{prefix}_VAF_mean"))
        if alt_cols:
            mean_exprs.append(pl.mean_horizontal(alt_cols).alias(f"{prefix}_ALT_DP_mean"))
        if ref_cols:
            mean_exprs.append(pl.mean_horizontal(ref_cols).alias(f"{prefix}_REF_DP_mean"))

    if mean_exprs:
        lf = lf.with_columns(mean_exprs)

    print("  Strelka repair: complete")
    return lf


def _streaming_join_one(df: pl.DataFrame, col_data: dict, caller_name: str) -> pl.DataFrame:
    """Join a single caller's column-oriented data and rename columns."""
    join_cols = ["CHROM", "POS", "REF", "ALT"]
    int_fields = {"DP", "AD_REF", "AD_ALT", "TAR", "TIR", "TOR", "AU", "CU", "GU", "TU", "POS"}
    float_fields = {"VAF_CALLER"}
    is_strelka = caller_name in CALLERS_STRELKA
    has_gt = caller_name in CALLERS_WITH_GT

    series_list = []
    data_col_names = []
    for cname, cvals in col_data.items():
        if cname in join_cols:
            series_list.append(pl.Series(cname, cvals, dtype=pl.Utf8 if cname != "POS" else pl.Int64))
        elif cname in int_fields:
            series_list.append(pl.Series(cname, cvals, dtype=pl.Int64))
            data_col_names.append(cname)
        elif cname in float_fields:
            series_list.append(pl.Series(cname, cvals, dtype=pl.Float64))
            data_col_names.append(cname)
        else:
            series_list.append(pl.Series(cname, cvals, dtype=pl.Utf8))
            data_col_names.append(cname)

    if not series_list or not data_col_names:
        return df

    caller_df = pl.DataFrame(series_list)
    rename_map = {c: f"{caller_name}_{c}" for c in data_col_names}
    caller_df = caller_df.select(join_cols + data_col_names).rename(rename_map)
    df = df.join(caller_df, on=join_cols, how="left")
    del caller_df, series_list

    if is_strelka:
        for suf, src in [("AD_REF", "TAR"), ("AD_ALT", "TIR")]:
            src_col = f"{caller_name}_{src}"
            tgt_col = f"{caller_name}_{suf}"
            if src_col in df.columns:
                df = df.with_columns(pl.col(src_col).alias(tgt_col))
    if not is_strelka and has_gt:
        for suf in ["AD_REF", "AD_ALT"]:
            col = f"{caller_name}_{suf}"
            if col not in df.columns:
                df = df.with_columns(pl.lit(None).alias(col))
    return df



def process_single_sample(row: dict, max_workers: int = 1, use_rust: bool = True,
                         large_sem=None, no_pileup: bool = False,
                         pileup_mode: str = "all",
                         bed_regions: list | None = None,
                         no_pre_norm: bool = False) -> dict:
    """Process one sample: parse rescue VCF + all caller VCFs + compute stats.

    Returns a dict with 'sample_id', 'df', and 'stats'.

    Args:
        row: Manifest row dict.
        max_workers: Threads for within-sample caller parsing.
        use_rust: Use Rust VCF parser.
        large_sem: Semaphore for large-sample exclusive access.
        no_pileup: Skip variant-wise BAM pileup (enabled by default).
        pileup_mode: "all" (default) or "filtered" (exclude NoConsensus).
        bed_regions: Optional merged BED regions for region-guided pileup (WES).
        no_pre_norm: Skip pre-decomposition caller VCF parsing.
    """
    sample_id = row["sample_id"]
    rescue_path = row["rescue_vcf_path"]
    base_dir = row["base_output_dir"]
    dir_name = row["dir_name"]
    vcf_prefix = row["vcf_prefix"]

    print(f"  [{sample_id}] Parsing rescue VCF (rust)...")
    rescue_df = _rust_parse_rescue(rescue_path)
    if rescue_df.is_empty():
        print(f"  [{sample_id}] WARNING: No variants in rescue VCF — sample SKIPPED")
        return {"sample_id": sample_id, "df": None, "stats": None}
    _mem(f"after rescue parse ({len(rescue_df)} vars)")

    # Auto-throttle: large samples process exclusively
    n_variants = len(rescue_df)
    is_large = n_variants > _LARGE_THRESHOLD and large_sem is not None
    if is_large:
        print(f"  [{sample_id}] Large sample ({n_variants} variants) — waiting for exclusive access...")
        try:
            acquired = large_sem.acquire(timeout=300)
        except TypeError:
            # Manager proxy may not support timeout — fall back to blocking acquire
            acquired = large_sem.acquire()
        if not acquired:
            print(f"  [{sample_id}] WARNING: Semaphore timeout — processing without exclusive access")
            is_large = False
        else:
            print(f"  [{sample_id}] Large sample acquired exclusive access")
    try:
        # Column pruning: split into slim (12 processing cols) and output (rest).
        slim_cols = [c for c in _SLIM_COLS if c in rescue_df.columns]
        output_cols = [c for c in rescue_df.columns if c not in _SLIM_COLS]
        rescue_slim = rescue_df.select(slim_cols)
        rescue_output = rescue_df.select(output_cols) if output_cols else None
        del rescue_df
        _mem("after column split")

        # Build target positions in chunks
        chunk_size = 100_000
        n_positions = len(rescue_slim)
        target_positions = set()
        for start in range(0, n_positions, chunk_size):
            end = min(start + chunk_size, n_positions)
            chunk = rescue_slim[start:end, ["CHROM", "POS", "REF", "ALT"]]
            for t in chunk.iter_rows():
                target_positions.add((t[0], t[1], t[2], t[3]))
            del chunk
        del start, end

        print(f"  [{sample_id}] {len(target_positions)} positions, streaming join...")
        base = os.path.join(base_dir, dir_name)

        # Streaming join: parse + join one caller at a time.
        df = rescue_slim.clone()
        for caller_name, cfg in CALLER_CONFIGS.items():
            # Pre-computed path from manifest (avoids per-caller glob at runtime)
            caller_vcf_col = f"caller_{caller_name.lower()}"
            caller_vcf_path = row.get(caller_vcf_col) or None
            name, cols = _parse_one_caller(caller_name, cfg, base, vcf_prefix,
                                           target_positions, vcf_path=caller_vcf_path)
            if cols:
                df = _streaming_join_one(df, cols, caller_name)
                del cols
            gc.collect()
            _mem(f"after caller: {caller_name}")

        del target_positions, rescue_slim
        gc.collect()
        _mem("after all callers joined")

        # Re-join output-only rescue columns before stats
        if rescue_output is not None:
            df = df.hstack(rescue_output)
            del rescue_output
            _mem("after hstack output cols")

        # ── Pre-decomposition caller VCF parsing (multi-allelic enrichment) ──
        # Parses pre-norm Mutect2 and DeepSomatic VCFs to recover full AD arrays,
        # per-allele AF, strand bias, and GT co-occurrence lost during vt decompose.
        if not no_pre_norm:
            for caller_name in ["DNA_mutect2", "RNA_mutect2", "DNA_deepsomatic", "RNA_deepsomatic"]:
                pre_norm_path = get_pre_norm_vcf_path(base_dir, dir_name, vcf_prefix, caller_name)
                if pre_norm_path:
                    try:
                        sample_suffix = CALLER_CONFIGS[caller_name].get("sample_suffix", "")
                        allele_registry = parse_pre_norm_multiallelic(
                            pre_norm_path, caller_name, sample_suffix,
                        )
                        if allele_registry is not None and not allele_registry.is_empty():
                            # Explode the allele registry: one row per (CHROM, POS, REF, individual ALT)
                            # for joining with the decomposed normalized records.
                            exploded_rows = []
                            for row in allele_registry.iter_rows(named=True):
                                chrom, pos, ref = row["CHROM"], row["POS"], row["REF"]
                                n_alts = row["n_original_alleles"]
                                for i in range(n_alts):
                                    exploded_rows.append({
                                        "CHROM": chrom,
                                        "POS": pos,
                                        "REF": ref,
                                        "ALT": row["original_alts"][i],
                                        f"pre_norm_{caller_name}_n_alleles": n_alts,
                                        f"pre_norm_{caller_name}_ad_ref": row["ad_ref"],
                                        f"pre_norm_{caller_name}_ad_alt": row["ad_alts"][i] if i < len(row["ad_alts"]) else None,
                                        f"pre_norm_{caller_name}_af": row["af_list"][i] if i < len(row["af_list"]) else None,
                                        f"pre_norm_{caller_name}_f1r2": row["f1r2_list"][i + 1] if row["f1r2_list"] and i + 1 < len(row["f1r2_list"]) else None,
                                        f"pre_norm_{caller_name}_f2r1": row["f2r1_list"][i + 1] if row["f2r1_list"] and i + 1 < len(row["f2r1_list"]) else None,
                                        f"pre_norm_{caller_name}_gt_alleles": row["gt_alleles"],
                                    })
                            if exploded_rows:
                                registry_df = pl.DataFrame(exploded_rows)
                                # Only select pre-norm columns + join keys, excluding
                                # "ALT" since it already exists in df (polars raises
                                # "duplicate output name" on join otherwise).
                                registry_cols = (
                                    ["CHROM", "POS", "REF", "ALT"]
                                    + [c for c in registry_df.columns
                                       if c not in ("CHROM", "POS", "REF", "ALT")]
                                )
                                df = df.join(
                                    registry_df.select(registry_cols),
                                    on=["CHROM", "POS", "REF", "ALT"],
                                    how="left",
                                )
                                print(f"  [{sample_id}] Pre-norm multi-allelic: {len(exploded_rows)} alleles across {len(allele_registry)} sites from {caller_name}")
                                del registry_df, exploded_rows
                        del allele_registry
                    except Exception as e:
                        print(f"  [{sample_id}] Pre-norm parsing failed for {caller_name}: {e}")
            _mem("after pre-norm multi-allelic enrichment")

        # Add sample metadata
        set_num = row.get("set_number", 0)
        if set_num is None:
            set_num = 0
        disease_val = row.get("disease", "")
        disease_norm = row.get("disease_normalized", "")
        df = df.with_columns([
            pl.lit(sample_id).alias("sample_id"),
            pl.lit(int(set_num)).cast(pl.Int64).alias("set_number"),
            pl.lit(str(disease_val)).alias("disease"),
            pl.lit(str(disease_norm)).alias("disease_normalized"),
        ])

        # Compute per-variant statistics (VAF, means)
        df = compute_all_per_variant(df)

        # Compute CxDy tiers via tiering engine
        df = compute_tiers_for_dataframe(df)

        # Compute RESCUED / CROSS_MODALITY backward compat columns from
        # modality_evidence_caller (produced by compute_tiers_for_dataframe).
        # These are derived columns for downstream chart functions that still
        # expect the legacy rescue VCF INFO field names.
        if "modality_evidence_caller" in df.columns:
            df = df.with_columns([
                pl.when(pl.col("modality_evidence_caller") == "cross_modality")
                .then(pl.lit("YES"))
                .otherwise(pl.lit("NO"))
                .alias("CROSS_MODALITY"),
                pl.when(pl.col("modality_evidence_caller").is_in(["cross_modality", "rna_rescued"]))
                .then(pl.lit("YES"))
                .otherwise(pl.lit("NO"))
                .alias("RESCUED"),
            ])

        # Variant-wise BAM pileup (enabled by default, --no-pileup to skip).
        # Computes per-position DP, strand bias, BQ, MQ from alignment BAMs.
        # Uses combined multi-BAM Rust function for efficiency (all BAMs in
        # one FFI call, shared position grouping, binary search inner loop).
        if not no_pileup:
            from .manifest_loader import get_manifest_bam_paths

            # Build positions list for all variants.
            # Per-variant filtering is handled by the unified filter pipeline,
            # not at the pileup level.
            cols_4 = ["CHROM", "POS", "REF", "ALT"]
            positions = [(row[0], row[1], row[2], row[3])
                         for row in df.select(cols_4).iter_rows()]

            n_pos = len(positions)
            print(f"  [{sample_id}] BAM pileup: {n_pos} positions")

            # Collect available BAM paths (manifest-first, fallback to glob)
            bam_paths = {
                bt: bp for bt, bp in get_manifest_bam_paths(
                    base_dir, dir_name, manifest_row=row).items()
                if bp
            }

            if bam_paths:
                # Pre-processing BAM integrity check — validate BGZF EOF marker
                # before hours of pileup computation. All-truncated samples skip
                # pileup entirely; partially truncated samples proceed with warnings.
                from .bam_stats import _check_bam_eof as check_bam_eof
                bam_ok = {}
                bam_bad = []
                for bt, bp in bam_paths.items():
                    ok, msg = check_bam_eof(bp)
                    bam_ok[bt] = ok
                    if not ok:
                        bam_bad.append((bt, msg))
                if bam_bad:
                    all_bad = len(bam_bad) == len(bam_paths)
                    for bt, msg in bam_bad:
                        if all_bad:
                            print(f"  [{sample_id}] ERROR: {msg}")
                        else:
                            print(f"  [{sample_id}] WARNING: {msg}")
                    if all_bad:
                        print(f"  [{sample_id}] All BAMs truncated — skipping pileup")
                        bam_paths = {}
                    else:
                        # Keep only valid BAMs
                        bam_paths = {bt: bp for bt, bp in bam_paths.items() if bam_ok[bt]}

            if bam_paths:
                try:
                    print(f"  [{sample_id}] BAM pileup: combined multi-BAM ({','.join(bam_paths.keys())})...")
                    pileup_results = pileup_variants_multi(bam_paths, positions, bed_regions)
                    for bt, pileup_df in pileup_results.items():
                        if pileup_df is not None and not pileup_df.is_empty():
                            join_key = ["CHROM", "POS", "REF", "ALT"]
                            rename = {c: f"BAM_{bt}_{c}" for c in pileup_df.columns
                                      if c not in join_key}
                            pileup_df = pileup_df.rename(rename)
                            df = df.join(
                                pileup_df.select(join_key + list(rename.values())),
                                on=join_key, how="left",
                            )
                            del pileup_df
                            _mem(f"after pileup {bt}")
                except BaseException as e:
                    # PanicException (from Rust panics) inherits from BaseException,
                    # not Exception. Catch BaseException to handle both panics and
                    # normal errors.
                    if isinstance(e, KeyboardInterrupt):
                        raise
                    print(f"  [{sample_id}] BAM pileup error: {e}")
                    # Fall back to per-BAM calls if multi-BAM fails
                    for bt, bam_path in bam_paths.items():
                        try:
                            pileup_df = pileup_variants(bam_path, positions)
                            if pileup_df is not None and not pileup_df.is_empty():
                                join_key = ["CHROM", "POS", "REF", "ALT"]
                                rename = {c: f"BAM_{bt}_{c}" for c in pileup_df.columns
                                          if c not in join_key}
                                pileup_df = pileup_df.rename(rename)
                                df = df.join(
                                    pileup_df.select(join_key + list(rename.values())),
                                    on=join_key, how="left",
                                )
                                del pileup_df
                        except BaseException as e2:
                            if isinstance(e2, KeyboardInterrupt):
                                raise
                            print(f"  [{sample_id}] BAM pileup {bt} fallback error: {e2}")
            del positions

        # Compute sample-level summary
        stats = sample_summary(df, sample_id)
        stats["set_number"] = row.get("set_number", 0)
        stats["disease"] = row.get("disease", "")

        _mem("after stats")
        print(f"  [{sample_id}] Done: {len(df)} variants")
        return {"sample_id": sample_id, "df": df, "stats": stats}
    finally:
        if is_large and large_sem is not None:
            large_sem.release()
            print(f"  [{sample_id}] Large sample released exclusive access")


# ── Multiprocessing worker (module-level, must be picklable) ──────────────

def _process_worker(args: tuple) -> dict:
    """Process one sample in a worker process. Writes parquet, returns only stats.

    This is the entry point for multiprocessing.Pool workers. It must be
    defined at module level so the 'fork' context can call it.

    Args:
        args: (row_dict, max_workers, use_rust, variant_dir_str, no_pileup, pileup_mode, bed_regions, no_pre_norm)

    Returns:
        {"sample_id": str, "stats": dict} — no DataFrame (too large to pickle).
        On error: {"sample_id": str, "error": str}
    """
    row, max_workers, use_rust, variant_dir_str, no_pileup, pileup_mode, bed_regions, no_pre_norm = args

    try:
        result = process_single_sample(row, max_workers=max_workers, use_rust=use_rust,
                                       large_sem=None, no_pileup=no_pileup,
                                       pileup_mode=pileup_mode, bed_regions=bed_regions,
                                       no_pre_norm=no_pre_norm)
    except BaseException:
        import traceback
        return {"sample_id": row["sample_id"], "error": traceback.format_exc()}

    sample_id = result["sample_id"]
    df = result["df"]
    stats = result["stats"]

    if df is not None:
        try:
            parquet_path = os.path.join(variant_dir_str, f"{sample_id}_variants.parquet")
            df.write_parquet(parquet_path)
        except Exception:
            import traceback
            return {"sample_id": sample_id, "error": f"parquet write failed:\n{traceback.format_exc()}"}
        finally:
            del df
            gc.collect()
            _malloc_trim()

    return {"sample_id": sample_id, "stats": stats}


def main():
    # Suppress cosmetic resource_tracker warning from multiprocessing.spawn.
    # The 6 POSIX named semaphores from mp.Pool's internal SimpleQueue objects
    # are cleaned up by the kernel on process exit regardless.
    import warnings
    warnings.filterwarnings('ignore', message='resource_tracker')

    parser = argparse.ArgumentParser(description="Seq2neo variant statistics")
    parser.add_argument("--manifest", required=True, help="Path to sample manifest TSV/Parquet")
    parser.add_argument("--output-dir", required=True, help="Output directory for statistics")
    parser.add_argument("--threads", type=int, default=6,
                        help="Threads for within-sample caller parsing (1-6, default: 6).")
    parser.add_argument("--sample-workers", type=int, default=1,
                        help="Parallel sample processing threads (default: 1). "
                             "Uses ThreadPoolExecutor (safe with htslib). Set to 4-8 for 32-core machines.")
    parser.add_argument("--max-samples", type=int, default=None, help="Limit number of samples")
    parser.add_argument("--set", type=int, default=None, help="Process only this set (1-4)")
    parser.add_argument("--sample-ids", nargs="*", default=None, help="Process specific sample IDs")
    parser.add_argument("--exclude-sample-ids", nargs="*", default=None,
                        help="Exclude specific sample IDs from processing")
    parser.add_argument("--no-validate", action="store_true", help="Skip rescue VCF validation")
    parser.add_argument("--tolerance", type=float, default=0.01, help="Validation tolerance")
    parser.add_argument("--pileup-mode", choices=["all", "filtered"], default="all",
                        help="Variant-wise BAM pileup mode: all variants (default) or exclude NoConsensus. "
                             "DEPRECATED: use --exclude-filters NoConsensus instead.")
    parser.add_argument("--include-filters", nargs="*", metavar="FILTER", default=None,
                        help="Only keep variants whose FILTER value is in this list (e.g. Somatic Germline). "
                             "Applied after all other filters.")
    parser.add_argument("--exclude-filters", nargs="*", metavar="FILTER", default=None,
                        help="Exclude variants whose FILTER value is in this list (e.g. NoConsensus Artifact). "
                             "Applied after all other filters.")
    parser.add_argument("--min-dna-callers", type=int, default=0,
                        help="Minimum N_DNA_CALLERS_SUPPORT (default: 0, no filter)")
    parser.add_argument("--min-rna-callers", type=int, default=0,
                        help="Minimum N_RNA_CALLERS_SUPPORT (default: 0, no filter)")
    parser.add_argument("--min-evidence-tier", choices=["cross_modality", "dna_confident", "rna_rescued",
                        "low_confidence"], default=None,
                        help="Minimum modality_evidence_caller tier. Variants below this tier are excluded. "
                             "Order: cross_modality > dna_confident > rna_rescued > low_confidence")
    parser.add_argument("--max-gnomad-af", type=float, default=None,
                        help="Maximum gnomAD AF threshold. Variants with GNOMAD_AF above this are excluded. "
                             "NULL GNOMAD_AF values pass through (no gnomAD entry = not excluded).")
    parser.add_argument("--exclude-multiallelic-conflict", action="store_true",
                        help="Exclude variants at multi-allelic sites with category conflicts "
                             "(flag_category_conflict). These are positions where different "
                             "callers assigned different biological FILTER categories.")
    parser.add_argument("--no-filter", action="store_true",
                        help="Disable the entire unified filter pipeline. All variants pass through unfiltered.")
    parser.add_argument("--no-hard-filter", action="store_true",
                        help="Disable Stage 0 hard exclusions (auto-drop of false variants). "
                             "Variants will still receive soft flags and confidence tiers.")
    parser.add_argument("--confidence-tier", choices=["high", "medium", "low"], default=None,
                        help="Filter output to only the specified confidence tier.")
    parser.add_argument("--no-pre-norm", action="store_true",
                        help="Skip pre-decomposition caller VCF parsing (for environments "
                             "where variant_calling/ VCFs are unavailable).")
    parser.add_argument("--no-bam", action="store_true", help="Skip all BAM processing")
    parser.add_argument("--resume", action="store_true", help="Skip variant processing and BAM stats — go straight to "
                        "aggregation + validation + visualizations using existing parquet files")
    parser.add_argument("--no-pileup", action="store_true", help="Skip variant-wise BAM pileup (whole-genome only)")
    parser.add_argument("--bed", type=str, default=None, help="Path to BED file for WES coverage denominator. "
                        "When provided, mean_coverage is recalculated using BED region total instead of "
                        "whole-genome reference lengths (which under-report coverage by ~100× for WES).")
    parser.add_argument("--bam-workers", type=int, default=8,
                        help="Threads for parallel BAM stats processing (default: 8). "
                             "Uses ThreadPoolExecutor (safe with htslib).")
    parser.add_argument("--wise", nargs="*", default=None, metavar="WISE",
                        help="Generate specific wise summaries and charts (space-separated). "
                             "Choices: set, disease, sample, tier, caller, chromosome, threshold. "
                             "Default: all wises.")
    parser.add_argument("--exclude-disease", nargs="*", default=None, metavar="DISEASE",
                        help="Exclude specific diseases (space-separated) from filtered parquet output. "
                             "Used for zero-shot experiments — variants from these diseases are "
                             "written to a separate variant_details_filtered/ directory.")
    parser.add_argument("--min-vaf", type=float, default=None,
                        help="Minimum VAF threshold for DNA evidence. Variants below this "
                             "threshold are excluded UNLESS >=2 RNA callers support with "
                             "RNA_DP>=10 (RNA caller support replaces RNA VAF which is "
                             "confounded by allele-specific expression). Applies to ALL "
                             "output with unified filter pipeline.")
    parser.add_argument("--min-dp", type=int, default=None,
                        help="Minimum DP threshold for filtered parquet output. Variants below this "
                             "threshold in BOTH modalities are excluded from variant_details_filtered/.")
    parser.add_argument("--verbose", action="store_true",
                        help="Show per-caller progress messages")
    parser.add_argument("--process-mode", choices=["thread", "spawn"], default=None,
                        help="Parallel execution mode: thread (ThreadPoolExecutor) or spawn "
                             "(multiprocessing.Pool with fresh Python process per worker). "
                             "Default: spawn when --sample-workers > 1, thread otherwise.")
    parser.add_argument("--max-tasks-per-child", type=int, default=1,
                        help="Max samples per worker process before restart (spawn mode only, default: 1).")
    parser.add_argument("--theme", choices=["default", "publishing"], default="default",
                        help="Chart theme: default (Altair built-in) or publishing (clean journal-ready style)")
    args = parser.parse_args()

    # ── --pileup-mode deprecation ──────────────────────────────────────────
    if args.pileup_mode == "filtered":
        import warnings
        warnings.warn(
            "--pileup-mode filtered is deprecated. Use --exclude-filters NoConsensus instead. "
            "Mapping \"filtered\" to --exclude-filters NoConsensus.",
            DeprecationWarning,
        )
        if args.exclude_filters is None:
            args.exclude_filters = ["NoConsensus"]

    # ── --exclude-disease deprecation ─────────────────────────────────────
    # Unified filter handles this. Keep the arg for backward compat:
    # map to --exclude-filters. The old variant_details_filtered/ section
    # is removed (unified filter applies to all output).
    if args.exclude_disease and not args.no_filter:
        print("NOTE: --exclude-disease is applied via unified filter as a FILTER-based exclusion. "
              "The old variant_details_filtered/ output directory is removed — filtered output is "
              "now the default variant_details/ parquet files.")

    # Load and filter manifest
    manifest = load_manifest(args.manifest)
    manifest = filter_complete(manifest)

    if args.set:
        manifest = manifest.filter(pl.col("set_number") == args.set)
    if args.sample_ids:
        manifest = manifest.filter(pl.col("sample_id").is_in(args.sample_ids))
    if args.max_samples:
        manifest = manifest.head(args.max_samples)
    if args.exclude_sample_ids:
        manifest = manifest.filter(~pl.col("sample_id").is_in(args.exclude_sample_ids))

    if manifest.is_empty():
        print("No samples to process.")
        sys.exit(0)

    use_rust = True  # Rust parser is the only backend (cyvcf2 removed)

    # Determine process mode (used by print below and execution logic below)
    process_mode = args.process_mode
    if process_mode is None:
        process_mode = "spawn" if args.sample_workers > 1 else "thread"

    print(f"Processing {len(manifest)} samples (parser=rust, caller_threads={args.threads}, sample_workers={args.sample_workers}, bam_workers={args.bam_workers}, process_mode={process_mode})")
    if len(manifest) > 20 and not args.no_validate:
        print("NOTE: >20 samples with validation enabled may be slow due to BAM pileup.")
        print("      Consider --no-validate for initial runs, then validate separately.")
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    rows = manifest.to_dicts()
    all_stats = []
    bam_stats_df = pl.DataFrame()

    # ── Read and merge BED regions early (shared between pileup and BAM stats) ──
    raw_bed_total = 0
    raw_bed_regions = None
    bed_regions = None   # merged regions for pileup (gap=500Kb)
    if args.bed:
        if not os.path.isfile(args.bed):
            print(f"WARNING: BED file not found: {args.bed}, using whole-genome mode")
        else:
            raw_bed_total, raw_bed_regions, bed_regions = read_and_merge_bed(args.bed, gap=500_000)
            print(f"BED raw total: {raw_bed_total:,} bp ({raw_bed_total / 1e6:.1f} Mbp), "
                  f"{len(raw_bed_regions)} raw intervals, "
                  f"{len(bed_regions)} merged regions (500Kb gap)")

    # ── Per-sample parquet directory (streaming, not memory-accumulated) ──
    variant_dir = output_dir / "variant_details"
    variant_dir.mkdir(parents=True, exist_ok=True)
    variant_dir_str = str(variant_dir)
    total_variants = 0

    # ── Launch BAM stats in background before variant processing ───────────
    # BAM stats (manifest-only, no variant data needed) runs concurrently
    # with variant processing, hiding BAM latency (~9 min for large samples)
    # behind variant processing (~2-3 min). Max wall-clock = max(variant, BAM)
    # instead of sum.
    bam_future = None
    bam_bg_executor = None
    if not args.no_bam and not args.resume:
        bam_bg_executor = ThreadPoolExecutor(max_workers=1)
        bam_future = bam_bg_executor.submit(
            compute_all_bam_stats, rows,
            max_workers=args.bam_workers,
            bed_total=raw_bed_total,
            bed_regions=raw_bed_regions,
        )
        print("BAM statistics running in background...")

    if args.resume:
        # Skip variant processing + BAM stats — parquet files already exist.
        print(f"Resuming from existing parquet files in {variant_dir_str}/")
        import glob as _glob

        # Count variants, respecting active filters if any
        if args.max_samples or args.set or args.sample_ids:
            sample_ids = {r["sample_id"] for r in rows}
            for pq in _glob.glob(os.path.join(variant_dir_str, "*_variants.parquet")):
                sid = os.path.basename(pq).replace("_variants.parquet", "")
                if sid in sample_ids:
                    try:
                        total_variants += pl.scan_parquet(pq).select(pl.len()).collect().item()
                    except Exception:
                        pass
        else:
            for pq in _glob.glob(os.path.join(variant_dir_str, "*_variants.parquet")):
                try:
                    total_variants += pl.scan_parquet(pq).select(pl.len()).collect().item()
                except Exception:
                    pass

        if total_variants == 0:
            print("No existing parquet files found. Run without --resume first.")
            sys.exit(1)
        print(f"Found {total_variants} variants in existing parquet files")

        # Reload per-sample stats so downstream TSVs + charts still generate
        stats_tsv = output_dir / "sample_summary.tsv"
        if stats_tsv.exists():
            all_stats = pl.read_csv(str(stats_tsv), separator="\t").to_dicts()
        else:
            # Fallback: check legacy CSV for backwards compatibility with old runs
            stats_csv = output_dir / "sample_summary.csv"
            if stats_csv.exists():
                all_stats = pl.read_csv(str(stats_csv)).to_dicts()
            else:
                # Recompute from parquet files
                print("  Recomputing sample_summary from parquet files...")
                import glob as _glob_resume
                for pq in sorted(_glob_resume.glob(str(variant_dir / "*_variants.parquet"))):
                    sid = os.path.basename(pq).replace("_variants.parquet", "")
                    df_tmp = pl.read_parquet(pq)
                    all_stats.append(sample_summary(df_tmp, sid))
                    del df_tmp
                print(f"  Recomputed {len(all_stats)} sample summaries")

        # Reload BAM stats so BAM charts still render
        bam_tsv = output_dir / "bam_stats.tsv"
        if bam_tsv.exists():
            bam_stats_df = pl.read_csv(str(bam_tsv), separator="\t")
        else:
            # Fallback: check legacy CSV for backwards compatibility
            bam_csv = output_dir / "bam_stats.csv"
            if bam_csv.exists():
                bam_stats_df = pl.read_csv(str(bam_csv))
            else:
                print("  No bam_stats found — BAM charts will be skipped")

        if not bam_stats_df.is_empty():
            from .bam_stats import ensure_bam_stats_columns
            bam_stats_df = ensure_bam_stats_columns(bam_stats_df)

        # Skip BAM stats (they were already computed)
        args.no_bam = True

    elif process_mode == "spawn" and args.sample_workers > 1:
        # ── Process-isolated parallel mode (spawn) ──────────────────────────
        # spawn creates a fresh Python interpreter per worker (fork+exec),
        # avoiding the fork+threads deadlock with polars' rayon pool.
        # maxtasksperchild=1 ensures each worker exits after one sample
        # → kernel reclaims all memory. No cross-process throttle needed —
        # process isolation alone keeps total RSS within limits.
        import multiprocessing as mp
        ctx = mp.get_context("spawn")

        worker_args = [
            (row, args.threads, use_rust, variant_dir_str, args.no_pileup, args.pileup_mode, bed_regions, args.no_pre_norm)
            for row in rows
        ]

        pool = ctx.Pool(
            processes=min(args.sample_workers, len(rows)),
            maxtasksperchild=args.max_tasks_per_child,
        )
        try:
            for result in pool.imap_unordered(_process_worker, worker_args):
                sid = result["sample_id"]
                if "error" in result:
                    print(f"[ERROR] {sid}: {result['error']}", flush=True)
                    continue
                if result["stats"] is not None:
                    all_stats.append(result["stats"])
                # Count variants from the written parquet file
                parquet_path = os.path.join(variant_dir_str, f"{sid}_variants.parquet")
                try:
                    n = pl.scan_parquet(parquet_path).select(pl.len()).collect().item()
                    total_variants += n
                except Exception:
                    pass
                print(f"[{len(all_stats)}/{len(manifest)}] {sid} - Done")
        finally:
            pool.close()   # signal no more tasks → workers exit gracefully
            pool.join()    # wait for workers → Finalize runs → sem_unlink

    elif args.sample_workers > 1:
        # ── Thread-based parallel mode ─────────────────────────────────────
        def _process_one(row, max_workers, use_rust):
            """Process one sample and write its variant details to parquet immediately."""
            nonlocal total_variants
            result = process_single_sample(
                row, max_workers=max_workers, use_rust=use_rust,
                large_sem=_THREAD_LARGE_SEM, no_pileup=args.no_pileup,
                pileup_mode=args.pileup_mode, bed_regions=bed_regions,
                no_pre_norm=args.no_pre_norm,
            )
            if result["df"] is not None:
                sid = result["sample_id"]
                parquet_path = str(variant_dir / f"{sid}_variants.parquet")
                result["df"].write_parquet(parquet_path)
                total_variants += len(result["df"])
                # Free memory immediately
                result["df"] = None
                gc.collect()
                _malloc_trim()
                _mem(f"after write+free [{sid}]")
            if result["stats"] is not None:
                all_stats.append(result["stats"])
            return result

        with ThreadPoolExecutor(max_workers=args.sample_workers) as executor:
            futures = {}
            for i, row in enumerate(rows):
                future = executor.submit(_process_one, row, args.threads, use_rust)
                futures[future] = (i, row["sample_id"])

            for future in as_completed(futures):
                i, sid = futures[future]
                try:
                    future.result()
                    print(f"[{i+1}/{len(manifest)}] {sid} - Done")
                except Exception as e:
                    import traceback
                    print(f"[{i+1}/{len(manifest)}] {sid} - ERROR: {e}")
                    traceback.print_exc()
    else:
        # ── Sequential mode ────────────────────────────────────────────────
        for i, row in enumerate(rows):
            sid = row["sample_id"]
            print(f"\n[{i+1}/{len(manifest)}] {sid}")
            try:
                result = process_single_sample(
                    row, max_workers=args.threads, use_rust=use_rust,
                    large_sem=_THREAD_LARGE_SEM, no_pileup=args.no_pileup,
                    pileup_mode=args.pileup_mode, bed_regions=bed_regions,
                )
                if result["df"] is not None:
                    parquet_path = str(variant_dir / f"{sid}_variants.parquet")
                    result["df"].write_parquet(parquet_path)
                    total_variants += len(result["df"])
                    del result["df"]
                    gc.collect()
                    _malloc_trim()
                    _mem(f"after write+free [{sid}]")
                if result["stats"] is not None:
                    all_stats.append(result["stats"])
            except Exception as e:
                import traceback
                print(f"  [{sid}] ERROR: {e}")
                traceback.print_exc()

    if total_variants == 0:
        print("No data processed.")
        sys.exit(1)

    # ── Lazy scan across all per-sample parquet files ──────────────────────
    # Use extra_columns='ignore' to handle schema drift when new columns are
    # added to compute_biological_flags (e.g., category_conflict_resolution)
    combined_df = pl.scan_parquet(str(variant_dir / "*_variants.parquet"), extra_columns="ignore")

    # ── Repair Strelka TAR/TIR inversion in old parquet files ────────────
    # Old parquet has AD_ALT == TAR (ref counts) instead of TIR (alt counts).
    # Auto-detects and skips if already correct.
    combined_df = _repair_strelka_columns(combined_df)

    # ── Add ML train/val/test partition column (Section 7.1) ──────────────
    # Partition by chromosome: chr1 → test, chr21/chr22 → val, rest → train.
    # This is a deterministic, reproducible split that avoids data leakage
    # between genomic regions.
    combined_df = combined_df.with_columns(
        pl.when(pl.col("CHROM") == "chr1").then(pl.lit("test"))
        .when(pl.col("CHROM").is_in(["chr21", "chr22"])).then(pl.lit("val"))
        .otherwise(pl.lit("train"))
        .alias("partition")
    )

    # ── Unified filter pipeline ────────────────────────────────────────────
    # Apply all CLI filters (--include-filters, --exclude-filters, --min-vaf,
    # --min-dp, --max-gnomad-af, --min-evidence-tier, --exclude-multiallelic-conflict)
    # at the LazyFrame level before any statistics are computed.
    # Build the filter expression once, evaluate it, and keep only passing rows.
    combined_df_cols = set(combined_df.collect_schema().names())
    filter_expr = build_unified_filter(args, combined_df_cols)

    # Compute cross-sample recurrence before filtering (needs all samples)
    if "CHROM" in combined_df_cols and "POS" in combined_df_cols:
        if "REF" in combined_df_cols and "ALT" in combined_df_cols:
            combined_df = combined_df.with_columns(
                pl.len().over(["CHROM", "POS", "REF", "ALT"]).alias("n_recurrent_samples")
            )
            print("Computed n_recurrent_samples (cross-sample recurrence)")

    if filter_expr is not None:
        n_before = combined_df.select(pl.len()).collect().item()
        combined_df = combined_df.filter(filter_expr)
        n_after = combined_df.select(pl.len()).collect().item()
        if n_before != n_after:
            print(f"Unified filter: {n_before} → {n_after} variants "
                  f"({n_before - n_after} excluded)")

    # Stage 0 hard exclusions (auto-drop false variants)
    if not args.no_hard_filter:
        from .statistics import build_hard_filter_expr
        hard_expr = build_hard_filter_expr(combined_df.collect_schema().names())
        if hard_expr is not None:
            n_before_hard = combined_df.select(pl.len()).collect().item()
            combined_df = combined_df.filter(~hard_expr)
            n_after_hard = combined_df.select(pl.len()).collect().item()
            if n_before_hard != n_after_hard:
                print(f"Hard filter (Stage 0): {n_before_hard} → {n_after_hard} variants "
                      f"({n_before_hard - n_after_hard} dropped)")

    # Stage 1+2: compute soft flags and confidence tiers
    from .statistics import compute_confidence_tier, compute_soft_flags
    combined_df_cols = set(combined_df.collect_schema().names())
    if "confidence_tier" not in combined_df_cols:
        # Collect, compute confidence/soft flags, then re-scan
        combined_eager = combined_df.collect()
        combined_eager = compute_confidence_tier(combined_eager)
        combined_eager = compute_soft_flags(combined_eager)
        # Write confidence tier summary
        from .statistics import write_tsv as _write_tsv
        ct_summary = combined_eager.group_by(["confidence_tier", "FILTER"]).agg(
            pl.len().alias("count")
        ).sort(["confidence_tier", "FILTER"])
        if not ct_summary.is_empty():
            _write_tsv(ct_summary, str(output_dir / "confidence_tier_summary.tsv"))
            print(f"Confidence tier summary: {output_dir / 'confidence_tier_summary.tsv'}")
        # Replace combined_df with the enriched version
        combined_df = pl.scan_parquet(variant_dir / "*_variants.parquet", extra_columns="ignore")
        # Re-apply filters on the new scan
        combined_df_cols = set(combined_df.collect_schema().names())
        new_filter = build_unified_filter(args, combined_df_cols)
        if new_filter is not None:
            combined_df = combined_df.filter(new_filter)
        # Add n_recurrent_samples again
        if all(c in combined_df_cols for c in ["CHROM", "POS", "REF", "ALT"]):
            combined_df = combined_df.with_columns(
                pl.len().over(["CHROM", "POS", "REF", "ALT"]).alias("n_recurrent_samples")
            )
        # Re-apply hard filter
        if not args.no_hard_filter:
            hard_expr2 = build_hard_filter_expr(combined_df.collect_schema().names())
            if hard_expr2 is not None:
                combined_df = combined_df.filter(~hard_expr2)
        # Join confidence_tier and soft_flags from the eager computation
        ct_lookup = combined_eager.select(["sample_id", "CHROM", "POS", "confidence_tier", "soft_flags"])
        # Note: This is a best-effort approach; for production, confidence_tier should be
        # computed at per-sample parquet write time and stored in the parquet.
        # For now, we use the eager version for statistics.
        combined_df = combined_eager.lazy()

    print(f"Variant details: {variant_dir}/ (lazy scan, {total_variants} variants across {len(rows)} samples)")

    # BAM statistics — collect background result or compute now
    if not args.no_bam and not args.resume:
        if bam_future is not None:
            try:
                bam_stats_df = bam_future.result(timeout=3600)
                print("BAM statistics completed (background)")
            except TimeoutError:
                print("WARNING: BAM stats timed out after 3600s - BAM charts will be skipped")
                bam_stats_df = pl.DataFrame()
            except Exception as e:
                print(f"WARNING: Background BAM stats failed: {e}")
                bam_stats_df = pl.DataFrame()
            finally:
                if bam_bg_executor is not None:
                    bam_bg_executor.shutdown(wait=False)
        else:
            print("Computing per-sample BAM statistics...")
            bam_stats_df = compute_all_bam_stats(rows, max_workers=args.bam_workers,
                                                  bed_total=raw_bed_total, bed_regions=raw_bed_regions)
        if not bam_stats_df.is_empty():
            write_tsv(bam_stats_df, str(output_dir / "bam_stats.tsv"))
            print(f"BAM stats: {output_dir / 'bam_stats.tsv'}")
    elif args.no_bam:
        print("Skipping BAM statistics (--no-bam)")
    # else: resume mode — BAM stats already loaded from TSV above

    # Sample summary
    if all_stats:
        sample_stats_df = pl.DataFrame(all_stats)
        write_tsv(sample_stats_df, str(output_dir / "sample_summary.tsv"))
        print(f"Sample summary: {output_dir / 'sample_summary.tsv'}")

        # Set summary
        set_summary_df = set_summary(sample_stats_df)
        if not set_summary_df.is_empty():
            write_tsv(set_summary_df, str(output_dir / "set_summary.tsv"))

        # Disease summary
        _mem("before disease_summary")
        disease_summary_df = disease_summary(combined_df)
        if not disease_summary_df.is_empty():
            write_tsv(disease_summary_df, str(output_dir / "disease_summary.tsv"))
        _mem("after disease_summary")

        # Tier summary
        tier_summary_df = compute_tier_summary(combined_df)
        if not tier_summary_df.is_empty():
            write_tsv(tier_summary_df, str(output_dir / "tier_summary.tsv"))
            print(f"Tier summary: {output_dir / 'tier_summary.tsv'}")
        _mem("after tier_summary")

        # Dataset summary (whole-dataset aggregates)
        ds_summary = dataset_summary(combined_df)
        if ds_summary:
            write_tsv(pl.DataFrame([ds_summary]), str(output_dir / "dataset_summary.tsv"))
            print(f"Dataset summary: {output_dir / 'dataset_summary.tsv'}")
        _mem("after dataset_summary")

        # Sample-tier summary (Level 4)
        sample_tier_df = sample_tier_summary(combined_df)
        if not sample_tier_df.is_empty():
            write_tsv(sample_tier_df, str(output_dir / "sample_tier_summary.tsv"))
            print(f"Sample-tier summary: {output_dir / 'sample_tier_summary.tsv'}")
        _mem("after sample_tier_summary")

    # Tier distribution (renamed from caller_overlap.tsv — it contains tier counts)
    from .statistics import caller_overlap_distribution, caller_overlap_matrix, caller_support_distribution
    tier_dist_df = caller_overlap_distribution(combined_df)
    if not tier_dist_df.is_empty():
        write_tsv(tier_dist_df, str(output_dir / "tier_distribution.tsv"))
        print(f"Tier distribution: {output_dir / 'tier_distribution.tsv'}")

    # Caller overlap matrix (pairwise co-occurrence)
    overlap_matrix_df = caller_overlap_matrix(combined_df)
    if not overlap_matrix_df.is_empty():
        write_tsv(overlap_matrix_df, str(output_dir / "caller_overlap_matrix.tsv"))
        print(f"Caller overlap matrix: {output_dir / 'caller_overlap_matrix.tsv'}")

    # Caller support distribution (histogram of N_SUPPORT_CALLERS)
    support_dist_df = caller_support_distribution(combined_df)
    if not support_dist_df.is_empty():
        write_tsv(support_dist_df, str(output_dir / "caller_support_distribution.tsv"))
        print(f"Caller support distribution: {output_dir / 'caller_support_distribution.tsv'}")
    _mem("after caller overlap/distribution")

    # Filter distribution
    from .statistics import filter_distribution, variant_type_distribution, vc_distribution
    filter_df = filter_distribution(combined_df)
    if not filter_df.is_empty():
        write_tsv(filter_df, str(output_dir / "filter_distribution.tsv"))

    # Variant type distribution
    vt_df = variant_type_distribution(combined_df)
    if not vt_df.is_empty():
        write_tsv(vt_df, str(output_dir / "variant_type_distribution.tsv"))
    _mem("after filter/vt distributions")

    # GT concordance
    concordance_data = gt_concordance(combined_df)
    if concordance_data:
        write_tsv(pl.DataFrame({"agreement_level": list(concordance_data.keys()),
                                 "count": list(concordance_data.values())}),
                   str(output_dir / "gt_concordance.tsv"))
    _mem("after gt_concordance")

    # ── Initialize chart list early (used by threshold analysis below) ──────
    figs = []

    # ── Wise-based summaries ──────────────────────────────────────────────────
    # Generate set, disease, sample, tier, caller, chromosome summaries
    # using the shared wise kernel. Output to stats/{wise}/ directories.
    # Controlled by --wise flag (default: all wises).
    wise_names = args.wise
    if wise_names is not None and len(wise_names) == 0:
        wise_names = None  # --wise with no args → all wises
    if wise_names is not None:
        wise_names = set(wise_names)

    wise_configs = [
        ("set", ["set_number"]),
        ("disease", ["disease_normalized"]),
        ("sample", ["sample_id"]),
        ("tier", ["final_tier"]),
        ("chromosome", ["CHROM"]),
        ("variant-category", ["FILTER"]),
    ]

    if wise_names is None or any(w in wise_names for w in ["set", "disease", "sample", "tier", "chromosome", "variant-category"]):
        print("Generating wise summaries...")
        for wise_name, group_cols in wise_configs:
            if wise_names is not None and wise_name not in wise_names:
                continue
            wise_dir = output_dir / "stats" / wise_name
            wise_dir.mkdir(parents=True, exist_ok=True)
            try:
                wise_df = compute_wise_summary(combined_df, group_cols)
                if not wise_df.is_empty():
                    if wise_name == "chromosome":
                        wise_df = _sort_chromosomes(wise_df, "CHROM")
                    write_tsv(wise_df, str(wise_dir / f"{wise_name}_summary.tsv"))
            except Exception as e:
                print(f"  WARNING: {wise_name}-wise summary failed: {e}")

    # Caller-wise summary (special aggregation — per-caller VAF/DP)
    if wise_names is None or "caller" in wise_names:
        caller_dir = output_dir / "stats" / "caller"
        caller_dir.mkdir(parents=True, exist_ok=True)
        try:
            caller_df = compute_caller_wise_summary(combined_df)
            if not caller_df.is_empty():
                write_tsv(caller_df, str(caller_dir / "caller_summary.tsv"))
        except Exception as e:
            print(f"  WARNING: caller-wise summary failed: {e}")

    # ── Threshold analysis ────────────────────────────────────────────────────
    if wise_names is None or "threshold" in wise_names:
        print("Generating threshold analysis...")
        threshold_dir = output_dir / "stats" / "threshold"
        threshold_dir.mkdir(parents=True, exist_ok=True)
        try:
            vaf_sweep_df = compute_vaf_threshold_sweep(combined_df)
            if not vaf_sweep_df.is_empty():
                write_tsv(vaf_sweep_df, str(threshold_dir / "vaf_threshold_sweep.tsv"))
                print(f"  VAF threshold sweep: {threshold_dir / 'vaf_threshold_sweep.tsv'}")
                figs.append(plot_vaf_threshold_sweep(vaf_sweep_df, str(output_dir)))
        except Exception as e:
            print(f"  WARNING: VAF threshold sweep failed: {e}")

        try:
            filter_matrix_df = compute_filter_effectiveness_matrix(combined_df)
            if not filter_matrix_df.is_empty():
                write_tsv(filter_matrix_df, str(threshold_dir / "filter_effectiveness.tsv"))
                print(f"  Filter effectiveness: {threshold_dir / 'filter_effectiveness.tsv'}")
                figs.append(plot_filter_effectiveness_heatmap(filter_matrix_df, str(output_dir)))
        except Exception as e:
            print(f"  WARNING: filter effectiveness matrix failed: {e}")

        try:
            dp_sweep_df = compute_dp_threshold_sweep(combined_df)
            if not dp_sweep_df.is_empty():
                write_tsv(dp_sweep_df, str(threshold_dir / "dp_threshold_sweep.tsv"))
                print(f"  DP threshold sweep: {threshold_dir / 'dp_threshold_sweep.tsv'}")
                figs.append(plot_dp_threshold_sweep(dp_sweep_df, str(output_dir)))
        except Exception as e:
            print(f"  WARNING: DP threshold sweep failed: {e}")

        # ── ML Threshold Guidance (Section 7) ────────────────────────────────
        try:
            cross_tab_df = compute_filter_vaf_dp_cross_tab(combined_df)
            if not cross_tab_df.is_empty():
                write_tsv(cross_tab_df, str(threshold_dir / "filter_vaf_dp_cross_tab.tsv"))
                print(f"  FILTER x VAF x DP cross-tab: {threshold_dir / 'filter_vaf_dp_cross_tab.tsv'}")
                figs.append(plot_filter_vaf_dp_heatmap(cross_tab_df, str(output_dir)))
        except Exception as e:
            print(f"  WARNING: FILTER x VAF x DP cross-tab failed: {e}")

        try:
            low_vaf_df = compute_low_vaf_rna_support(combined_df)
            if not low_vaf_df.is_empty():
                write_tsv(low_vaf_df, str(threshold_dir / "low_vaf_rna_support.tsv"))
                print(f"  Low VAF RNA support: {threshold_dir / 'low_vaf_rna_support.tsv'}")
                figs.append(plot_low_vaf_rna_support(low_vaf_df, str(output_dir)))
        except Exception as e:
            print(f"  WARNING: Low VAF RNA support failed: {e}")

        # Partition summary (7.4)
        try:
            from .statistics import _ensure_eager
            part_df = _ensure_eager(combined_df)
            if "partition" in part_df.columns:
                part_summary = part_df.group_by("partition").agg([
                    pl.len().alias("n_variants"),
                    pl.col("sample_id").n_unique().alias("n_samples") if "sample_id" in part_df.columns else pl.lit(0).alias("n_samples"),
                ]).sort("partition")
                write_tsv(part_summary, str(threshold_dir / "partition_summary.tsv"))
                print(f"  Partition summary: {threshold_dir / 'partition_summary.tsv'}")

                # Disease x partition (7.5)
                if "disease_normalized" in part_df.columns:
                    disease_part = part_df.group_by(["disease_normalized", "partition"]).agg([
                        pl.len().alias("n_variants"),
                    ]).sort(["disease_normalized", "partition"])
                    write_tsv(disease_part, str(threshold_dir / "disease_partition_summary.tsv"))
                    print(f"  Disease x partition: {threshold_dir / 'disease_partition_summary.tsv'}")
                del part_df
        except Exception as e:
            print(f"  WARNING: Partition summary failed: {e}")

    # ── FP Cross-Tabulation (Section 8) ──────────────────────────────────────
    if wise_names is None or "threshold" in wise_names:
        try:
            fp_ct_df = compute_fp_cross_tab(combined_df)
            if not fp_ct_df.is_empty():
                fp_dir = output_dir / "stats" / "threshold"
                fp_dir.mkdir(parents=True, exist_ok=True)
                write_tsv(fp_ct_df, str(fp_dir / "fp_cross_tab.tsv"))
                print(f"  FP cross-tab: {fp_dir / 'fp_cross_tab.tsv'}")
                figs.append(plot_fp_cross_tab_heatmap(fp_ct_df, str(output_dir)))
        except Exception as e:
            print(f"  WARNING: FP cross-tabulation failed: {e}")

    # ── Somatic Modality Sub-Classification (Section 9) ──────────────────────
    if wise_names is None or "tier" in wise_names:
        try:
            modality_df = compute_somatic_modality(combined_df)
            if not modality_df.is_empty():
                modality_dir = output_dir / "stats" / "tier"
                modality_dir.mkdir(parents=True, exist_ok=True)
                write_tsv(modality_df, str(modality_dir / "somatic_modality.tsv"))
                print(f"  Somatic modality: {modality_dir / 'somatic_modality.tsv'}")
                figs.append(plot_somatic_modality_pie(modality_df, str(output_dir)))
                figs.append(plot_somatic_modality_bars(modality_df, str(output_dir)))
        except Exception as e:
            print(f"  WARNING: Somatic modality sub-classification failed: {e}")

    # ── Rescue Analytics (Section 10) ─────────────────────────────────────────
    try:
        rescue_dir = output_dir / "stats" / "rescue"
        rescue_dir.mkdir(parents=True, exist_ok=True)
        rescue_plot_dir = output_dir / "plots" / "rescue"
        rescue_plot_dir.mkdir(parents=True, exist_ok=True)

        # 48: Rescue breakdown by set
        breakdown_df = compute_rescue_breakdown(combined_df)
        if not breakdown_df.is_empty():
            write_tsv(breakdown_df, str(rescue_dir / "rescue_breakdown.tsv"))
            figs.append(plot_rescue_breakdown(breakdown_df, str(rescue_plot_dir)))
            figs.append(plot_rescue_rate_trend(breakdown_df, str(rescue_plot_dir)))
            print(f"  Rescue breakdown: {rescue_dir / 'rescue_breakdown.tsv'}")

        # 49: Rescue by FILTER
        rescue_filter_df = compute_rescue_by_filter(combined_df)
        if not rescue_filter_df.is_empty():
            write_tsv(rescue_filter_df, str(rescue_dir / "rescue_by_filter.tsv"))
            figs.append(plot_rescue_by_filter(rescue_filter_df, str(rescue_plot_dir)))

        # 50: Rescue cross-tab (RESCUED × FILTER × set)
        rescue_ct_df = compute_rescue_cross_tab(combined_df)
        if not rescue_ct_df.is_empty():
            write_tsv(rescue_ct_df, str(rescue_dir / "rescue_cross_tab.tsv"))
            figs.append(plot_rescue_cross_tab_heatmap(rescue_ct_df, str(rescue_plot_dir)))

        # 51-52: Rescue VAF/DP distributions
        rescue_vaf_dp = compute_rescue_vaf_dp(combined_df)
        if not rescue_vaf_dp.is_empty():
            write_tsv(rescue_vaf_dp, str(rescue_dir / "rescue_vaf_dp.tsv"))
        figs.append(plot_rescue_vaf_boxplot(combined_df, str(rescue_plot_dir)))
        figs.append(plot_rescue_dp_boxplot(combined_df, str(rescue_plot_dir)))

        # 53: Per-sample rescue distribution
        sample_rescue_df = sample_rescue_summary(combined_df)
        if not sample_rescue_df.is_empty():
            write_tsv(sample_rescue_df, str(rescue_dir / "sample_rescue_summary.tsv"))
            figs.append(plot_rescue_sample_distribution(sample_rescue_df, str(rescue_plot_dir)))

        # 54: Rescue by tier
        rescue_tier_df = compute_rescue_by_tier(combined_df)
        if not rescue_tier_df.is_empty():
            write_tsv(rescue_tier_df, str(rescue_dir / "rescue_by_tier.tsv"))
            figs.append(plot_rescue_by_tier(rescue_tier_df, str(rescue_plot_dir)))

        # 56: Rescue by caller support
        rescue_caller_df = compute_rescue_by_caller_support(combined_df)
        if not rescue_caller_df.is_empty():
            write_tsv(rescue_caller_df, str(rescue_dir / "rescue_by_caller_support.tsv"))
            figs.append(plot_rescue_caller_support(rescue_caller_df, str(rescue_plot_dir)))

        print(f"  Rescue analytics: {rescue_dir}/ ({len(list(rescue_dir.glob('*.tsv')))} TSVs)")
    except Exception as e:
        print(f"  WARNING: Rescue analytics failed: {e}")

    # ── Multi-Allelic Visualization (Section 2) ───────────────────────────────
    try:
        multi_allelic_dir = output_dir / "plots" / "multi_allelic"
        multi_allelic_dir.mkdir(parents=True, exist_ok=True)
        figs.append(plot_multiallelic_classification(combined_df, str(multi_allelic_dir)))
        figs.append(plot_allele_balance_scatter(combined_df, str(multi_allelic_dir)))
        figs.append(plot_vaf_sum_histogram(combined_df, str(multi_allelic_dir)))
        figs.append(plot_category_conflict_summary(combined_df, str(multi_allelic_dir)))
    except Exception as e:
        print(f"  WARNING: Multi-allelic visualization failed: {e}")

    # ── BAM Visualization (Section 5) ─────────────────────────────────────────
    try:
        bam_plot_dir = output_dir / "plots" / "bam"
        bam_plot_dir.mkdir(parents=True, exist_ok=True)
        figs.append(plot_bam_metrics_sample_wise(bam_stats_df, str(bam_plot_dir)))
        figs.append(plot_bam_coverage_distribution(bam_stats_df, str(bam_plot_dir)))
    except Exception as e:
        print(f"  WARNING: BAM visualization failed: {e}")

    gc.collect()
    _malloc_trim()
    _mem("after aggregation cleanup")

    # Validation — read per-sample parquet files one at a time to keep memory low
    if not args.no_validate:
        print("Running rescue VCF validation...")
        import glob as _glob
        report_rows = []
        bam_rows = []
        sid_to_set = {}
        for parquet_path in sorted(_glob.glob(str(variant_dir / "*_variants.parquet"))):
            sid = os.path.basename(parquet_path).replace("_variants.parquet", "")
            df = pl.read_parquet(parquet_path)
            report_rows.extend(validate_sample(df, sid, args.tolerance))
            bam_rows.append(validate_bam_one(df, sid))
            if "set_number" in df.columns:
                sn = df["set_number"].drop_nulls().unique()
                if len(sn) > 0:
                    sid_to_set[sid] = str(sn[0])
            del df
        report = pl.DataFrame(report_rows) if report_rows else pl.DataFrame()
        # Join set_number to validation report for heatmap faceting
        if not report.is_empty() and sid_to_set:
            set_df = pl.DataFrame({"sample_id": list(sid_to_set.keys()),
                                   "set_number": list(sid_to_set.values())})
            report = report.join(set_df, on="sample_id", how="left")
        if not report.is_empty():
            write_tsv(report, str(output_dir / "rescue_validation_report.tsv"))
            summary = validation_summary(report)
            if not summary.is_empty():
                write_tsv(summary, str(output_dir / "rescue_validation_summary.tsv"))
                print(f"  Validation report: {output_dir / 'rescue_validation_report.tsv'}")
        bam_report = pl.DataFrame(bam_rows) if bam_rows else pl.DataFrame()
        if not bam_report.is_empty():
            write_tsv(bam_report, str(output_dir / "bam_validation.tsv"))
            print(f"  BAM validation: {output_dir / 'bam_validation.tsv'}")
    else:
        report = None

    # ═══════════════════════════════════════════════════════════════════════════
    # Visualizations — per-wise chart generation
    # ═══════════════════════════════════════════════════════════════════════════
    print("Generating visualizations...")

    # Activate chart theme if requested (Section 10)
    if args.theme == "publishing":
        import altair as alt
        alt.themes.enable("publishing")
        print("  Using publishing theme for charts")

    # Determine which wises to generate (from --wise flag)
    all_wise_names = ["set", "disease", "sample", "tier", "caller", "chromosome", "variant-category"]
    if args.wise is not None:
        active_wises = [w for w in args.wise if w in all_wise_names] if args.wise else all_wise_names
    else:
        active_wises = all_wise_names

    # Wise chart registry: wise_name → [(chart_fn, kwargs), ...]
    # Each chart function receives the wise-specific plot directory as output_dir.
    _WISE_CHART_REGISTRY = {
        "set": [
            (plot_vc_distribution, {"group_col": "set_number"}),
            (plot_caller_overlap, {"group_col": "set_number"}),
            (plot_vaf_distribution, {"color_col": "set_number"}),
            (plot_dna_vs_rna_vaf, {"color_col": "set_number"}),
            (plot_dna_vs_rna_dp, {"group_col": "set_number"}),
            (plot_gt_concordance, {"group_col": "set_number"}),
            (plot_cosmic_gnomad_annotation, {"group_col": "set_number"}),
            (plot_variant_type_distribution, {"group_col": "set_number"}),
            (plot_ti_tv_ratio, {"group_col": "set_number"}),
            (plot_cross_modality, {"group_col": "set_number"}),
            (plot_ref_alt_dp_scatter, {"group_col": "set_number"}),
            (plot_filter_distribution, {"group_col": "set_number"}),
            (plot_redi_evidence, {"group_col": "set_number"}),
            (plot_caller_agreement_matrix, {}),
            (plot_tier_quality_distribution, {}),
            (plot_caller_concordance_vs_vaf, {"color_col": "set_number"}),
            (plot_dna_vs_rna_per_caller, {"color_col": "set_number"}),
        ],
        "disease": [
            (plot_vc_distribution, {"group_col": "disease_normalized"}),
            (plot_vaf_distribution, {"color_col": "disease_normalized"}),
            (plot_dna_vs_rna_vaf, {"color_col": "disease_normalized"}),
            (plot_dna_vs_rna_dp, {"group_col": "disease_normalized"}),
            (plot_gt_concordance, {"group_col": "disease_normalized"}),
            (plot_cosmic_gnomad_annotation, {"group_col": "disease_normalized"}),
            (plot_variant_type_distribution, {"group_col": "disease_normalized"}),
            (plot_ti_tv_ratio, {"group_col": "disease_normalized"}),
            (plot_cross_modality, {"group_col": "disease_normalized"}),
            (plot_ref_alt_dp_scatter, {"group_col": "disease_normalized"}),
            (plot_filter_distribution, {"group_col": "disease_normalized"}),
            (plot_redi_evidence, {"group_col": "disease_normalized"}),
            (plot_caller_agreement_matrix, {}),
            (plot_tier_quality_distribution, {}),
            (plot_caller_concordance_vs_vaf, {"color_col": "disease_normalized"}),
            (plot_dna_vs_rna_per_caller, {"color_col": "disease_normalized"}),
        ],
        "sample": [
            (plot_vc_distribution, {"group_col": "sample_id", "facet_col": "set_number"}),
            (plot_ti_tv_ratio, {"group_col": "sample_id", "facet_col": "set_number"}),
            (plot_cross_modality, {"group_col": "sample_id", "facet_col": "set_number"}),
            (plot_cosmic_gnomad_annotation, {"group_col": "sample_id", "facet_col": "set_number"}),
            (plot_variant_type_distribution, {"group_col": "sample_id", "facet_col": "set_number"}),
        ],
        "tier": [
            (plot_vc_distribution, {"group_col": "final_tier"}),
            (plot_variant_type_distribution, {"group_col": "final_tier"}),
            (plot_ti_tv_ratio, {"group_col": "final_tier"}),
            (plot_cross_modality, {"group_col": "final_tier"}),
            (plot_filter_distribution, {"group_col": "final_tier"}),
            (plot_redi_evidence, {"group_col": "final_tier"}),
            (plot_cosmic_gnomad_annotation, {"group_col": "final_tier"}),
            (plot_vaf_boxplot_per_tier, {}),
            (plot_dp_boxplot_per_tier, {}),
            (plot_gt_concordance_per_tier, {}),
            (plot_tiered_caller_overlap, {}),
            (plot_tiered_variant_types, {}),
            (plot_per_tier_vaf_boxplot, {}),
            (plot_per_tier_dp_boxplot, {}),
            (plot_database_enrichment_by_tier, {}),
            (plot_caller_concordance_vs_vaf, {"color_col": "final_tier"}),
        ],
        "caller": [
            (plot_vaf_distribution, {}),
            (plot_dp_distribution, {}),
            (plot_caller_agreement_matrix, {}),
            (plot_dna_vs_rna_per_caller, {}),
        ],
        "chromosome": [
            (plot_chromosome_density, {}),
            (plot_vc_distribution, {"group_col": "CHROM"}),
            (plot_variant_type_distribution, {"group_col": "CHROM"}),
            (plot_ti_tv_ratio, {"group_col": "CHROM"}),
            (plot_cross_modality, {"group_col": "CHROM"}),
            (plot_filter_distribution, {"group_col": "CHROM"}),
            (plot_cosmic_gnomad_annotation, {"group_col": "CHROM"}),
        ],
        "variant-category": [
            (plot_vc_distribution, {"group_col": "FILTER"}),
            (plot_variant_type_distribution, {"group_col": "FILTER"}),
            (plot_ti_tv_ratio, {"group_col": "FILTER"}),
            (plot_cross_modality, {"group_col": "FILTER"}),
            (plot_redi_evidence, {"group_col": "FILTER"}),
            (plot_cosmic_gnomad_annotation, {"group_col": "FILTER"}),
            (plot_dna_vs_rna_per_caller, {"color_col": "FILTER"}),
            # NOTE: plot_vaf_distribution and plot_caller_concordance_vs_vaf excluded —
            # they use color_col which triggers .facet() on mark_boxplot (composite mark),
            # which Altair cannot facet ("data must be specified at the top level").
            (plot_caller_agreement_matrix, {}),
            (plot_tier_quality_distribution, {}),
        ],
    }

    # Global charts (not per-wise — plotted once into top-level plots/)
    global_charts = []
    global_charts.append(plot_chromosome_density(combined_df, str(output_dir)))
    global_charts.append(plot_tier_quality_distribution(combined_df, str(output_dir)))
    global_charts.append(plot_caller_agreement_matrix(combined_df, str(output_dir)))

    # Per-wise chart generation
    for wise_name in active_wises:
        wise_plot_dir = output_dir / "plots" / wise_name
        wise_plot_dir.mkdir(parents=True, exist_ok=True)
        chart_entries = _WISE_CHART_REGISTRY.get(wise_name, [])
        for chart_fn, kwargs in chart_entries:
            try:
                result = chart_fn(combined_df, str(wise_plot_dir), **kwargs)
                if result is not None:
                    figs.append(result)
            except Exception as e:
                print(f"  WARNING: {chart_fn.__name__} ({wise_name}) failed: {e}")

    # Threshold charts (global, pre-computed data)
    figs.append(plot_caller_concordance_vs_vaf(combined_df, str(output_dir)))
    figs.append(plot_database_enrichment_by_tier(combined_df, str(output_dir)))

    # BAM charts (sample-wise only)
    if not args.no_bam and not bam_stats_df.is_empty():
        figs.append(plot_bam_metrics_bars(bam_stats_df, str(output_dir)))
    if not args.no_pileup:
        figs.append(plot_bam_coverage_violin(combined_df, str(output_dir)))
        figs.append(plot_bam_dp_distribution(combined_df, str(output_dir)))

    # Per-sample charts (use sample_stats_df, not combined_df)
    if all_stats:
        sample_stats_df = pl.DataFrame(all_stats)
        figs.append(plot_per_sample_distribution(sample_stats_df, str(output_dir)))
        sample_tier_df = sample_tier_summary(combined_df)
        if not sample_tier_df.is_empty():
            figs.append(plot_per_sample_tier_distribution(sample_tier_df, str(output_dir)))
            figs.append(plot_per_tier_cross_sample_vaf(sample_tier_df, str(output_dir)))

    # Validation heatmap (global)
    if not args.no_validate and report is not None:
        figs.append(plot_validation_heatmap(report, str(output_dir)))

    # Master dashboard
    generate_dashboard(figs, str(output_dir))

    # ── variant_details_filtered/ removed ─────────────────────────────────
    # The separate variant_details_filtered/ output directory is removed.
    # The unified filter pipeline (build_unified_filter) applies all CLI
    # filters (--include-filters, --exclude-filters, --min-vaf, --min-dp,
    # --max-gnomad-af, --min-evidence-tier, --exclude-multiallelic-conflict,
    # --exclude-disease) directly to combined_df before any statistics or
    # charts, so ALL output (TSVs, charts, summaries) reflects the filtered
    # view. No separate filtered parquet copy is needed.

    print(f"\nAll outputs written to: {output_dir}")
    print("Done.")


if __name__ == "__main__":
    main()
