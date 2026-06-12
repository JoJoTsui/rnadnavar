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

from .bam_stats import compute_all_bam_stats
from .bam_validation import validate_bam_vs_caller as validate_bam_one
from .caller_parser import _parse_one_caller, CALLERS_STRELKA, CALLERS_WITH_GT
from .manifest_loader import CALLER_CONFIGS, filter_complete, load_manifest
from .rescue_parser import parse_rescue_vcf as _py_parse_rescue
from .rust_vcf import parse_rescue_vcf as _rust_parse_rescue
from .rescue_validator import validate_all_samples, validate_sample, validation_summary
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


def _streaming_join_one(df: pl.DataFrame, col_data: dict, caller_name: str) -> pl.DataFrame:
    """Join a single caller's column-oriented data and rename columns."""
    from .caller_parser import CALLERS_STRELKA, CALLERS_WITH_GT
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
        for suf, src in [("AD_REF", "TOR"), ("AD_ALT", "TAR")]:
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



def process_single_sample(row: dict, max_workers: int = 1, use_rust: bool = True, large_sem=None) -> dict:
    """Process one sample: parse rescue VCF + all caller VCFs + compute stats.

    Returns a dict with 'sample_id', 'df', and 'stats'.

    Args:
        row: Manifest row dict.
        max_workers: Threads for within-sample caller parsing.
        use_rust: Use Rust VCF parser.
        large_sem: Semaphore for large-sample exclusive access (threading.Semaphore
                   or multiprocessing.Manager.Semaphore proxy). If None, large
                   samples are not throttled.
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
            name, cols = _parse_one_caller(caller_name, cfg, base, vcf_prefix, target_positions)
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
        args: (row_dict, max_workers, use_rust, variant_dir_str)

    Returns:
        {"sample_id": str, "stats": dict} — no DataFrame (too large to pickle).
        On error: {"sample_id": str, "error": str}
    """
    row, max_workers, use_rust, variant_dir_str = args

    try:
        result = process_single_sample(row, max_workers=max_workers, use_rust=use_rust, large_sem=None)
    except Exception:
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
    parser.add_argument("--verbose", action="store_true",
                        help="Show per-caller progress messages")
    parser.add_argument("--process-mode", choices=["thread", "spawn"], default=None,
                        help="Parallel execution mode: thread (ThreadPoolExecutor) or spawn "
                             "(multiprocessing.Pool with fresh Python process per worker). "
                             "Default: spawn when --sample-workers > 1, thread otherwise.")
    parser.add_argument("--max-tasks-per-child", type=int, default=1,
                        help="Max samples per worker process before restart (spawn mode only, default: 1).")
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

    # Determine process mode (used by print below and execution logic below)
    process_mode = args.process_mode
    if process_mode is None:
        process_mode = "spawn" if args.sample_workers > 1 else "thread"

    print(f"Processing {len(manifest)} samples (parser={args.parser}, caller_threads={args.threads}, sample_workers={args.sample_workers}, bam_workers={args.bam_workers}, process_mode={process_mode})")
    if len(manifest) > 20 and not args.no_validate:
        print("NOTE: >20 samples with validation enabled may be slow due to BAM pileup.")
        print("      Consider --no-validate for initial runs, then validate separately.")
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    rows = manifest.to_dicts()
    all_stats = []
    bam_stats_df = pl.DataFrame()

    # ── Per-sample parquet directory (streaming, not memory-accumulated) ──
    variant_dir = output_dir / "variant_details"
    variant_dir.mkdir(parents=True, exist_ok=True)
    variant_dir_str = str(variant_dir)
    total_variants = 0

    if process_mode == "spawn" and args.sample_workers > 1:
        # ── Process-isolated parallel mode (spawn) ──────────────────────────
        # spawn creates a fresh Python interpreter per worker (fork+exec),
        # avoiding the fork+threads deadlock with polars' rayon pool.
        # maxtasksperchild=1 ensures each worker exits after one sample
        # → kernel reclaims all memory. No cross-process throttle needed —
        # process isolation alone keeps total RSS within limits.
        import multiprocessing as mp
        ctx = mp.get_context("spawn")

        worker_args = [
            (row, args.threads, use_rust, variant_dir_str)
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
                large_sem=_THREAD_LARGE_SEM,
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
                    large_sem=_THREAD_LARGE_SEM,
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
    combined_df = pl.scan_parquet(str(variant_dir / "*_variants.parquet"))
    print(f"Variant details: {variant_dir}/ (lazy scan, {total_variants} variants across {len(rows)} samples)")

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
        _mem("before disease_summary")
        disease_summary_df = disease_summary(combined_df)
        if not disease_summary_df.is_empty():
            disease_summary_df.write_csv(str(output_dir / "disease_summary.csv"))
        _mem("after disease_summary")

        # Tier summary
        tier_summary_df = compute_tier_summary(combined_df)
        if not tier_summary_df.is_empty():
            tier_summary_df.write_csv(str(output_dir / "tier_summary.csv"))
            print(f"Tier summary: {output_dir / 'tier_summary.csv'}")
        _mem("after tier_summary")

        # Dataset summary (whole-dataset aggregates)
        ds_summary = dataset_summary(combined_df)
        if ds_summary:
            pl.DataFrame([ds_summary]).write_csv(str(output_dir / "dataset_summary.csv"))
            print(f"Dataset summary: {output_dir / 'dataset_summary.csv'}")
        _mem("after dataset_summary")

        # Sample-tier summary (Level 4)
        sample_tier_df = sample_tier_summary(combined_df)
        if not sample_tier_df.is_empty():
            sample_tier_df.write_csv(str(output_dir / "sample_tier_summary.csv"))
            print(f"Sample-tier summary: {output_dir / 'sample_tier_summary.csv'}")
        _mem("after sample_tier_summary")

    # Caller overlap
    from .statistics import caller_overlap_distribution
    overlap_df = caller_overlap_distribution(combined_df)
    if not overlap_df.is_empty():
        overlap_df.write_csv(str(output_dir / "caller_overlap.csv"))
    _mem("after caller_overlap")

    # Filter distribution
    from .statistics import filter_distribution, variant_type_distribution, vc_distribution
    filter_df = filter_distribution(combined_df)
    if not filter_df.is_empty():
        filter_df.write_csv(str(output_dir / "filter_distribution.csv"))

    # Variant type distribution
    vt_df = variant_type_distribution(combined_df)
    if not vt_df.is_empty():
        vt_df.write_csv(str(output_dir / "variant_type_distribution.csv"))
    _mem("after filter/vt distributions")

    # GT concordance
    concordance_data = gt_concordance(combined_df)
    if concordance_data:
        pl.DataFrame({"agreement_level": list(concordance_data.keys()), "count": list(concordance_data.values())}).write_csv(
            str(output_dir / "gt_concordance.csv")
        )
    _mem("after gt_concordance")
    gc.collect()
    _malloc_trim()
    _mem("after aggregation cleanup")

    # Validation — read per-sample parquet files one at a time to keep memory low
    if not args.no_validate:
        print("Running rescue VCF validation...")
        import glob as _glob
        report_rows = []
        bam_rows = []
        for parquet_path in sorted(_glob.glob(str(variant_dir / "*_variants.parquet"))):
            sid = os.path.basename(parquet_path).replace("_variants.parquet", "")
            df = pl.read_parquet(parquet_path)
            report_rows.extend(validate_sample(df, sid, args.tolerance))
            bam_rows.append(validate_bam_one(df, sid))
            del df
        report = pl.DataFrame(report_rows) if report_rows else pl.DataFrame()
        if not report.is_empty():
            report.write_csv(str(output_dir / "rescue_validation_report.csv"))
            summary = validation_summary(report)
            if not summary.is_empty():
                summary.write_csv(str(output_dir / "rescue_validation_summary.csv"))
                print(f"  Validation report: {output_dir / 'rescue_validation_report.csv'}")
        bam_report = pl.DataFrame(bam_rows) if bam_rows else pl.DataFrame()
        if not bam_report.is_empty():
            bam_report.write_csv(str(output_dir / "bam_validation.csv"))
            print(f"  BAM validation: {output_dir / 'bam_validation.csv'}")
    else:
        report = None

    # Visualizations
    print("Generating visualizations...")
    figs = []

    # Pass lazy frame directly to chart functions.
    # Each chart function calls _eager(df) internally, which triggers
    # polars scan + collect. polars' query optimizer ensures that only
    # the columns needed by each chart are read from the parquet files.
    # Scatter charts sample internally to 5-10K rows.

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
