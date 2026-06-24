"""Per-sample per-modality BAM statistics.

Computes read-level statistics from alignment BAM files for both DNA and RNA
modalities. Uses Rust stats_core (noodles-bam) exclusively — no pysam fallback.
Parallelizes across samples via ThreadPoolExecutor.
"""

import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any

import polars as pl

# Standard BGZF EOF marker: an empty gzip block (28 bytes) that terminates
# every valid BGZF-compressed file (BAM, VCF.gz, etc.). Files truncated during
# transfer or storage lack this marker.
_BGZF_EOF_MARKER = bytes([
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00,
    0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00,
])


def _check_bam_eof(bam_path: str) -> tuple[bool, str]:
    """Check whether a BAM file has a valid BGZF EOF marker.

    Reads the last 28 bytes of the file and compares against the standard
    BGZF empty-block footer. This is an O(1) check — no parsing required.

    Returns:
        (is_valid, reason) — True if EOF is valid, False with explanation if not.
    """
    if not bam_path or not os.path.isfile(bam_path):
        return False, f"file not found: {bam_path}"
    try:
        file_size = os.path.getsize(bam_path)
        if file_size < 28:
            return False, f"file too small to be valid BAM ({file_size} bytes): {bam_path}"
        with open(bam_path, "rb") as f:
            f.seek(-28, os.SEEK_END)
            tail = f.read(28)
        if tail == _BGZF_EOF_MARKER:
            return True, ""
        # Also accepts BGZF EOF blocks where the deflate stored-block header
        # uses ISIZE=0x00000000 instead of 0x0003 — the first 16 bytes
        # (gzip header + BGZF extra subfield) are identical in both variants.
        if tail[:16] == _BGZF_EOF_MARKER[:16]:
            return True, ""
        return False, f"no valid BGZF EOF marker (file may be truncated): {bam_path}"
    except OSError as e:
        return False, f"cannot read BAM file: {bam_path}: {e}"


def read_and_merge_bed(bed_path: str, gap: int = 100_000) -> tuple[int, list[tuple[str, int, int]], list[tuple[str, int, int]]]:
    """Read BED file, merge adjacent intervals, return raw/merged totals + regions.

    Reads a BED file (0-based, 3+ column format: chrom, start, end, ...),
    sorts intervals by (chromosome, start), and merges intervals on the same
    chromosome when they are within `gap` bp of each other.

    Merging within gap reduces ~50K raw target intervals to ~300 contiguous
    regions for WES — critical for efficient BAM pileup region-guided queries.

    For BAM stats coverage, raw (unmerged) intervals MUST be used as the
    on-target denominator. Using merged intervals inflates bed_total (e.g.,
    170 Mbp → 2,780 Mbp with gap=500Kb), making on-target coverage
    indistinguishable from whole-genome.

    Returns:
        (raw_bed_total, raw_bed_regions, merged_bed_regions)
        - raw_bed_total: sum of all unmerged interval lengths (coverage denominator)
        - raw_bed_regions: list of (chrom, start, end) before merging (on-target counting)
        - merged_bed_regions: list of (chrom, start, end) after merging (pileup queries)
    """
    intervals = []
    with open(bed_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            try:
                start = int(parts[1])
                end = int(parts[2])
                if end > start:
                    intervals.append((parts[0], start, end))
            except (ValueError, IndexError):
                continue

    if not intervals:
        return 0, [], []

    # Raw total: sum of all unmerged intervals (correct coverage denominator)
    raw_total = sum(end - start for _, start, end in intervals)

    # Sort by chromosome then start position
    intervals.sort(key=lambda x: (x[0], x[1]))

    # Merge adjacent intervals within gap (for efficient pileup region queries)
    merged = []
    for chrom, start, end in intervals:
        if (merged and merged[-1][0] == chrom
            and start - merged[-1][2] <= gap):
            merged[-1] = (chrom, merged[-1][1], max(merged[-1][2], end))
        else:
            merged.append((chrom, start, end))

    return raw_total, intervals, merged


try:
    import stats_core
    if not hasattr(stats_core, 'bam_stats'):
        raise ImportError("stats_core.bam_stats not found")
except ImportError:
    raise ImportError(
        "Rust BAM stats backend (stats_core.bam_stats) is required but could not be loaded. "
        "Build the Rust module: cd bin/vcf_stats/seq2neo/stats_core && ./build_rust.sh"
    )


# Three BAM types per sample
BAM_TYPES = {
    "DN": {"suffix": "DN", "label": "DNA Normal"},
    "DT": {"suffix": "DT", "label": "DNA Tumor"},
    "RT": {"suffix": "RT", "label": "RNA Tumor"},
}


def _locate_bam_file(base_dir: str, dir_name: str, bam_type: str) -> str | None:
    """Locate the BAM file for a given BAM type (DN, DT, or RT).

    DN BAMs: preprocessing/mapped/{prefix}DN/{prefix}DN.sorted.bam
    DT BAMs: preprocessing/mapped/{prefix}DT/{prefix}DT.sorted.bam
    RT BAMs: preprocessing/mapped/{prefix}RT/{prefix}RT.bam
             OR vcf_realignment/preprocessing/mapped/{prefix}RT/*.bam
    """
    import glob

    suffix = BAM_TYPES[bam_type]["suffix"]

    if bam_type in ("DN", "DT"):
        # DNA normal and tumor: preprocessing/mapped/{prefix}DN/{prefix}DN.sorted.bam
        search_paths = [
            os.path.join(base_dir, dir_name, "preprocessing", "mapped", f"*{suffix}", "*.sorted.bam"),
            os.path.join(base_dir, dir_name, "preprocessing", "mapped", f"*{suffix}", "*.bam"),
        ]
    else:  # RT
        # RNA tumor: preprocessing/mapped/{prefix}RT/{prefix}RT.bam
        # or vcf_realignment/preprocessing/mapped/{prefix}RT/*.bam
        search_paths = [
            os.path.join(base_dir, dir_name, "preprocessing", "mapped", f"*{suffix}", "*.bam"),
            os.path.join(base_dir, dir_name, "vcf_realignment", "preprocessing", "mapped", f"*{suffix}", "*.bam"),
        ]

    for pattern in search_paths:
        files = glob.glob(pattern, recursive=True)
        if files:
            return files[0]

    return None


# ---------------------------------------------------------------------------
# Section 4: BAM Alignment Stats (Rust-only, expanded metrics from Rust)
# ---------------------------------------------------------------------------


def _compute_bam_stats_rust(bam_path: str,
                            bed_regions: list[tuple[str, int, int]] | None = None) -> dict[str, Any] | None:
    """Compute BAM statistics using Rust stats_core (noodles-bam).

    Passes max_reads=0 to read the entire BAM file without sampling.
    When bed_regions is provided, passes them to Rust for on-target coverage
    calculation (WES mode).
    """
    try:
        if bed_regions:
            # Separate into three parallel lists for FFI
            bed_chroms = [r[0] for r in bed_regions]
            bed_starts = [r[1] for r in bed_regions]
            bed_ends = [r[2] for r in bed_regions]
            raw = stats_core.bam_stats_bed(bam_path, 0, bed_chroms, bed_starts, bed_ends)
        else:
            raw = stats_core.bam_stats(bam_path, 0)  # max_reads=0 → no limit
        return {
            "total_reads": raw["total_reads"],
            "mapped_reads": raw["mapped_reads"],
            "mapping_rate_pct": round(raw["mapping_rate"], 2),
            "mean_coverage": round(raw["mean_coverage"], 4),
            "mean_insert_size": round(raw["mean_insert_size"], 1),
            "mean_mapq": round(raw["mean_mapq"], 1),
            "duplication_rate_pct": round(raw["duplication_rate_pct"], 2),
            "properly_paired_pct": round(raw["properly_paired_pct"], 2),
            "insert_size_stddev": round(raw["insert_size_stddev"], 1),
        }
    except Exception as e:
        print(f"  [BAM STATS] Rust error on {bam_path}: {e}")
        return None


def compute_bam_stats(bam_path: str, bed_total: int = 0,
                     bed_regions: list[tuple[str, int, int]] | None = None) -> dict[str, Any] | None:
    """Compute statistics from a BAM file using the Rust stats_core backend.

    Args:
        bam_path: Path to the BAM file.
        bed_total: Total length of BED regions (for WES coverage denominator).
                   When > 0, used as coverage denominator instead of whole-genome
                   reference lengths. Default 0 = whole-genome.
        bed_regions: BED regions [(chrom, start, end), ...] for on-target
                     coverage. When provided, only bases overlapping BED regions
                     are counted toward coverage. Default None = whole-genome.

    Returns dict with: total_reads, mapped_reads, mapping_rate_pct,
    mean_coverage, mean_insert_size, mean_mapq,
    duplication_rate_pct, properly_paired_pct, insert_size_stddev,
    cov_1x_pct, cov_10x_pct, cov_20x_pct, cov_50x_pct, cov_100x_pct.
    Returns None if BAM is unreadable.
    """
    if not bam_path or not os.path.isfile(bam_path):
        return None

    has_bed = bool(bed_regions)
    if has_bed:
        print(f"  [BAM STATS] Using Rust bam_stats_bed ({len(bed_regions)} regions, bed_total={bed_total})")
    result = _compute_bam_stats_rust(bam_path, bed_regions)
    if result is None:
        print(f"  [BAM STATS] Rust failed for {bam_path}")
        return None

    # Coverage bins via dedicated Rust function
    if bed_regions:
        try:
            cov_bins = stats_core.coverage_bins(bam_path, bed_regions)
            if cov_bins is not None:
                result.update(cov_bins)
            else:
                result["cov_1x_pct"] = None
                result["cov_10x_pct"] = None
                result["cov_20x_pct"] = None
                result["cov_50x_pct"] = None
                result["cov_100x_pct"] = None
        except TypeError as e:
            print(f"  [BAM STATS] coverage_bins type error for {bam_path}: {e}")
            result["cov_1x_pct"] = None
            result["cov_10x_pct"] = None
            result["cov_20x_pct"] = None
            result["cov_50x_pct"] = None
            result["cov_100x_pct"] = None
        except RuntimeError as e:
            print(f"  [BAM STATS] coverage_bins Rust error for {bam_path}: {e}")
            result["cov_1x_pct"] = None
            result["cov_10x_pct"] = None
            result["cov_20x_pct"] = None
            result["cov_50x_pct"] = None
            result["cov_100x_pct"] = None
        except OSError as e:
            print(f"  [BAM STATS] coverage_bins file error for {bam_path}: {e}")
            result["cov_1x_pct"] = None
            result["cov_10x_pct"] = None
            result["cov_20x_pct"] = None
            result["cov_50x_pct"] = None
            result["cov_100x_pct"] = None
    else:
        result["cov_1x_pct"] = None
        result["cov_10x_pct"] = None
        result["cov_20x_pct"] = None
        result["cov_50x_pct"] = None
        result["cov_100x_pct"] = None

    # Diagnostic: show path taken + coverage value
    mode = "BED" if bed_regions else "WG"
    cov = result.get("mean_coverage", "N/A")
    print(f"  [BAM STATS] mode={mode} mean_coverage={cov}")
    return result


def compute_sample_bam_stats(
    base_output_dir: str,
    dir_name: str,
    sample_id: str,
    set_number: int,
    bed_total: int = 0,
    bed_regions: list[tuple[str, int, int]] | None = None,
) -> list[dict[str, Any]]:
    """Compute BAM statistics for a single sample (DNA + RNA modalities).

    Returns a list of dicts, one per modality, for aggregation into a DataFrame.
    Each dict contains sample metadata and BAM stats. Modalities with missing
    BAM files produce rows with null stats.
    """
    results = []

    for bam_type in ["DN", "DT", "RT"]:
        bam_path = _locate_bam_file(base_output_dir, dir_name, bam_type)
        stats = compute_bam_stats(bam_path, bed_total, bed_regions) if bam_path else None

        row = {
            "sample_id": sample_id,
            "set_number": set_number,
            "bam_type": bam_type,
            "bam_label": BAM_TYPES[bam_type]["label"],
            "bam_path": bam_path or "",
            "has_bam": bam_path is not None,
        }
        if stats:
            row.update(stats)
        else:
            row.update({
                "total_reads": None,
                "mapped_reads": None,
                "mapping_rate_pct": None,
                "mean_coverage": None,
                "mean_insert_size": None,
                "mean_mapq": None,
                "duplication_rate_pct": None,
                "properly_paired_pct": None,
                "insert_size_stddev": None,
                "cov_1x_pct": None,
                "cov_10x_pct": None,
                "cov_20x_pct": None,
                "cov_50x_pct": None,
                "cov_100x_pct": None,
            })

        results.append(row)

    return results


# Expected columns for backward compatibility when reloading from TSV.
_BAM_STATS_COLUMNS = [
    "sample_id", "set_number", "bam_type", "bam_label", "bam_path", "has_bam",
    "total_reads", "mapped_reads", "mapping_rate_pct", "mean_coverage",
    "mean_insert_size", "mean_mapq",
    "duplication_rate_pct", "properly_paired_pct", "insert_size_stddev",
    "cov_1x_pct", "cov_10x_pct", "cov_20x_pct", "cov_50x_pct", "cov_100x_pct",
]


def ensure_bam_stats_columns(df: pl.DataFrame) -> pl.DataFrame:
    """Ensure all expected bam_stats columns exist, null-filling missing ones.

    Provides backward compatibility when reloading bam_stats TSV from --resume
    that was written by an older version lacking the expanded metrics.
    """
    for col in _BAM_STATS_COLUMNS:
        if col not in df.columns:
            df = df.with_columns(pl.lit(None).alias(col))
    return df


def compute_all_bam_stats(manifest_rows: list[dict], max_workers: int = 8, bed_total: int = 0,
                         bed_regions: list[tuple[str, int, int]] | None = None) -> pl.DataFrame:
    """Compute BAM statistics for all samples in the manifest.

    Args:
        manifest_rows: List of manifest row dicts with sample_id, base_output_dir,
            dir_name, set_number.
        max_workers: Number of threads for parallel sample processing (default: 8).
            Uses ThreadPoolExecutor (threads, safe with htslib).
        bed_total: Total length of BED regions for WES coverage denominator.
                   When > 0, passed through to compute_bam_stats for coverage
                   recalculation. Default 0 = whole-genome.
        bed_regions: Merged BED regions [(chrom, start, end), ...] for on-target
                     coverage calculation. Passed to Rust for accurate WES
                     coverage. Default None = whole-genome.

    Returns:
        polars DataFrame with one row per sample per modality.
    """
    all_rows = []

    if max_workers > 1 and len(manifest_rows) > 1:
        # Parallel BAM processing via ThreadPoolExecutor (threads safe with htslib)
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = {}
            for row in manifest_rows:
                future = executor.submit(
                    compute_sample_bam_stats,
                    base_output_dir=row["base_output_dir"],
                    dir_name=row["dir_name"],
                    sample_id=row["sample_id"],
                    set_number=row["set_number"],
                    bed_total=bed_total,
                    bed_regions=bed_regions,
                )
                futures[future] = row["sample_id"]

            for future in as_completed(futures):
                sid = futures[future]
                try:
                    sample_results = future.result()
                    all_rows.extend(sample_results)
                    ok_flags = []
                    for r in sample_results:
                        ok_flags.append(f"{r['bam_type']}={'OK' if r['has_bam'] else 'missing'}")
                    print(f"  [{sid}] BAM stats: {', '.join(ok_flags)}")
                except Exception as e:
                    print(f"  [{sid}] BAM stats ERROR: {e}")
    else:
        # Sequential processing
        for row in manifest_rows:
            sample_results = compute_sample_bam_stats(
                base_output_dir=row["base_output_dir"],
                dir_name=row["dir_name"],
                sample_id=row["sample_id"],
                set_number=row["set_number"],
                bed_total=bed_total,
                bed_regions=bed_regions,
            )
            all_rows.extend(sample_results)
            ok_flags = []
            for r in sample_results:
                ok_flags.append(f"{r['bam_type']}={'OK' if r['has_bam'] else 'missing'}")
            print(f"  [{row['sample_id']}] BAM stats: {', '.join(ok_flags)}")

    if not all_rows:
        df = pl.DataFrame()
        # Ensure column schema even for empty result
        df = ensure_bam_stats_columns(df)
        return df

    # Sort by sample_id then bam_type for deterministic output ordering.
    # Without this, as_completed() produces non-deterministic row order.
    all_rows.sort(key=lambda r: (r["sample_id"], str(r["bam_type"])))

    df = pl.DataFrame(all_rows)
    df = ensure_bam_stats_columns(df)
    return df
