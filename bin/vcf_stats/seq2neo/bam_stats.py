"""Per-sample per-modality BAM statistics.

Computes read-level statistics from alignment BAM files for both DNA and RNA
modalities. Uses Rust stats_core (noodles-bam) when available, falls back to
pysam. Parallelizes across samples via ThreadPoolExecutor.
"""

import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any

import polars as pl


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


def _get_bam_ref_total(bam_path: str) -> int:
    """Get the total reference length from a BAM header (pysam, header-only).

    Fast — reads only the header, not the alignment records.
    Used to post-hoc recalculate mean_coverage with a custom denominator
    (e.g., WES BED region total) without needing the Rust module rebuilt.
    """
    try:
        import pysam
        bam = pysam.AlignmentFile(bam_path, "rb")
        total = sum(bam.lengths) if bam.lengths else 0
        bam.close()
        return total
    except Exception:
        return 0

try:
    import stats_core
    HAS_RUST_BAM = hasattr(stats_core, 'bam_stats')
except ImportError:
    HAS_RUST_BAM = False

try:
    import pysam
    HAS_PYSAM = True
except ImportError:
    HAS_PYSAM = False


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


def _compute_bam_stats_pysam(bam_path: str, bed_total: int = 0,
                             bed_regions: list[tuple[str, int, int]] | None = None) -> dict[str, Any] | None:
    """Compute BAM statistics using pysam (Python fallback)."""
    if not HAS_PYSAM:
        return None

    try:
        bam = pysam.AlignmentFile(bam_path, "rb")

        total_reads = 0
        mapped_reads = 0
        total_mapq = 0
        total_insert = 0
        insert_count = 0
        total_length = 0
        on_target_bases = 0

        # Reference lengths for coverage estimation — use BED total for WES,
        # otherwise use BAM header reference lengths (whole-genome).
        ref_lengths = bed_total if bed_total > 0 else (sum(bam.lengths) if bam.lengths else 1)

        # Build per-chromosome BED interval index for on-target check
        bed_index: dict[str, list[tuple[int, int]]] = {}
        if bed_regions:
            for chrom, start, end in bed_regions:
                bed_index.setdefault(chrom, []).append((start, end))

        for read in bam.fetch():
            total_reads += 1
            if not read.is_unmapped:
                mapped_reads += 1
                total_mapq += read.mapping_quality
                read_len = read.query_length or 0
                total_length += read_len
                # Insert size: only count properly paired reads (TLEN can be
                # arbitrarily large for supplementary/improper pairs)
                if (read.is_proper_pair and not read.is_supplementary
                    and not read.is_secondary
                    and read.template_length and read.template_length > 0):
                    total_insert += read.template_length
                    insert_count += 1

                # On-target check for WES coverage accuracy
                if bed_index:
                    chrom = read.reference_name
                    if chrom in bed_index:
                        read_start = read.reference_start
                        read_end = read.reference_end or (read_start + read_len)
                        for int_start, int_end in bed_index[chrom]:
                            if read_start < int_end and read_end > int_start:
                                overlap_start = max(read_start, int_start)
                                overlap_end = min(read_end, int_end)
                                if overlap_end > overlap_start:
                                    on_target_bases += overlap_end - overlap_start
                                break  # one interval is enough per read

        bam.close()

        if total_reads == 0:
            return None

        mapping_rate = mapped_reads / total_reads * 100 if total_reads > 0 else 0
        mean_mapq = total_mapq / mapped_reads if mapped_reads > 0 else 0
        mean_insert = total_insert / insert_count if insert_count > 0 else 0
        coverage_bases = on_target_bases if bed_regions else total_length
        mean_coverage = coverage_bases / ref_lengths if coverage_bases > 0 else 0

        return {
            "total_reads": total_reads,
            "mapped_reads": mapped_reads,
            "mapping_rate_pct": round(mapping_rate, 2),
            "mean_coverage": round(mean_coverage, 4),
            "mean_insert_size": round(mean_insert, 1),
            "mean_mapq": round(mean_mapq, 1),
        }
    except Exception as e:
        print(f"  [BAM STATS] Error reading {bam_path}: {e}")
        return None


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
        }
    except Exception as e:
        print(f"  [BAM STATS] Rust error on {bam_path}: {e}")
        return None


def compute_bam_stats(bam_path: str, bed_total: int = 0,
                     bed_regions: list[tuple[str, int, int]] | None = None) -> dict[str, Any] | None:
    """Compute basic statistics from a BAM file.

    Uses Rust stats_core when available (faster), falls back to pysam.

    Args:
        bam_path: Path to the BAM file.
        bed_total: Total length of BED regions (for WES coverage denominator).
                   When > 0, used as coverage denominator instead of whole-genome
                   reference lengths. Default 0 = whole-genome.
        bed_regions: Merged BED regions [(chrom, start, end), ...] for on-target
                     coverage. When provided, only bases overlapping BED regions
                     are counted toward coverage. Default None = whole-genome.

    Returns dict with: total_reads, mapped_reads, mapping_rate_pct,
    mean_coverage, mean_insert_size, mean_mapq. Returns None if BAM is unreadable.
    """
    if not bam_path or not os.path.isfile(bam_path):
        return None

    if HAS_RUST_BAM:
        has_bed = bool(bed_regions)
        if has_bed:
            print(f"  [BAM STATS] Using Rust bam_stats_bed ({len(bed_regions)} regions, bed_total={bed_total})")
        result = _compute_bam_stats_rust(bam_path, bed_regions)
        if result is not None:
            # When bed_regions provided, Rust already computes on-target coverage.
            # No post-hoc recalculation needed.
            if not bed_regions and bed_total > 0:
                # Legacy path: no BED regions but BED total provided.
                # Recalculate using WES denominator (kept for backwards compat).
                bam_ref = _get_bam_ref_total(bam_path)
                if bam_ref > 0 and result.get("mean_coverage"):
                    result["mean_coverage"] = round(
                        result["mean_coverage"] * (bam_ref / bed_total), 4
                    )
            # Diagnostic: show path taken + coverage value
            mode = "BED" if bed_regions else "WG"
            cov = result.get("mean_coverage", "N/A")
            print(f"  [BAM STATS] mode={mode} mean_coverage={cov}")
            return result
        # Fall through to pysam on Rust failure
        print(f"  [BAM STATS] Rust failed for {bam_path}, falling back to pysam")

    result = _compute_bam_stats_pysam(bam_path, bed_total, bed_regions)
    if result:
        mode = "BED" if bed_regions else "WG"
        cov = result.get("mean_coverage", "N/A")
        print(f"  [BAM STATS] mode={mode} (pysam) mean_coverage={cov}")
    return result
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
            })

        results.append(row)

    return results


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
                     coverage calculation. Passed to Rust/pysam for accurate WES
                     coverage. Default None = whole-genome.

    Returns:
        polars DataFrame with one row per sample per modality.
    """
    all_rows = []

    if max_workers > 1 and len(manifest_rows) > 1:
        # Parallel BAM processing via ThreadPoolExecutor (threads, safe with htslib
        # even when using pysam fallback — threads share the htslib state safely)
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
        return pl.DataFrame()

    # Sort by sample_id then bam_type for deterministic output ordering.
    # Without this, as_completed() produces non-deterministic row order.
    all_rows.sort(key=lambda r: (r["sample_id"], str(r["bam_type"])))

    return pl.DataFrame(all_rows)
