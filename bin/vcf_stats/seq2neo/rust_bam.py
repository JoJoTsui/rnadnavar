"""Python wrapper for Rust BAM pileup via stats_core.

Uses stats_core (Rust) exclusively. No pysam fallback — if Rust fails,
the error propagates so the root cause can be debugged and fixed.
"""

import os
from concurrent.futures import ThreadPoolExecutor, as_completed

import polars as pl

try:
    import stats_core
    HAS_RUST_BAM = hasattr(stats_core, 'pileup_variants')
    HAS_RUST_MULTI = hasattr(stats_core, 'pileup_variants_multi')
except ImportError:
    HAS_RUST_BAM = False
    HAS_RUST_MULTI = False


def _require_rust():
    """Raise ImportError if the Rust stats_core module is not available."""
    if not HAS_RUST_BAM:
        raise ImportError(
            "stats_core.pileup_variants not found. "
            "Build the Rust module: cd bin/vcf_stats/seq2neo/stats_core && ./build_rust.sh"
        )


def pileup_variants(bam_path: str, positions: list[tuple[str, int, str, str]]) -> pl.DataFrame | None:
    """Perform BAM pileup at specific variant positions using Rust stats_core.

    Args:
        bam_path: Path to the BAM file.
        positions: List of (chrom, pos, ref_base, alt_base) tuples.

    Returns:
        DataFrame with columns: CHROM, POS, DP, REF_DP, ALT_DP,
        F1R2_ref, F2R1_ref, F1R2_alt, F2R1_alt, mean_BQ, mean_MQ.
        Returns None if BAM file doesn't exist.

    Raises:
        ImportError: stats_core module not built.
        RuntimeError: Rust pileup failed.
    """
    if not os.path.isfile(bam_path):
        return None

    _require_rust()

    chroms = [p[0] for p in positions]
    poss = [p[1] for p in positions]
    refs = [p[2] for p in positions]
    alts = [p[3] for p in positions]

    try:
        result = stats_core.pileup_variants(bam_path, chroms, poss, refs, alts)
    except BaseException as e:
        if isinstance(e, KeyboardInterrupt):
            raise
        raise RuntimeError(f"Rust BAM pileup failed for {bam_path}: {e}") from e

    if result:
        # Add coordinate columns (Rust returns only metric columns)
        result["CHROM"] = chroms
        result["POS"] = poss
        result["REF"] = refs
        result["ALT"] = alts
        return pl.DataFrame(result)
    return None


def pileup_variants_multi(
    bam_paths: dict[str, str],
    positions: list[tuple[str, int, str, str]],
    bed_regions: list[tuple[str, int, int]] | None = None,
) -> dict[str, pl.DataFrame]:
    """Perform combined BAM pileup across multiple BAMs in a single FFI call.

    Opens all BAMs once, shares position grouping and window iteration,
    and uses binary search to match reads to positions within each region.

    Args:
        bam_paths: Dict mapping BAM label (e.g., "DN", "DT", "RT") to BAM path.
                   Only BAMs with existing files are included.
        positions: List of (chrom, pos, ref_base, alt_base) tuples.
        bed_regions: Optional merged BED regions [(chrom, start, end), ...] for
                     region-guided queries (WES mode).

    Returns:
        Dict mapping BAM label to DataFrame with pileup columns. Returns empty
        dict if no BAMs or no positions provided.

    Raises:
        ImportError: stats_core module not built.
        RuntimeError: Rust multi-BAM pileup failed.
    """
    if not bam_paths or not positions:
        return {}

    _require_rust()

    chroms = [p[0] for p in positions]
    poss = [p[1] for p in positions]
    refs = [p[2] for p in positions]
    alts = [p[3] for p in positions]

    # Preferred path: Rust multi-BAM in a single FFI call
    if HAS_RUST_MULTI:
        bam_labels = list(bam_paths.keys())
        bam_file_paths = [bam_paths[l] for l in bam_labels]

        if bed_regions:
            bed_chroms = [r[0] for r in bed_regions]
            bed_starts = [r[1] for r in bed_regions]
            bed_ends = [r[2] for r in bed_regions]
        else:
            bed_chroms, bed_starts, bed_ends = [], [], []

        try:
            multi_result = stats_core.pileup_variants_multi(
                bam_file_paths, bam_labels,
                chroms, poss, refs, alts,
                bed_chroms, bed_starts, bed_ends,
            )
        except BaseException as e:
            if isinstance(e, KeyboardInterrupt):
                raise
            raise RuntimeError(f"Rust multi-BAM pileup failed: {e}") from e

        if multi_result:
            output = {}
            for label in bam_labels:
                if label in multi_result:
                    result_dict = multi_result[label]
                    result_dict["CHROM"] = chroms
                    result_dict["POS"] = poss
                    result_dict["REF"] = refs
                    result_dict["ALT"] = alts
                    output[label] = pl.DataFrame(result_dict)
            return output

        # Multi-BAM returned empty — fall through to per-BAM Rust
        print(f"  [BAM pileup] Rust multi-BAM returned empty, trying per-BAM ({len(bam_paths)} BAMs, {len(positions)} positions)")

    # Fallback: per-BAM Rust calls (used when HAS_RUST_MULTI is False or multi returned empty)
    print(f"  [BAM pileup] Using per-BAM Rust calls ({len(bam_paths)} BAMs, {len(positions)} positions)")
    result = {}
    if len(bam_paths) > 1:
        with ThreadPoolExecutor(max_workers=min(len(bam_paths), 4)) as executor:
            futures = {
                executor.submit(pileup_variants, path, positions): label
                for label, path in bam_paths.items()
            }
            for future in as_completed(futures):
                label = futures[future]
                try:
                    df = future.result()
                except Exception as e:
                    print(f"  [BAM pileup] {label}: error: {e}")
                    continue
                if df is not None and not df.is_empty():
                    result[label] = df
                    print(f"  [BAM pileup] {label}: {len(df)} rows, DP non-null: {(df['DP'].is_not_null()).sum()}")
                else:
                    print(f"  [BAM pileup] {label}: no data returned")
    else:
        for label, path in bam_paths.items():
            df = pileup_variants(path, positions)
            if df is not None and not df.is_empty():
                result[label] = df
                print(f"  [BAM pileup] {label}: {len(df)} rows, DP non-null: {(df['DP'].is_not_null()).sum()}")
            else:
                print(f"  [BAM pileup] {label}: no data returned")
    return result
