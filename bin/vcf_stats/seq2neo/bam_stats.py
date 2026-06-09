"""Per-sample per-modality BAM statistics.

Computes read-level statistics from alignment BAM files for both DNA and RNA
modalities. Uses pysam for BAM parsing. Falls back gracefully when BAM files
are missing or pysam is unavailable.
"""

import os
from pathlib import Path
from typing import Any

import polars as pl

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


def compute_bam_stats(bam_path: str) -> dict[str, Any] | None:
    """Compute basic statistics from a BAM file.

    Returns dict with: total_reads, mapped_reads, unmapped_reads, mapping_rate,
    mean_coverage, mean_insert_size, mean_mapq. Returns None if BAM is unreadable.
    """
    if not HAS_PYSAM:
        return None

    if not bam_path or not os.path.isfile(bam_path):
        return None

    try:
        bam = pysam.AlignmentFile(bam_path, "rb")

        total_reads = 0
        mapped_reads = 0
        total_mapq = 0
        total_insert = 0
        insert_count = 0
        total_length = 0

        # Reference lengths for coverage estimation
        ref_lengths = sum(bam.lengths) if bam.lengths else 1

        for read in bam.fetch():
            total_reads += 1
            if not read.is_unmapped:
                mapped_reads += 1
                total_mapq += read.mapping_quality
                total_length += read.query_length or 0
                if read.template_length and read.template_length > 0:
                    total_insert += read.template_length
                    insert_count += 1

            # Sample: stop after 1M reads for performance
            if total_reads >= 1_000_000:
                break

        bam.close()

        if total_reads == 0:
            return None

        mapping_rate = mapped_reads / total_reads * 100 if total_reads > 0 else 0
        mean_mapq = total_mapq / mapped_reads if mapped_reads > 0 else 0
        mean_insert = total_insert / insert_count if insert_count > 0 else 0
        mean_coverage = total_length / ref_lengths if total_length > 0 else 0

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


def compute_sample_bam_stats(
    base_output_dir: str,
    dir_name: str,
    sample_id: str,
    set_number: int,
) -> list[dict[str, Any]]:
    """Compute BAM statistics for a single sample (DNA + RNA modalities).

    Returns a list of dicts, one per modality, for aggregation into a DataFrame.
    Each dict contains sample metadata and BAM stats. Modalities with missing
    BAM files produce rows with null stats.
    """
    results = []

    for bam_type in ["DN", "DT", "RT"]:
        bam_path = _locate_bam_file(base_output_dir, dir_name, bam_type)
        stats = compute_bam_stats(bam_path) if bam_path else None

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


def compute_all_bam_stats(manifest_rows: list[dict]) -> pl.DataFrame:
    """Compute BAM statistics for all samples in the manifest.

    Args:
        manifest_rows: List of manifest row dicts with sample_id, base_output_dir,
            dir_name, set_number.

    Returns:
        polars DataFrame with one row per sample per modality.
    """
    all_rows = []
    for row in manifest_rows:
        sample_results = compute_sample_bam_stats(
            base_output_dir=row["base_output_dir"],
            dir_name=row["dir_name"],
            sample_id=row["sample_id"],
            set_number=row["set_number"],
        )
        all_rows.extend(sample_results)
        ok_flags = []
        for r in sample_results:
            ok_flags.append(f"{r['bam_type']}={'OK' if r['has_bam'] else 'missing'}")
        print(f"  [{row['sample_id']}] BAM stats: {', '.join(ok_flags)}")

    if not all_rows:
        return pl.DataFrame()

    return pl.DataFrame(all_rows)
