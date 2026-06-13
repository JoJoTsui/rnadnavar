"""Python wrapper for Rust BAM pileup with pysam fallback.

Uses stats_core (Rust) when available, falls back to pysam when not.
"""

import os
import polars as pl

try:
    import stats_core
    HAS_RUST_BAM = hasattr(stats_core, 'pileup_variants')
except ImportError:
    HAS_RUST_BAM = False

try:
    import pysam
    HAS_PYSAM = True
except ImportError:
    HAS_PYSAM = False


def pileup_variants(bam_path: str, positions: list[tuple[str, int, str, str]]) -> pl.DataFrame | None:
    """Perform BAM pileup at specific variant positions.

    Args:
        bam_path: Path to the BAM file.
        positions: List of (chrom, pos, ref_base, alt_base) tuples.

    Returns:
        DataFrame with columns: CHROM, POS, DP, REF_DP, ALT_DP,
        F1R2_ref, F2R1_ref, F1R2_alt, F2R1_alt, mean_BQ, mean_MQ.
        Returns None if BAM is unavailable or parsing fails.
    """
    if not os.path.isfile(bam_path):
        return None

    if HAS_RUST_BAM:
        try:
            chroms = [p[0] for p in positions]
            poss = [p[1] for p in positions]
            refs = [p[2] for p in positions]
            alts = [p[3] for p in positions]
            result = stats_core.pileup_variants(bam_path, chroms, poss, refs, alts)
            if result:
                # Add coordinate columns (Rust returns only metric columns)
                result["CHROM"] = chroms
                result["POS"] = poss
                result["REF"] = refs
                result["ALT"] = alts
                return pl.DataFrame(result)
        except Exception as e:
            print(f"  [WARNING] Rust BAM pileup failed: {e}, falling back to pysam")

    # Python fallback using pysam
    if HAS_PYSAM:
        return _pileup_pysam(bam_path, positions)

    return None


def _pileup_pysam(bam_path: str, positions: list[tuple[str, int, str, str]]) -> pl.DataFrame:
    """Pysam-based pileup at variant positions (fallback)."""
    import numpy as np
    bam = pysam.AlignmentFile(bam_path, "rb")

    rows = []
    for chrom, pos, ref_base, alt_base in positions:
        row = {
            "CHROM": chrom, "POS": pos, "REF": ref_base, "ALT": alt_base,
            "DP": None, "REF_DP": None, "ALT_DP": None,
            "F1R2_ref": None, "F2R1_ref": None,
            "F1R2_alt": None, "F2R1_alt": None,
            "mean_BQ": None, "mean_MQ": None,
        }

        try:
            reads = list(bam.fetch(chrom, pos - 1, pos))
            if not reads:
                rows.append(row)
                continue

            dp = 0
            ref_dp = 0
            alt_dp = 0
            f1r2_ref = 0
            f2r1_ref = 0
            f1r2_alt = 0
            f2r1_alt = 0
            bq_sum = 0.0
            bq_count = 0
            mq_sum = 0.0

            for read in reads:
                if read.is_unmapped or read.is_duplicate:
                    continue
                try:
                    pos_in_read = pos - read.reference_start - 1
                    if pos_in_read < 0 or pos_in_read >= len(read.query_sequence):
                        continue
                    base = read.query_sequence[pos_in_read]
                    q = read.query_qualities[pos_in_read] if read.query_qualities else None
                    dp += 1
                    if base == ref_base:
                        ref_dp += 1
                        if read.is_reverse:
                            f2r1_ref += 1
                        else:
                            f1r2_ref += 1
                    elif base == alt_base:
                        alt_dp += 1
                        if read.is_reverse:
                            f2r1_alt += 1
                        else:
                            f1r2_alt += 1
                    if q is not None:
                        bq_sum += q
                        bq_count += 1
                    if read.mapping_quality is not None:
                        mq_sum += read.mapping_quality
                except Exception:
                    continue

            row["DP"] = dp
            row["REF_DP"] = ref_dp
            row["ALT_DP"] = alt_dp
            row["F1R2_ref"] = f1r2_ref
            row["F2R1_ref"] = f2r1_ref
            row["F1R2_alt"] = f1r2_alt
            row["F2R1_alt"] = f2r1_alt
            row["mean_BQ"] = round(bq_sum / bq_count, 1) if bq_count > 0 else None
            row["mean_MQ"] = round(mq_sum / dp, 1) if dp > 0 else None
        except Exception:
            pass

        rows.append(row)

    bam.close()
    if not rows:
        return pl.DataFrame()
    return pl.DataFrame(rows)
