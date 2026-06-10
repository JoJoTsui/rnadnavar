"""Parse individual caller VCFs to extract ground-truth FORMAT fields.

Each caller has different FORMAT fields — see CALLER_CONFIGS in manifest_loader.
This module handles caller-specific extraction:

- Mutect2/DeepSomatic: GT, AD (REF+ALT), DP, plus pre-computed AF/VAF
- Strelka: DP (tier1), TAR/TIR/TOR (no GT, no AD)

Design:
  1. Build a set of (chrom, pos) target positions from the rescue VCF.
  2. Scan each caller VCF sequentially, extracting FORMAT fields only for
     positions in the target set.
  3. Early termination: stop scanning when all target positions are found.
  4. Missing callers → null-filled columns.
"""

import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Any

import numpy as np
import polars as pl
from cyvcf2 import VCF

from .manifest_loader import CALLER_CONFIGS, _find_vcf_file

try:
    import stats_core
    HAS_RUST_CALLER = hasattr(stats_core, 'parse_caller_vcf')
except ImportError:
    HAS_RUST_CALLER = False

# ── Callers that have GT/AD ───────────────────────────────────────────────
CALLERS_WITH_GT = {"DNA_mutect2", "RNA_mutect2", "DNA_deepsomatic", "RNA_deepsomatic"}
CALLERS_WITH_AD = {"DNA_mutect2", "RNA_mutect2", "DNA_deepsomatic", "RNA_deepsomatic"}
CALLERS_WITH_PRECOMPUTED_VAF = {
    "DNA_mutect2": "AF",
    "RNA_mutect2": "AF",
    "DNA_deepsomatic": "VAF",
    "RNA_deepsomatic": "VAF",
}
# Strelka uses TAR/TOR instead of AD
CALLERS_STRELKA = {"DNA_strelka", "RNA_strelka"}


def _find_sample_index(vcf: VCF, sample_suffix: str) -> int:
    """Find the sample index in a VCF by suffix matching."""
    for i, name in enumerate(vcf.samples):
        if name.endswith(sample_suffix) or name == sample_suffix:
            return i
    return -1


def parse_single_caller(
    vcf_path: str,
    target_positions: set[tuple[str, int]],
    sample_suffix: str,
    caller_name: str,
) -> dict[str, list]:
    """Parse a single caller VCF for FORMAT fields at target positions.

    Args:
        vcf_path: Path to the caller VCF file.
        target_positions: Set of (CHROM, POS) tuples to extract.
        sample_suffix: Sample suffix to match (DT, RT, TUMOR).
        caller_name: Caller key from CALLER_CONFIGS.

    Returns:
        Dict mapping column name -> list of values (same length as target_positions
        but in the order of scanning, not in target order).
    """
    result: dict[str, list] = {
        "CHROM": [],
        "POS": [],
    }

    is_strelka = caller_name in CALLERS_STRELKA
    has_gt = caller_name in CALLERS_WITH_GT
    has_ad = caller_name in CALLERS_WITH_AD
    has_mutect2_fields = has_gt and "mutect2" in caller_name.lower()
    precomp_vaf_key = CALLERS_WITH_PRECOMPUTED_VAF.get(caller_name)

    # Initialize result keys based on caller type
    result["DP"] = []
    if has_ad:
        result["AD_REF"] = []
        result["AD_ALT"] = []
    if has_gt:
        result["GT"] = []
    if precomp_vaf_key:
        result["VAF_CALLER"] = []
    if is_strelka:
        result["TAR"] = []
        result["TIR"] = []
        result["TOR"] = []
        result["AU"] = []  # A-allele count tier1
        result["CU"] = []  # C-allele count tier1
        result["GU"] = []  # G-allele count tier1
        result["TU"] = []  # T-allele count tier1
    if has_mutect2_fields:
        result["SB"] = []  # Strand bias (4 values)
        result["FAD"] = []  # Fragment allele depth

    if not os.path.isfile(vcf_path):
        return result

    remaining = set(target_positions)
    if not remaining:
        return result

    try:
        reader = VCF(vcf_path)
        sample_idx = _find_sample_index(reader, sample_suffix)
        if sample_idx == -1:
            reader.close()
            return result

        for record in reader:
            key = (record.CHROM, record.POS)
            if key not in remaining:
                continue

            result["CHROM"].append(record.CHROM)
            result["POS"].append(record.POS)

            # DP
            dp = None
            if "DP" in record.FORMAT:
                try:
                    dp_val = record.format("DP")[sample_idx]
                    dp = int(dp_val[0]) if dp_val is not None and len(dp_val) > 0 else None
                except Exception:
                    dp = None
            result["DP"].append(dp)

            # AD (Mutect2 + DeepSomatic)
            if has_ad:
                ad_ref, ad_alt = None, None
                if "AD" in record.FORMAT:
                    try:
                        ad_val = record.format("AD")[sample_idx]
                        if ad_val is not None and len(ad_val) >= 2:
                            ad_ref = int(ad_val[0])
                            ad_alt = int(ad_val[1])
                    except Exception:
                        pass
                result["AD_REF"].append(ad_ref)
                result["AD_ALT"].append(ad_alt)

            # GT
            if has_gt:
                gt = None
                if "GT" in record.FORMAT:
                    try:
                        gt_val = record.format("GT")[sample_idx]
                        if gt_val is not None and len(gt_val) > 0:
                            gt = str(gt_val[0])
                    except Exception:
                        pass
                result["GT"].append(gt)

            # Pre-computed VAF/AF
            if precomp_vaf_key and precomp_vaf_key in record.FORMAT:
                try:
                    vaf_val = record.format(precomp_vaf_key)[sample_idx]
                    vaf = float(vaf_val[0]) if vaf_val is not None and len(vaf_val) > 0 else None
                except Exception:
                    vaf = None
                result["VAF_CALLER"].append(vaf)

            # Strelka-specific: TAR, TIR, TOR, AU/CU/GU/TU
            if is_strelka:
                for tar_field in ["TAR", "TIR", "TOR"]:
                    value = None
                    if tar_field in record.FORMAT:
                        try:
                            f_val = record.format(tar_field)[sample_idx]
                            value = int(f_val[0]) if f_val is not None and len(f_val) > 0 else None
                        except Exception:
                            pass
                    result[tar_field].append(value)
                # Per-allele counts (tier1 only)
                for allele_field in ["AU", "CU", "GU", "TU"]:
                    value = None
                    if allele_field in record.FORMAT:
                        try:
                            f_val = record.format(allele_field)[sample_idx]
                            value = int(f_val[0]) if f_val is not None and len(f_val) > 0 else None
                        except Exception:
                            pass
                    result[allele_field].append(value)

            # Mutect2-specific: SB (strand bias), FAD (fragment allele depth)
            if has_mutect2_fields:
                for field_name in ["SB", "FAD"]:
                    value = None
                    if field_name in record.FORMAT:
                        try:
                            f_val = record.format(field_name)[sample_idx]
                            if f_val is not None and len(f_val) > 0:
                                if field_name == "SB":
                                    # SB has 4 values: F1R2, F2R1 for ref and alt
                                    value = ",".join(str(v) for v in f_val) if len(f_val) > 1 else str(f_val[0])
                                else:
                                    value = str(f_val[0])
                        except Exception:
                            pass
                    result[field_name].append(value)

            remaining.discard(key)
            if not remaining:
                break

        reader.close()
    except Exception as e:
        print(f"  [ERROR] parse_single_caller({caller_name}, {vcf_path}): {e}")
        import traceback
        traceback.print_exc()

    return result


def build_caller_results_lookup(
    caller_result: dict[str, list],
) -> dict[tuple[str, int, str, str], dict[str, Any]]:
    """Convert caller result lists to a position-keyed lookup dict.

    Uses (CHROM, POS, REF, ALT) 4-tuple keys for correct matching at
    multiallelic sites. Falls back to (CHROM, POS) if REF/ALT not present.
    """
    lookup = {}
    chroms = caller_result.get("CHROM", [])
    poss = caller_result.get("POS", [])
    refs = caller_result.get("REF", [])
    alts = caller_result.get("ALT", [])
    keys = [k for k in caller_result.keys() if k not in ("CHROM", "POS", "REF", "ALT")]

    for i in range(len(chroms)):
        if refs and alts:
            key = (chroms[i], poss[i], refs[i], alts[i])
        else:
            # Fallback for cyvcf2 results (2-tuple)
            key = (chroms[i], poss[i], "", "")
        lookup[key] = {k: caller_result[k][i] for k in keys}
    return lookup


def _parse_one_caller(
    caller_name: str,
    cfg: dict,
    base: str,
    vcf_prefix: str,
    target_positions: set[tuple],
) -> tuple[str, dict]:
    """Parse a single caller VCF. Returns (caller_name, lookup_dict).

    Uses Rust stats_core.parse_caller_vcf when available (faster, GIL-released).
    Falls back to cyvcf2 when Rust is unavailable.
    """
    subdir = cfg["subdir"].format(prefix=vcf_prefix)
    vcf_path = _find_vcf_file(base, subdir, cfg["pattern"])

    if vcf_path is None:
        print(f"    [{caller_name}] VCF not found, skipping")
        return (caller_name, {})

    print(f"    [{caller_name}] scanning {vcf_path}...")

    if HAS_RUST_CALLER and target_positions:
        # Use Rust parser with 4-tuple target positions
        try:
            sample = next(iter(target_positions))
            if len(sample) >= 4:
                chroms = [t[0] for t in target_positions]
                poss = [t[1] for t in target_positions]
                refs = [t[2] for t in target_positions]
                alts = [t[3] for t in target_positions]
                result = stats_core.parse_caller_vcf(
                    vcf_path, chroms, poss, refs, alts,
                    cfg["sample_suffix"], caller_name,
                )
                lookup = build_caller_results_lookup(result)
                print(f"    [{caller_name}] found {len(lookup)} variants (Rust)")
                return (caller_name, lookup)
        except Exception as e:
            print(f"    [{caller_name}] Rust parser failed ({e}), falling back to cyvcf2")

    # Python fallback: use 2-tuple positions for cyvcf2
    pos2 = {(t[0], t[1]) for t in target_positions}
    result = parse_single_caller(
        vcf_path, pos2, cfg["sample_suffix"], caller_name
    )
    lookup = build_caller_results_lookup(result)
    print(f"    [{caller_name}] found {len(lookup)} variants")
    return (caller_name, lookup)


def parse_all_callers(
    base_output_dir: str,
    dir_name: str,
    vcf_prefix: str,
    target_positions: set[tuple],
    max_workers: int = 1,
) -> dict[str, dict[str, Any]]:
    """Parse all 6 caller VCFs and return position-keyed results.

    Args:
        base_output_dir: Sample's base output directory.
        dir_name: Sample's directory name.
        vcf_prefix: VCF prefix for this sample.
        target_positions: Set of (CHROM, POS, REF, ALT) 4-tuples to extract.
        max_workers: Number of threads for parallel caller parsing.

    Returns:
        Nested dict: {caller_name: {(chrom, pos): {field: value}}}
    """
    base = os.path.join(base_output_dir, dir_name)
    caller_data = {}

    if max_workers <= 1:
        # Sequential
        for caller_name, cfg in CALLER_CONFIGS.items():
            name, lookup = _parse_one_caller(caller_name, cfg, base, vcf_prefix, target_positions)
            caller_data[name] = lookup
    else:
        # Thread-parallel (threads are safe with cyvcf2/htslib; processes are not)
        with ThreadPoolExecutor(max_workers=min(max_workers, len(CALLER_CONFIGS))) as executor:
            futures = {
                executor.submit(
                    _parse_one_caller, caller_name, cfg, base, vcf_prefix, target_positions
                ): caller_name
                for caller_name, cfg in CALLER_CONFIGS.items()
            }
            for future in as_completed(futures):
                name, lookup = future.result()
                caller_data[name] = lookup

    return caller_data


def join_caller_columns(
    rescue_df: pl.DataFrame,
    caller_data: dict[str, dict[str, Any]],
) -> pl.DataFrame:
    """Join caller FORMAT columns onto the rescue DataFrame using polars join.

    Uses (CHROM, POS, REF, ALT) 4-column matching for correct behavior at
    multiallelic sites. Each caller's lookup dict is converted to a temporary
    DataFrame and left-joined. Missing callers get null-filled columns.

    For each caller, adds columns: {caller}_DP, {caller}_AD_REF, {caller}_AD_ALT,
    {caller}_GT (where available), {caller}_VAF_CALLER (where available), etc.
    """
    df = rescue_df.clone()
    join_cols = ["CHROM", "POS", "REF", "ALT"]

    # Ensure the rescue DF has the join columns
    for col in join_cols:
        if col not in df.columns:
            raise KeyError(f"Rescue DataFrame missing join column: {col}")

    for caller_name in CALLER_CONFIGS:
        lookup = caller_data.get(caller_name, {})
        is_strelka = caller_name in CALLERS_STRELKA
        has_gt = caller_name in CALLERS_WITH_GT

        if not lookup:
            # Missing caller — null-fill all expected columns
            null_fields = ["DP", "AD_REF", "AD_ALT"]
            if has_gt:
                null_fields.extend(["GT", "VAF_CALLER"])
            if is_strelka:
                null_fields.extend(["TAR", "TIR", "TOR", "AU", "CU", "GU", "TU"])
            if has_gt and "mutect2" in caller_name.lower():
                null_fields.extend(["SB", "FAD"])
            for field in null_fields:
                df = df.with_columns(pl.lit(None).alias(f"{caller_name}_{field}"))
            continue

        # Build caller DataFrame from lookup dict
        rows = []
        for key, fields in lookup.items():
            row = {"CHROM": key[0], "POS": key[1]}
            if len(key) >= 4:
                row["REF"] = key[2]
                row["ALT"] = key[3]
            row.update(fields)
            rows.append(row)

        caller_df = pl.DataFrame(rows) if rows else pl.DataFrame(schema={c: pl.Utf8 for c in join_cols})

        # Select columns to join (all except join columns)
        data_cols = [c for c in caller_df.columns if c not in join_cols]

        if not data_cols:
            continue

        # Rename data columns with caller prefix for join
        rename_map = {c: f"{caller_name}_{c}" for c in data_cols}
        caller_df = caller_df.select(list(join_cols) + data_cols).rename(rename_map)

        # Left join on all 4 coordinate columns
        df = df.join(caller_df, on=join_cols, how="left")

        # Ensure Strelka has AD_REF/AD_ALT columns (may come as TAR/TOR)
        if is_strelka:
            tor_col = f"{caller_name}_TOR"
            tar_col = f"{caller_name}_TAR"
            ad_ref_col = f"{caller_name}_AD_REF"
            ad_alt_col = f"{caller_name}_AD_ALT"
            if tor_col in df.columns:
                df = df.with_columns(pl.col(tor_col).alias(ad_ref_col))
            else:
                df = df.with_columns(pl.lit(None).alias(ad_ref_col))
            if tar_col in df.columns:
                df = df.with_columns(pl.col(tar_col).alias(ad_alt_col))
            else:
                df = df.with_columns(pl.lit(None).alias(ad_alt_col))

        # Ensure Mutect2/DeepSomatic have AD_REF/AD_ALT
        if not is_strelka and has_gt:
            for suf in ["AD_REF", "AD_ALT"]:
                col = f"{caller_name}_{suf}"
                if col not in df.columns:
                    df = df.with_columns(pl.lit(None).alias(col))

    return df
