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

            # Strelka-specific: TAR, TIR, TOR
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
) -> dict[tuple[str, int], dict[str, Any]]:
    """Convert caller result lists to a position-keyed lookup dict."""
    lookup = {}
    chroms = caller_result.get("CHROM", [])
    poss = caller_result.get("POS", [])
    keys = [k for k in caller_result.keys() if k not in ("CHROM", "POS")]

    for i in range(len(chroms)):
        key = (chroms[i], poss[i])
        lookup[key] = {k: caller_result[k][i] for k in keys}
    return lookup


def _parse_one_caller(
    caller_name: str,
    cfg: dict,
    base: str,
    vcf_prefix: str,
    target_positions: set[tuple[str, int]],
) -> tuple[str, dict]:
    """Parse a single caller VCF. Returns (caller_name, lookup_dict)."""
    subdir = cfg["subdir"].format(prefix=vcf_prefix)
    vcf_path = _find_vcf_file(base, subdir, cfg["pattern"])

    if vcf_path is None:
        print(f"    [{caller_name}] VCF not found, skipping")
        return (caller_name, {})

    print(f"    [{caller_name}] scanning {vcf_path}...")
    result = parse_single_caller(
        vcf_path, target_positions, cfg["sample_suffix"], caller_name
    )
    lookup = build_caller_results_lookup(result)
    print(f"    [{caller_name}] found {len(lookup)} variants")
    return (caller_name, lookup)


def parse_all_callers(
    base_output_dir: str,
    dir_name: str,
    vcf_prefix: str,
    target_positions: set[tuple[str, int]],
    max_workers: int = 1,
) -> dict[str, dict[str, Any]]:
    """Parse all 6 caller VCFs and return position-keyed results.

    Args:
        base_output_dir: Sample's base output directory.
        dir_name: Sample's directory name.
        vcf_prefix: VCF prefix for this sample.
        target_positions: Set of (CHROM, POS) to extract.
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
    """Join caller FORMAT columns onto the rescue DataFrame by (CHROM, POS).

    For each caller, adds columns: {caller}_DP, {caller}_AD_REF, {caller}_AD_ALT,
    {caller}_GT (where available), {caller}_VAF (computed), etc.

    Strelka: AD_REF = TOR[0], AD_ALT = TAR[0]
    """
    df = rescue_df.clone()
    chroms = df["CHROM"].to_list()
    poss = df["POS"].to_list()
    n = len(chroms)

    for caller_name in CALLER_CONFIGS:
        lookup = caller_data.get(caller_name, {})
        is_strelka = caller_name in CALLERS_STRELKA
        has_gt = caller_name in CALLERS_WITH_GT

        # DP
        dp_vals = []
        for i in range(n):
            key = (chroms[i], poss[i])
            entry = lookup.get(key, {})
            dp_vals.append(entry.get("DP"))
        df = df.with_columns(pl.Series(f"{caller_name}_DP", dp_vals, dtype=pl.Int64))

        # AD (or Strelka proxy)
        if is_strelka:
            ad_ref_vals = []
            ad_alt_vals = []
            for i in range(n):
                key = (chroms[i], poss[i])
                entry = lookup.get(key, {})
                ad_ref_vals.append(entry.get("TOR"))  # TOR[0] = other reads tier1
                ad_alt_vals.append(entry.get("TAR"))  # TAR[0] = alt reads tier1
            df = df.with_columns([
                pl.Series(f"{caller_name}_AD_REF", ad_ref_vals, dtype=pl.Int64).alias(f"{caller_name}_AD_REF"),
                pl.Series(f"{caller_name}_AD_ALT", ad_alt_vals, dtype=pl.Int64).alias(f"{caller_name}_AD_ALT"),
            ])
        else:
            ad_ref_vals = []
            ad_alt_vals = []
            for i in range(n):
                key = (chroms[i], poss[i])
                entry = lookup.get(key, {})
                ad_ref_vals.append(entry.get("AD_REF"))
                ad_alt_vals.append(entry.get("AD_ALT"))
            df = df.with_columns([
                pl.Series(f"{caller_name}_AD_REF", ad_ref_vals, dtype=pl.Int64).alias(f"{caller_name}_AD_REF"),
                pl.Series(f"{caller_name}_AD_ALT", ad_alt_vals, dtype=pl.Int64).alias(f"{caller_name}_AD_ALT"),
            ])

        # GT
        if has_gt:
            gt_vals = []
            for i in range(n):
                key = (chroms[i], poss[i])
                entry = lookup.get(key, {})
                gt_vals.append(entry.get("GT"))
            df = df.with_columns(pl.Series(f"{caller_name}_GT", gt_vals, dtype=pl.Utf8))

        # Pre-computed VAF/AF from caller
        precomp_key = CALLERS_WITH_PRECOMPUTED_VAF.get(caller_name)
        if precomp_key:
            vaf_caller_vals = []
            for i in range(n):
                key = (chroms[i], poss[i])
                entry = lookup.get(key, {})
                vaf_caller_vals.append(entry.get("VAF_CALLER"))
            df = df.with_columns(
                pl.Series(f"{caller_name}_VAF_CALLER", vaf_caller_vals, dtype=pl.Float64)
            )

    return df
