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
    vcf_path: str | None = None,
) -> tuple[str, dict]:
    """Parse a single caller VCF. Returns (caller_name, lookup_dict).

    Uses Rust stats_core.parse_caller_vcf when available (faster, GIL-released).
    Falls back to cyvcf2 when Rust is unavailable.

    Args:
        vcf_path: Pre-resolved VCF path from manifest. When provided, skips glob.
    """
    if vcf_path and os.path.isfile(vcf_path):
        pass  # Use manifest-provided path
    else:
        subdir = cfg["subdir"].format(prefix=vcf_prefix)
        vcf_path = _find_vcf_file(base, subdir, cfg["pattern"])

    if vcf_path is None:
        print(f"    [{caller_name}] VCF not found, skipping")
        return (caller_name, {})

    print(f"    [{caller_name}] scanning {vcf_path}...")

    if HAS_RUST_CALLER and target_positions:
        # Use Rust parser with 4-tuple target positions.
        # Returns column-oriented data directly — no intermediate lookup dict.
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
                n_found = len(result.get("CHROM", []))
                print(f"    [{caller_name}] found {n_found} variants (Rust, column-oriented)")
                return (caller_name, result)  # column-oriented: {col: [vals]}
        except Exception as e:
            print(f"    [{caller_name}] Rust parser failed ({e}), falling back to cyvcf2")

    # Python fallback: use 2-tuple positions for cyvcf2
    pos2 = {(t[0], t[1]) for t in target_positions}
    result = parse_single_caller(
        vcf_path, pos2, cfg["sample_suffix"], caller_name
    )
    # Convert to column-oriented via build_caller_results_lookup (needed for cyvcf2 compat)
    lookup = build_caller_results_lookup(result)
    # Convert lookup back to column-oriented for consistent return type
    cols = _lookup_to_columns(lookup)
    print(f"    [{caller_name}] found {len(cols.get('CHROM', []))} variants")
    return (caller_name, cols)


def _lookup_to_columns(lookup: dict) -> dict[str, list]:
    """Convert a row-oriented lookup dict to column-oriented data.

    Used only for the cyvcf2 fallback path — Rust parser returns columns directly.
    """
    if not lookup:
        return {}
    cols: dict[str, list] = {"CHROM": [], "POS": [], "REF": [], "ALT": []}
    field_names: set[str] = set()
    for key, fields in lookup.items():
        cols["CHROM"].append(key[0])
        cols["POS"].append(key[1])
        cols["REF"].append(key[2] if len(key) >= 4 else "")
        cols["ALT"].append(key[3] if len(key) >= 4 else "")
        for fname in fields:
            field_names.add(fname)
    # Add field columns
    for fname in field_names:
        col = []
        for key, fields in lookup.items():
            col.append(fields.get(fname))
        cols[fname] = col
    return cols


def parse_all_callers(
    base_output_dir: str,
    dir_name: str,
    vcf_prefix: str,
    target_positions: set[tuple],
    max_workers: int = 1,
) -> dict[str, dict[str, list]]:
    """Parse all 6 caller VCFs and return column-oriented results.

    Returns:
        Dict: {caller_name: {column_name: [values]}} — column-oriented.
    """
    base = os.path.join(base_output_dir, dir_name)
    caller_data: dict[str, dict[str, list]] = {}

    # With Rust parser, caller VCFs are fast (~0.3-1.2s each, 4s total).
    # ThreadPoolExecutor adds overhead and creates 6+ threads — avoid when
    # already inside an outer ThreadPoolExecutor (sample_workers > 1).
    use_parallel = max_workers > 1 and not HAS_RUST_CALLER

    if not use_parallel:
        for caller_name, cfg in CALLER_CONFIGS.items():
            name, cols = _parse_one_caller(caller_name, cfg, base, vcf_prefix, target_positions)
            caller_data[name] = cols
    else:
        with ThreadPoolExecutor(max_workers=min(max_workers, len(CALLER_CONFIGS))) as executor:
            futures = {
                executor.submit(
                    _parse_one_caller, caller_name, cfg, base, vcf_prefix, target_positions
                ): caller_name
                for caller_name, cfg in CALLER_CONFIGS.items()
            }
            for future in as_completed(futures):
                name, cols = future.result()
                caller_data[name] = cols

    return caller_data


def join_caller_columns(
    rescue_df: pl.DataFrame,
    caller_data: dict[str, dict[str, list]],
) -> pl.DataFrame:
    """Join caller FORMAT columns onto the rescue DataFrame using polars join.

    Receives column-oriented data directly (no row-oriented lookup dict),
    eliminating ~5 GB of Python small-object overhead per 1.4M-variant sample.

    Uses (CHROM, POS, REF, ALT) 4-column matching. Missing callers get null-filled columns.
    """
    df = rescue_df.clone()
    join_cols = ["CHROM", "POS", "REF", "ALT"]

    for col in join_cols:
        if col not in df.columns:
            raise KeyError(f"Rescue DataFrame missing join column: {col}")

    for caller_name in CALLER_CONFIGS:
        col_data = caller_data.get(caller_name, {})
        is_strelka = caller_name in CALLERS_STRELKA
        has_gt = caller_name in CALLERS_WITH_GT

        if not col_data or not col_data.get("CHROM"):
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

        # Build caller DataFrame from column-oriented data with explicit dtypes
        int_fields = {"DP", "AD_REF", "AD_ALT", "TAR", "TIR", "TOR", "AU", "CU", "GU", "TU", "POS"}
        float_fields = {"VAF_CALLER"}
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

        if not series_list:
            continue

        caller_df = pl.DataFrame(series_list)

        if not data_col_names:
            continue

        # Rename data columns with caller prefix for join
        rename_map = {c: f"{caller_name}_{c}" for c in data_col_names}
        caller_df = caller_df.select(join_cols + data_col_names).rename(rename_map)

        # Left join on all 4 coordinate columns
        df = df.join(caller_df, on=join_cols, how="left")

        # Free intermediate data immediately
        del caller_df, series_list, col_data

        # Ensure Strelka has AD_REF/AD_ALT columns
        if is_strelka:
            for suf, src in [("AD_REF", "TAR"), ("AD_ALT", "TIR")]:
                src_col = f"{caller_name}_{src}"
                tgt_col = f"{caller_name}_{suf}"
                if src_col in df.columns:
                    df = df.with_columns(pl.col(src_col).alias(tgt_col))
                elif tgt_col not in df.columns:
                    df = df.with_columns(pl.lit(None).alias(tgt_col))

        if not is_strelka and has_gt:
            for suf in ["AD_REF", "AD_ALT"]:
                col = f"{caller_name}_{suf}"
                if col not in df.columns:
                    df = df.with_columns(pl.lit(None).alias(col))

        # Free caller_data entry to help GC
        caller_data[caller_name] = {}

    return df


def parse_pre_norm_multiallelic(
    vcf_path: str,
    caller_name: str,
    sample_suffix: str,
) -> pl.DataFrame | None:
    """Parse multi-allelic records from a pre-normalization caller VCF.

    Extracts per-allele AD, AF/VAF, F1R2/F2R1 (Mutect2 only), and GT for
    records where ALT contains a comma (multi-allelic sites).

    Args:
        vcf_path: Path to the pre-normalization VCF file.
        caller_name: Caller key from CALLER_CONFIGS.
        sample_suffix: Sample suffix to match (DT, RT).

    Returns:
        DataFrame with columns: CHROM, POS, REF, ALT_original, n_original_alleles,
        original_alts, ad_ref, ad_alts, af_list, f1r2_list, f2r1_list, gt_str, gt_alleles.
        Returns None if the file is missing or no multi-allelic records are found.
    """
    import os as _os
    from cyvcf2 import VCF as _VCF

    if not _os.path.isfile(vcf_path):
        return None

    is_mutect2 = "mutect2" in caller_name.lower()
    is_deepsomatic = "deepsomatic" in caller_name.lower()

    af_field = "AF" if is_mutect2 else ("VAF" if is_deepsomatic else None)
    if af_field is None:
        # Not a supported caller for multi-allelic parsing
        return None

    rows: list[dict] = []
    try:
        reader = _VCF(vcf_path)
        sample_idx = None
        for i, name in enumerate(reader.samples):
            if name.endswith(sample_suffix) or name == sample_suffix:
                sample_idx = i
                break
        if sample_idx is None:
            reader.close()
            return None

        for record in reader:
            alts = record.ALT
            if len(alts) <= 1:
                continue  # Skip non-multi-allelic

            alt_str = ",".join(str(a) for a in alts)

            # AD array (full, all alleles)
            ad_ref = None
            ad_alts = None
            if "AD" in record.FORMAT:
                try:
                    ad_val = record.format("AD")[sample_idx]
                    if ad_val is not None and len(ad_val) >= 2:
                        ad_ref = int(ad_val[0])
                        ad_alts = [int(x) for x in ad_val[1:]]
                except Exception:
                    pass

            # Per-allele AF/VAF
            af_list = None
            if af_field in record.FORMAT:
                try:
                    af_val = record.format(af_field)[sample_idx]
                    if af_val is not None:
                        af_list = [float(x) for x in af_val]
                except Exception:
                    pass

            # F1R2 and F2R1 (Mutect2 only)
            f1r2_list = None
            f2r1_list = None
            if is_mutect2:
                for field, target in [("F1R2", "f1r2"), ("F2R1", "f2r1")]:
                    if field in record.FORMAT:
                        try:
                            f_val = record.format(field)[sample_idx]
                            if f_val is not None and len(f_val) >= 2:
                                if target == "f1r2":
                                    f1r2_list = [int(x) for x in f_val[1:]]
                                else:
                                    f2r1_list = [int(x) for x in f_val[1:]]
                        except Exception:
                            pass

            # GT as allele indices
            gt_str = None
            gt_alleles = None
            if "GT" in record.FORMAT:
                try:
                    gt_val = record.format("GT")[sample_idx]
                    if gt_val is not None and len(gt_val) > 0:
                        gt_str = "/".join(str(v) for v in gt_val)
                        gt_alleles = [int(v) if v is not None else -1 for v in gt_val]
                except Exception:
                    pass

            rows.append({
                "CHROM": record.CHROM,
                "POS": record.POS,
                "REF": record.REF,
                "ALT_original": alt_str,
                "n_original_alleles": len(alts),
                "original_alts": [str(a) for a in alts],
                "ad_ref": ad_ref,
                "ad_alts": ad_alts,
                "af_list": af_list,
                "f1r2_list": f1r2_list,
                "f2r1_list": f2r1_list,
                "gt_str": gt_str,
                "gt_alleles": gt_alleles,
            })

        reader.close()
    except Exception:
        return None

    if not rows:
        return None

    return pl.DataFrame(rows)
