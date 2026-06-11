"""Python wrapper for Rust VCF parser with cyvcf2 fallback.

Uses stats_core (Rust) when available, falls back to rescue_parser.py (cyvcf2)
when the Rust module is not built.

The Rust parser now supports two modes:
- parse_rescue: legacy row-oriented (PyDict per record) — kept for backward compat
- parse_rescue_columns: column-oriented ({col: [vals]}) — eliminates 7M PyDicts
"""

import polars as pl

from .rescue_parser import (
    ALL_RESCUE_FIELDS,
    RESCUE_FLAG_FIELDS,
    RESCUE_FLOAT_FIELDS,
    RESCUE_INT_FIELDS,
    _add_derived_columns_polars,
)

try:
    import stats_core
    HAS_RUST = True
    HAS_RUST_COLUMNS = hasattr(stats_core, 'parse_rescue_columns')
except ImportError:
    HAS_RUST = False
    HAS_RUST_COLUMNS = False


def _fill_missing_fields(df: pl.DataFrame) -> pl.DataFrame:
    """Add null-filled columns for any ALL_RESCUE_FIELDS not present in the DataFrame.

    Uses correct dtypes: Int64, Float64, Boolean, or Utf8 based on field type.
    """
    missing = [f for f in ALL_RESCUE_FIELDS if f not in df.columns]
    if missing:
        exprs = []
        for f in missing:
            if f in RESCUE_INT_FIELDS:
                exprs.append(pl.lit(None, dtype=pl.Int64).alias(f))
            elif f in RESCUE_FLOAT_FIELDS:
                exprs.append(pl.lit(None, dtype=pl.Float64).alias(f))
            elif f in RESCUE_FLAG_FIELDS:
                exprs.append(pl.lit(None, dtype=pl.Boolean).alias(f))
            else:
                exprs.append(pl.lit(None, dtype=pl.Utf8).alias(f))
        df = df.with_columns(exprs)
    return df


def _cast_columns_fallback(df: pl.DataFrame) -> pl.DataFrame:
    """Cast string columns to correct types (used only for the row-oriented path).

    The legacy row-oriented Rust parser and cyvcf2 fallback both return all
    values as strings. This function casts them to the types expected by
    downstream statistics code.
    """
    # Fill missing fields
    df = _fill_missing_fields(df)
    for field in RESCUE_INT_FIELDS:
        if field in df.columns:
            df = df.with_columns(
                pl.col(field).cast(pl.Utf8, strict=False).cast(pl.Int64, strict=False)
            )
    for field in RESCUE_FLOAT_FIELDS:
        if field in df.columns:
            df = df.with_columns(
                pl.col(field).cast(pl.Utf8, strict=False).cast(pl.Float64, strict=False)
            )
    for field in RESCUE_FLAG_FIELDS:
        if field in df.columns:
            df = df.with_columns(
                pl.when(pl.col(field).cast(pl.Utf8) == "true")
                .then(True)
                .otherwise(False)
                .alias(field)
            )
    return df


def parse_rescue_vcf_columns(vcf_path: str) -> pl.DataFrame:
    """Parse a rescue VCF using the Rust column-oriented parser.

    Returns column-oriented data directly from Rust, builds polars DataFrame
    from pl.Series (no intermediate PyDicts), and casts columns to correct types.

    This path eliminates the 7M-PyDict bottleneck and the _cast_columns cascade.
    variant_type and ti_tv are computed in Rust during the parse.
    """
    if not HAS_RUST_COLUMNS:
        # Fall back to row-oriented Rust or cyvcf2
        return _parse_rescue_vcf_fallback(vcf_path)

    try:
        cols = stats_core.parse_rescue_columns(vcf_path)
    except Exception as e:
        print(f"  [WARNING] Rust column-oriented parser failed: {e}, falling back")
        return _parse_rescue_vcf_fallback(vcf_path)

    if not cols:
        return pl.DataFrame()

    # Build polars DataFrame from column-oriented data with explicit dtypes
    series_list = []
    join_keys = {"CHROM", "POS", "REF", "ALT"}

    for col_name, values in cols.items():
        if col_name in join_keys:
            if col_name == "POS":
                series_list.append(pl.Series(col_name, values, dtype=pl.Int64))
            else:
                series_list.append(pl.Series(col_name, values, dtype=pl.Utf8))
        elif col_name in RESCUE_INT_FIELDS:
            # Rust returns strings for INFO; cast to Int64
            s = pl.Series(col_name, values, dtype=pl.Utf8)
            series_list.append(s.cast(pl.Int64, strict=False))
        elif col_name in RESCUE_FLOAT_FIELDS:
            s = pl.Series(col_name, values, dtype=pl.Utf8)
            series_list.append(s.cast(pl.Float64, strict=False))
        elif col_name in RESCUE_FLAG_FIELDS:
            series_list.append(
                pl.Series(col_name, [v == "true" if v is not None else False for v in values], dtype=pl.Boolean)
            )
        elif col_name in ("variant_type", "ti_tv"):
            # Derived columns from Rust — variant_type is str, ti_tv is bool
            if col_name == "variant_type":
                series_list.append(pl.Series(col_name, values, dtype=pl.Utf8))
            else:
                # ti_tv: Rust returns Option<bool>, Python gets None/True/False
                series_list.append(pl.Series(col_name, values, dtype=pl.Boolean))
        else:
            # Other INFO string fields
            series_list.append(pl.Series(col_name, values, dtype=pl.Utf8))

    df = pl.DataFrame(series_list)

    # Ensure all expected fields are present (some may be missing from header)
    df = _fill_missing_fields(df)

    return df


def _parse_rescue_vcf_fallback(vcf_path: str) -> pl.DataFrame:
    """Fallback: row-oriented Rust parser or cyvcf2 parser.

    Used when the column-oriented Rust parser is unavailable or fails.
    Uses polars-native derived column expressions instead of .to_list().
    """
    # Try legacy row-oriented Rust parser
    if HAS_RUST:
        try:
            records = stats_core.parse_rescue(vcf_path)
            if records:
                df = pl.DataFrame(records)
                df = _cast_columns_fallback(df)
                df = _add_derived_columns_polars(df)
                return df
        except Exception as e:
            print(f"  [WARNING] Rust VCF parser failed: {e}, falling back to cyvcf2")

    # Python fallback (cyvcf2 already computes variant_type/ti_tv per-row)
    from vcf_stats.seq2neo.rescue_parser import parse_rescue_vcf as _py_parse
    return _py_parse(vcf_path)


def parse_rescue_vcf(vcf_path: str) -> pl.DataFrame:
    """Parse a rescue VCF file, returning a polars DataFrame.

    Uses the Rust column-oriented parser when available (fastest, lowest memory).
    Falls back to row-oriented Rust parser, then to cyvcf2.
    """
    # Prefer column-oriented Rust parser
    if HAS_RUST_COLUMNS:
        try:
            return parse_rescue_vcf_columns(vcf_path)
        except Exception:
            pass  # fall through to fallback
    return _parse_rescue_vcf_fallback(vcf_path)
