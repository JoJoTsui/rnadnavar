"""Python wrapper for Rust VCF parser with cyvcf2 fallback.

Uses stats_core (Rust) when available, falls back to rescue_parser.py (cyvcf2)
when the Rust module is not built.
"""

import polars as pl

from .rescue_parser import (
    ALL_RESCUE_FIELDS,
    RESCUE_FLAG_FIELDS,
    RESCUE_FLOAT_FIELDS,
    RESCUE_INT_FIELDS,
)

try:
    import stats_core
    HAS_RUST = True
except ImportError:
    HAS_RUST = False


def _cast_columns(df: pl.DataFrame) -> pl.DataFrame:
    """Cast INFO columns to their correct types (matching rescue_parser.py output).

    The Rust VCF parser returns all values as strings. This function casts
    columns to the types expected by downstream statistics code.
    Also fills missing fields that are in ALL_RESCUE_FIELDS but not discovered
    by the Rust parser (fields present in VCF records but not declared in header).
    """
    # Fill missing fields that Rust parser didn't discover from header
    for field in ALL_RESCUE_FIELDS:
        if field not in df.columns:
            df = df.with_columns(pl.lit(None).alias(field))
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
            # Flag fields from Rust are strings ("true" or "" from INFO presence/absence)
            # Utf8View → Boolean direct cast not supported; cast via Utf8 first
            df = df.with_columns(
                pl.when(pl.col(field).cast(pl.Utf8) == "true")
                .then(True)
                .otherwise(False)
                .alias(field)
            )

    return df


def parse_rescue_vcf(vcf_path: str) -> pl.DataFrame:
    """Parse a rescue VCF file, returning a polars DataFrame.

    Uses Rust (noodles-vcf) when available, falls back to cyvcf2.
    Column types are cast to match the cyvcf2 parser output.
    """
    if HAS_RUST:
        try:
            records = stats_core.parse_rescue(vcf_path)
            if records:
                df = pl.DataFrame(records)
                df = _cast_columns(df)
                return df
        except Exception as e:
            print(f"  [WARNING] Rust VCF parser failed: {e}, falling back to cyvcf2")

    # Python fallback
    from vcf_stats.seq2neo.rescue_parser import parse_rescue_vcf as _py_parse
    return _py_parse(vcf_path)
