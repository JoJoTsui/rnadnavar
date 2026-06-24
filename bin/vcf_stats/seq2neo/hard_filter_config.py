"""Hard filter conditions as a Python config module.

Defines variant quality flag conditions that are observed during the
statistics stage (not dropped). Each condition produces a boolean flag
column and contributes to the hard_filter_flags and n_hard_flags summary
columns.

These flags are for OBSERVATION in the statistics stage.
Actual exclusion happens downstream in dataset preparation.

Pattern matches bin/common/tier_config.py — Python data with lambda-based
expression builders, version-controlled alongside the code that uses it.

Usage:
    from vcf_stats.seq2neo.hard_filter_config import (
        HARD_FILTER_CONDITIONS,
        build_hard_filter_flag_exprs,
        build_hard_filter_summary_exprs,
    )
"""

import polars as pl
from typing import Any

# ── Hard filter conditions ──────────────────────────────────────────────────
#
# Each condition has:
#   name            - Short kebab-case identifier
#   description     - Human-readable explanation
#   severity        - high / medium / low (used for chart coloring)
#   flag_column     - Boolean column name: flag_hard_<name>
#   required_columns - Columns this condition depends on
#   build_expr      - fn(available_columns) -> pl.Expr | None
#
# build_expr receives the set of available columns and returns None if
# required columns are missing, allowing graceful degradation when
# parquet schemas differ across runs.

HARD_FILTER_CONDITIONS: list[dict[str, Any]] = [
    {
        "name": "no_caller_support",
        "description": "No DNA/RNA caller supports this variant (N_SUPPORT_CALLERS == 0)",
        "severity": "high",
        "flag_column": "flag_hard_no_support",
        "required_columns": ["N_SUPPORT_CALLERS"],
        "build_expr": lambda cols: (
            pl.col("N_SUPPORT_CALLERS") == 0
            if "N_SUPPORT_CALLERS" in cols else None
        ),
    },
    {
        "name": "vaf_overflow",
        "description": "Sum of VAF across alleles > 1.1 (physics violation)",
        "severity": "high",
        "flag_column": "flag_hard_vaf_overflow",
        "required_columns": ["flag_vaf_overflow"],
        "build_expr": lambda cols: (
            pl.col("flag_vaf_overflow")
            if "flag_vaf_overflow" in cols else None
        ),
    },
    {
        "name": "noise_allele",
        "description": "Multiallelic site with one real allele + sequencing noise",
        "severity": "high",
        "flag_column": "flag_hard_noise_allele",
        "required_columns": ["multiallelic_class"],
        "build_expr": lambda cols: (
            pl.col("multiallelic_class") == "noise"
            if "multiallelic_class" in cols else None
        ),
    },
    {
        "name": "no_coverage",
        "description": "No modality has minimum coverage (DNA_DP < 5 AND RNA_DP < 5)",
        "severity": "medium",
        "flag_column": "flag_hard_no_coverage",
        "required_columns": ["DNA_DP_mean", "RNA_DP_mean"],
        "build_expr": lambda cols: (
            (pl.col("DNA_DP_mean") < 5) & (pl.col("RNA_DP_mean") < 5)
            if "DNA_DP_mean" in cols and "RNA_DP_mean" in cols else None
        ),
    },
    {
        "name": "no_alt_evidence",
        "description": "No alt allele evidence in either modality (DNA_ALT_DP < 2 AND RNA_ALT_DP < 2)",
        "severity": "medium",
        "flag_column": "flag_hard_no_alt_evidence",
        "required_columns": ["DNA_ALT_DP_mean", "RNA_ALT_DP_mean"],
        "build_expr": lambda cols: (
            (pl.col("DNA_ALT_DP_mean") < 2) & (pl.col("RNA_ALT_DP_mean") < 2)
            if "DNA_ALT_DP_mean" in cols and "RNA_ALT_DP_mean" in cols else None
        ),
    },
    {
        "name": "germline_low_vaf",
        "description": "Called Germline but VAF is suspiciously low",
        "severity": "medium",
        "flag_column": "flag_hard_germline_low_vaf",
        "required_columns": ["flag_germline_low_vaf"],
        "build_expr": lambda cols: (
            pl.col("flag_germline_low_vaf")
            if "flag_germline_low_vaf" in cols else None
        ),
    },
    {
        "name": "somatic_loh",
        "description": "Called Somatic with loss-of-heterozygosity pattern",
        "severity": "medium",
        "flag_column": "flag_hard_somatic_loh",
        "required_columns": ["flag_somatic_loh"],
        "build_expr": lambda cols: (
            pl.col("flag_somatic_loh")
            if "flag_somatic_loh" in cols else None
        ),
    },
    {
        "name": "reference_with_signal",
        "description": "Called Reference but has substantive alt signal (VAF > 0.10, ALT_DP ≥ 3)",
        "severity": "low",
        "flag_column": "flag_hard_reference_signal",
        "required_columns": ["FILTER", "DNA_VAF_mean", "DNA_ALT_DP_mean"],
        "build_expr": lambda cols: (
            (pl.col("FILTER") == "Reference")
            & (pl.col("DNA_VAF_mean") > 0.10)
            & (pl.col("DNA_ALT_DP_mean") >= 3)
            if all(c in cols for c in ["FILTER", "DNA_VAF_mean", "DNA_ALT_DP_mean"])
            else None
        ),
    },
]


def load_hard_filter_config() -> list[dict[str, Any]]:
    """Return the hard filter conditions list."""
    return HARD_FILTER_CONDITIONS


def build_hard_filter_flag_exprs(columns: set[str] | list[str]) -> list[pl.Expr]:
    """Build flag column expressions for all applicable conditions.

    Returns a list of pl.Expr that evaluate to boolean values (True when
    the condition is met). Only includes conditions whose required columns
    are present in the input column set.

    Args:
        columns: Available column names from the parquet schema.

    Returns:
        List of boolean pl.Expr, each aliased to its flag_column name.
    """
    cols = set(columns)
    exprs: list[pl.Expr] = []
    for cond in HARD_FILTER_CONDITIONS:
        expr = cond["build_expr"](cols)
        if expr is not None:
            exprs.append(expr.alias(cond["flag_column"]))
    return exprs


def build_hard_filter_summary_exprs(columns: set[str] | list[str]) -> list[pl.Expr]:
    """Build summary column expressions: hard_filter_flags and n_hard_flags.

    hard_filter_flags: comma-joined string of active flag column names.
    n_hard_flags: integer count of active hard filter flags.

    Args:
        columns: Available column names (must include the flag columns
                 produced by build_hard_filter_flag_exprs).

    Returns:
        List of two pl.Expr: [hard_filter_flags, n_hard_flags].
    """
    flag_cols = [
        cond["flag_column"]
        for cond in HARD_FILTER_CONDITIONS
        if cond["flag_column"] in columns
    ]

    # Build comma-joined flags string: pl.col("a").cast(str) + "," + ...
    flags_str_expr: pl.Expr | None = None
    for fc in flag_cols:
        col_as_flag = pl.when(pl.col(fc)).then(pl.lit(fc.replace("flag_hard_", ""))).otherwise(pl.lit(""))
        if flags_str_expr is None:
            flags_str_expr = col_as_flag
        else:
            flags_str_expr = (
                pl.when(pl.col(fc))
                .then(
                    flags_str_expr + pl.lit(",") + pl.lit(fc.replace("flag_hard_", ""))
                )
                .otherwise(flags_str_expr)
            )

    if flags_str_expr is None:
        flags_str_expr = pl.lit("")

    # Strip leading/trailing commas (from the conditional concatenation)
    # Simpler approach: build with concat_str and filter empty strings
    flag_name_exprs = [
        pl.when(pl.col(fc)).then(pl.lit(fc.replace("flag_hard_", ""))).otherwise(pl.lit(None))
        for fc in flag_cols
    ]
    flags_str_expr = pl.concat_str(flag_name_exprs, separator=",", ignore_nulls=True)

    # Count of active flags: sum of boolean columns (cast to int)
    n_flags_expr: pl.Expr = pl.lit(0)
    for fc in flag_cols:
        n_flags_expr = n_flags_expr + pl.col(fc).cast(pl.Int32)

    return [
        flags_str_expr.alias("hard_filter_flags"),
        n_flags_expr.alias("n_hard_flags"),
    ]


def get_hard_filter_flag_columns() -> list[str]:
    """Return the list of all flag column names."""
    return [cond["flag_column"] for cond in HARD_FILTER_CONDITIONS]
