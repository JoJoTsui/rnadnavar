"""Python wrapper for Rust VCF parser with cyvcf2 fallback.

Uses stats_core (Rust) when available, falls back to rescue_parser.py (cyvcf2)
when the Rust module is not built.
"""

import polars as pl

try:
    import stats_core
    HAS_RUST = True
except ImportError:
    HAS_RUST = False


def parse_rescue_vcf(vcf_path: str) -> pl.DataFrame:
    """Parse a rescue VCF file, returning a polars DataFrame.

    Uses Rust (noodles-vcf) when available, falls back to cyvcf2.
    """
    if HAS_RUST:
        try:
            records = stats_core.parse_rescue(vcf_path)
            if records:
                # Convert list of dicts to polars DataFrame
                return pl.DataFrame(records)
        except Exception as e:
            print(f"  [WARNING] Rust VCF parser failed: {e}, falling back to cyvcf2")

    # Python fallback
    from vcf_stats.seq2neo.rescue_parser import parse_rescue_vcf as _py_parse
    return _py_parse(vcf_path)
