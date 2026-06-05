"""Seq2neo variant statistics package.

Integrated into bin/vcf_stats/ — provides tools for parsing seq2neo pipeline
output VCFs, computing statistics, validating rescue VCF correctness, and
generating interactive visualizations.

Modules:
    manifest_loader  — Load sample manifest and construct VCF paths
    rescue_parser    — Parse rescue VCF for metadata/annotations
    caller_parser    — Parse caller VCFs for DP/AD/GT (ground truth)
    statistics       — Compute variant statistics using polars
    rescue_validator — Cross-validate rescue VCF INFO against caller ground truth
    visualizer       — Plotly express dashboard + static exports
    cli              — CLI entry point
"""

from .manifest_loader import (
    filter_complete,
    get_all_caller_vcf_paths,
    get_vcf_prefix,
    load_manifest,
)
from .rescue_parser import parse_rescue_vcf, rescue_info_fields

__all__ = [
    "load_manifest",
    "get_vcf_prefix",
    "filter_complete",
    "get_all_caller_vcf_paths",
    "parse_rescue_vcf",
    "rescue_info_fields",
]
