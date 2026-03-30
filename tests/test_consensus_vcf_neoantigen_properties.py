"""
Property-based tests for AD_BY_CALLER and AF_BY_CALLER Presence.

# Feature: neoantigen-workflow, Property 6: AD_BY_CALLER and AF_BY_CALLER Presence

Validates: Requirements 2.10, 2.11, 2.12

For any consensus VCF produced when run_consensus_vcf.py is invoked with --neoantigen,
every variant record must contain both AD_BY_CALLER and AF_BY_CALLER INFO fields.
Conversely, for any consensus VCF produced without --neoantigen, neither field must
appear in the header or records.
"""

import os
import sys
import tempfile

import pytest
from cyvcf2 import VCF
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

# Ensure bin/ is importable
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin"))

from vcf_utils.io_utils import create_output_header, write_union_vcf  # noqa: E402

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

CHROMOSOMES = ["chr1", "chr2", "chr7", "chr17", "chrX"]

# ---------------------------------------------------------------------------
# Hypothesis strategies
# ---------------------------------------------------------------------------

chrom_st = st.sampled_from(CHROMOSOMES)
pos_st = st.integers(min_value=1, max_value=250_000_000)
depth_st = st.integers(min_value=1, max_value=500)
vaf_st = st.floats(min_value=0.0, max_value=1.0, allow_nan=False, allow_infinity=False)

snp_pairs_st = st.sampled_from(
    [
        ("A", "T"),
        ("A", "C"),
        ("A", "G"),
        ("T", "A"),
        ("T", "C"),
        ("T", "G"),
        ("C", "A"),
        ("C", "T"),
        ("C", "G"),
        ("G", "A"),
        ("G", "T"),
        ("G", "C"),
    ]
)

callers_st = st.lists(
    st.sampled_from(["mutect2", "strelka", "deepsomatic", "sage"]),
    min_size=1,
    max_size=4,
    unique=True,
)

# ---------------------------------------------------------------------------
# Minimal VCF template header builder
# ---------------------------------------------------------------------------


def _build_template_vcf_content(chrom: str = "chr1") -> str:
    """Build a minimal VCF content string to use as a template header."""
    return (
        "##fileformat=VCFv4.2\n"
        "##FILTER=<ID=PASS,Description=\"All filters passed\">\n"
        f"##contig=<ID={chrom},length=248956422>\n"
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
        "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">\n"
        "##FORMAT=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">\n"
        "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">\n"
        "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR\n"
    )


def _make_template_vcf(chrom: str = "chr1") -> "VCF":
    """Write a minimal VCF to a temp file and return a cyvcf2.VCF handle."""
    content = _build_template_vcf_content(chrom)
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".vcf", delete=False, prefix="test_template_"
    ) as f:
        f.write(content)
        path = f.name
    return VCF(path), path


# ---------------------------------------------------------------------------
# Variant data structure builder
# ---------------------------------------------------------------------------


def _build_variant_data(
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    callers: list,
    ref_depth: int,
    alt_depth: int,
    vaf: float,
) -> dict:
    """
    Build a minimal variant_data dict matching what aggregate_variants() returns,
    suitable for passing to write_union_vcf().

    Each caller gets an AD entry (ref,alt) and a VAF entry.
    """
    total_dp = ref_depth + alt_depth
    ad_str = f"{ref_depth},{alt_depth}"

    genotypes = {}
    for caller in callers:
        genotypes[caller] = {
            "GT": "0/1",
            "DP": total_dp,
            "VAF": vaf,
            "AD": ad_str,
        }

    gt_by_caller = ["0/1"] * len(callers)
    dp_by_caller = [total_dp] * len(callers)
    vaf_by_caller = [vaf] * len(callers)

    key = f"{chrom.replace('chr', '')}:{pos}:{ref}:{alt}"

    return {
        key: {
            "CHROM": chrom,
            "POS": pos,
            "REF": ref,
            "ALT": alt,
            "ids": [],
            "qualities": [50.0],
            "callers": list(callers),
            "filters_original": ["PASS"] * len(callers),
            "filters_normalized": ["Somatic"] * len(callers),
            "filters_category": ["Somatic"] * len(callers),
            "passes_consensus": True,
            "rescued": False,
            "genotypes": genotypes,
            "gt_aggregated": {
                "consensus_gt": "0/1",
                "gt_by_caller": gt_by_caller,
                "dp_values": dp_by_caller,
                "dp_mean": float(total_dp),
                "dp_min": total_dp,
                "dp_max": total_dp,
                "dp_by_caller": dp_by_caller,
                "vaf_values": vaf_by_caller,
                "vaf_mean": vaf,
                "vaf_min": vaf,
                "vaf_max": vaf,
                "vaf_by_caller": vaf_by_caller,
            },
        }
    }


# ---------------------------------------------------------------------------
# Helper: run write_union_vcf and return output path
# ---------------------------------------------------------------------------


def _run_write_union_vcf(
    variant_data: dict,
    callers: list,
    neoantigen: bool,
    chrom: str = "chr1",
) -> str:
    """
    Write variant_data to a temp VCF using write_union_vcf() and return the output path.
    """
    template_vcf, template_path = _make_template_vcf(chrom)

    out_fd, out_path = tempfile.mkstemp(suffix=".vcf", prefix="test_consensus_out_")
    os.close(out_fd)

    try:
        write_union_vcf(
            variant_data=variant_data,
            template_header=template_vcf,
            sample_name="TUMOR",
            out_file=out_path,
            output_format="vcf",
            all_callers=callers,
            modality_map=None,
            snv_threshold=1,
            indel_threshold=1,
            include_non_canonical=True,
            neoantigen=neoantigen,
        )
    finally:
        template_vcf.close()
        try:
            os.unlink(template_path)
        except OSError:
            pass

    return out_path


# ---------------------------------------------------------------------------
# Assertion helpers
# ---------------------------------------------------------------------------


def _get_info_field_ids(vcf_path: str) -> set:
    """Return the set of INFO field IDs declared in the VCF header."""
    vcf = VCF(vcf_path)
    ids = {h["ID"] for h in vcf.header_iter() if h["HeaderType"] == "INFO"}
    vcf.close()
    return ids


def _get_record_info_keys(vcf_path: str) -> list:
    """Return a list of sets, one per record, of INFO keys present in that record."""
    vcf = VCF(vcf_path)
    # Collect all INFO field IDs declared in the header
    all_info_ids = [h["ID"] for h in vcf.header_iter() if h["HeaderType"] == "INFO"]
    result = []
    for variant in vcf:
        # Check which declared INFO fields are actually set on this record
        present = set()
        for field_id in all_info_ids:
            try:
                val = variant.INFO.get(field_id)
                if val is not None:
                    present.add(field_id)
            except Exception:
                pass
        result.append(present)
    vcf.close()
    return result


# ---------------------------------------------------------------------------
# Property 6: AD_BY_CALLER and AF_BY_CALLER Presence — neoantigen=True
# ---------------------------------------------------------------------------


@given(
    chrom=chrom_st,
    pos=pos_st,
    bases=snp_pairs_st,
    callers=callers_st,
    ref_depth=depth_st,
    alt_depth=depth_st,
    vaf=vaf_st,
)
@settings(max_examples=100, suppress_health_check=[HealthCheck.function_scoped_fixture])
def test_p6_neoantigen_true_header_declares_ad_af_by_caller(
    chrom, pos, bases, callers, ref_depth, alt_depth, vaf
):
    """
    **Validates: Requirements 2.10, 2.11, 2.12**

    Property 6: AD_BY_CALLER and AF_BY_CALLER Presence — neoantigen=True, header check.
    When write_union_vcf() is called with neoantigen=True, the output VCF header
    MUST declare both AD_BY_CALLER and AF_BY_CALLER as INFO fields.

    # Feature: neoantigen-workflow, Property 6: AD_BY_CALLER and AF_BY_CALLER Presence
    """
    ref, alt = bases
    variant_data = _build_variant_data(
        chrom, pos, ref, alt, callers, ref_depth, alt_depth, vaf
    )
    out_path = _run_write_union_vcf(variant_data, callers, neoantigen=True, chrom=chrom)
    try:
        info_ids = _get_info_field_ids(out_path)
        assert "AD_BY_CALLER" in info_ids, (
            f"AD_BY_CALLER not declared in header INFO fields. "
            f"Declared INFO fields: {info_ids}"
        )
        assert "AF_BY_CALLER" in info_ids, (
            f"AF_BY_CALLER not declared in header INFO fields. "
            f"Declared INFO fields: {info_ids}"
        )
    finally:
        try:
            os.unlink(out_path)
        except OSError:
            pass


@given(
    chrom=chrom_st,
    pos=pos_st,
    bases=snp_pairs_st,
    callers=callers_st,
    ref_depth=depth_st,
    alt_depth=depth_st,
    vaf=vaf_st,
)
@settings(max_examples=100, suppress_health_check=[HealthCheck.function_scoped_fixture])
def test_p6_neoantigen_true_records_contain_ad_af_by_caller(
    chrom, pos, bases, callers, ref_depth, alt_depth, vaf
):
    """
    **Validates: Requirements 2.10, 2.11, 2.12**

    Property 6: AD_BY_CALLER and AF_BY_CALLER Presence — neoantigen=True, record check.
    When write_union_vcf() is called with neoantigen=True, every variant record
    MUST contain both AD_BY_CALLER and AF_BY_CALLER in its INFO column.

    # Feature: neoantigen-workflow, Property 6: AD_BY_CALLER and AF_BY_CALLER Presence
    """
    ref, alt = bases
    variant_data = _build_variant_data(
        chrom, pos, ref, alt, callers, ref_depth, alt_depth, vaf
    )
    out_path = _run_write_union_vcf(variant_data, callers, neoantigen=True, chrom=chrom)
    try:
        record_info_keys = _get_record_info_keys(out_path)
        assert len(record_info_keys) > 0, "No variant records written to output VCF"
        for i, keys in enumerate(record_info_keys):
            assert "AD_BY_CALLER" in keys, (
                f"Record {i}: AD_BY_CALLER missing from INFO. "
                f"Present INFO keys: {keys}"
            )
            assert "AF_BY_CALLER" in keys, (
                f"Record {i}: AF_BY_CALLER missing from INFO. "
                f"Present INFO keys: {keys}"
            )
    finally:
        try:
            os.unlink(out_path)
        except OSError:
            pass


# ---------------------------------------------------------------------------
# Property 6: AD_BY_CALLER and AF_BY_CALLER Presence — neoantigen=False
# ---------------------------------------------------------------------------


@given(
    chrom=chrom_st,
    pos=pos_st,
    bases=snp_pairs_st,
    callers=callers_st,
    ref_depth=depth_st,
    alt_depth=depth_st,
    vaf=vaf_st,
)
@settings(max_examples=100, suppress_health_check=[HealthCheck.function_scoped_fixture])
def test_p6_neoantigen_false_header_excludes_ad_af_by_caller(
    chrom, pos, bases, callers, ref_depth, alt_depth, vaf
):
    """
    **Validates: Requirements 2.10, 2.11, 2.12**

    Property 6: AD_BY_CALLER and AF_BY_CALLER Presence — neoantigen=False, header check.
    When write_union_vcf() is called with neoantigen=False, the output VCF header
    MUST NOT declare AD_BY_CALLER or AF_BY_CALLER as INFO fields.

    # Feature: neoantigen-workflow, Property 6: AD_BY_CALLER and AF_BY_CALLER Presence
    """
    ref, alt = bases
    variant_data = _build_variant_data(
        chrom, pos, ref, alt, callers, ref_depth, alt_depth, vaf
    )
    out_path = _run_write_union_vcf(
        variant_data, callers, neoantigen=False, chrom=chrom
    )
    try:
        info_ids = _get_info_field_ids(out_path)
        assert "AD_BY_CALLER" not in info_ids, (
            f"AD_BY_CALLER unexpectedly declared in header when neoantigen=False. "
            f"Declared INFO fields: {info_ids}"
        )
        assert "AF_BY_CALLER" not in info_ids, (
            f"AF_BY_CALLER unexpectedly declared in header when neoantigen=False. "
            f"Declared INFO fields: {info_ids}"
        )
    finally:
        try:
            os.unlink(out_path)
        except OSError:
            pass


@given(
    chrom=chrom_st,
    pos=pos_st,
    bases=snp_pairs_st,
    callers=callers_st,
    ref_depth=depth_st,
    alt_depth=depth_st,
    vaf=vaf_st,
)
@settings(max_examples=100, suppress_health_check=[HealthCheck.function_scoped_fixture])
def test_p6_neoantigen_false_records_exclude_ad_af_by_caller(
    chrom, pos, bases, callers, ref_depth, alt_depth, vaf
):
    """
    **Validates: Requirements 2.10, 2.11, 2.12**

    Property 6: AD_BY_CALLER and AF_BY_CALLER Presence — neoantigen=False, record check.
    When write_union_vcf() is called with neoantigen=False, no variant record
    MUST contain AD_BY_CALLER or AF_BY_CALLER in its INFO column.

    # Feature: neoantigen-workflow, Property 6: AD_BY_CALLER and AF_BY_CALLER Presence
    """
    ref, alt = bases
    variant_data = _build_variant_data(
        chrom, pos, ref, alt, callers, ref_depth, alt_depth, vaf
    )
    out_path = _run_write_union_vcf(
        variant_data, callers, neoantigen=False, chrom=chrom
    )
    try:
        record_info_keys = _get_record_info_keys(out_path)
        for i, keys in enumerate(record_info_keys):
            assert "AD_BY_CALLER" not in keys, (
                f"Record {i}: AD_BY_CALLER unexpectedly present in INFO when neoantigen=False. "
                f"Present INFO keys: {keys}"
            )
            assert "AF_BY_CALLER" not in keys, (
                f"Record {i}: AF_BY_CALLER unexpectedly present in INFO when neoantigen=False. "
                f"Present INFO keys: {keys}"
            )
    finally:
        try:
            os.unlink(out_path)
        except OSError:
            pass
