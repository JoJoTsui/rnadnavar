"""
Property-based tests for FORMAT Harmonization Completeness.

# Feature: neoantigen-workflow, Property 2: FORMAT Harmonization Completeness

Validates: Requirements 2.2, 2.5

For any per-caller normalized VCF processed by FORMAT_HARMONIZER, the output VCF
must contain GT, AD, AF, DP, and GQ as declared FORMAT fields in both the header
and every variant record's FORMAT column.
"""

import os
import sys
import tempfile

import pytest
from cyvcf2 import VCF
from hypothesis import given, settings, HealthCheck
from hypothesis import strategies as st

# Ensure bin/ is importable
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin"))
from harmonize_vcf_format import harmonize  # noqa: E402

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

REQUIRED_FORMAT_FIELDS = {"GT", "AD", "AF", "DP", "GQ"}

CHROMOSOMES = ["chr1", "chr2", "chr7", "chr17", "chrX"]

# ---------------------------------------------------------------------------
# VCF template builders
# ---------------------------------------------------------------------------


def _format_header_lines(fields: list[tuple[str, str, str, str]]) -> str:
    """Build FORMAT header lines from (ID, Number, Type, Description) tuples."""
    lines = []
    for fid, number, ftype, desc in fields:
        lines.append(f'##FORMAT=<ID={fid},Number={number},Type={ftype},Description="{desc}">')
    return "\n".join(lines) + "\n"


def _vcf_header_base(fmt_fields: list[tuple[str, str, str, str]], sample_name: str = "TUMOR") -> str:
    """Build a complete VCF header with ##fileformat first, then FORMAT lines, then column header."""
    fmt_lines = _format_header_lines(fmt_fields)
    return (
        "##fileformat=VCFv4.2\n"
        "##FILTER=<ID=PASS,Description=\"All filters passed\">\n"
        "##contig=<ID=chr1,length=248956422>\n"
        + fmt_lines
        + f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}\n"
    )


def _variant_line(chrom, pos, ref, alt, fmt_keys, fmt_values) -> str:
    """Build a single VCF variant line."""
    fmt_str = ":".join(fmt_keys)
    val_str = ":".join(str(v) for v in fmt_values)
    return f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\t{fmt_str}\t{val_str}\n"


# ---------------------------------------------------------------------------
# Per-caller VCF builders
# ---------------------------------------------------------------------------


def build_mutect2_vcf(chrom, pos, ref, alt, ref_depth, alt_depth, gq) -> str:
    total = ref_depth + alt_depth
    af = round(alt_depth / total, 4) if total > 0 else 0.0
    header = _vcf_header_base([
        ("GT", "1", "String", "Genotype"),
        ("AD", "R", "Integer", "Allelic depths for the ref and alt alleles"),
        ("AF", "A", "Float", "Allele fractions of alternate alleles in the tumor"),
        ("DP", "1", "Integer", "Approximate read depth"),
        ("GQ", "1", "Integer", "Genotype Quality"),
        ("F1R2", "R", "Integer", "Count of reads in F1R2 pair orientation"),
        ("F2R1", "R", "Integer", "Count of reads in F2R1 pair orientation"),
    ])
    variant = _variant_line(
        chrom, pos, ref, alt,
        ["GT", "AD", "AF", "DP", "GQ", "F1R2", "F2R1"],
        ["0/1", f"{ref_depth},{alt_depth}", str(af), str(total), str(gq),
         f"{ref_depth},{0}", f"{0},{alt_depth}"],
    )
    return header + variant


def build_strelka_snv_vcf(chrom, pos, ref, alt, au, cu, gu, tu, gq) -> str:
    """Strelka2 SNV VCF with AU/CU/GU/TU tier counts (each is tier1,tier2)."""
    header = _vcf_header_base([
        ("GT", "1", "String", "Genotype"),
        ("AU", "2", "Integer", "Number of 'A' alleles used in tiers 1,2"),
        ("CU", "2", "Integer", "Number of 'C' alleles used in tiers 1,2"),
        ("GU", "2", "Integer", "Number of 'G' alleles used in tiers 1,2"),
        ("TU", "2", "Integer", "Number of 'T' alleles used in tiers 1,2"),
        ("DP", "1", "Integer", "Read depth for tier1"),
        ("GQ", "1", "Integer", "Genotype Quality"),
    ])
    total = au + cu + gu + tu
    variant = _variant_line(
        chrom, pos, ref, alt,
        ["GT", "AU", "CU", "GU", "TU", "DP", "GQ"],
        ["0/1", f"{au},0", f"{cu},0", f"{gu},0", f"{tu},0", str(total), str(gq)],
    )
    return header + variant


def build_strelka_indel_vcf(chrom, pos, ref, alt, tar, tir, gq) -> str:
    """Strelka2 Indel VCF with TAR/TIR counts (each is tier1,tier2)."""
    header = _vcf_header_base([
        ("GT", "1", "String", "Genotype"),
        ("TAR", "2", "Integer", "Reads strongly supporting the ref allele tier 1,2"),
        ("TIR", "2", "Integer", "Reads strongly supporting the indel allele tier 1,2"),
        ("DP", "1", "Integer", "Read depth for tier1"),
        ("GQ", "1", "Integer", "Genotype Quality"),
    ])
    total = tar + tir
    variant = _variant_line(
        chrom, pos, ref, alt,
        ["GT", "TAR", "TIR", "DP", "GQ"],
        ["0/1", f"{tar},0", f"{tir},0", str(total), str(gq)],
    )
    return header + variant


def build_sage_vcf(chrom, pos, ref, alt, ref_depth, alt_depth, gq) -> str:
    """SAGE VCF with AD/DP/VF (VF is the allele frequency field in SAGE)."""
    total = ref_depth + alt_depth
    vf = round(alt_depth / total, 4) if total > 0 else 0.0
    header = _vcf_header_base([
        ("GT", "1", "String", "Genotype"),
        ("AD", "R", "Integer", "Allelic depths for the ref and alt alleles"),
        ("DP", "1", "Integer", "Approximate read depth"),
        ("VF", "1", "Float", "Variant allele frequency"),
        ("GQ", "1", "Integer", "Genotype Quality"),
    ])
    variant = _variant_line(
        chrom, pos, ref, alt,
        ["GT", "AD", "DP", "VF", "GQ"],
        ["0/1", f"{ref_depth},{alt_depth}", str(total), str(vf), str(gq)],
    )
    return header + variant


def build_deepsomatic_vcf(chrom, pos, ref, alt, ref_depth, alt_depth, gq) -> str:
    """DeepSomatic VCF with AD/DP/VAF (VAF is the allele frequency field)."""
    total = ref_depth + alt_depth
    vaf = round(alt_depth / total, 4) if total > 0 else 0.0
    header = _vcf_header_base([
        ("GT", "1", "String", "Genotype"),
        ("AD", "R", "Integer", "Allelic depths for the ref and alt alleles"),
        ("DP", "1", "Integer", "Approximate read depth"),
        ("VAF", "A", "Float", "Variant allele fractions"),
        ("GQ", "1", "Integer", "Genotype Quality"),
    ])
    variant = _variant_line(
        chrom, pos, ref, alt,
        ["GT", "AD", "DP", "VAF", "GQ"],
        ["0/1", f"{ref_depth},{alt_depth}", str(total), str(vaf), str(gq)],
    )
    return header + variant


# ---------------------------------------------------------------------------
# Hypothesis strategies
# ---------------------------------------------------------------------------

chrom_st = st.sampled_from(CHROMOSOMES)
pos_st = st.integers(min_value=1, max_value=250_000_000)
depth_st = st.integers(min_value=0, max_value=500)
gq_st = st.integers(min_value=0, max_value=99)

# SNP bases
snp_pairs_st = st.sampled_from([
    ("A", "T"), ("A", "C"), ("A", "G"),
    ("T", "A"), ("T", "C"), ("T", "G"),
    ("C", "A"), ("C", "T"), ("C", "G"),
    ("G", "A"), ("G", "T"), ("G", "C"),
])

# Indel pairs (ref longer than alt or vice versa)
indel_pairs_st = st.sampled_from([
    ("AT", "A"), ("ACG", "A"), ("A", "AT"), ("A", "ACG"),
    ("ATCG", "A"), ("A", "ATCG"),
])


# ---------------------------------------------------------------------------
# Helper: run harmonizer and return output VCF path
# ---------------------------------------------------------------------------


def run_harmonizer(vcf_content: str, caller: str) -> str:
    """Write vcf_content to a temp file, run harmonize(), return output path."""
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".vcf", delete=False, prefix=f"test_{caller}_in_"
    ) as f:
        f.write(vcf_content)
        input_path = f.name

    output_path = input_path.replace("_in_", "_out_").replace(".vcf", ".harmonized.vcf")

    try:
        harmonize(input_path, caller, output_path)
    finally:
        try:
            os.unlink(input_path)
        except OSError:
            pass

    return output_path


# ---------------------------------------------------------------------------
# Assertion helpers
# ---------------------------------------------------------------------------


def assert_format_header_completeness(output_vcf_path: str, caller: str):
    """Assert GT, AD, AF, DP, GQ are declared in the FORMAT header."""
    vcf = VCF(output_vcf_path)
    declared_format_ids = {h["ID"] for h in vcf.header_iter() if h["HeaderType"] == "FORMAT"}
    vcf.close()

    missing = REQUIRED_FORMAT_FIELDS - declared_format_ids
    assert not missing, (
        f"[{caller}] Missing FORMAT header declarations: {missing}. "
        f"Declared: {declared_format_ids}"
    )


def assert_format_record_completeness(output_vcf_path: str, caller: str):
    """Assert every variant record's FORMAT column contains AD, AF, DP, GQ."""
    # GT is always present as the first field in cyvcf2 records
    required_in_records = {"AD", "AF", "DP", "GQ"}

    vcf = VCF(output_vcf_path)
    for i, variant in enumerate(vcf):
        record_fmt = set(variant.FORMAT)
        missing = required_in_records - record_fmt
        assert not missing, (
            f"[{caller}] Record {i} at {variant.CHROM}:{variant.POS} "
            f"missing FORMAT fields: {missing}. "
            f"Record FORMAT: {record_fmt}"
        )
    vcf.close()


# ---------------------------------------------------------------------------
# Property tests
# ---------------------------------------------------------------------------


@given(
    chrom=chrom_st,
    pos=pos_st,
    bases=snp_pairs_st,
    ref_depth=depth_st,
    alt_depth=depth_st,
    gq=gq_st,
)
@settings(max_examples=100, suppress_health_check=[HealthCheck.function_scoped_fixture])
def test_p2_mutect2_format_completeness(chrom, pos, bases, ref_depth, alt_depth, gq):
    """
    **Validates: Requirements 2.2, 2.5**

    Property 2: FORMAT Harmonization Completeness — mutect2 caller.
    For any mutect2 VCF, the harmonized output must declare GT/AD/AF/DP/GQ
    in the FORMAT header and include AD/AF/DP/GQ in every record's FORMAT column.
    """
    ref, alt = bases
    vcf_content = build_mutect2_vcf(chrom, pos, ref, alt, ref_depth, alt_depth, gq)
    output_path = run_harmonizer(vcf_content, "mutect2")
    try:
        assert_format_header_completeness(output_path, "mutect2")
        assert_format_record_completeness(output_path, "mutect2")
    finally:
        try:
            os.unlink(output_path)
        except OSError:
            pass


@given(
    chrom=chrom_st,
    pos=pos_st,
    bases=snp_pairs_st,
    au=depth_st,
    cu=depth_st,
    gu=depth_st,
    tu=depth_st,
    gq=gq_st,
)
@settings(max_examples=100, suppress_health_check=[HealthCheck.function_scoped_fixture])
def test_p2_strelka_snv_format_completeness(chrom, pos, bases, au, cu, gu, tu, gq):
    """
    **Validates: Requirements 2.2, 2.5**

    Property 2: FORMAT Harmonization Completeness — strelka_snv caller.
    For any strelka SNV VCF, the harmonized output must declare GT/AD/AF/DP/GQ
    in the FORMAT header and include AD/AF/DP/GQ in every record's FORMAT column.
    """
    ref, alt = bases
    vcf_content = build_strelka_snv_vcf(chrom, pos, ref, alt, au, cu, gu, tu, gq)
    output_path = run_harmonizer(vcf_content, "strelka_snv")
    try:
        assert_format_header_completeness(output_path, "strelka_snv")
        assert_format_record_completeness(output_path, "strelka_snv")
    finally:
        try:
            os.unlink(output_path)
        except OSError:
            pass


@given(
    chrom=chrom_st,
    pos=pos_st,
    indel=indel_pairs_st,
    tar=depth_st,
    tir=depth_st,
    gq=gq_st,
)
@settings(max_examples=100, suppress_health_check=[HealthCheck.function_scoped_fixture])
def test_p2_strelka_indel_format_completeness(chrom, pos, indel, tar, tir, gq):
    """
    **Validates: Requirements 2.2, 2.5**

    Property 2: FORMAT Harmonization Completeness — strelka_indel caller.
    For any strelka Indel VCF, the harmonized output must declare GT/AD/AF/DP/GQ
    in the FORMAT header and include AD/AF/DP/GQ in every record's FORMAT column.
    """
    ref, alt = indel
    vcf_content = build_strelka_indel_vcf(chrom, pos, ref, alt, tar, tir, gq)
    output_path = run_harmonizer(vcf_content, "strelka_indel")
    try:
        assert_format_header_completeness(output_path, "strelka_indel")
        assert_format_record_completeness(output_path, "strelka_indel")
    finally:
        try:
            os.unlink(output_path)
        except OSError:
            pass


@given(
    chrom=chrom_st,
    pos=pos_st,
    bases=snp_pairs_st,
    ref_depth=depth_st,
    alt_depth=depth_st,
    gq=gq_st,
)
@settings(max_examples=100, suppress_health_check=[HealthCheck.function_scoped_fixture])
def test_p2_sage_format_completeness(chrom, pos, bases, ref_depth, alt_depth, gq):
    """
    **Validates: Requirements 2.2, 2.5**

    Property 2: FORMAT Harmonization Completeness — sage caller.
    For any SAGE VCF, the harmonized output must declare GT/AD/AF/DP/GQ
    in the FORMAT header and include AD/AF/DP/GQ in every record's FORMAT column.
    """
    ref, alt = bases
    vcf_content = build_sage_vcf(chrom, pos, ref, alt, ref_depth, alt_depth, gq)
    output_path = run_harmonizer(vcf_content, "sage")
    try:
        assert_format_header_completeness(output_path, "sage")
        assert_format_record_completeness(output_path, "sage")
    finally:
        try:
            os.unlink(output_path)
        except OSError:
            pass


@given(
    chrom=chrom_st,
    pos=pos_st,
    bases=snp_pairs_st,
    ref_depth=depth_st,
    alt_depth=depth_st,
    gq=gq_st,
)
@settings(max_examples=100, suppress_health_check=[HealthCheck.function_scoped_fixture])
def test_p2_deepsomatic_format_completeness(chrom, pos, bases, ref_depth, alt_depth, gq):
    """
    **Validates: Requirements 2.2, 2.5**

    Property 2: FORMAT Harmonization Completeness — deepsomatic caller.
    For any DeepSomatic VCF, the harmonized output must declare GT/AD/AF/DP/GQ
    in the FORMAT header and include AD/AF/DP/GQ in every record's FORMAT column.
    """
    ref, alt = bases
    vcf_content = build_deepsomatic_vcf(chrom, pos, ref, alt, ref_depth, alt_depth, gq)
    output_path = run_harmonizer(vcf_content, "deepsomatic")
    try:
        assert_format_header_completeness(output_path, "deepsomatic")
        assert_format_record_completeness(output_path, "deepsomatic")
    finally:
        try:
            os.unlink(output_path)
        except OSError:
            pass


# ===========================================================================
# Property 3: FORMAT Harmonization Correctness
# Feature: neoantigen-workflow, Property 3: FORMAT Harmonization Correctness
# Validates: Requirements 2.3
# ===========================================================================

FLOAT_TOLERANCE = 1e-4

# ---------------------------------------------------------------------------
# Strelka2 SNV correctness
# ---------------------------------------------------------------------------

@given(
    chrom=chrom_st,
    pos=pos_st,
    bases=snp_pairs_st,
    au=depth_st,
    cu=depth_st,
    gu=depth_st,
    tu=depth_st,
    gq=gq_st,
)
@settings(max_examples=100, suppress_health_check=[HealthCheck.function_scoped_fixture])
def test_p3_strelka_snv_correctness(chrom, pos, bases, au, cu, gu, tu, gq):
    """
    **Validates: Requirements 2.3**

    Property 3: FORMAT Harmonization Correctness — strelka SNV.
    Given AU/CU/GU/TU tier-0 counts and known REF/ALT bases:
    - AD[0] == REF_BASE_U[0] (e.g., if REF=A, AD[0] == AU)
    - AD[1] == ALT_BASE_U[0] (e.g., if ALT=T, AD[1] == TU)
    - AF ≈ alt_count / (ref_count + alt_count) within float tolerance 1e-4
    - DP == AU + CU + GU + TU
    """
    # Feature: neoantigen-workflow, Property 3: FORMAT Harmonization Correctness
    ref, alt = bases
    base_counts = {"A": au, "C": cu, "G": gu, "T": tu}
    ref_count = base_counts[ref.upper()]
    alt_count = base_counts[alt.upper()]
    expected_dp = au + cu + gu + tu

    vcf_content = build_strelka_snv_vcf(chrom, pos, ref, alt, au, cu, gu, tu, gq)
    output_path = run_harmonizer(vcf_content, "strelka_snv")
    try:
        vcf = VCF(output_path)
        for variant in vcf:
            ad = variant.format("AD")
            af = variant.format("AF")
            dp = variant.format("DP")

            assert ad is not None, "AD field missing from strelka SNV output"
            assert af is not None, "AF field missing from strelka SNV output"
            assert dp is not None, "DP field missing from strelka SNV output"

            ad_vals = list(ad[0])
            af_val = float(af[0][0])
            dp_val = int(dp[0][0])

            assert ad_vals[0] == ref_count, (
                f"AD[0] (ref) expected {ref_count} (REF={ref}), got {ad_vals[0]}"
            )
            assert ad_vals[1] == alt_count, (
                f"AD[1] (alt) expected {alt_count} (ALT={alt}), got {ad_vals[1]}"
            )

            total = ref_count + alt_count
            if total > 0:
                expected_af = alt_count / total
                assert abs(af_val - expected_af) < FLOAT_TOLERANCE, (
                    f"AF expected ≈{expected_af:.6f}, got {af_val:.6f}"
                )

            assert dp_val == expected_dp, (
                f"DP expected {expected_dp} (AU+CU+GU+TU), got {dp_val}"
            )
        vcf.close()
    finally:
        try:
            os.unlink(output_path)
        except OSError:
            pass


# ---------------------------------------------------------------------------
# Strelka2 Indel correctness
# ---------------------------------------------------------------------------

@given(
    chrom=chrom_st,
    pos=pos_st,
    indel=indel_pairs_st,
    tar=depth_st,
    tir=depth_st,
    gq=gq_st,
)
@settings(max_examples=100, suppress_health_check=[HealthCheck.function_scoped_fixture])
def test_p3_strelka_indel_correctness(chrom, pos, indel, tar, tir, gq):
    """
    **Validates: Requirements 2.3**

    Property 3: FORMAT Harmonization Correctness — strelka Indel.
    Given TAR/TIR counts:
    - AD[0] == TAR, AD[1] == TIR
    - AF ≈ TIR / (TAR + TIR) when TAR+TIR > 0
    - DP == TAR + TIR
    """
    # Feature: neoantigen-workflow, Property 3: FORMAT Harmonization Correctness
    ref, alt = indel
    expected_dp = tar + tir

    vcf_content = build_strelka_indel_vcf(chrom, pos, ref, alt, tar, tir, gq)
    output_path = run_harmonizer(vcf_content, "strelka_indel")
    try:
        vcf = VCF(output_path)
        for variant in vcf:
            ad = variant.format("AD")
            af = variant.format("AF")
            dp = variant.format("DP")

            assert ad is not None, "AD field missing from strelka Indel output"
            assert af is not None, "AF field missing from strelka Indel output"
            assert dp is not None, "DP field missing from strelka Indel output"

            ad_vals = list(ad[0])
            af_val = float(af[0][0])
            dp_val = int(dp[0][0])

            assert ad_vals[0] == tar, (
                f"AD[0] (TAR) expected {tar}, got {ad_vals[0]}"
            )
            assert ad_vals[1] == tir, (
                f"AD[1] (TIR) expected {tir}, got {ad_vals[1]}"
            )

            total = tar + tir
            if total > 0:
                expected_af = tir / total
                assert abs(af_val - expected_af) < FLOAT_TOLERANCE, (
                    f"AF expected ≈{expected_af:.6f}, got {af_val:.6f}"
                )

            assert dp_val == expected_dp, (
                f"DP expected {expected_dp} (TAR+TIR), got {dp_val}"
            )
        vcf.close()
    finally:
        try:
            os.unlink(output_path)
        except OSError:
            pass


# ---------------------------------------------------------------------------
# SAGE correctness
# ---------------------------------------------------------------------------

@given(
    chrom=chrom_st,
    pos=pos_st,
    bases=snp_pairs_st,
    ref_depth=depth_st,
    alt_depth=depth_st,
    gq=gq_st,
)
@settings(max_examples=100, suppress_health_check=[HealthCheck.function_scoped_fixture])
def test_p3_sage_correctness(chrom, pos, bases, ref_depth, alt_depth, gq):
    """
    **Validates: Requirements 2.3**

    Property 3: FORMAT Harmonization Correctness — SAGE.
    Given AD (ref_depth, alt_depth) and VF value:
    - Output AF ≈ VF (the VF value is preserved as AF)
    """
    # Feature: neoantigen-workflow, Property 3: FORMAT Harmonization Correctness
    ref, alt = bases
    total = ref_depth + alt_depth
    expected_vf = round(alt_depth / total, 4) if total > 0 else 0.0

    vcf_content = build_sage_vcf(chrom, pos, ref, alt, ref_depth, alt_depth, gq)
    output_path = run_harmonizer(vcf_content, "sage")
    try:
        vcf = VCF(output_path)
        for variant in vcf:
            af = variant.format("AF")
            assert af is not None, "AF field missing from SAGE output"
            af_val = float(af[0][0])
            assert abs(af_val - expected_vf) < FLOAT_TOLERANCE, (
                f"SAGE AF expected ≈{expected_vf:.6f} (from VF), got {af_val:.6f}"
            )
        vcf.close()
    finally:
        try:
            os.unlink(output_path)
        except OSError:
            pass


# ---------------------------------------------------------------------------
# DeepSomatic correctness
# ---------------------------------------------------------------------------

@given(
    chrom=chrom_st,
    pos=pos_st,
    bases=snp_pairs_st,
    ref_depth=depth_st,
    alt_depth=depth_st,
    gq=gq_st,
)
@settings(max_examples=100, suppress_health_check=[HealthCheck.function_scoped_fixture])
def test_p3_deepsomatic_correctness(chrom, pos, bases, ref_depth, alt_depth, gq):
    """
    **Validates: Requirements 2.3**

    Property 3: FORMAT Harmonization Correctness — DeepSomatic.
    Given AD (ref_depth, alt_depth) and VAF value:
    - Output AF ≈ VAF (the VAF value is preserved as AF)
    """
    # Feature: neoantigen-workflow, Property 3: FORMAT Harmonization Correctness
    ref, alt = bases
    total = ref_depth + alt_depth
    expected_vaf = round(alt_depth / total, 4) if total > 0 else 0.0

    vcf_content = build_deepsomatic_vcf(chrom, pos, ref, alt, ref_depth, alt_depth, gq)
    output_path = run_harmonizer(vcf_content, "deepsomatic")
    try:
        vcf = VCF(output_path)
        for variant in vcf:
            af = variant.format("AF")
            assert af is not None, "AF field missing from DeepSomatic output"
            af_val = float(af[0][0])
            assert abs(af_val - expected_vaf) < FLOAT_TOLERANCE, (
                f"DeepSomatic AF expected ≈{expected_vaf:.6f} (from VAF), got {af_val:.6f}"
            )
        vcf.close()
    finally:
        try:
            os.unlink(output_path)
        except OSError:
            pass


# ===========================================================================
# Property 4: FORMAT Harmonization Idempotency
# Feature: neoantigen-workflow, Property 4: FORMAT Harmonization Idempotency
# Validates: Requirements 2.2, 2.3
# ===========================================================================


def _read_vcf_format_values(vcf_path: str) -> list[dict]:
    """
    Read all variant records from a VCF and return a list of dicts with
    AD, AF, DP, GQ values for each record (for comparison).
    """
    records = []
    vcf = VCF(vcf_path)
    for variant in vcf:
        ad_fmt = variant.format("AD")
        af_fmt = variant.format("AF")
        dp_fmt = variant.format("DP")
        gq_fmt = variant.format("GQ")

        ad_val = list(ad_fmt[0]) if ad_fmt is not None else None
        af_val = float(af_fmt[0][0]) if af_fmt is not None else None
        dp_val = int(dp_fmt[0][0]) if dp_fmt is not None else None
        gq_val = int(gq_fmt[0][0]) if gq_fmt is not None else None

        records.append({
            "chrom": variant.CHROM,
            "pos": variant.POS,
            "ref": variant.REF,
            "alt": variant.ALT,
            "AD": ad_val,
            "AF": af_val,
            "DP": dp_val,
            "GQ": gq_val,
        })
    vcf.close()
    return records


def _assert_records_identical(records1: list[dict], records2: list[dict], caller: str):
    """Assert that two lists of variant record dicts are identical."""
    assert len(records1) == len(records2), (
        f"[{caller}] Idempotency: record count differs: {len(records1)} vs {len(records2)}"
    )
    for i, (r1, r2) in enumerate(zip(records1, records2)):
        assert r1["chrom"] == r2["chrom"], f"[{caller}] Record {i}: CHROM differs"
        assert r1["pos"] == r2["pos"], f"[{caller}] Record {i}: POS differs"

        # AD comparison
        if r1["AD"] is not None and r2["AD"] is not None:
            assert list(r1["AD"]) == list(r2["AD"]), (
                f"[{caller}] Record {i}: AD differs: {r1['AD']} vs {r2['AD']}"
            )
        else:
            assert r1["AD"] == r2["AD"], (
                f"[{caller}] Record {i}: AD differs (None mismatch): {r1['AD']} vs {r2['AD']}"
            )

        # AF comparison (float tolerance)
        if r1["AF"] is not None and r2["AF"] is not None:
            assert abs(r1["AF"] - r2["AF"]) < FLOAT_TOLERANCE, (
                f"[{caller}] Record {i}: AF differs: {r1['AF']} vs {r2['AF']}"
            )
        else:
            assert r1["AF"] == r2["AF"], (
                f"[{caller}] Record {i}: AF differs (None mismatch): {r1['AF']} vs {r2['AF']}"
            )

        # DP comparison
        assert r1["DP"] == r2["DP"], (
            f"[{caller}] Record {i}: DP differs: {r1['DP']} vs {r2['DP']}"
        )

        # GQ comparison
        assert r1["GQ"] == r2["GQ"], (
            f"[{caller}] Record {i}: GQ differs: {r1['GQ']} vs {r2['GQ']}"
        )


@given(
    chrom=chrom_st,
    pos=pos_st,
    bases=snp_pairs_st,
    ref_depth=depth_st,
    alt_depth=depth_st,
    gq=gq_st,
)
@settings(max_examples=100, suppress_health_check=[HealthCheck.function_scoped_fixture])
def test_p4_mutect2_idempotency(chrom, pos, bases, ref_depth, alt_depth, gq):
    """
    **Validates: Requirements 2.2, 2.3**

    Property 4: FORMAT Harmonization Idempotency — mutect2 caller.
    Running harmonize() twice on a mutect2 VCF must produce identical
    AD/AF/DP/GQ values in the second output as in the first output.

    # Feature: neoantigen-workflow, Property 4: FORMAT Harmonization Idempotency
    """
    ref, alt = bases
    vcf_content = build_mutect2_vcf(chrom, pos, ref, alt, ref_depth, alt_depth, gq)

    output1 = run_harmonizer(vcf_content, "mutect2")
    output2 = output1.replace(".harmonized.vcf", ".harmonized2.vcf")
    try:
        harmonize(output1, "mutect2", output2)

        records1 = _read_vcf_format_values(output1)
        records2 = _read_vcf_format_values(output2)
        _assert_records_identical(records1, records2, "mutect2")
    finally:
        for path in [output1, output2]:
            if path:
                try:
                    os.unlink(path)
                except OSError:
                    pass


@given(
    chrom=chrom_st,
    pos=pos_st,
    bases=snp_pairs_st,
    au=depth_st,
    cu=depth_st,
    gu=depth_st,
    tu=depth_st,
    gq=gq_st,
)
@settings(max_examples=100, suppress_health_check=[HealthCheck.function_scoped_fixture])
def test_p4_strelka_snv_idempotency(chrom, pos, bases, au, cu, gu, tu, gq):
    """
    **Validates: Requirements 2.2, 2.3**

    Property 4: FORMAT Harmonization Idempotency — strelka_snv caller.
    Running harmonize() twice on a strelka SNV VCF must produce identical
    AD/AF/DP/GQ values in the second output as in the first output.

    # Feature: neoantigen-workflow, Property 4: FORMAT Harmonization Idempotency
    """
    ref, alt = bases
    vcf_content = build_strelka_snv_vcf(chrom, pos, ref, alt, au, cu, gu, tu, gq)

    output1 = run_harmonizer(vcf_content, "strelka_snv")
    output2 = output1.replace(".harmonized.vcf", ".harmonized2.vcf")
    try:
        harmonize(output1, "strelka_snv", output2)

        records1 = _read_vcf_format_values(output1)
        records2 = _read_vcf_format_values(output2)
        _assert_records_identical(records1, records2, "strelka_snv")
    finally:
        for path in [output1, output2]:
            try:
                os.unlink(path)
            except OSError:
                pass


@given(
    chrom=chrom_st,
    pos=pos_st,
    indel=indel_pairs_st,
    tar=depth_st,
    tir=depth_st,
    gq=gq_st,
)
@settings(max_examples=100, suppress_health_check=[HealthCheck.function_scoped_fixture])
def test_p4_strelka_indel_idempotency(chrom, pos, indel, tar, tir, gq):
    """
    **Validates: Requirements 2.2, 2.3**

    Property 4: FORMAT Harmonization Idempotency — strelka_indel caller.
    Running harmonize() twice on a strelka Indel VCF must produce identical
    AD/AF/DP/GQ values in the second output as in the first output.

    # Feature: neoantigen-workflow, Property 4: FORMAT Harmonization Idempotency
    """
    ref, alt = indel
    vcf_content = build_strelka_indel_vcf(chrom, pos, ref, alt, tar, tir, gq)

    output1 = run_harmonizer(vcf_content, "strelka_indel")
    output2 = output1.replace(".harmonized.vcf", ".harmonized2.vcf")
    try:
        harmonize(output1, "strelka_indel", output2)

        records1 = _read_vcf_format_values(output1)
        records2 = _read_vcf_format_values(output2)
        _assert_records_identical(records1, records2, "strelka_indel")
    finally:
        for path in [output1, output2]:
            try:
                os.unlink(path)
            except OSError:
                pass


@given(
    chrom=chrom_st,
    pos=pos_st,
    bases=snp_pairs_st,
    ref_depth=depth_st,
    alt_depth=depth_st,
    gq=gq_st,
)
@settings(max_examples=100, suppress_health_check=[HealthCheck.function_scoped_fixture])
def test_p4_sage_idempotency(chrom, pos, bases, ref_depth, alt_depth, gq):
    """
    **Validates: Requirements 2.2, 2.3**

    Property 4: FORMAT Harmonization Idempotency — sage caller.
    Running harmonize() twice on a SAGE VCF must produce identical
    AD/AF/DP/GQ values in the second output as in the first output.

    # Feature: neoantigen-workflow, Property 4: FORMAT Harmonization Idempotency
    """
    ref, alt = bases
    vcf_content = build_sage_vcf(chrom, pos, ref, alt, ref_depth, alt_depth, gq)

    output1 = run_harmonizer(vcf_content, "sage")
    output2 = output1.replace(".harmonized.vcf", ".harmonized2.vcf")
    try:
        harmonize(output1, "sage", output2)

        records1 = _read_vcf_format_values(output1)
        records2 = _read_vcf_format_values(output2)
        _assert_records_identical(records1, records2, "sage")
    finally:
        for path in [output1, output2]:
            try:
                os.unlink(path)
            except OSError:
                pass


@given(
    chrom=chrom_st,
    pos=pos_st,
    bases=snp_pairs_st,
    ref_depth=depth_st,
    alt_depth=depth_st,
    gq=gq_st,
)
@settings(max_examples=100, suppress_health_check=[HealthCheck.function_scoped_fixture])
def test_p4_deepsomatic_idempotency(chrom, pos, bases, ref_depth, alt_depth, gq):
    """
    **Validates: Requirements 2.2, 2.3**

    Property 4: FORMAT Harmonization Idempotency — deepsomatic caller.
    Running harmonize() twice on a DeepSomatic VCF must produce identical
    AD/AF/DP/GQ values in the second output as in the first output.

    # Feature: neoantigen-workflow, Property 4: FORMAT Harmonization Idempotency
    """
    ref, alt = bases
    vcf_content = build_deepsomatic_vcf(chrom, pos, ref, alt, ref_depth, alt_depth, gq)

    output1 = run_harmonizer(vcf_content, "deepsomatic")
    output2 = output1.replace(".harmonized.vcf", ".harmonized2.vcf")
    try:
        harmonize(output1, "deepsomatic", output2)

        records1 = _read_vcf_format_values(output1)
        records2 = _read_vcf_format_values(output2)
        _assert_records_identical(records1, records2, "deepsomatic")
    finally:
        for path in [output1, output2]:
            try:
                os.unlink(path)
            except OSError:
                pass


# ===========================================================================
# Property 5: FORMAT Field Preservation
# Feature: neoantigen-workflow, Property 5: FORMAT Field Preservation
# Validates: Requirements 2.4
# ===========================================================================


@given(
    chrom=chrom_st,
    pos=pos_st,
    bases=snp_pairs_st,
    ref_depth=depth_st,
    alt_depth=depth_st,
    gq=gq_st,
)
@settings(max_examples=100, suppress_health_check=[HealthCheck.function_scoped_fixture])
def test_p5_strelka_snv_non_target_fields_preserved(chrom, pos, bases, ref_depth, alt_depth, gq):
    """
    **Validates: Requirements 2.4**

    Property 5: FORMAT Field Preservation — strelka SNV caller.
    For any strelka SNV VCF, the non-target FORMAT fields AU, CU, GU, TU
    must be present and have the same tier1 values in the harmonized output
    as in the input VCF.

    # Feature: neoantigen-workflow, Property 5: FORMAT Field Preservation
    """
    ref, alt = bases
    # Use ref_depth/alt_depth as AU/CU/GU/TU values for simplicity
    au, cu, gu, tu = ref_depth, alt_depth, gq, (ref_depth + alt_depth) % 501
    vcf_content = build_strelka_snv_vcf(chrom, pos, ref, alt, au, cu, gu, tu, gq)
    output_path = run_harmonizer(vcf_content, "strelka_snv")
    try:
        vcf = VCF(output_path)
        # Verify non-target fields are declared in header
        declared_ids = {h["ID"] for h in vcf.header_iter() if h["HeaderType"] == "FORMAT"}
        for field in ("AU", "CU", "GU", "TU"):
            assert field in declared_ids, (
                f"[strelka_snv] Non-target field {field} missing from FORMAT header. "
                f"Declared: {declared_ids}"
            )
        # Verify tier1 values are preserved in each record
        expected = {"AU": au, "CU": cu, "GU": gu, "TU": tu}
        for variant in vcf:
            for field, expected_val in expected.items():
                fmt = variant.format(field)
                assert fmt is not None, (
                    f"[strelka_snv] Non-target field {field} missing from record "
                    f"{variant.CHROM}:{variant.POS}"
                )
                # Number=2 fields: stored as (tier1, tier2); verify tier1
                tier1 = int(fmt[0][0])
                assert tier1 == expected_val, (
                    f"[strelka_snv] Field {field} tier1 expected {expected_val}, got {tier1} "
                    f"at {variant.CHROM}:{variant.POS}"
                )
        vcf.close()
    finally:
        try:
            os.unlink(output_path)
        except OSError:
            pass


@given(
    chrom=chrom_st,
    pos=pos_st,
    indel=indel_pairs_st,
    tar=depth_st,
    tir=depth_st,
    gq=gq_st,
)
@settings(max_examples=100, suppress_health_check=[HealthCheck.function_scoped_fixture])
def test_p5_strelka_indel_non_target_fields_preserved(chrom, pos, indel, tar, tir, gq):
    """
    **Validates: Requirements 2.4**

    Property 5: FORMAT Field Preservation — strelka indel caller.
    For any strelka Indel VCF, the non-target FORMAT fields TAR and TIR
    must be present and have the same tier1 values in the harmonized output
    as in the input VCF.

    # Feature: neoantigen-workflow, Property 5: FORMAT Field Preservation
    """
    ref, alt = indel
    vcf_content = build_strelka_indel_vcf(chrom, pos, ref, alt, tar, tir, gq)
    output_path = run_harmonizer(vcf_content, "strelka_indel")
    try:
        vcf = VCF(output_path)
        # Verify non-target fields are declared in header
        declared_ids = {h["ID"] for h in vcf.header_iter() if h["HeaderType"] == "FORMAT"}
        for field in ("TAR", "TIR"):
            assert field in declared_ids, (
                f"[strelka_indel] Non-target field {field} missing from FORMAT header. "
                f"Declared: {declared_ids}"
            )
        # Verify tier1 values are preserved in each record
        expected = {"TAR": tar, "TIR": tir}
        for variant in vcf:
            for field, expected_val in expected.items():
                fmt = variant.format(field)
                assert fmt is not None, (
                    f"[strelka_indel] Non-target field {field} missing from record "
                    f"{variant.CHROM}:{variant.POS}"
                )
                # Number=2 fields: stored as (tier1, tier2); verify tier1
                tier1 = int(fmt[0][0])
                assert tier1 == expected_val, (
                    f"[strelka_indel] Field {field} tier1 expected {expected_val}, got {tier1} "
                    f"at {variant.CHROM}:{variant.POS}"
                )
        vcf.close()
    finally:
        try:
            os.unlink(output_path)
        except OSError:
            pass


@given(
    chrom=chrom_st,
    pos=pos_st,
    bases=snp_pairs_st,
    ref_depth=depth_st,
    alt_depth=depth_st,
    gq=gq_st,
)
@settings(max_examples=100, suppress_health_check=[HealthCheck.function_scoped_fixture])
def test_p5_sage_non_target_fields_preserved(chrom, pos, bases, ref_depth, alt_depth, gq):
    """
    **Validates: Requirements 2.4**

    Property 5: FORMAT Field Preservation — SAGE caller.
    For any SAGE VCF, the non-target FORMAT field VF must be present and
    have the same value in the harmonized output as in the input VCF.

    # Feature: neoantigen-workflow, Property 5: FORMAT Field Preservation
    """
    ref, alt = bases
    total = ref_depth + alt_depth
    expected_vf = round(alt_depth / total, 4) if total > 0 else 0.0

    vcf_content = build_sage_vcf(chrom, pos, ref, alt, ref_depth, alt_depth, gq)
    output_path = run_harmonizer(vcf_content, "sage")
    try:
        vcf = VCF(output_path)
        # Verify VF is declared in header
        declared_ids = {h["ID"] for h in vcf.header_iter() if h["HeaderType"] == "FORMAT"}
        assert "VF" in declared_ids, (
            f"[sage] Non-target field VF missing from FORMAT header. "
            f"Declared: {declared_ids}"
        )
        # Verify VF value is preserved in each record
        for variant in vcf:
            vf_fmt = variant.format("VF")
            assert vf_fmt is not None, (
                f"[sage] Non-target field VF missing from record "
                f"{variant.CHROM}:{variant.POS}"
            )
            vf_val = float(vf_fmt[0][0])
            assert abs(vf_val - expected_vf) < 1e-4, (
                f"[sage] VF expected {expected_vf:.6f}, got {vf_val:.6f} "
                f"at {variant.CHROM}:{variant.POS}"
            )
        vcf.close()
    finally:
        try:
            os.unlink(output_path)
        except OSError:
            pass


@given(
    chrom=chrom_st,
    pos=pos_st,
    bases=snp_pairs_st,
    ref_depth=depth_st,
    alt_depth=depth_st,
    gq=gq_st,
)
@settings(max_examples=100, suppress_health_check=[HealthCheck.function_scoped_fixture])
def test_p5_deepsomatic_non_target_fields_preserved(chrom, pos, bases, ref_depth, alt_depth, gq):
    """
    **Validates: Requirements 2.4**

    Property 5: FORMAT Field Preservation — DeepSomatic caller.
    For any DeepSomatic VCF, the non-target FORMAT field VAF must be present
    and have the same value in the harmonized output as in the input VCF.

    # Feature: neoantigen-workflow, Property 5: FORMAT Field Preservation
    """
    ref, alt = bases
    total = ref_depth + alt_depth
    expected_vaf = round(alt_depth / total, 4) if total > 0 else 0.0

    vcf_content = build_deepsomatic_vcf(chrom, pos, ref, alt, ref_depth, alt_depth, gq)
    output_path = run_harmonizer(vcf_content, "deepsomatic")
    try:
        vcf = VCF(output_path)
        # Verify VAF is declared in header
        declared_ids = {h["ID"] for h in vcf.header_iter() if h["HeaderType"] == "FORMAT"}
        assert "VAF" in declared_ids, (
            f"[deepsomatic] Non-target field VAF missing from FORMAT header. "
            f"Declared: {declared_ids}"
        )
        # Verify VAF value is preserved in each record
        for variant in vcf:
            vaf_fmt = variant.format("VAF")
            assert vaf_fmt is not None, (
                f"[deepsomatic] Non-target field VAF missing from record "
                f"{variant.CHROM}:{variant.POS}"
            )
            vaf_val = float(vaf_fmt[0][0])
            assert abs(vaf_val - expected_vaf) < 1e-4, (
                f"[deepsomatic] VAF expected {expected_vaf:.6f}, got {vaf_val:.6f} "
                f"at {variant.CHROM}:{variant.POS}"
            )
        vcf.close()
    finally:
        try:
            os.unlink(output_path)
        except OSError:
            pass


@given(
    chrom=chrom_st,
    pos=pos_st,
    bases=snp_pairs_st,
    ref_depth=depth_st,
    alt_depth=depth_st,
    gq=gq_st,
)
@settings(max_examples=100, suppress_health_check=[HealthCheck.function_scoped_fixture])
def test_p5_mutect2_non_target_fields_preserved(chrom, pos, bases, ref_depth, alt_depth, gq):
    """
    **Validates: Requirements 2.4**

    Property 5: FORMAT Field Preservation — Mutect2 caller.
    For any Mutect2 VCF, the non-target FORMAT fields F1R2 and F2R1 must be
    present and have the same values in the harmonized output as in the input VCF.

    # Feature: neoantigen-workflow, Property 5: FORMAT Field Preservation
    """
    ref, alt = bases
    # F1R2 = [ref_depth, 0], F2R1 = [0, alt_depth] as set in build_mutect2_vcf
    expected_f1r2 = [ref_depth, 0]
    expected_f2r1 = [0, alt_depth]

    vcf_content = build_mutect2_vcf(chrom, pos, ref, alt, ref_depth, alt_depth, gq)
    output_path = run_harmonizer(vcf_content, "mutect2")
    try:
        vcf = VCF(output_path)
        # Verify F1R2 and F2R1 are declared in header
        declared_ids = {h["ID"] for h in vcf.header_iter() if h["HeaderType"] == "FORMAT"}
        for field in ("F1R2", "F2R1"):
            assert field in declared_ids, (
                f"[mutect2] Non-target field {field} missing from FORMAT header. "
                f"Declared: {declared_ids}"
            )
        # Verify values are preserved in each record
        for variant in vcf:
            f1r2_fmt = variant.format("F1R2")
            f2r1_fmt = variant.format("F2R1")
            assert f1r2_fmt is not None, (
                f"[mutect2] Non-target field F1R2 missing from record "
                f"{variant.CHROM}:{variant.POS}"
            )
            assert f2r1_fmt is not None, (
                f"[mutect2] Non-target field F2R1 missing from record "
                f"{variant.CHROM}:{variant.POS}"
            )
            f1r2_vals = list(f1r2_fmt[0])
            f2r1_vals = list(f2r1_fmt[0])
            assert f1r2_vals == expected_f1r2, (
                f"[mutect2] F1R2 expected {expected_f1r2}, got {f1r2_vals} "
                f"at {variant.CHROM}:{variant.POS}"
            )
            assert f2r1_vals == expected_f2r1, (
                f"[mutect2] F2R1 expected {expected_f2r1}, got {f2r1_vals} "
                f"at {variant.CHROM}:{variant.POS}"
            )
        vcf.close()
    finally:
        try:
            os.unlink(output_path)
        except OSError:
            pass
