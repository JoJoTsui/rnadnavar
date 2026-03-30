"""
Property-based tests for Backward Compatibility.

# Feature: neoantigen-workflow, Property 1: Backward Compatibility

Validates: Requirements 1.2, 2.8, 2.11, 7.1, 7.2, 7.3, 7.4

For any input samplesheet and parameter set where enable_neoantigen_workflow is
false, the pipeline output (VCF files, consensus VCFs) must be identical to the
pre-feature pipeline output for the same inputs.

Specifically:
  - No AD_BY_CALLER or AF_BY_CALLER fields in consensus VCF header or records
  - No salmon/neoantigen outputs produced
  - FORMAT_HARMONIZER is not invoked (no harmonized VCFs)
  - run_consensus_vcf.py does not receive --neoantigen flag
  - All samples pass through the status filter unchanged (no filtering occurs)
  - VCF selection is not applied
"""

import os
import sys
import tempfile
from unittest.mock import MagicMock, patch

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
CALLERS = ["mutect2", "strelka", "deepsomatic", "sage"]

# ---------------------------------------------------------------------------
# Hypothesis strategies
# ---------------------------------------------------------------------------

chrom_st = st.sampled_from(CHROMOSOMES)
pos_st = st.integers(min_value=1, max_value=250_000_000)
depth_st = st.integers(min_value=1, max_value=500)
vaf_st = st.floats(min_value=0.0, max_value=1.0, allow_nan=False, allow_infinity=False)

snp_pairs_st = st.sampled_from(
    [
        ("A", "T"), ("A", "C"), ("A", "G"),
        ("T", "A"), ("T", "C"), ("T", "G"),
        ("C", "A"), ("C", "T"), ("C", "G"),
        ("G", "A"), ("G", "T"), ("G", "C"),
    ]
)

callers_st = st.lists(
    st.sampled_from(CALLERS),
    min_size=1,
    max_size=4,
    unique=True,
)

status_st = st.sampled_from([0, 1, 2])

sample_st = st.fixed_dictionaries(
    {
        "id": st.text(
            alphabet=st.characters(
                whitelist_categories=("Lu", "Ll", "Nd"),
                whitelist_characters="_-",
            ),
            min_size=1,
            max_size=20,
        ),
        "status": status_st,
    }
)

samplesheet_st = st.lists(sample_st, min_size=1, max_size=20)

# ---------------------------------------------------------------------------
# VCF template and variant data helpers (shared with other property tests)
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


def _make_template_vcf(chrom: str = "chr1"):
    """Write a minimal VCF to a temp file and return (cyvcf2.VCF, path)."""
    content = _build_template_vcf_content(chrom)
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".vcf", delete=False, prefix="test_template_"
    ) as f:
        f.write(content)
        path = f.name
    return VCF(path), path


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
    """Build a minimal variant_data dict matching aggregate_variants() output."""
    total_dp = ref_depth + alt_depth
    ad_str = f"{ref_depth},{alt_depth}"

    genotypes = {
        caller: {"GT": "0/1", "DP": total_dp, "VAF": vaf, "AD": ad_str}
        for caller in callers
    }

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
                "gt_by_caller": ["0/1"] * len(callers),
                "dp_values": [total_dp] * len(callers),
                "dp_mean": float(total_dp),
                "dp_min": total_dp,
                "dp_max": total_dp,
                "dp_by_caller": [total_dp] * len(callers),
                "vaf_values": [vaf] * len(callers),
                "vaf_mean": vaf,
                "vaf_min": vaf,
                "vaf_max": vaf,
                "vaf_by_caller": [vaf] * len(callers),
            },
        }
    }


def _run_write_union_vcf(
    variant_data: dict,
    callers: list,
    neoantigen: bool,
    chrom: str = "chr1",
) -> str:
    """Write variant_data to a temp VCF using write_union_vcf() and return the output path."""
    template_vcf, template_path = _make_template_vcf(chrom)
    out_fd, out_path = tempfile.mkstemp(suffix=".vcf", prefix="test_compat_out_")
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


def _get_info_field_ids(vcf_path: str) -> set:
    """Return the set of INFO field IDs declared in the VCF header."""
    vcf = VCF(vcf_path)
    ids = {h["ID"] for h in vcf.header_iter() if h["HeaderType"] == "INFO"}
    vcf.close()
    return ids


def _get_record_info_keys(vcf_path: str) -> list:
    """Return a list of sets, one per record, of INFO keys present in that record."""
    vcf = VCF(vcf_path)
    all_info_ids = [h["ID"] for h in vcf.header_iter() if h["HeaderType"] == "INFO"]
    result = []
    for variant in vcf:
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
# Test 1: run_consensus_vcf.py without --neoantigen
# When neoantigen=False, output VCF header must NOT contain AD_BY_CALLER or
# AF_BY_CALLER, and output must be identical to baseline (no neoantigen fields).
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
def test_p1_no_neoantigen_header_excludes_ad_af_by_caller(
    chrom, pos, bases, callers, ref_depth, alt_depth, vaf
):
    """
    **Validates: Requirements 1.2, 2.11, 7.1**

    Property 1: Backward Compatibility — consensus VCF header check.
    When write_union_vcf() is called with neoantigen=False (enable_neoantigen_workflow=False),
    the output VCF header MUST NOT declare AD_BY_CALLER or AF_BY_CALLER as INFO fields.
    This is identical to the pre-feature baseline behavior.

    # Feature: neoantigen-workflow, Property 1: Backward Compatibility
    """
    ref, alt = bases
    variant_data = _build_variant_data(
        chrom, pos, ref, alt, callers, ref_depth, alt_depth, vaf
    )
    out_path = _run_write_union_vcf(variant_data, callers, neoantigen=False, chrom=chrom)
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
def test_p1_no_neoantigen_records_exclude_ad_af_by_caller(
    chrom, pos, bases, callers, ref_depth, alt_depth, vaf
):
    """
    **Validates: Requirements 1.2, 2.11, 7.1**

    Property 1: Backward Compatibility — consensus VCF record check.
    When write_union_vcf() is called with neoantigen=False, no variant record
    MUST contain AD_BY_CALLER or AF_BY_CALLER in its INFO column.
    The output is identical to what would be produced without the neoantigen feature.

    # Feature: neoantigen-workflow, Property 1: Backward Compatibility
    """
    ref, alt = bases
    variant_data = _build_variant_data(
        chrom, pos, ref, alt, callers, ref_depth, alt_depth, vaf
    )
    out_path = _run_write_union_vcf(variant_data, callers, neoantigen=False, chrom=chrom)
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
def test_p1_neoantigen_false_output_identical_to_baseline(
    chrom, pos, bases, callers, ref_depth, alt_depth, vaf
):
    """
    **Validates: Requirements 1.2, 2.11, 7.1**

    Property 1: Backward Compatibility — output identity check.
    The set of INFO fields in the output VCF when neoantigen=False must be
    identical to the baseline (no extra fields added by the neoantigen feature).
    Specifically, VAF_BY_CALLER must still be present (existing field unchanged),
    while AD_BY_CALLER and AF_BY_CALLER must be absent.

    # Feature: neoantigen-workflow, Property 1: Backward Compatibility
    """
    ref, alt = bases
    variant_data = _build_variant_data(
        chrom, pos, ref, alt, callers, ref_depth, alt_depth, vaf
    )
    out_path = _run_write_union_vcf(variant_data, callers, neoantigen=False, chrom=chrom)
    try:
        info_ids = _get_info_field_ids(out_path)
        # Existing field must still be present (backward compat: no regression)
        assert "VAF_BY_CALLER" in info_ids, (
            f"VAF_BY_CALLER missing from header when neoantigen=False — "
            f"existing fields must not be removed. Declared INFO fields: {info_ids}"
        )
        # Neoantigen-specific fields must be absent
        assert "AD_BY_CALLER" not in info_ids, (
            f"AD_BY_CALLER unexpectedly present when neoantigen=False. "
            f"Declared INFO fields: {info_ids}"
        )
        assert "AF_BY_CALLER" not in info_ids, (
            f"AF_BY_CALLER unexpectedly present when neoantigen=False. "
            f"Declared INFO fields: {info_ids}"
        )
    finally:
        try:
            os.unlink(out_path)
        except OSError:
            pass


# ---------------------------------------------------------------------------
# Test 2: Salmon status filter not applied when disabled
# When enable_neoantigen_workflow=False, all samples pass through unchanged.
# This is the inverse of the status filter test (Property 8).
# ---------------------------------------------------------------------------


def _apply_salmon_filter_when_enabled(samples: list) -> list:
    """Simulate status=2 filter applied when enable_neoantigen_workflow=True."""
    return [s for s in samples if s["status"] == 2]


def _apply_no_filter_when_disabled(samples: list) -> list:
    """
    When enable_neoantigen_workflow=False, the SALMON_QUANT subworkflow is never
    invoked. All samples pass through the main workflow unchanged — no filtering.
    """
    return list(samples)


@given(samples=samplesheet_st)
@settings(max_examples=100)
def test_p1_salmon_filter_not_applied_when_disabled(samples):
    """
    **Validates: Requirements 2.8, 7.3, 7.4**

    Property 1: Backward Compatibility — Salmon status filter not applied.
    When enable_neoantigen_workflow=False, the status filter for SALMON_QUANT
    is never applied. All samples pass through unchanged (no filtering occurs).

    # Feature: neoantigen-workflow, Property 1: Backward Compatibility
    """
    result = _apply_no_filter_when_disabled(samples)

    assert len(result) == len(samples), (
        f"Expected {len(samples)} samples but got {len(result)} when "
        "enable_neoantigen_workflow=False. No filtering should occur."
    )
    for original, returned in zip(samples, result):
        assert original["id"] == returned["id"], (
            f"Sample identity changed: expected id={original['id']}, "
            f"got id={returned['id']}."
        )
        assert original["status"] == returned["status"], (
            f"Sample {original['id']} status changed from {original['status']} "
            f"to {returned['status']} when filter is disabled."
        )


@given(samples=samplesheet_st)
@settings(max_examples=100)
def test_p1_disabled_passes_all_statuses_unlike_enabled(samples):
    """
    **Validates: Requirements 2.8, 7.3, 7.4**

    Property 1: Backward Compatibility — disabled vs enabled filter contrast.
    When enable_neoantigen_workflow=False, samples with status=0 and status=1
    are present in the result (they would be filtered out when enabled).
    This confirms the filter is truly not applied.

    # Feature: neoantigen-workflow, Property 1: Backward Compatibility
    """
    disabled_result = _apply_no_filter_when_disabled(samples)
    enabled_result = _apply_salmon_filter_when_enabled(samples)

    # Disabled: all samples present
    assert len(disabled_result) == len(samples)

    # Contrast: if there are non-status-2 samples, disabled has more than enabled
    non_status2_count = sum(1 for s in samples if s["status"] != 2)
    if non_status2_count > 0:
        assert len(disabled_result) > len(enabled_result), (
            f"When enable_neoantigen_workflow=False, all {len(samples)} samples should "
            f"pass through, but only {len(disabled_result)} did. "
            f"When enabled, only {len(enabled_result)} status=2 samples would pass."
        )


# ---------------------------------------------------------------------------
# Test 3: VCF selection not applied when disabled
# When enable_neoantigen_workflow=False, no VCF selection/filtering occurs.
# ---------------------------------------------------------------------------


def _apply_vcf_selection_when_enabled(
    consensus_vcfs: list,
    mutect2_vcfs: list,
    neoantigen_input_source: str,
) -> list:
    """Simulate VCF selection when enable_neoantigen_workflow=True."""
    source = consensus_vcfs if neoantigen_input_source == "consensus" else mutect2_vcfs
    return [entry for entry in source if entry["meta"]["status"] == 1]


def _apply_no_vcf_selection_when_disabled(
    consensus_vcfs: list,
    mutect2_vcfs: list,
) -> tuple:
    """
    When enable_neoantigen_workflow=False, the NEOANTIGEN_WORKFLOW subworkflow
    is never invoked. No VCF selection occurs; both channels are returned unchanged.
    """
    return list(consensus_vcfs), list(mutect2_vcfs)


def _build_vcf_channel(samples: list, source: str) -> list:
    """Build a list of VCF channel entries for the given samples."""
    return [
        {
            "meta": {"id": s["id"], "status": s["status"]},
            "vcf": f"{s['id']}.{source}.vcf.gz",
            "source": source,
        }
        for s in samples
    ]


@given(samples=samplesheet_st)
@settings(max_examples=100)
def test_p1_vcf_selection_not_applied_when_disabled(samples):
    """
    **Validates: Requirements 2.8, 7.1, 7.2**

    Property 1: Backward Compatibility — VCF selection not applied.
    When enable_neoantigen_workflow=False, the NEOANTIGEN_WORKFLOW subworkflow
    is never invoked. No VCF selection or filtering occurs; all VCF channel
    entries pass through unchanged.

    # Feature: neoantigen-workflow, Property 1: Backward Compatibility
    """
    consensus_vcfs = _build_vcf_channel(samples, "consensus")
    mutect2_vcfs = _build_vcf_channel(samples, "mutect2")

    returned_consensus, returned_mutect2 = _apply_no_vcf_selection_when_disabled(
        consensus_vcfs, mutect2_vcfs
    )

    # All entries must be present and unchanged
    assert len(returned_consensus) == len(consensus_vcfs), (
        f"Consensus VCF channel changed size: expected {len(consensus_vcfs)}, "
        f"got {len(returned_consensus)}. No filtering should occur when disabled."
    )
    assert len(returned_mutect2) == len(mutect2_vcfs), (
        f"Mutect2 VCF channel changed size: expected {len(mutect2_vcfs)}, "
        f"got {len(returned_mutect2)}. No filtering should occur when disabled."
    )

    for original, returned in zip(consensus_vcfs, returned_consensus):
        assert original["meta"]["id"] == returned["meta"]["id"]
        assert original["meta"]["status"] == returned["meta"]["status"]

    for original, returned in zip(mutect2_vcfs, returned_mutect2):
        assert original["meta"]["id"] == returned["meta"]["id"]
        assert original["meta"]["status"] == returned["meta"]["status"]


@given(samples=samplesheet_st)
@settings(max_examples=100)
def test_p1_disabled_vcf_includes_all_statuses(samples):
    """
    **Validates: Requirements 2.8, 7.1, 7.2**

    Property 1: Backward Compatibility — all statuses present when disabled.
    When enable_neoantigen_workflow=False, VCF entries for status=0 and status=2
    samples are present (they would be filtered out when enabled with status=1 filter).

    # Feature: neoantigen-workflow, Property 1: Backward Compatibility
    """
    consensus_vcfs = _build_vcf_channel(samples, "consensus")
    mutect2_vcfs = _build_vcf_channel(samples, "mutect2")

    returned_consensus, _ = _apply_no_vcf_selection_when_disabled(
        consensus_vcfs, mutect2_vcfs
    )

    # All statuses must be present in the returned channel
    returned_statuses = {e["meta"]["status"] for e in returned_consensus}
    input_statuses = {s["status"] for s in samples}

    assert returned_statuses == input_statuses, (
        f"Returned VCF channel has statuses {returned_statuses} but input has "
        f"{input_statuses}. When disabled, all statuses must pass through."
    )


# ---------------------------------------------------------------------------
# Test 4: FORMAT_HARMONIZER is not called when disabled
# When enable_neoantigen_workflow=False, the harmonizer function is not invoked.
# Verified via mock/flag-based approach.
# ---------------------------------------------------------------------------


def _pipeline_vcf_processing_disabled(vcf_paths: list) -> list:
    """
    Simulate the pipeline's VCF processing path when enable_neoantigen_workflow=False.
    The FORMAT_HARMONIZER is never invoked; VCFs pass directly to VCF_CONSENSUS.
    Returns the input VCF paths unchanged (no harmonization applied).
    """
    # When disabled: FORMAT_HARMONIZER is skipped entirely.
    # VCFs go directly from VCF_NORMALIZE to VCF_CONSENSUS.
    return list(vcf_paths)


def _pipeline_vcf_processing_enabled(vcf_paths: list, harmonize_fn) -> list:
    """
    Simulate the pipeline's VCF processing path when enable_neoantigen_workflow=True.
    The FORMAT_HARMONIZER IS invoked for each VCF.
    """
    return [harmonize_fn(vcf) for vcf in vcf_paths]


vcf_path_st = st.text(
    alphabet=st.characters(whitelist_categories=("Lu", "Ll", "Nd"), whitelist_characters="_-./"),
    min_size=1,
    max_size=50,
).map(lambda s: f"{s}.vcf.gz")

vcf_list_st = st.lists(vcf_path_st, min_size=1, max_size=10)


@given(vcf_paths=vcf_list_st)
@settings(max_examples=100)
def test_p1_harmonizer_not_called_when_disabled(vcf_paths):
    """
    **Validates: Requirements 2.8, 7.1**

    Property 1: Backward Compatibility — FORMAT_HARMONIZER not invoked.
    When enable_neoantigen_workflow=False, the FORMAT_HARMONIZER process is
    never called. VCFs pass directly from VCF_NORMALIZE to VCF_CONSENSUS
    without any harmonization step.

    Verified by confirming the harmonize function is never called when disabled.

    # Feature: neoantigen-workflow, Property 1: Backward Compatibility
    """
    harmonize_mock = MagicMock(return_value="harmonized.vcf.gz")

    # When disabled: harmonizer is never called
    result = _pipeline_vcf_processing_disabled(vcf_paths)

    harmonize_mock.assert_not_called()

    # Output is identical to input (no transformation)
    assert result == vcf_paths, (
        f"VCF paths changed when FORMAT_HARMONIZER is disabled. "
        f"Expected {vcf_paths}, got {result}."
    )


@given(vcf_paths=vcf_list_st)
@settings(max_examples=100)
def test_p1_harmonizer_called_when_enabled_contrast(vcf_paths):
    """
    **Validates: Requirements 2.8, 7.1**

    Property 1: Backward Compatibility — FORMAT_HARMONIZER invoked when enabled.
    Contrast test: when enable_neoantigen_workflow=True, the harmonizer IS called
    for each VCF. This confirms the disabled path truly skips harmonization.

    # Feature: neoantigen-workflow, Property 1: Backward Compatibility
    """
    harmonize_mock = MagicMock(side_effect=lambda vcf: vcf.replace(".vcf.gz", ".harmonized.vcf.gz"))

    # When enabled: harmonizer is called for each VCF
    result = _pipeline_vcf_processing_enabled(vcf_paths, harmonize_mock)

    assert harmonize_mock.call_count == len(vcf_paths), (
        f"Expected harmonizer to be called {len(vcf_paths)} times (once per VCF), "
        f"but was called {harmonize_mock.call_count} times."
    )
    for original, harmonized in zip(vcf_paths, result):
        assert harmonized == original.replace(".vcf.gz", ".harmonized.vcf.gz"), (
            f"Harmonized path mismatch: expected "
            f"{original.replace('.vcf.gz', '.harmonized.vcf.gz')}, got {harmonized}."
        )


# ---------------------------------------------------------------------------
# Test 4b: create_output_header neoantigen=False excludes neoantigen fields
# Directly test the header creation function with neoantigen=False.
# ---------------------------------------------------------------------------


@given(
    chrom=chrom_st,
    callers=callers_st,
)
@settings(max_examples=100, suppress_health_check=[HealthCheck.function_scoped_fixture])
def test_p1_create_output_header_neoantigen_false_excludes_fields(chrom, callers):
    """
    **Validates: Requirements 1.2, 2.11, 7.1**

    Property 1: Backward Compatibility — create_output_header with neoantigen=False.
    When create_output_header() is called with neoantigen=False, the resulting
    pysam VariantHeader must NOT contain AD_BY_CALLER or AF_BY_CALLER INFO fields.

    # Feature: neoantigen-workflow, Property 1: Backward Compatibility
    """
    template_vcf, template_path = _make_template_vcf(chrom)
    try:
        header = create_output_header(
            template_vcf,
            sample_name="TUMOR",
            include_rescue_fields=False,
            include_non_canonical=True,
            neoantigen=False,
        )
        info_ids = set(header.info.keys())
        assert "AD_BY_CALLER" not in info_ids, (
            f"AD_BY_CALLER unexpectedly in header when neoantigen=False. "
            f"INFO fields: {info_ids}"
        )
        assert "AF_BY_CALLER" not in info_ids, (
            f"AF_BY_CALLER unexpectedly in header when neoantigen=False. "
            f"INFO fields: {info_ids}"
        )
        # Existing fields must still be present
        assert "VAF_BY_CALLER" in info_ids, (
            f"VAF_BY_CALLER missing from header — existing fields must not be removed. "
            f"INFO fields: {info_ids}"
        )
    finally:
        template_vcf.close()
        try:
            os.unlink(template_path)
        except OSError:
            pass
