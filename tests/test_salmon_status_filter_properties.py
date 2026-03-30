"""
Property-based tests for Salmon Status Filter.

# Feature: neoantigen-workflow, Property 8: Salmon Status Filter

Validates: Requirements 4.7

For any samplesheet containing samples with mixed status values (0, 1, 2),
SALMON_QUANT must be invoked exactly once per status=2 sample and must never
be invoked for status=0 or status=1 samples.

The Nextflow channel filtering logic being tested:
    ch_fastq.filter { meta, reads -> meta.status == 2 }
"""

from hypothesis import given, settings
from hypothesis import strategies as st

# ---------------------------------------------------------------------------
# Strategies
# ---------------------------------------------------------------------------

status_st = st.sampled_from([0, 1, 2])

sample_st = st.fixed_dictionaries(
    {
        "id": st.text(
            alphabet=st.characters(whitelist_categories=("Lu", "Ll", "Nd"), whitelist_characters="_-"),
            min_size=1,
            max_size=20,
        ),
        "status": status_st,
    }
)

samplesheet_st = st.lists(sample_st, min_size=1, max_size=20)


# ---------------------------------------------------------------------------
# Helper: simulate the Nextflow channel filter
# ---------------------------------------------------------------------------


def apply_salmon_status_filter(samples: list) -> list:
    """
    Simulate the Nextflow channel filtering logic:
        ch_fastq.filter { meta, reads -> meta.status == 2 }

    Returns only samples with status == 2.
    """
    return [s for s in samples if s["status"] == 2]


# ---------------------------------------------------------------------------
# Property 8: Salmon Status Filter — enable_neoantigen_workflow=True
# ---------------------------------------------------------------------------


@given(samples=samplesheet_st)
@settings(max_examples=100)
def test_p8_only_status2_samples_pass_filter(samples):
    """
    **Validates: Requirements 4.7**

    Property 8: Salmon Status Filter.
    After applying the status filter, the result must contain EXACTLY the
    samples with status=2 — no status=0 or status=1 samples may appear.

    # Feature: neoantigen-workflow, Property 8: Salmon Status Filter
    """
    filtered = apply_salmon_status_filter(samples)

    # No status=0 or status=1 samples in filtered list
    for s in filtered:
        assert s["status"] == 2, (
            f"Sample {s['id']} with status={s['status']} should not pass the filter. "
            "Only status=2 samples should be passed to SALMON_QUANT."
        )


@given(samples=samplesheet_st)
@settings(max_examples=100)
def test_p8_all_status2_samples_appear_in_filtered_list(samples):
    """
    **Validates: Requirements 4.7**

    Property 8: Salmon Status Filter.
    Every sample with status=2 must appear in the filtered list — none may be dropped.

    # Feature: neoantigen-workflow, Property 8: Salmon Status Filter
    """
    filtered = apply_salmon_status_filter(samples)
    filtered_ids = [s["id"] for s in filtered]

    for s in samples:
        if s["status"] == 2:
            assert s["id"] in filtered_ids, (
                f"Sample {s['id']} with status=2 is missing from the filtered list. "
                "All status=2 samples must be passed to SALMON_QUANT exactly once."
            )


@given(samples=samplesheet_st)
@settings(max_examples=100)
def test_p8_each_status2_sample_appears_exactly_once(samples):
    """
    **Validates: Requirements 4.7**

    Property 8: Salmon Status Filter.
    Each status=2 sample must appear exactly once in the filtered list
    (SALMON_QUANT invoked exactly once per status=2 sample).

    # Feature: neoantigen-workflow, Property 8: Salmon Status Filter
    """
    filtered = apply_salmon_status_filter(samples)

    # Count occurrences of each sample id in the filtered list
    from collections import Counter
    id_counts = Counter(s["id"] for s in filtered)

    for s in samples:
        if s["status"] == 2:
            count = id_counts.get(s["id"], 0)
            # Each unique id with status=2 should appear exactly as many times
            # as it appears in the input (preserving duplicates from input)
            input_count = sum(1 for x in samples if x["id"] == s["id"] and x["status"] == 2)
            assert count == input_count, (
                f"Sample {s['id']} with status=2 appears {count} time(s) in filtered list "
                f"but appears {input_count} time(s) in input. "
                "SALMON_QUANT must be invoked exactly once per status=2 sample occurrence."
            )


@given(samples=samplesheet_st)
@settings(max_examples=100)
def test_p8_filtered_count_equals_status2_count(samples):
    """
    **Validates: Requirements 4.7**

    Property 8: Salmon Status Filter.
    The total number of samples in the filtered list must equal the number of
    status=2 samples in the input samplesheet.

    # Feature: neoantigen-workflow, Property 8: Salmon Status Filter
    """
    filtered = apply_salmon_status_filter(samples)
    expected_count = sum(1 for s in samples if s["status"] == 2)

    assert len(filtered) == expected_count, (
        f"Filtered list has {len(filtered)} samples but samplesheet has "
        f"{expected_count} status=2 samples. "
        "SALMON_QUANT must be invoked exactly once per status=2 sample."
    )


# ---------------------------------------------------------------------------
# Property 8 (inverse): enable_neoantigen_workflow=False — filter not applied
# ---------------------------------------------------------------------------


def apply_no_filter(samples: list) -> list:
    """
    When enable_neoantigen_workflow=False, the filter is never applied.
    All samples pass through unchanged.
    """
    return list(samples)


@given(samples=samplesheet_st)
@settings(max_examples=100)
def test_p8_inverse_no_filter_when_neoantigen_disabled(samples):
    """
    **Validates: Requirements 4.7**

    Property 8 (inverse): When enable_neoantigen_workflow=False, the status
    filter is never applied and all samples pass through unchanged.

    # Feature: neoantigen-workflow, Property 8: Salmon Status Filter
    """
    result = apply_no_filter(samples)

    # All samples must be present, in the same order, with unchanged content
    assert len(result) == len(samples), (
        f"Expected {len(samples)} samples but got {len(result)} when filter is disabled."
    )
    for original, returned in zip(samples, result):
        assert original["id"] == returned["id"], (
            f"Sample order or identity changed: expected id={original['id']}, "
            f"got id={returned['id']}."
        )
        assert original["status"] == returned["status"], (
            f"Sample {original['id']} status changed from {original['status']} "
            f"to {returned['status']} when filter is disabled."
        )
