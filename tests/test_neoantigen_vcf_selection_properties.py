"""
Property-based tests for Neoantigen VCF Selection Correctness.

# Feature: neoantigen-workflow, Property 7: Neoantigen VCF Selection Correctness

Validates: Requirements 3.1, 3.2

For any pipeline run with enable_neoantigen_workflow=True, the neoantigen output
VCF must contain only samples with status=1, and the VCF source must match
neoantigen_input_source: the consensus VCF when 'consensus', the Mutect2-filtered
VCF when 'mutect2'.

The Nextflow channel filtering logic being tested:

    # When neoantigen_input_source == 'consensus':
    selected = [s for s in consensus_vcfs if s['meta']['status'] == 1]

    # When neoantigen_input_source == 'mutect2':
    selected = [s for s in mutect2_vcfs if s['meta']['status'] == 1]
"""

from hypothesis import given, settings
from hypothesis import strategies as st

# ---------------------------------------------------------------------------
# Strategies
# ---------------------------------------------------------------------------

status_st = st.sampled_from([0, 1, 2])

sample_id_st = st.text(
    alphabet=st.characters(whitelist_categories=("Lu", "Ll", "Nd"), whitelist_characters="_-"),
    min_size=1,
    max_size=20,
)

sample_st = st.fixed_dictionaries(
    {
        "id": sample_id_st,
        "status": status_st,
    }
)

sample_set_st = st.lists(sample_st, min_size=1, max_size=20)

neoantigen_input_source_st = st.sampled_from(["consensus", "mutect2"])


# ---------------------------------------------------------------------------
# Helper: build VCF channel entries from a sample list
# ---------------------------------------------------------------------------


def _build_vcf_channel(samples: list, source: str) -> list:
    """
    Build a list of VCF channel entries (simulating a Nextflow channel) for the
    given samples. Each entry is a dict with 'meta' and 'source' keys.

    'source' indicates which VCF type this entry belongs to ('consensus' or 'mutect2').
    """
    return [
        {
            "meta": {"id": s["id"], "status": s["status"]},
            "vcf": f"{s['id']}.{source}.vcf.gz",
            "tbi": f"{s['id']}.{source}.vcf.gz.tbi",
            "source": source,
        }
        for s in samples
    ]


# ---------------------------------------------------------------------------
# Helper: simulate the Nextflow VCF selection logic
# ---------------------------------------------------------------------------


def apply_vcf_selection(
    consensus_vcfs: list,
    mutect2_vcfs: list,
    neoantigen_input_source: str,
) -> list:
    """
    Simulate the Nextflow channel filtering logic from NEOANTIGEN_WORKFLOW:

        ch_source_vcf = params.neoantigen_input_source == 'consensus'
                        ? ch_consensus_vcf : ch_mutect2_vcf
        ch_neoantigen_vcf = ch_source_vcf.filter { meta, vcf, tbi -> meta.status == 1 }

    Returns only status=1 entries from the appropriate source channel.
    """
    if neoantigen_input_source == "consensus":
        source_vcfs = consensus_vcfs
    else:
        source_vcfs = mutect2_vcfs

    return [entry for entry in source_vcfs if entry["meta"]["status"] == 1]


# ---------------------------------------------------------------------------
# Property 7: Neoantigen VCF Selection Correctness
# ---------------------------------------------------------------------------


@given(samples=sample_set_st, neoantigen_input_source=neoantigen_input_source_st)
@settings(max_examples=100)
def test_p7_only_status1_samples_selected(samples, neoantigen_input_source):
    """
    **Validates: Requirements 3.1, 3.2**

    Property 7: Neoantigen VCF Selection Correctness.
    After applying the VCF selection logic, the result must contain ONLY
    samples with status=1 — no status=0 or status=2 samples may appear.

    # Feature: neoantigen-workflow, Property 7: Neoantigen VCF Selection Correctness
    """
    consensus_vcfs = _build_vcf_channel(samples, "consensus")
    mutect2_vcfs = _build_vcf_channel(samples, "mutect2")

    selected = apply_vcf_selection(consensus_vcfs, mutect2_vcfs, neoantigen_input_source)

    for entry in selected:
        assert entry["meta"]["status"] == 1, (
            f"Sample {entry['meta']['id']} with status={entry['meta']['status']} "
            f"should not be selected. Only status=1 samples should appear in the "
            f"neoantigen VCF output."
        )


@given(samples=sample_set_st, neoantigen_input_source=neoantigen_input_source_st)
@settings(max_examples=100)
def test_p7_all_status1_samples_are_selected(samples, neoantigen_input_source):
    """
    **Validates: Requirements 3.1, 3.2**

    Property 7: Neoantigen VCF Selection Correctness.
    Every sample with status=1 must appear in the selected output — none may be dropped.

    # Feature: neoantigen-workflow, Property 7: Neoantigen VCF Selection Correctness
    """
    consensus_vcfs = _build_vcf_channel(samples, "consensus")
    mutect2_vcfs = _build_vcf_channel(samples, "mutect2")

    selected = apply_vcf_selection(consensus_vcfs, mutect2_vcfs, neoantigen_input_source)
    selected_ids = [entry["meta"]["id"] for entry in selected]

    for s in samples:
        if s["status"] == 1:
            assert s["id"] in selected_ids, (
                f"Sample {s['id']} with status=1 is missing from the selected output. "
                f"All status=1 samples must be included in the neoantigen VCF."
            )


@given(samples=sample_set_st, neoantigen_input_source=neoantigen_input_source_st)
@settings(max_examples=100)
def test_p7_selected_vcfs_come_from_correct_source(samples, neoantigen_input_source):
    """
    **Validates: Requirements 3.1, 3.2**

    Property 7: Neoantigen VCF Selection Correctness.
    All selected VCF entries must originate from the source matching
    neoantigen_input_source: 'consensus' entries when source='consensus',
    'mutect2' entries when source='mutect2'.

    # Feature: neoantigen-workflow, Property 7: Neoantigen VCF Selection Correctness
    """
    consensus_vcfs = _build_vcf_channel(samples, "consensus")
    mutect2_vcfs = _build_vcf_channel(samples, "mutect2")

    selected = apply_vcf_selection(consensus_vcfs, mutect2_vcfs, neoantigen_input_source)

    for entry in selected:
        assert entry["source"] == neoantigen_input_source, (
            f"Sample {entry['meta']['id']} VCF comes from source='{entry['source']}' "
            f"but neoantigen_input_source='{neoantigen_input_source}'. "
            f"Selected VCFs must come from the configured source."
        )


@given(samples=sample_set_st)
@settings(max_examples=100)
def test_p7_consensus_source_does_not_include_mutect2_vcfs(samples):
    """
    **Validates: Requirements 3.1, 3.2**

    Property 7: Neoantigen VCF Selection Correctness — source isolation.
    When neoantigen_input_source='consensus', no mutect2 VCF entries appear
    in the selected output.

    # Feature: neoantigen-workflow, Property 7: Neoantigen VCF Selection Correctness
    """
    consensus_vcfs = _build_vcf_channel(samples, "consensus")
    mutect2_vcfs = _build_vcf_channel(samples, "mutect2")

    selected = apply_vcf_selection(consensus_vcfs, mutect2_vcfs, "consensus")

    for entry in selected:
        assert entry["source"] != "mutect2", (
            f"Sample {entry['meta']['id']} has a mutect2 VCF in the selected output "
            f"when neoantigen_input_source='consensus'. Only consensus VCFs should be selected."
        )


@given(samples=sample_set_st)
@settings(max_examples=100)
def test_p7_mutect2_source_does_not_include_consensus_vcfs(samples):
    """
    **Validates: Requirements 3.1, 3.2**

    Property 7: Neoantigen VCF Selection Correctness — source isolation.
    When neoantigen_input_source='mutect2', no consensus VCF entries appear
    in the selected output.

    # Feature: neoantigen-workflow, Property 7: Neoantigen VCF Selection Correctness
    """
    consensus_vcfs = _build_vcf_channel(samples, "consensus")
    mutect2_vcfs = _build_vcf_channel(samples, "mutect2")

    selected = apply_vcf_selection(consensus_vcfs, mutect2_vcfs, "mutect2")

    for entry in selected:
        assert entry["source"] != "consensus", (
            f"Sample {entry['meta']['id']} has a consensus VCF in the selected output "
            f"when neoantigen_input_source='mutect2'. Only mutect2 VCFs should be selected."
        )


@given(samples=sample_set_st, neoantigen_input_source=neoantigen_input_source_st)
@settings(max_examples=100)
def test_p7_selected_count_equals_status1_count(samples, neoantigen_input_source):
    """
    **Validates: Requirements 3.1, 3.2**

    Property 7: Neoantigen VCF Selection Correctness.
    The total number of selected entries must equal the number of status=1
    samples in the input set.

    # Feature: neoantigen-workflow, Property 7: Neoantigen VCF Selection Correctness
    """
    consensus_vcfs = _build_vcf_channel(samples, "consensus")
    mutect2_vcfs = _build_vcf_channel(samples, "mutect2")

    selected = apply_vcf_selection(consensus_vcfs, mutect2_vcfs, neoantigen_input_source)
    expected_count = sum(1 for s in samples if s["status"] == 1)

    assert len(selected) == expected_count, (
        f"Selected {len(selected)} VCF entries but expected {expected_count} "
        f"(number of status=1 samples). "
        f"neoantigen_input_source='{neoantigen_input_source}'."
    )


@given(samples=sample_set_st, neoantigen_input_source=neoantigen_input_source_st)
@settings(max_examples=100)
def test_p7_empty_result_when_no_status1_samples(samples, neoantigen_input_source):
    """
    **Validates: Requirements 3.1, 3.2**

    Property 7: Neoantigen VCF Selection Correctness — empty channel case.
    When no samples have status=1, the selected output must be empty
    (matching the log.warn behavior in NEOANTIGEN_WORKFLOW).

    # Feature: neoantigen-workflow, Property 7: Neoantigen VCF Selection Correctness
    """
    # Override all statuses to non-1 values
    non_status1_samples = [{"id": s["id"], "status": 0 if s["status"] != 2 else 2} for s in samples]

    consensus_vcfs = _build_vcf_channel(non_status1_samples, "consensus")
    mutect2_vcfs = _build_vcf_channel(non_status1_samples, "mutect2")

    selected = apply_vcf_selection(consensus_vcfs, mutect2_vcfs, neoantigen_input_source)

    assert len(selected) == 0, (
        f"Expected empty selection when no status=1 samples present, "
        f"but got {len(selected)} entries. "
        f"neoantigen_input_source='{neoantigen_input_source}'."
    )
