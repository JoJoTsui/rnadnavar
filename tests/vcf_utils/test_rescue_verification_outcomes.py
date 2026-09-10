"""Verification may withhold positive RNA-dependent calls, never erase vetoes."""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "bin"))
from vcf_utils.classification import compute_unified_classification_rescue
from vcf_utils.variant_statistics import compute_rescue_statistics, print_statistics


def classify(callers, labels, status="inconclusive", eligible=None, rescue_min_rna=2):
    modalities = {caller: caller.split("_")[0] for caller in callers}
    data = {"callers": callers, "filters_normalized": labels, "is_snv": True,
            "dna_verification_status": status}
    if eligible is not None:
        data["support_callers"] = eligible
    result = compute_unified_classification_rescue(
        data, modalities,
        rescue_config={"rescue_promotion_min_rna_callers": rescue_min_rna},
        with_rationale=True,
    )
    return result, data


@pytest.mark.parametrize("label", ["Artifact", "Germline", "Reference", "RNAedit"])
@pytest.mark.parametrize("status", ["inconclusive", "rejected", "confirmed"])
def test_verification_preserves_negative_dna_classification(label, status):
    (result, rationale), data = classify(
        ["DNA_consensus"], [label], status)
    assert result == label
    assert "rule:dna_verification" not in rationale
    assert not data["rescue_promoted"]
    assert not data["rescued"]


def test_dna_artifact_veto_precedes_verification():
    (result, rationale), _ = classify(
        ["DNA_consensus", "RNA_consensus"], ["Artifact", "Somatic"])
    assert result == "Artifact"
    assert "dna_artifact_veto" in rationale


def test_verification_withholds_rna_only_somatic():
    (result, rationale), data = classify(["RNA_consensus"], ["Somatic"])
    assert result == data["final_classification"] == "NoConsensus"
    assert "rule:rna_without_dna_gate|class:NoConsensus" in rationale
    assert not data["rescue_promoted"]


@pytest.mark.parametrize("status,expected", [("confirmed", "Somatic"), (None, "Somatic"),
                                              ("inconclusive", "NoConsensus")])
def test_one_dna_vote_does_not_bypass_promotion_gate(status, expected):
    (result, _), data = classify(
        ["DNA_mutect2", "RNA_strelka"], ["Somatic", "Somatic"], status,
        rescue_min_rna=1,
    )
    assert result == expected
    assert data["rescue_promoted"] == (expected == "Somatic")
    assert data["rescued"] == (expected == "Somatic")


def test_dna_consensus_somatic_is_independent_of_verification():
    (result, _), _ = classify(["DNA_consensus", "RNA_consensus"], ["Somatic", "Somatic"])
    assert result == "Somatic"


def test_dna_votes_must_be_eligible_for_verification_bypass():
    callers = ["DNA_mutect2", "DNA_strelka", "RNA_consensus"]
    labels = ["Somatic"] * 3
    assert classify(callers, labels)[0][0] == "NoConsensus"
    assert classify(callers, labels, eligible=["DNA_mutect2"])[0][0] == "NoConsensus"


def test_rescue_outcomes_are_bounded_and_reconcile(capsys):
    records = {
        "retained": {"final_classification": "Somatic", "is_snv": True, "rescued": True},
        "new": {"final_classification": "Somatic", "is_snv": False, "rescued": True},
        "lost": {"final_classification": "Artifact", "is_snv": True, "rescued": True},
    }
    dna = {key: {"filter_normalized": "Somatic"} for key in ["retained", "lost"]}
    stats = compute_rescue_statistics(records, dna, {})
    assert stats["cross_modality"] == 0  # Promotions need not be consensus overlap.
    assert stats["rescued"] == 2
    assert stats["rescued_union_fraction"] == pytest.approx(2 / 3)
    assert stats["dna_somatic_baseline"] == 2
    assert stats["output_somatic"] == 2
    assert stats["somatic_retained_from_dna"] == 1
    assert stats["somatic_new_vs_dna"] == 1
    assert stats["somatic_lost_from_dna"] == 1
    assert "rescue_rate" not in stats
    print_statistics(stats, "rescue")
    output = capsys.readouterr().out
    assert "Rescued / union records: 2/3 (66.67%)" in output
    assert "not truth-based accuracy" in output


def test_empty_union_and_required_final_labels():
    assert compute_rescue_statistics({}, {}, {})["rescued_union_fraction"] == 0
    with pytest.raises(ValueError, match="finalized"):
        compute_rescue_statistics({"x": {"is_snv": True}}, {}, {})
    with pytest.raises(TypeError, match="labeled"):
        compute_rescue_statistics({}, set(), set())
