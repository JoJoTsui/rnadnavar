import pytest

from vcf_utils.refined_rescue_policy import biological_veto, decide, nominates, rna_supports


@pytest.mark.parametrize("count,expected", [(None, False), (0, False), (2, False), (3, True), (4, True)])
def test_rna_pass_requires_eligible_alt_floor(count, expected):
    assert rna_supports("mutect2", "PASS", count) == expected
    assert not rna_supports("mutect2", "weak_evidence", count)


def test_filtered_dna_nomination_is_not_somatic_vote():
    assert nominates("mutect2", "weak_evidence", 1)
    assert nominates("strelka", "LowEVS", 1)
    assert nominates("deepsomatic", "PASS", 1)
    assert not nominates("deepsomatic", "RefCall", 10)
    assert not nominates("unknown", "PASS", 10)


@pytest.mark.parametrize("value", [None, ".", -1, 0, 0.5, float("nan"), float("inf")])
def test_invalid_or_zero_alt_cannot_nominate(value):
    assert not nominates("mutect2", "PASS", value)


def test_duplicate_rna_rows_and_unknown_callers_do_not_create_votes():
    for callers in (["mutect2", "mutect2"], ["mutect2", "unknown"], ["RNA_consensus"]):
        assert not decide("A", "G", "Somatic", {"strelka"}, callers, {})[0]
    assert decide("A", "G", "Somatic", {"strelka"}, {"mutect2", "deepsomatic"}, {})[0]


@pytest.mark.parametrize("value", ["nan", "inf", "oops", "-0.01", "1.01", "0.0001,."])
def test_available_malformed_population_af_fails_closed(value):
    assert biological_veto({"gnomAD_AF": value}) == "invalid_population_af"


def test_population_threshold_and_editing_guard():
    assert biological_veto({"GNOMAD_AF": "0.001"}) is None
    assert biological_veto({"GNOMAD_AF": "0.0011"}) == "common_population_af"
    info = {"REDI_ACCESSION": "record", "REDI_CANONICAL": "YES", "N_DNA_CALLERS_SOMATIC": "0"}
    assert biological_veto(info) == "canonical_editing_without_dna_somatic"
    info["N_DNA_CALLERS_SOMATIC"] = "1"
    assert biological_veto(info) is None


@pytest.mark.parametrize("ref,alt,label", [("A", "AT", "Somatic"), ("A", "*", "Somatic"),
                                             ("A", "C,G", "Somatic"), ("A", "G", "Artifact")])
def test_only_somatic_snp_candidates(ref, alt, label):
    assert not decide(ref, alt, label, {"strelka"}, {"mutect2", "deepsomatic"}, {})[0]
