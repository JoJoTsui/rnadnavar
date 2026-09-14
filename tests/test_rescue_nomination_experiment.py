"""Tests for the isolated rescue investigation, not a production policy switch."""
import gzip
from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "examples/seqc2/scripts"))
from test_rescue_nomination_gate import biological_veto, nominated, tumor_alt
from investigate_verified_rescue import scan

KEY = ("chr1", 10, "A", "C")


def candidate(filt="RefCall", ad="39,3"):
    return {"filter": filt, "alt_index": 1,
            "samples": {"TUMOR": {"AD": ad}, "NORMAL": {"AD": "0,99"}}}


def test_reference_observation_is_not_variant_nomination():
    panel = {"deepsomatic": {KEY: [candidate()]}, "mutect2": {}, "strelka": {}}
    assert not nominated(panel, KEY)
    panel["mutect2"][KEY] = [candidate("weak_evidence", "20,1")]
    assert nominated(panel, KEY)  # Deliberately not a PASS/Somatic-vote requirement.


@pytest.mark.parametrize("ad", ["20,0", ".", "20,."])
def test_normal_alt_reads_cannot_supply_missing_tumor_evidence(ad):
    panel = {"deepsomatic": {KEY: [candidate("PASS", ad)]}, "mutect2": {}, "strelka": {}}
    assert not nominated(panel, KEY)


def test_sample_ambiguity_fails_closed_and_multiallelic_ad_uses_correct_alt():
    row = candidate("PASS", "20,0,4")
    row["alt_index"] = 2
    assert tumor_alt(row, "mutect2", "G") == 4
    row["samples"]["OTHER"] = {"AD": "0,99"}
    assert tumor_alt(row, "mutect2", "G") is None


def test_strelka_uses_tumor_tier_one_not_normal_or_tier_two():
    row = {"samples": {"WES_LL_N_1": {"CU": "20,30"}, "WES_LL_T_1": {"CU": "0,8"}}}
    assert tumor_alt(row, "strelka", "C") == 0


def test_biological_veto_threshold_missingness_and_editing():
    assert not biological_veto([{"info": {}}])
    assert not biological_veto([{"info": {"GNOMAD_AF": "."}}])
    assert not biological_veto([{"info": {"GNOMAD_AF": "0.001"}}])
    assert biological_veto([{"info": {"GNOMAD_AF": "0,0.0995979"}}])
    editing = {"REDI_ACCESSION": "EDHSAAAK8017", "REDI_CANONICAL": "YES", "N_DNA_CALLERS_SOMATIC": "0"}
    assert biological_veto([{"info": editing}])
    assert not biological_veto([{"info": dict(editing, N_DNA_CALLERS_SOMATIC="1")}])
    with pytest.raises(ValueError, match="Malformed"):
        biological_veto([{"info": {"GNOMAD_AF": "bad"}}])


def test_extractor_preserves_sample_fields_and_alt_index(tmp_path):
    path = tmp_path / "input.vcf.gz"
    with gzip.open(path, "wt") as handle:
        handle.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\tTUMOR\n")
        handle.write("chr1\t10\t.\tA\tG,C\t0\tRefCall\tFLAG;DP=30\tGT:AD\t0/0:20,0,0\t0/2:20,0,3\n")
    row = scan(path, {KEY})[KEY][0]
    assert row["alt_index"] == 2
    assert row["info"] == {"FLAG": True, "DP": "30"}
    assert tumor_alt(row, "deepsomatic", "C") == 3
    assert scan(path, {("chr2", 10, "A", "C")}) == {}
