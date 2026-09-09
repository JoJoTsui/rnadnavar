import pytest

from vcf_utils.caller_evidence import (
    bind_evidence, decode_evidence, encode_evidence, evidence_from_record, merge_evidence, read_info,
)


def raw(**kwargs):
    return evidence_from_record({}, "mutect2", "1:10:A:T", genotype={"GT": "0|1", "AD": "9,3", "DP": 12}, **kwargs)


def test_round_trip_is_vcf_safe_and_missing_normal_explicit():
    entries = raw(modality="DNA", sample_id="tumor")
    encoded = encode_evidence(entries)
    assert not any(c in encoded for c in ';,= \t\n')
    assert decode_evidence(encoded) == entries
    normal = next(e for e in entries if e["sample_role"] == "normal")
    assert normal["available"] is False
    assert normal["measurements"] == {}


def test_conflicts_keep_measurements_and_merge_source_provenance():
    first = raw(source="one")
    second = raw(source="two")
    conflict = raw(source="three")
    next(e for e in conflict if e["sample_role"] == "tumor")["measurements"]["DP"] = 15
    merged = merge_evidence(first, second, conflict)
    assert merged == merge_evidence(conflict, second, first)
    tumors = [e for e in merged if e["sample_role"] == "tumor"]
    assert len(tumors) == 2
    assert next(e for e in tumors if e["measurements"]["DP"] == 12)["sources"] == ["one", "two"]


def test_legacy_phased_gt_and_ad_survive_and_missing_sources_enrich():
    legacy = {"GT_BY_CALLER": "DNA_mutect2:0|1|DNA_strelka:0/1",
              "AD_BY_CALLER": "DNA_mutect2:9,3|DNA_strelka:8,4"}
    entries = evidence_from_record(legacy, "DNA_consensus", "1:10:A:T")
    tumor = {e["caller"]: e["measurements"] for e in entries if e["sample_role"] == "tumor"}
    assert tumor["mutect2"] == {"GT": "0|1", "AD": "9,3"}
    assert tumor["strelka"]["AD"] == "8,4"
    legacy["CALLER_EVIDENCE"] = encode_evidence(raw(modality="DNA"))
    enriched = evidence_from_record(legacy, "DNA_consensus", "1:10:A:T")
    assert {e["caller"] for e in enriched} == {"mutect2", "strelka"}


def test_unknown_binding_does_not_relabel_known_modality_or_round():
    initial = raw()
    entries = evidence_from_record({"CALLER_EVIDENCE": encode_evidence(initial)},
                                   "RNA_consensus", "1:10:A:T", alignment_round="realigned")
    assert {e["modality"] for e in entries} == {"RNA"}
    assert {e["alignment_round"] for e in entries} == {"realigned"}
    rebound = evidence_from_record({"CALLER_EVIDENCE": encode_evidence(entries)},
                                   "DNA_consensus", "1:10:A:T", alignment_round="original")
    assert rebound == entries


def test_read_info_preserves_types_and_invalid_provenance():
    class Variant:
        INFO = [("AD", (10, 2)), ("FLAG", True), ("AF", 0.2)]
    assert read_info(Variant()) == {"AD": [10, 2], "FLAG": True, "AF": 0.2}
    entries = evidence_from_record({}, "strelka", "1:10:A:T",
                                   genotype={"DP": None, "INVALID_FIELDS": ["DP"], "raw": float("nan")})
    tumor = next(e for e in decode_evidence(encode_evidence(entries)) if e["sample_role"] == "tumor")
    assert tumor["measurements"]["INVALID_FIELDS"] == ["DP"]
    assert tumor["measurements"]["raw"] == "nan"


def test_invalid_schema_not_silently_dropped():
    with pytest.raises(ValueError):
        decode_evidence('{"version":"paired-v99","entries":[]}')


def test_round_identity_and_mixed_writer_binding():
    entries = merge_evidence(raw(modality="DNA"), raw(modality="RNA"))
    bound = bind_evidence(entries, alignment_round={"DNA": "first", "RNA": "realignment"})
    assert {(e["modality"], e["alignment_round"]) for e in bound} == {
        ("DNA", "first"), ("RNA", "realignment")}
    first_rna = raw(modality="RNA", alignment_round="first")
    assert len(merge_evidence(bound, first_rna)) == 6
    assert {e["alignment_round"] for e in entries} == {"unknown"}


def test_actual_vcf_round_trip(tmp_path):
    cyvcf2 = pytest.importorskip("cyvcf2")
    entries = raw(modality="DNA", normal_genotype={"GT": "0/0", "AD": "20,0"})
    path = tmp_path / "evidence.vcf"
    path.write_text(
        '##fileformat=VCFv4.2\n'
        '##contig=<ID=1>\n'
        '##INFO=<ID=CALLER_EVIDENCE,Number=1,Type=String,Description="Evidence">\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
        '1\t10\t.\tA\tT\t.\tPASS\tCALLER_EVIDENCE=' + encode_evidence(entries) + '\n')
    vcf = cyvcf2.VCF(str(path))
    info = read_info(next(vcf))
    assert decode_evidence(info["CALLER_EVIDENCE"]) == entries
    vcf.close()


def test_legacy_numeric_format_does_not_duplicate_canonical_measurement():
    entries = evidence_from_record({}, "mutect2", "1:10:A:T", modality="DNA",
                                   alignment_round="first", genotype={"VAF": 0.25001, "DP": 12})
    info = {"CALLER_EVIDENCE": encode_evidence(entries),
            "VAF_BY_CALLER": "DNA_mutect2:0.2500", "DP_BY_CALLER": "DNA_mutect2:12.0"}
    parsed = evidence_from_record(info, "DNA_consensus", "1:10:A:T", alignment_round="first")
    assert parsed == entries
    conflicting_round = evidence_from_record(info, "DNA_consensus", "1:10:A:T", alignment_round="realignment")
    assert len(conflicting_round) > len(entries)
    info["VAF_BY_CALLER"] = "DNA_mutect2:0.2600"
    assert len(evidence_from_record(info, "DNA_consensus", "1:10:A:T")) > len(entries)


def test_ambiguous_legacy_raw_text_preserved_without_invented_genotype():
    text = "DNA_mutect2:0|1||RNA_mutect2:0|1"
    entries = evidence_from_record({"GT_BY_CALLER": text}, "DNA_consensus", "1:10:A:T")
    assert len(entries) == 1
    assert entries[0]["available"] is False
    assert entries[0]["measurements"] == {"LEGACY_RAW_FIELDS": {"GT_BY_CALLER": text}}
