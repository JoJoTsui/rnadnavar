from copy import deepcopy

import pytest

from vcf_utils.refined_native_policy import evaluate, trace
from vcf_utils.classification import compute_unified_classification_consensus


def candidate():
    return {
        "REF": "AC", "ALT": "A", "is_snv": False,
        "callers": ["deepsomatic", "mutect2"],
        "filters_original": ["PASS", "contamination;weak_evidence"],
        "filters_normalized": ["Artifact", "Artifact"],
        "qualities_by_caller": {"deepsomatic": 20},
        "native_evidence": {
            "deepsomatic": {"tumor_AD": [20, 3]},
            "mutect2": {"tumor_AD": [20, 2], "normal_AD": [10, 0],
                        "normal_DP": [10], "GERMQ": 20, "TLOD": 1},
        },
    }


def test_three_branches_and_missing_ecnt_independence():
    data = candidate()
    assert evaluate(data) == "indel_high_confidence"
    data["qualities_by_caller"]["deepsomatic"] = 10
    assert evaluate(data) is None
    data["native_evidence"]["mutect2"]["tumor_SB"] = [10, 10, 1, 1]
    assert evaluate(data) == "indel_moderate_strand"
    data["qualities_by_caller"]["deepsomatic"] = 1
    data["filters_original"] = ["RefCall", "PASS"]
    data["native_evidence"]["mutect2"].update(tumor_AD=[20, 3], ECNT=1)
    assert evaluate(data) == "indel_reciprocal_single_event"
    data["native_evidence"]["mutect2"]["ECNT"] = 2
    assert evaluate(data) is None


@pytest.mark.parametrize("field,value", [
    ("normal_DP", None), ("normal_DP", [float("nan")]),
    ("normal_DP", [9]), ("normal_DP", [10.5]),
    ("normal_AD", [10, 1]), ("normal_AD", [10, -1]),
    ("GERMQ", float("inf")), ("GERMQ", [20, 30]),
    ("TLOD", 0), ("tumor_AD", [20, 1]),
])
def test_required_native_measurements_fail_closed(field, value):
    data = candidate()
    data["native_evidence"]["mutect2"][field] = value
    assert evaluate(data) is None


def test_snp_soft_tokens_not_exact_string_veto():
    data = candidate()
    data.update(is_snv=True, REF="C", ALT="T")
    data["filters_original"] = ["RefCall", "weak_evidence;contamination"]
    data["native_evidence"]["mutect2"].update(TLOD=12, GERMQ=60)
    assert evaluate(data) == "snv_soft_corroborated"
    data["filters_original"][1] += ";orientation"
    assert evaluate(data) is None
    data["filters_original"][0] = "PASS"
    assert evaluate(data) == "snv_native_pass"


def test_classifier_opt_in_vetoes_and_legacy_unchanged():
    data = candidate()
    old = compute_unified_classification_consensus(data, 2, 2)
    data.update(refined_native_enabled=True, refined_native_branch=evaluate(data))
    data["refined_native_trace"] = trace(data, data["refined_native_branch"])
    label, rationale = compute_unified_classification_consensus(data, 2, 2, with_rationale=True)
    assert label == "Somatic" and "indel_high_confidence" in rationale
    assert "normal_DP:10.0" in rationale
    for veto in ({"GNOMAD_AF": 0.01}, {"dna_verification_status": "rejected"}):
        blocked = deepcopy(data)
        blocked.update(veto)
        assert compute_unified_classification_consensus(blocked, 2, 2) != "Somatic"
    data["refined_native_enabled"] = False
    assert compute_unified_classification_consensus(data, 2, 2) == old


def test_unadmitted_indel_retains_threshold_but_snp_does_not_fallback():
    data = candidate()
    data["native_evidence"] = {}
    data["filters_normalized"] = ["Somatic", "Somatic"]
    data["support_callers"] = set(data["callers"])
    data["passes_consensus"] = True
    data.update(refined_native_enabled=True, refined_native_branch=None,
                refined_native_trace=trace(data, "not_admitted"))
    assert compute_unified_classification_consensus(data, 2, 2) == "Somatic"
    data["is_snv"] = True
    assert compute_unified_classification_consensus(data, 2, 2) == "NoConsensus"


@pytest.mark.parametrize("alt", ["A,T", "<DEL>", "*", "."])
def test_non_biallelic_sequence_indels_fail_closed(alt):
    data = candidate()
    data["ALT"] = alt
    assert evaluate(data) is None


@pytest.mark.parametrize("reverse", [False, True])
def test_reader_uses_header_roles_and_raw_dp(tmp_path, reverse):
    from vcf_utils.aggregation import read_variants_from_vcf
    samples = [("case-17", "0/1:20:17,3:8,9,1,2"),
               ("control-42", "0/0:.:10,0:5,5,0,0")]
    if reverse:
        samples.reverse()
    path = tmp_path / "case.mutect2.vcf"
    path.write_text(
        '##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n'
        '##normal_sample=control-42\n'
        '##INFO=<ID=ECNT,Number=1,Type=Integer,Description="Event count">\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">\n'
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depth">\n'
        '##FORMAT=<ID=SB,Number=4,Type=Integer,Description="Strand counts">\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'
        + '\t'.join(name for name, _ in samples) + '\n'
        + 'chr1\t10\t.\tAC\tA\t50\tPASS\tECNT=1\tGT:DP:AD:SB\t'
        + '\t'.join(value for _, value in samples) + '\n'
    )
    result = next(iter(read_variants_from_vcf(path, "mutect2", refined_native=True).values()))
    evidence = result["native_evidence"]
    assert evidence["tumor_AD"] == [17, 3]
    assert evidence["normal_AD"] == [10, 0]
    assert evidence["tumor_SB"] == [8, 9, 1, 2]
    data = candidate()
    data["native_evidence"]["mutect2"].update(evidence)
    assert evaluate(data) is None  # DP must not be derived from normal AD.
    legacy = next(iter(read_variants_from_vcf(path, "mutect2").values()))
    assert "tumor_SB" not in legacy["native_evidence"]
