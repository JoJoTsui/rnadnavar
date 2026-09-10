from vcf_utils.classification import compute_unified_classification_consensus


def test_native_evidence_snv_promotes_with_rationale():
    data = {"is_snv": True, "native_evidence_pass": True, "callers": [], "filters_normalized": []}
    label, rationale = compute_unified_classification_consensus(data, 2, 2, with_rationale=True)
    assert label == "Somatic"
    assert "native_evidence_snv" in rationale


def test_native_evidence_does_not_promote_indel():
    data = {"is_snv": False, "native_evidence_pass": True, "callers": [], "filters_normalized": []}
    label = compute_unified_classification_consensus(data, 2, 2)
    assert label == "NoConsensus"
