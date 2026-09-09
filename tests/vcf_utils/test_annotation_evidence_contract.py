from vcf_utils.variant_classifier import VariantClassifier


def _info(**overrides):
    base = {
        "FILTER": "NoConsensus",
        "N_DNA_CALLERS": 3,
        "N_RNA_CALLERS": 0,
        "FILTERS_NORMALIZED": "DNA_mutect2:Somatic|DNA_strelka:Somatic|DNA_deepsomatic:Reference",
        "AD_BY_CALLER": "DNA_mutect2:36,2|DNA_strelka:36,2|DNA_deepsomatic:36,2",
    }
    base.update(overrides)
    return base


def test_annotation_does_not_resurrect_below_alt_floor():
    result = VariantClassifier().classify_variant_from_info(_info(FILTER="Somatic", CLASSIFICATION_RATIONALE="rule:insufficient_modality_support|class:NoConsensus"))
    assert result.classification == "NoConsensus"


def test_annotation_uses_eligible_somatic_support():
    result = VariantClassifier().classify_variant_from_info(
        _info(AD_BY_CALLER="DNA_mutect2:30,4|DNA_strelka:30,4|DNA_deepsomatic:30,0")
    )
    assert result.classification == "Somatic"


def test_rejected_dna_verification_blocks_annotation_promotion():
    result = VariantClassifier().classify_variant_from_info(
        _info(
            FILTER="NoConsensus",
            DNA_VERIFICATION="rejected",
            FILTERS_NORMALIZED=(
                "DNA_mutect2:Somatic|DNA_strelka:Somatic|DNA_deepsomatic:Somatic|"
                "RNA_mutect2:Somatic"
            ),
            AD_BY_CALLER=(
                "DNA_mutect2:30,4|DNA_strelka:30,4|DNA_deepsomatic:30,4|"
                "RNA_mutect2:30,4"
            ),
            COSMIC_CNT=20,
        )
    )
    assert result.classification == "NoConsensus"


def test_malformed_or_missing_ad_is_not_new_evidence_for_stale_somatic():
    result = VariantClassifier().classify_variant_from_info(
        _info(
            FILTER="Somatic",
            CLASSIFICATION_RATIONALE="rule:insufficient_modality_support|class:NoConsensus",
            AD_BY_CALLER="DNA_mutect2:missing|DNA_strelka:30|DNA_deepsomatic:bad,alt",
        )
    )
    assert result.classification == "NoConsensus"


def test_opt_in_baseline_retention_respects_af_and_verification_vetoes():
    from vcf_utils.classification import compute_unified_classification_consensus

    base = {
        "callers": ["deepsomatic"],
        "filters_normalized": ["Somatic"],
        "support_callers": {"deepsomatic"},
        "preserve_baseline": True,
        "is_snv": True,
    }
    assert compute_unified_classification_consensus(base, 2, 2) == "Somatic"
    assert compute_unified_classification_consensus(
        {**base, "gnomad_af": 0.01}, 2, 2
    ) != "Somatic"
    assert compute_unified_classification_consensus(
        {**base, "dna_verification_status": "inconclusive"}, 2, 2
    ) != "Somatic"
