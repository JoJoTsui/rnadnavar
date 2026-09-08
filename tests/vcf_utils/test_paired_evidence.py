from vcf_utils.classification import compute_unified_classification_rescue
from vcf_utils.aggregation import (
    aggregate_genotypes,
    aggregate_variants,
    resolve_normal_sample_index,
    resolve_tumor_sample_index,
)


def test_paired_sample_resolution_is_explicit_and_complementary():
    samples = ["WES_LL_N_1", "WES_LL_T_1"]
    assert resolve_tumor_sample_index(samples, "mutect2") == 1
    assert resolve_normal_sample_index(samples, "mutect2") == 0


def test_tumor_only_vcf_has_unavailable_normal():
    assert resolve_normal_sample_index(["tumor"], "deepsomatic") is None


def test_ambiguous_multi_sample_vcf_does_not_guess_normal():
    assert resolve_normal_sample_index(["A", "B", "C"], "unknown") is None


def test_aggregate_genotypes_keeps_caller_af_distinct_from_ad():
    result = aggregate_genotypes(
        {
            "mutect2": {
                "GT": "0/1",
                "DP": 20,
                "AD": "18,2",
                "VAF": 0.071,
                "GQ": 40,
            }
        },
        ["mutect2"],
    )
    assert result["dp_by_caller"] == [20]
    assert result["vaf_by_caller"] == [0.071]
    assert result["alt_count_by_caller"] == [2]


def test_multiallelic_alt_support_uses_all_alternates():
    result = aggregate_genotypes(
        {"caller": {"GT": "1/2", "ALT_INDICES": [1, 2], "DP": 30, "AD": "10,3,7", "VAF": None, "GQ": 50}},
        ["caller"],
    )
    assert result["alt_count_by_caller"] == [7]


def test_unverified_rna_only_nomination_is_not_somatic():
    data = {
        "callers": ["RNA_consensus"],
        "filters_normalized": ["Somatic"],
        "caller_modality_map": {"RNA_consensus": "RNA"},
        "is_snv": True,
        "support_callers": {"RNA_consensus"},
        "dna_verification_status": "inconclusive",
    }
    assert compute_unified_classification_rescue(
        data, {"RNA_consensus": "RNA"}, snv_threshold=2, indel_threshold=2
    ) == "NoConsensus"


def test_source_evidence_survives_aggregation():
    variants = {
        "chr1:10:A:G": {
            "CHROM": "chr1", "POS": 10, "REF": "A", "ALT": "G",
            "is_snv": True, "caller": "DNA_consensus",
            "filter_original": "Somatic", "filter_normalized": "Somatic",
            "filter_category": "Somatic", "quality": 50.0,
            "genotype": {"GT": "0/1", "DP": 10, "AD": "8,2", "VAF": 0.2},
            "normal_genotype": None, "source_evidence": {"GT_BY_CALLER": "DNA_mutect2:0/1"},
            "id": None, "classification": "Somatic",
        }
    }
    out = aggregate_variants([("DNA_consensus", variants, "DNA")], 1, 1, min_alt_support=0)
    assert out["chr1:10:A:G"]["source_evidence"]["GT_BY_CALLER"] == "DNA_mutect2:0/1"
