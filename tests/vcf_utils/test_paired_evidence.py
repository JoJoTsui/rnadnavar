from vcf_utils.aggregation import (
    aggregate_genotypes,
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
