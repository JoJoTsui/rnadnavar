"""Unavailable measurements must not become support or numeric summaries."""
import pytest
from vcf_utils.aggregation import (
    _counts_toward_support, aggregate_genotypes, sanitize_genotype,
    tumor_alt_count_from_genotype,
)


@pytest.mark.parametrize('value', [float('nan'), float('inf'), -1, 1.5, '.', None])
def test_invalid_or_unavailable_vaf_is_not_aggregated(value):
    result = aggregate_genotypes({'caller': {'VAF': value}}, ['caller'])
    assert result['vaf_values'] == []
    assert result['vaf_mean'] is None


@pytest.mark.parametrize('value', [-2147483648, -1, float('inf'), 2.5, '.', None])
def test_invalid_depth_is_not_aggregated(value):
    result = aggregate_genotypes({'caller': {'DP': value}}, ['caller'])
    assert result['dp_values'] == []


def test_measured_zero_remains_distinct_from_unavailable():
    result = sanitize_genotype({'DP': 0, 'AD': '0,0', 'VAF': 0})
    assert result['DP'] == result['VAF'] == 0
    assert result['AD'] == '0,0'
    assert result['INVALID_FIELDS'] == []


@pytest.mark.parametrize('ad', [None, '.', '10,.', '10,-2147483648', '10,-1', '10', '10,2'])
def test_positive_alt_floor_requires_valid_sufficient_ad(ad):
    assert not _counts_toward_support({'classification': 'Somatic', 'genotype': {'AD': ad}}, 3)


def test_explicit_disabled_floor_allows_unavailable_ad():
    assert _counts_toward_support({'classification': 'Somatic', 'genotype': {}}, 0)
    assert not _counts_toward_support({'classification': 'Artifact', 'genotype': {}}, 0)


def test_sanitization_retains_invalid_field_provenance():
    result = sanitize_genotype({'DP': -1, 'AD': '1,-2', 'VAF': float('nan')})
    assert result['INVALID_FIELDS'] == ['AD', 'DP', 'VAF']
    assert all(result[k] is None for k in ('AD', 'DP', 'VAF'))


def test_missing_allele_index_is_not_another_alleles_support():
    assert tumor_alt_count_from_genotype({'AD': '10,4', 'ALT_INDICES': [2]}) is None
