from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'examples/seqc2/scripts'))
from audit_frozen_baseline_conflicts import af_comparison, partition_label


def test_af_requires_exact_matching_observation_and_valid_scalar():
    assert af_comparison({'GNOMAD_AF': '.1'}, [{'af': .1}])['af_status'] == 'exact_allele_confirmed'
    assert af_comparison({'GNOMAD_AF': '.1'}, [])['af_status'] == 'not_confirmed'
    assert af_comparison({'gnomad_af': '.'}, [])['af_status'] == 'unavailable'
    for bad in ('nan', '-1', '1.2', '.1,.2'):
        assert af_comparison({'GNOMAD_AF': bad}, [])['af_status'] == 'invalid_or_non_scalar'


def test_unmatched_is_not_false_positive():
    key = ('chr1', 1, 'A', 'C')
    assert partition_label(key, {'FP': set()}) == ['unscored_or_representation_unresolved']
    assert partition_label(key, {'FP': {key}, 'TP': set()}) == ['FP']
