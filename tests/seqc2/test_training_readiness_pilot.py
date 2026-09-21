import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location('readiness', ROOT/'examples/seq2neo/scripts/evaluate_three_class_training_readiness.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_selection_bounded_order_independent():
    sites = [('chr1', i, 'A', 'T') for i in range(1, 101)]
    forward, reverse = [], []
    for site in sites:
        module.choose(forward, site, 8)
    for site in reversed(sites):
        module.choose(reverse, site, 8)
    assert len(forward) == 8
    assert sorted(forward) == sorted(reverse)


def test_reference_detection_limit_is_not_old_pilot_threshold():
    def counts(n):
        return dict(depth=n, ref=n, alt=0, other=0)
    assert module.assess('Reference', 'A', 'T', counts(60), counts(60))[0] == 'WITHHELD'
    assert module.assess('Reference', 'A', 'T', counts(298), counts(299))[0] == 'WITHHELD'
    assert module.assess('Reference', 'A', 'T', counts(299), counts(299))[0] == 'SUPPORTED'


def test_missing_evidence_and_indels_abstain():
    assert module.assess('Reference', 'A', 'T', {}, {})[0] == 'WITHHELD'
    assert module.assess('Germline', 'A', 'AT', {}, {})[1] == 'requires_haplotype_aware_validation'
