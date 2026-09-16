"""Evidence boundaries for the exploratory assay, not production-policy tests."""
import copy
import importlib.util
from pathlib import Path
import sys

SCRIPTS = Path(__file__).resolve().parents[2] / 'examples/seqc2/scripts'
sys.path.insert(0, str(SCRIPTS))
SPEC = importlib.util.spec_from_file_location('indel_assay', SCRIPTS / 'test_corroborated_indel_policy.py')
ASSAY = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(ASSAY)


def evidence():
    ds = {'filter': 'PASS', 'qual': '20', 'samples': {'WES_LL_T_1': {'AD': '7,3'}}}
    m2 = {'filter': 'weak_evidence;contamination', 'info': 'GERMQ=20;TLOD=4',
          'samples': {'WES_LL_T_1': {'AD': '8,2'}, 'WES_LL_N_1': {'AD': '10,0', 'DP': '10'}}}
    return ds, m2


def test_corroboration_accepts_boundary_and_filter_order():
    ds, m2 = evidence()
    assert ASSAY.qualifies(ds, m2)
    m2['filter'] = 'contamination;weak_evidence'
    assert ASSAY.qualifies(ds, m2)


def test_normal_evidence_required_and_never_substitutes_for_tumor():
    ds, m2 = evidence()
    for change in ({'AD': '9,1', 'DP': '10'}, {'AD': '9,0', 'DP': '9'}, {}):
        changed = copy.deepcopy(m2)
        changed['samples']['WES_LL_N_1'] = change
        assert not ASSAY.qualifies(ds, changed)
    m2['samples']['WES_LL_T_1']['AD'] = '9,1'
    assert not ASSAY.qualifies(ds, m2)


def test_hard_rejection_missing_and_nonfinite_confidence_fail():
    ds, m2 = evidence()
    m2['filter'] = 'contamination;germline'
    assert not ASSAY.qualifies(ds, m2)
    ds, m2 = evidence()
    for quality in ('nan', 'inf', '.', '19.9'):
        ds['qual'] = quality
        assert not ASSAY.qualifies(ds, m2)
