import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location('gate', ROOT / 'examples/seq2neo/scripts/validate_native_gated_policy.py')
gate = importlib.util.module_from_spec(spec)
spec.loader.exec_module(gate)


def metrics(p, f1, tp, fp, fn):
    return {'precision': p, 'f1': f1, 'tp': tp, 'fp': fp, 'fn': fn}


def test_release_gate_accepts_candidate_beating_deepsomatic():
    slices = []
    for cohort in ("SEQC2_WES_LL", "SEQC2_WES_IL", "SEQC2_WGS_IL", "HG008_WGS"):
        for variant_type in ("snp", "indel", "records"):
            slices.append({
                'cohort': cohort, 'variant_type': variant_type,
                'candidate': metrics(.968, .804, 570, 19, 259),
                'deepsomatic': metrics(.967, .798, 563, 19, 266),
                'valid_input_failures': 0,
            })
    result = gate.evaluate({'slices': slices})
    assert result['status'] == 'passed'


def test_release_gate_rejects_precision_loss_and_failures():
    result = gate.evaluate({'slices': [{
        'cohort': 'HG008', 'variant_type': 'indel',
        'candidate': metrics(.90, .80, 10, 5, 2),
        'deepsomatic': metrics(.95, .80, 10, 2, 2),
        'valid_input_failures': 1,
    }]})
    assert result['status'] == 'failed'
    assert 'HG008:indel' in result['reasons'][0]


def test_release_gate_rejects_each_dominance_violation():
    base = metrics(.95, .80, 10, 2, 5)
    cases = [
        metrics(.94, .80, 10, 2, 5),  # precision
        metrics(.95, .79, 10, 2, 5),  # f1
        metrics(.95, .80, 9, 2, 5),   # tp
        metrics(.95, .80, 10, 3, 5),  # fp
        metrics(.95, .80, 10, 2, 6),  # fn
    ]
    for candidate in cases:
        result = gate.evaluate({'slices': [{'candidate': candidate, 'deepsomatic': base}]})
        assert result['status'] == 'failed'


def test_release_gate_rejects_missing_or_nonfinite_metrics():
    result = gate.evaluate({'slices': [{'candidate': {'precision': 'nan'}, 'deepsomatic': metrics(.9, .8, 1, 1, 1)}]})
    assert result['status'] == 'failed'
