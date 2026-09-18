"""The archived execution gate must never mistake partial results for approval."""
import importlib.util
import json
from pathlib import Path
import sys

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'examples/seqc2/scripts'))
spec = importlib.util.spec_from_file_location('record_separated', ROOT / 'examples/seqc2/scripts/record_separated_three_class.py')
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


@pytest.mark.parametrize('failure', ['running', 'structural', 'parity', 'collision'])
def test_fail_closed_before_creating_archive(tmp_path, monkeypatch, failure):
    source = tmp_path / 'inputs' / 'seqc2_wes_ll'
    source.mkdir(parents=True)
    stage = dict(somatic_parity=True, adapter={'somatic_membership_mismatches': 0},
                 structural={'issues': {}, 'sources_unchanged': True},
                 negative_evidence={'status': 'complete_not_training_approved'})
    report = dict(status='complete_candidate_validation_not_training_approved',
                  sources_unchanged=True, code_unchanged=True,
                  stages={s: stage for s in ('consensus', 'first', 'realignment')})
    collision = dict(status='pass_targeted_screen_not_accuracy', supported_known_somatic_collisions=0)
    if failure == 'running': report['status'] = 'running'
    if failure == 'structural': stage['structural']['issues'] = {'contradiction': 1}
    if failure == 'parity': stage['somatic_parity'] = False
    if failure == 'collision': collision['supported_known_somatic_collisions'] = 1
    (source / 'validation.json').write_text(json.dumps(report))
    (source / 'negative_collision_check.json').write_text(json.dumps(collision))
    output = tmp_path / 'archive'
    monkeypatch.setattr(sys, 'argv', ['record', '--root', str(source.parent), '--outdir', str(output)])
    with pytest.raises(ValueError): m.main()
    assert not output.exists()
