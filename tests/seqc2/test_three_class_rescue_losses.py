"""Characterization tests: document the provisional policy's loss mechanisms."""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'examples/seqc2/scripts'))
from audit_three_class_rescue_losses import explain


def record(label, info='.'):
    return f'chr1\t1\t.\tA\tT\t.\t{label}\t{info}'


def test_legacy_label_alone_can_veto_baseline():
    result = explain(record('Somatic'), record('Reference'))
    assert result['label'] == 'NoConsensus'
    assert result['causes'] == ['inherited_label_Reference']


def test_population_veto_and_verification_are_distinct():
    result = explain(record('Somatic'), record('NoConsensus', 'GNOMAD_AF=0.1'))
    assert result['causes'] == ['common_population_af']
    result = explain(record('Somatic', 'DNA_VERIFICATION=rejected'), record('Somatic'))
    assert result['causes'] == ['verification_rejected']
    assert result['label'] == 'NoConsensus'


def test_missing_rescue_and_low_af_do_not_erase_baseline():
    assert explain(record('Somatic'), None)['label'] == 'Somatic'
    assert explain(record('Somatic'), record('NoConsensus', 'GNOMAD_AF=0.00001'))['label'] == 'Somatic'


def test_overlapping_causes_are_not_double_counted_as_separate_records():
    result = explain(record('Somatic'), record('Germline', 'GNOMAD_AF=0.1'))
    assert result['causes'] == ['common_population_af', 'inherited_label_Germline']
