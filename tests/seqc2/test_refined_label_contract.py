import importlib.util
from pathlib import Path
import sys
from urllib.parse import quote

SCRIPTS = Path(__file__).resolve().parents[2] / 'examples/seqc2/scripts'
sys.path.insert(0, str(SCRIPTS))
spec = importlib.util.spec_from_file_location('label_contract', SCRIPTS / 'audit_refined_label_contract.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_flags_and_snapshot(tmp_path):
    path = tmp_path / 'test.vcf'
    original = 'chr1\t1\t.\tA\tC\t.\tSomatic\t.'
    fields = ('UNIFIED_FILTER=Somatic;CLASSIFICATION_RATIONALE=rule:test|class:Somatic;'
              'GATE_POLICY=test;RESCUE_PROMOTED=YES;PASSES_CONSENSUS_DNA=NO;'
              'GATE_RNA_ELIGIBLE=mutect2|strelka;GATE_DNA_NOMINATORS=mutect2;'
              'GATE_SOURCE_RECORD=' + quote(original, safe=''))
    path.write_text(original[:-1] + fields + '\n')
    result = module.audit(path)
    assert result['status'] == 'structural_pass_not_training_approval'
    path.write_text(original[:-1] + fields.replace('UNIFIED_FILTER=Somatic', 'UNIFIED_FILTER=Artifact') + '\n')
    assert module.audit(path)['issues'] == {'unified_filter_mismatch': 1}


def test_missing_metadata_fails(tmp_path):
    path = tmp_path / 'test.vcf'
    path.write_text('chr1\t1\t.\tA\tC\t.\tPASS\t.\n')
    assert set(module.audit(path)['issues']) == {'invalid_label', 'missing_rationale', 'unified_filter_mismatch'}
