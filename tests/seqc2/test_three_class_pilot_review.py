"""Exact pilot review and cache-preserving execution tuning."""
import importlib.util
import json
from pathlib import Path
import sys

import pysam
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'examples/seq2neo/scripts'))
spec = importlib.util.spec_from_file_location('pilot_review', ROOT/'examples/seq2neo/scripts/review_three_class_pilot.py')
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def vcf(path, rows):
    header = pysam.VariantHeader()
    header.contigs.add('chr1', length=100)
    for name in ('Somatic','Germline','Reference','NoConsensus'):
        header.filters.add(name, None, None, name)
    for name in ('CLASSIFICATION_RATIONALE','TRAINING_ELIGIBLE','UNIFIED_FILTER','THREE_CLASS_POLICY',
                 'THREE_CLASS_BASELINE_FILTER','THREE_CLASS_NATIVE_FILTER'):
        header.info.add(name, 1, 'String', name)
    with pysam.VariantFile(str(path), 'wz', header=header) as writer:
        for pos, cls, info in rows:
            record = header.new_record(contig='chr1', start=pos-1, alleles=('A','G'))
            record.filter.add(cls)
            for k,v in info.items(): record.info[k] = v
            writer.write(record)
    return path


def fixtures(tmp_path, change=None):
    base = vcf(tmp_path/'base.vcf.gz', [(1,'Somatic',{}),(2,'NoConsensus',{})])
    native = vcf(tmp_path/'native.vcf.gz', [(1,'Reference',{'CLASSIFICATION_RATIONALE':'three_class_policy:native_three_class_v1'}),
        (2,'Germline',{'CLASSIFICATION_RATIONALE':'three_class_policy:native_three_class_v1'}),
        (3,'NoConsensus',{})])
    rows = []
    for pos,cls,b,n in ((1,'Somatic','Somatic','Reference'),(2,'Germline','NoConsensus','Germline'),
                        (3,'NoConsensus','MISSING','NoConsensus')):
        info = dict(TRAINING_ELIGIBLE='NO',UNIFIED_FILTER=cls,THREE_CLASS_POLICY=m.cohort.THREE_CLASS_POLICY,
                    THREE_CLASS_BASELINE_FILTER=b,THREE_CLASS_NATIVE_FILTER=n)
        rows.append((pos, cls, info))
    if change == 'lost_somatic': rows[0] = (1,'Reference',rows[0][2])
    if change == 'invented_reference': rows[2] = (3,'Reference',rows[2][2])
    if change == 'missing': rows.pop()
    if change == 'duplicate': rows.append(rows[-1])
    if change == 'training': rows[0][2]['TRAINING_ELIGIBLE'] = 'YES'
    if change == 'provenance': rows[0][2]['THREE_CLASS_NATIVE_FILTER'] = 'MISSING'
    final = vcf(tmp_path/'final.vcf.gz', rows)
    return base, native, final


def test_exact_review_preserves_somatic_and_native_negatives(tmp_path):
    out = m.review_stage(*fixtures(tmp_path), tmp_path/'review.sqlite')
    assert out['counts'] == {'Somatic':1,'Germline':1,'NoConsensus':1}
    assert out['records'] == 3


@pytest.mark.parametrize('change', ['lost_somatic','invented_reference','missing','duplicate','training','provenance'])
def test_review_rejects_class_and_provenance_errors(tmp_path, change):
    with pytest.raises(ValueError):
        m.review_stage(*fixtures(tmp_path, change), tmp_path/'review.sqlite')


def test_execution_tuning_keeps_validated_driver_and_policy_frozen():
    cfg = json.loads((ROOT/'examples/seq2neo/config/separated_three_class_v2_cohort.json').read_text())
    assert cfg['workers'] == 3
    assert cfg['memory_gib_per_worker'] == 16
    assert cfg['reserve_memory_gib'] == 8
    # Preserve the completed pilot identity; scheduling does not enter it.
    assert m.cohort.digest(ROOT/'examples/seq2neo/scripts/run_refined_cohort.py') == '4b3603a195c66d127fcc72345c5b407795f753e9c021cad0c6c448ce8e7bbaa7'
    summary = json.loads((ROOT/cfg['validation_summary']).read_text())
    for rel, sha in summary['validated_code'].items():
        assert m.cohort.digest(ROOT/rel) == sha
    assert (cfg['workers']*cfg['memory_gib_per_worker'] + cfg['reserve_memory_gib']) == 56


def test_resource_check_retains_cgroup_ceiling(tmp_path, monkeypatch):
    exists, read_text = Path.exists, Path.read_text
    files = {'/sys/fs/cgroup/memory.max':str(78*1024**3), '/sys/fs/cgroup/cpu.max':'4600000 100000'}
    monkeypatch.setattr(Path, 'exists', lambda p: str(p) in files if str(p).startswith('/sys/fs/cgroup/') else exists(p))
    monkeypatch.setattr(Path, 'read_text', lambda p, *a, **kw: files[str(p)] if str(p) in files else read_text(p, *a, **kw))
    cfg = dict(workers=3, memory_gib_per_worker=16, reserve_memory_gib=8,
               minimum_free_gib=0, disk_expansion_factor=1, output_root=str(tmp_path), work_root=str(tmp_path))
    assert m.cohort.resource_plan(cfg, 0)['workers'] == 3
    cfg['workers'] = 5
    with pytest.raises(ValueError, match='cgroup memory'):
        m.cohort.resource_plan(cfg, 0)
