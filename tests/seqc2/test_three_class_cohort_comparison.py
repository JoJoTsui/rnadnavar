import importlib.util
import json
from pathlib import Path
import subprocess

import pytest

ROOT=Path(__file__).resolve().parents[2]
spec=importlib.util.spec_from_file_location('cohort_comparison',ROOT/'examples/seq2neo/scripts/compare_three_class_somatic_sets.py')
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def test_query_retains_full_alt_and_rejects_duplicates(monkeypatch):
    row=b'chr1\t1\tA\tG,T\n'
    def query(cmd, **kwargs):
        assert cmd[2:4]==['-i','FILTER="Somatic"']
        assert kwargs['check'] is True
        return subprocess.CompletedProcess(cmd,0,stdout=row)
    monkeypatch.setattr(m.subprocess,'run',query)
    assert m.alleles(Path('a.vcf.gz'))=={row}
    monkeypatch.setattr(m.subprocess,'run',lambda *a,**kw:subprocess.CompletedProcess(a[0],0,stdout=row+row))
    with pytest.raises(ValueError,match='Duplicate'):m.alleles(Path('a.vcf.gz'))


@pytest.fixture
def cohort(tmp_path,monkeypatch):
    old,new=tmp_path/'old',tmp_path/'new'
    manifest=tmp_path/'manifest.tsv'
    ids=[f'sample{i}' for i in range(66)]
    manifest.write_text('sample_id\n'+'\n'.join(ids)+'\n')
    for root,policy in ((old,'seqc2_refined_v2+seqc2_refined_gate_v1'),(new,'separated_three_class_v2')):
        for sid in ids:
            folder=root/sid/'attempt001';folder.mkdir(parents=True)
            state=dict(status='candidate_complete_not_training_approved',identity={'policy':policy},
                       output=str(folder),outputs={},final_artifacts={})
            for stage,prior in (('consensus','refined.vcf.gz'),('rescue','refined.rescue.vcf.gz')):
                path=folder/(prior if root==old else f'three_class.{stage}.vcf.gz')
                path.write_bytes(b'chr1\t1\tA\tG,T\n')
                state['outputs'][str(path)]='declared_hash'
                state['final_artifacts'][stage]=str(path)
            (root/sid/'state.json').write_text(json.dumps(state))
    monkeypatch.setattr(m,'alleles',lambda p:set(p.read_bytes().splitlines(keepends=True)))
    return old,new,manifest,tmp_path/'reports'/'parity.json'


def test_all_66_exact_parity_and_no_overwrite(cohort):
    result=m.compare(*cohort)
    assert result['status']=='exact_somatic_membership_parity_pass'
    assert len(result['samples'])==66 and not result['training_approved']
    assert all(v['added']==v['removed']==0 for s in result['samples'] for v in s['stages'].values())
    with pytest.raises(ValueError,match='fresh'):m.compare(*cohort)


def test_equal_counts_but_different_alleles_fail(cohort):
    path=cohort[1]/'sample0'/'attempt001'/'three_class.rescue.vcf.gz'
    path.write_bytes(b'chr1\t2\tA\tG,T\n')
    with pytest.raises(ValueError,match='membership changed'):m.compare(*cohort)
    result=json.loads(cohort[3].read_text())
    assert result['status']=='failed'
    outcome=result['samples'][0]['stages']['rescue']
    assert outcome['added']==outcome['removed']==1


def test_wrong_policy_fails(cohort):
    path=cohort[1]/'sample0'/'state.json'
    state=json.loads(path.read_text());state['identity']['policy']='obsolete'
    path.write_text(json.dumps(state))
    with pytest.raises(ValueError,match='Wrong new'):m.compare(*cohort)
