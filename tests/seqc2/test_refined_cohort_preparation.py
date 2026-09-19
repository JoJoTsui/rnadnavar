"""Isolated cohort preparation and synthetic execution; never uses real samples."""
import csv
import importlib.util
import json
from pathlib import Path
import shutil
import sys
import subprocess

import pysam
import pytest

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location('refined_cohort', ROOT / 'examples/seq2neo/scripts/run_refined_cohort.py')
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


@pytest.fixture
def cohort(tmp_path):
    source = tmp_path / 'sources' / 'vcf_realignment'
    source.mkdir(parents=True)
    fasta = tmp_path / 'reference.fa'
    fasta.write_text('>chr1\n' + 'A' * 100 + '\n')
    pysam.faidx(str(fasta))
    header = ('##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100>\n'
              '##FILTER=<ID=Somatic,Description="Somatic">\n'
              '##INFO=<ID=GNOMAD_AF,Number=1,Type=Float,Description="AF">\n'
              '##INFO=<ID=REDI_CANONICAL,Number=1,Type=String,Description="Editing">\n'
              '##INFO=<ID=TLOD,Number=A,Type=Float,Description="Score">\n'
              '##INFO=<ID=GERMQ,Number=1,Type=Integer,Description="Score">\n'
              '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
              '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="AD">\n'
              '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="DP">\n')
    header += '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Native confidence">\n'
    header += ''.join(f'##FORMAT=<ID={b}U,Number=2,Type=Integer,Description="Counts">\n' for b in 'ACGT')
    columns = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO'
    row = dict(sample_id='sample1', base_output_dir=str(source.parent))
    for mod in ('dna', 'rna'):
        for caller in m.CALLERS:
            plain = source / f'{mod}.{caller}.vcf'
            plain.write_text(header + columns + '\tFORMAT\tTUMOR\n'
                             'chr1\t10\t.\tA\tG\t50\tPASS\tTLOD=20;GERMQ=80\tGT:AD:DP:AU:CU:GU:TU\t0/1:20,10:30:20,20:0,0:10,10:0,0\n')
            pysam.tabix_compress(str(plain), str(plain)+'.gz')
            pysam.tabix_index(str(plain)+'.gz', preset='vcf')
            row[f'caller_{mod}_{caller}'] = str(plain)+'.gz'
    rescue = source / 'rescue.vcf'
    rescue.write_text(header + columns + '\nchr1\t10\t.\tA\tG\t50\tSomatic\tGNOMAD_AF=0;REDI_CANONICAL=NO\n')
    pysam.tabix_compress(str(rescue), str(rescue)+'.gz')
    pysam.tabix_index(str(rescue)+'.gz', preset='vcf')
    row['rescue_vcf_path'] = str(rescue)+'.gz'
    manifest = tmp_path / 'manifest.tsv'
    with manifest.open('w') as out:
        writer = csv.DictWriter(out, fieldnames=list(row), delimiter='\t')
        writer.writeheader(); writer.writerow(row)
    return dict(policy=m.POLICY, manifest=str(manifest), fasta=str(fasta), expected_samples=1,
                pilot_samples=['sample1'], output_root=str(tmp_path/'output'), work_root=str(tmp_path/'work'),
                workers=1, memory_gib_per_worker=4, reserve_memory_gib=1, minimum_free_gib=0, disk_expansion_factor=1)


def test_prepare_is_read_only_and_requires_opt_in(cohort, monkeypatch):
    plan = m.prepare(cohort)
    assert len(plan['samples']) == 1
    assert not Path(cohort['output_root']).exists()
    cfg = Path(cohort['manifest']).with_suffix('.json')
    cfg.write_text(json.dumps(cohort))
    monkeypatch.setattr(sys, 'argv', ['runner', '--config', str(cfg)])
    monkeypatch.setattr(m, 'run_sample', lambda *a: pytest.fail('Preparation must not generate'))
    m.main()
    assert not Path(cohort['work_root']).exists()


@pytest.mark.parametrize('change,match', [
    ({'expected_samples':2}, 'count'), ({'policy':'old'}, 'policy'),
    ({'pilot_samples':['absent']}, 'Pilot'), ({'workers':0}, 'Positive'),
])
def test_invalid_config_fails(cohort, change, match):
    cohort.update(change)
    with pytest.raises(ValueError, match=match): m.prepare(cohort)


def test_overlap_and_missing_index_fail(cohort):
    source = Path(cohort['manifest']).parent/'sources'
    bad = dict(cohort, output_root=str(source/'new'))
    with pytest.raises(ValueError, match='overlaps'): m.prepare(bad)
    index = next(source.rglob('*.tbi'))
    index.unlink()
    with pytest.raises(ValueError, match='Missing index'): m.prepare(cohort)


def test_reference_and_native_schema_fail(cohort):
    source = next((Path(cohort['manifest']).parent/'sources').rglob('dna.mutect2.vcf.gz'))
    with pytest.raises(ValueError, match='dictionary'): m.inspect_vcf(source, {'chr1':99}, 'mutect2')
    rescue = source.parent/'rescue.vcf.gz'
    with pytest.raises(ValueError, match='sample count'): m.inspect_vcf(rescue, {'chr1':100}, 'mutect2')


@pytest.mark.skipif(not shutil.which('bcftools') or not shutil.which('prlimit'), reason='Requires bcftools and prlimit')
def test_synthetic_execution_integrity_resume_and_manifest(cohort):
    plan = m.prepare(cohort)
    sample = plan['samples'][0]
    identity = dict(policy=m.POLICY, manifest=plan['manifest_sha256'], reference=m.digest(cohort['fasta']),
                    reference_fai=m.digest(cohort['fasta']+'.fai'), code=m.code_hashes())
    sid, status = m.run_sample(sample, cohort, identity)
    assert status == 'candidate_complete_not_training_approved'
    root = Path(cohort['output_root'])
    state = json.loads((root/sid/'state.json').read_text())
    assert all(m.digest(p) == h for p,h in state['source_hashes'].items())
    assert not any('nextflow' in cmd for cmd in state['commands'])
    assert m.run_sample(sample, cohort, identity)[1] == 'cached_candidate_complete'
    m.write_candidate_manifest(root, [sample])
    assert 'refined.rescue.vcf.gz' in (root/'candidate_manifest.tsv').read_text()
    Path(next(iter(state['outputs']))).write_bytes(b'corrupt synthetic output')
    assert not m.completed(root/sid, identity, state['source_hashes'])


@pytest.fixture
def three_class_cohort(cohort):
    folder=Path(cohort['manifest']).parent
    evidence=folder/'synthetic_validation_evidence.json'
    evidence.write_text('{}\n')
    # Synthetic gate is for the one-record test only; never authorizes real data.
    gate=dict(policy=m.THREE_CLASS_POLICY,candidate_execution_checks_pass=True,
              biological_training_approved=False,
              validated_code={p:m.digest(m.REPO/p) for p in m.required_validation_code()},
              archived_sha256={evidence.name:m.digest(evidence)},datasets={})
    for ds in ('seqc2_wes_ll','seqc2_wgs_il','hg008_wgs'):
        gate['datasets'][ds]=dict(stages={s:{'somatic_parity':True} for s in ('consensus','first','realignment')},
            negative_collision_coverage={'status':'rescue_negatives_subset_of_challenged_consensus_negatives'})
    path=folder/'synthetic_validation.json'
    path.write_text(json.dumps(gate))
    cohort.update(policy=m.THREE_CLASS_POLICY,candidate_only=True,validation_summary=str(path))
    return cohort


@pytest.mark.parametrize('problem', ['approval','code','evidence','coverage','opt_in'])
def test_three_class_gate_fails_closed(three_class_cohort, problem):
    cfg=three_class_cohort
    path=Path(cfg['validation_summary']); report=json.loads(path.read_text())
    if problem=='approval':report['biological_training_approved']=True
    if problem=='code':report['validated_code']['bin/apply_three_class_labels.py']='0'*64
    if problem=='evidence':report['archived_sha256']['synthetic_validation_evidence.json']='0'*64
    if problem=='coverage':report['datasets']['hg008_wgs']['stages']['realignment']['somatic_parity']=False
    if problem=='opt_in':cfg['candidate_only']=False
    path.write_text(json.dumps(report))
    with pytest.raises(ValueError):m.prepare(cfg)
    assert not Path(cfg['output_root']).exists()


@pytest.mark.skipif(not shutil.which('bcftools') or not shutil.which('prlimit'), reason='Requires bcftools and prlimit')
def test_three_class_synthetic_execution_and_no_baseline_fallback(three_class_cohort):
    cfg=three_class_cohort
    plan=m.prepare(cfg)
    identity=dict(policy=cfg['policy'],manifest=plan['manifest_sha256'],reference=m.digest(cfg['fasta']),
                  reference_fai=m.digest(cfg['fasta']+'.fai'),code=m.code_hashes(),
                  validation_gate=plan['validation_gate']['sha256'])
    sample=plan['samples'][0]
    sid,status=m.run_sample(sample,cfg,identity)
    assert status=='candidate_complete_not_training_approved'
    root=Path(cfg['output_root']); state_path=root/sid/'state.json'
    state=json.loads(state_path.read_text())
    assert len([c for c in state['commands'] if '--experimental-three-class' in c])==1
    assert not any('--experimental-three-class' in c for c in state['commands'] if any('apply_refined_rescue.py' in a for a in c))
    assert not any('nextflow' in c for c in state['commands'])
    for p in state['final_artifacts'].values():
        with pysam.VariantFile(p) as reader:
            rows=list(reader)
            assert len(rows)==1 and set(rows[0].filter)=={'Somatic'}
            assert rows[0].info['THREE_CLASS_POLICY']==m.THREE_CLASS_POLICY
            assert rows[0].info['TRAINING_ELIGIBLE']=='NO'
    # Scheduling-only tuning must reuse byte-verified pilot results rather than
    # recompute native nominations or change the validated policy identity.
    cfg['workers'] = 3
    assert m.run_sample(sample,cfg,identity)[1]=='cached_candidate_complete'
    m.write_candidate_manifest(root,[sample])
    with (root/'candidate_manifest.tsv').open() as handle: row=next(csv.DictReader(handle,delimiter='\t'))
    assert row['truth_vcf']==row['training_label_vcf']==''
    assert row['candidate_vcf'].endswith('three_class.rescue.vcf.gz')
    state['final_artifacts']['rescue']=str(Path(state['output'])/'refined.rescue.vcf.gz')
    state_path.write_text(json.dumps(state))
    assert not m.completed(root/sid,identity,state['source_hashes'])
    with pytest.raises(ValueError,match='Invalid separated'):m.write_candidate_manifest(root,[sample])


@pytest.mark.parametrize('args', [['--execute'],['--pilot','--approve-pilot'],['--prepare','--pilot']])
def test_three_class_wrapper_rejects_unsafe_phase_arguments(args):
    result=subprocess.run(['bash',str(ROOT/'examples/seq2neo/run_three_class_cohort.sh'),*args],
                          capture_output=True,text=True)
    assert result.returncode==2
    assert 'Log:' not in result.stdout  # Reject before creating shared logs/output.


def test_three_class_config_uses_new_shared_namespaces():
    cfg=json.loads((ROOT/'examples/seq2neo/config/separated_three_class_v2_cohort.json').read_text())
    old=json.loads((ROOT/'examples/seq2neo/config/refined_native_v2_cohort.json').read_text())
    shared='/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/'
    assert cfg['output_root'].startswith(shared+'rnadnavar/examples/seq2neo/')
    assert cfg['work_root'].startswith(shared+'nf_work/')
    assert cfg['output_root']!=old['output_root'] and cfg['work_root']!=old['work_root']
    assert cfg['manifest']==old['manifest'] and cfg['fasta']==old['fasta']
    assert cfg['expected_samples']==66 and cfg['candidate_only'] is True
    assert not Path(cfg['validation_summary']).is_absolute()
