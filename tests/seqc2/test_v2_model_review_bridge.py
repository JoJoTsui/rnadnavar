import importlib.util
import json
from pathlib import Path
import csv
import pysam
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

ROOT=Path(__file__).resolve().parents[2]
spec=importlib.util.spec_from_file_location('bridge',ROOT/'examples/seq2neo/scripts/build_v2_model_review_bridge.py')
bridge=importlib.util.module_from_spec(spec)
spec.loader.exec_module(bridge)
check_spec=importlib.util.spec_from_file_location('check_bridge',ROOT/'examples/seq2neo/scripts/check_v2_cohort_handoff.py')
check_bridge=importlib.util.module_from_spec(check_spec)
check_spec.loader.exec_module(check_bridge)


@pytest.mark.parametrize('cohort',[False,True])
def test_exact_review_bridge_and_no_overwrite(tmp_path,cohort):
    source=tmp_path/'handoff';source.mkdir()
    header=pysam.VariantHeader();header.contigs.add('chr1',length=100)
    template=tmp_path/'template.vcf.gz'
    with pysam.VariantFile(str(template),'wz',header=header):pass
    ids = {f'SAMPLE_{i}' for i in range(63)} if cohort else bridge.IDS
    stem = 'cohort63' if cohort else 'three_samples'
    bam=tmp_path/'reads.bam';bam.write_bytes(b'fixture');Path(str(bam)+'.bai').write_bytes(b'fixture')
    split=tmp_path/'split.tsv'
    split.write_text('sample_id\tsample_pool\n'+''.join(f'{sid}\t'+('reserved' if sid=='SAMPLE_0' else 'train_pool')+'\n' for sid in sorted(ids)))
    manifest=source/(stem+'.review.tsv')
    with manifest.open('w') as f:
        writer=csv.DictWriter(f,fieldnames=['sample_id','rescue_vcf_path','bam_dt','bam_dn','bam_rt'],delimiter='\t');writer.writeheader()
        for sid in sorted(ids):
            writer.writerow(dict(sample_id=sid,rescue_vcf_path=str(template),bam_dt=str(bam),bam_dn=str(bam),bam_rt=str(bam)))
    rows=[];counts={}
    for sid in sorted(ids):
        for pos,label in enumerate(['Reference','Germline','Somatic'],1):
            rows.append(dict(sample_id=sid,CHROM='chr1',POS=pos,REF='A',ALT='T',FILTER=label,
                             THREE_CLASS_POLICY='separated_three_class_v2',training_eligible=False))
            counts[sid+'|'+label]=1
    parquet=source/(stem+'.review.parquet');pq.write_table(pa.Table.from_pylist(rows),parquet)
    audit=dict(status='review_ready_not_training_approved',training_approved=False,
               outputs={str(p):bridge.sha(p) for p in [manifest,parquet]},three_sample_counts=counts,selected_counts=counts,quarantined_counts={})
    (source/'audit.json').write_text(json.dumps(audit))
    out=tmp_path/'out'
    report=bridge.run(source,out,cohort,split)
    assert report['counts']==counts and report['training_approved'] is False
    model=json.loads((out/'samples.review.json').read_text())
    assert {r['sample_id'] for r in model['samples']}==ids
    if cohort:
        assert report['pool_counts']=={'reserved':1,'train_pool':62}
        assert len(json.loads((out/'reserved.review.json').read_text())['samples'])==1
        digest=bridge.sha(out/'report.json')
        checked=check_bridge.check(out,digest)
        assert checked['total_records']==189 and checked['training_approved'] is False
        with pytest.raises(ValueError,match='Report hash mismatch'):
            check_bridge.check(out,'0'*64)
        bam.write_bytes(b'changed size')
        with pytest.raises(ValueError,match='Alignment path/size changed'):
            check_bridge.check(out,digest)
        bam.write_bytes(b'fixture')
        subset=out/'train_pool.review.json'
        original=subset.read_bytes();subset.write_bytes(original+b' ')
        with pytest.raises(ValueError,match='Output hash mismatch'):
            check_bridge.check(out,digest)
        subset.write_bytes(original)
    for row in model['samples']:
        with pysam.VariantFile(row['vcf']) as vcf:
            assert all(r.info['TRAINING_ELIGIBLE']=='NO' for r in vcf.fetch('chr1'))
    with pytest.raises(FileExistsError):bridge.run(source,out,cohort,split)
    if cohort:
        Path(str(bam)+'.bai').unlink()
        with pytest.raises(ValueError,match='Missing alignment/index'):
            bridge.run(source,tmp_path/'missing',cohort,split)
