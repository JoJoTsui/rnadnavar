#!/usr/bin/env python3
"""Export selected v2 review rows to minimal model VCFs, without approving training."""
import argparse
from collections import Counter
from contextlib import ExitStack
import csv
import hashlib
import json
from pathlib import Path

import pyarrow.parquet as pq
import pysam

IDS = {'PRJNA298376_4007','PRJNA298376_4060','PRJNA298376_4072'}


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda:handle.read(1024*1024),b''):
            h.update(chunk)
    return h.hexdigest()


def run(handoff, out, cohort=False, split_manifest=None):
    audit_path = handoff/'audit.json'
    audit = json.loads(audit_path.read_text())
    stem = 'cohort63' if cohort else 'three_samples'
    parquet, manifest = handoff/(stem+'.review.parquet'), handoff/(stem+'.review.tsv')
    if audit['status'] != 'review_ready_not_training_approved' or audit['training_approved'] is not False:
        raise ValueError('Require the unapproved v2 review handoff')
    hashes = {str(p):sha(p) for p in (parquet,manifest,audit_path,Path(__file__))}
    if any(hashes[str(p)] != audit['outputs'][str(p)] for p in (parquet,manifest)):
        raise ValueError('Changed handoff')
    with manifest.open() as handle:
        rows = list(csv.DictReader(handle,delimiter='\t'))
    ids = {r['sample_id'] for r in rows}
    expected = {k.split('|')[0] for k in audit['selected_counts']} if cohort else IDS
    if len(rows) != (63 if cohort else 3) or ids != expected or len(ids)!=len(rows):
        raise ValueError('Require exact full three-sample IDs')
    pools = {}
    alignment_inventory = []
    if cohort:
        excluded = {'PRJNA298330_4032','PRJNA298376_4081','PRJNA298376_4255'}
        if ids & excluded or split_manifest is None:
            raise ValueError('Excluded sample present or missing split manifest')
        hashes[str(split_manifest)] = sha(split_manifest)
        with split_manifest.open() as handle:
            for row in csv.DictReader(handle,delimiter='\t'):
                sid = row['sample_id']
                if sid in pools: raise ValueError('Duplicate split sample')
                pools[sid] = row['sample_pool']
        if not ids <= pools.keys() or any(pools[sid] not in {'train_pool','reserved'} for sid in ids):
            raise ValueError('Missing or unsupported sample split')
        for row in rows:
            for modality in ('dn','dt','rt'):
                path = Path(row['bam_'+modality])
                indexes = [Path(str(path)+'.bai'),path.with_suffix('.bai'),Path(str(path)+'.crai'),path.with_suffix('.crai')]
                index = next((p for p in indexes if p.is_file()),None)
                if not path.is_file() or index is None: raise ValueError('Missing alignment/index: '+str(path))
                alignment_inventory.append(dict(sample_id=row['sample_id'],modality=modality,path=str(path),size=path.stat().st_size,index=str(index)))
    out.mkdir(parents=True,exist_ok=False)
    report = dict(status='running',training_approved=False,sources=hashes,
                  scope='Minimal model-input review VCF; annotation-rich evidence remains in source Parquet',
                  class_order={'Reference':0,'Germline':1,'Somatic':2})
    counts = Counter()
    allele_hashes = {sid:hashlib.sha256() for sid in ids}
    model_rows = []
    try:
        with ExitStack() as stack:
            writers = {}
            for row in rows:
                header = pysam.VariantHeader()
                with pysam.VariantFile(row['rescue_vcf_path']) as source:
                    for name, contig in source.header.contigs.items():
                        header.contigs.add(name,length=contig.length)
                for label in ('Reference','Germline','Somatic'):
                    header.filters.add(label,None,None,'Workflow-derived review class; not approved truth')
                header.info.add('TRAINING_ELIGIBLE',1,'String','NO: review-only artifact')
                header.info.add('THREE_CLASS_POLICY',1,'String','Source policy')
                path = out/(row['sample_id']+'.review.vcf.gz')
                writers[row['sample_id']] = stack.enter_context(pysam.VariantFile(str(path),'wz',header=header))
                model_rows.append(dict(sample_id=row['sample_id'],tumor_dna=row['bam_dt'],
                                       normal_dna=row['bam_dn'],tumor_rna=row['bam_rt'],vcf=str(path)))
                if cohort: model_rows[-1]['sample_pool'] = pools[row['sample_id']]
            for batch in pq.ParquetFile(parquet).iter_batches(batch_size=32768,
                    columns=['sample_id','CHROM','POS','REF','ALT','FILTER','THREE_CLASS_POLICY','training_eligible'],use_threads=False):
                for row in batch.to_pylist():
                    sid,label = row['sample_id'],row['FILTER']
                    if sid not in ids or label not in report['class_order'] or row['training_eligible'] is not False or row['THREE_CLASS_POLICY']!='separated_three_class_v2':
                        raise ValueError('Invalid source row')
                    writer = writers[sid]
                    rec = writer.new_record(contig=row['CHROM'],start=row['POS']-1,
                        alleles=(row['REF'],*row['ALT'].split(',')),filter=label,
                        info={'TRAINING_ELIGIBLE':'NO','THREE_CLASS_POLICY':'separated_three_class_v2'})
                    writer.write(rec)
                    counts[sid+'|'+label]+=1
                    allele_hashes[sid].update(f"{row['CHROM']}\t{row['POS']}\t{row['REF']}\t{row['ALT']}\t{label}\n".encode())
        if dict(counts) != audit['selected_counts' if cohort else 'three_sample_counts']:
            raise ValueError('VCF export class counts disagree')
        for row in model_rows:
            pysam.tabix_index(row['vcf'],preset='vcf')
            h = hashlib.sha256()
            with pysam.VariantFile(row['vcf']) as source:
                for rec in source:
                    h.update(f"{rec.contig}\t{rec.pos}\t{rec.ref}\t{','.join(rec.alts)}\t{';'.join(rec.filter)}\n".encode())
            if h.hexdigest()!=allele_hashes[row['sample_id']].hexdigest():
                raise ValueError('Full allele/class roundtrip mismatch')
        # Extra review status is explicit; the inspected model loader ignores it.
        with (out/'samples.review.json').open('x') as handle:
            json.dump(dict(schema_version=1,training_approved=False,samples=model_rows),handle,indent=2)
        if cohort:
            for pool in ('train_pool','reserved'):
                with (out/(pool+'.review.json')).open('x') as handle:
                    json.dump(dict(schema_version=1,training_approved=False,samples=[r for r in model_rows if r['sample_pool']==pool]),handle,indent=2)
            report.update(alignment_inventory=alignment_inventory,
                          pool_counts=dict(Counter(pools[sid] for sid in ids)),
                          excluded_samples=sorted(excluded),quarantined_counts=audit['quarantined_counts'])
        if any(sha(Path(p))!=value for p,value in hashes.items()):
            raise ValueError('Sources changed during export')
        report.update(status='review_bridge_verified_not_training_approved',counts=dict(counts),
                      allele_class_sha256={sid:h.hexdigest() for sid,h in allele_hashes.items()},
                      outputs={str(p):sha(p) for p in out.iterdir() if p.is_file()})
    except Exception as exc:
        report.update(status='failed',error=str(exc))
        raise
    finally:
        with (out/'report.json').open('x') as handle:
            json.dump(report,handle,indent=2)
    return report


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--handoff',type=Path,required=True)
    p.add_argument('--outdir',type=Path,required=True)
    p.add_argument('--cohort',action='store_true',help='Export all 63 reviewed samples, preserving sample pools')
    p.add_argument('--split-manifest',type=Path)
    a=p.parse_args()
    print(run(a.handoff,a.outdir,a.cohort,a.split_manifest)['status'])
