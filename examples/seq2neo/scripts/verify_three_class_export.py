#!/usr/bin/env python3
"""Read every exported Parquet row; verify counts, sample coverage and eligibility."""
import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path

import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq


def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda:handle.read(1024*1024),b''):h.update(chunk)
    return h.hexdigest()


def verify(summary_path, manifest_path, report_path):
    if report_path.exists():raise ValueError('Use a fresh verification report')
    summary=json.loads(summary_path.read_text())
    parquet=Path(summary['parquet'])
    with manifest_path.open() as handle:manifest=list(csv.DictReader(handle,delimiter='\t'))
    expected={s['sample_id']:s for s in summary['samples_detail']}
    if (len(expected)!=summary['samples'] or len(manifest)!=len(expected)
            or len({r['sample_id'] for r in manifest})!=len(manifest)
            or {r['sample_id'] for r in manifest}!=set(expected)):
        raise ValueError('Manifest/sample coverage mismatch')
    for row in manifest:
        sample=expected[row['sample_id']]
        if (row['training_label_vcf'] or row['label_qc_verdict']!='NOT_APPROVED'
                or row['candidate_policy']!='separated_three_class_v2'
                or row['rescue_vcf_path']!=sample['new_rescue']
                or row['consensus_vcf_path']!=sample['new_consensus']
                or row['variant_parquet_path']!=str(parquet)):
            raise ValueError('Manifest path/approval mismatch')
    original_sha=sha(parquet)
    if original_sha!=summary['parquet_sha256']:raise ValueError('Parquet checksum mismatch')
    counts=Counter();types=Counter();per_sample={sid:Counter() for sid in expected}
    columns=['sample_id','FILTER','variant_type','training_eligible','THREE_CLASS_POLICY',
             'TRAINING_ELIGIBLE','POS']
    reader=pq.ParquetFile(parquet)
    total=0
    for batch in reader.iter_batches(batch_size=65536,columns=columns,use_threads=False):
        table=pa.Table.from_batches([batch])
        if any(table[name].null_count for name in columns):raise ValueError('Missing class/eligibility data')
        if pc.any(table['training_eligible']).as_py():raise ValueError('Training approval must remain false')
        if (pc.any(pc.not_equal(table['TRAINING_ELIGIBLE'],'NO')).as_py()
                or pc.any(pc.not_equal(table['THREE_CLASS_POLICY'],'separated_three_class_v2')).as_py()
                or pc.any(pc.less_equal(table['POS'],0)).as_py()):
            raise ValueError('Wrong policy/eligibility/coordinates')
        for group in table.group_by(['sample_id','FILTER','variant_type'],use_threads=False).aggregate([('POS','count')]).to_pylist():
            sid,cls,kind,n=(group[k] for k in ('sample_id','FILTER','variant_type','POS_count'))
            if sid not in expected or cls not in {'Somatic','Germline','Reference'}:
                raise ValueError('Unexpected sample/class')
            counts[cls]+=n;types[cls+':'+kind]+=n;per_sample[sid][cls]+=n
        total+=batch.num_rows
    if (total!=summary['candidate_rows'] or dict(counts)!=summary['class_counts']
            or dict(types)!=summary['class_type_counts']
            or any(dict(per_sample[sid])!=expected[sid]['exported_counts'] for sid in expected)):
        raise ValueError('Parquet/sample class-count mismatch')
    if sha(parquet)!=original_sha:raise ValueError('Parquet changed during verification')
    report=dict(status='all_parquet_rows_verified',training_approved=False,samples=len(expected),
                rows=total,class_counts=dict(counts),class_type_counts=dict(types),
                parquet=str(parquet),parquet_sha256=original_sha,
                summary=str(summary_path),summary_sha256=sha(summary_path),
                manifest=str(manifest_path),manifest_sha256=sha(manifest_path),
                verifier_sha256=sha(Path(__file__)))
    report_path.parent.mkdir(parents=True,exist_ok=True)
    with report_path.open('x') as handle:json.dump(report,handle,indent=2);handle.write('\n')
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('summary','manifest','report'):parser.add_argument('--'+name,type=Path,required=True)
    args=parser.parse_args()
    print(json.dumps(verify(args.summary,args.manifest,args.report),indent=2))
