#!/usr/bin/env python3
"""CPU-only smoke probe of a named EvoSomatic copy; never builds a cache or trains."""
import argparse
from collections import Counter
import json
from pathlib import Path
import sys

import pysam
from build_v2_model_review_bridge import sha, IDS


def run(model_root, bridge, out):
    if out.exists():
        raise ValueError('Use a fresh probe report')
    report_path=bridge/'report.json'
    verified=json.loads(report_path.read_text())
    if verified['status']!='review_bridge_verified_not_training_approved':
        raise ValueError('Bridge is not verified')
    for p,h in verified['outputs'].items():
        if sha(Path(p))!=h:raise ValueError('Changed bridge artifact')
    files=[model_root/'evosomatic/sample_manifest.py',model_root/'evosomatic/phase1_v2.py',
           model_root/'scripts/build_phase1_v2.py',Path(__file__),report_path]
    hashes={str(p):sha(p) for p in files}
    sys.dont_write_bytecode=True
    sys.path.insert(0,str(model_root))
    from evosomatic.sample_manifest import load_sample_manifest
    from evosomatic.phase1_v2 import parse_labeled_vcf
    samples=load_sample_manifest(str(bridge/'samples.review.json'))
    if {s.sample_id for s in samples}!=IDS or len(samples)!=3:
        raise ValueError('Loader sample identity mismatch')
    result=dict(training_approved=False,model_root=str(model_root),sources=hashes,
                scope='Named local copy only; no claim this is the training owner execution version',samples=[])
    for sample in samples:
        with pysam.VariantFile(sample.vcf) as vcf:
            contigs=set(vcf.header.contigs)
            input_counts=Counter(next(iter(r.filter)) for r in vcf)
        # Reproduce the inspected local split_contigs function; report rather
        # than silently changing its non-autosomal training behavior.
        records=parse_labeled_vcf(sample.vcf,sample.sample_id,
            train_chroms=contigs-{'chr1','chr21','chr22'},valid_chroms={'chr21','chr22'},test_chroms={'chr1'})
        counts=Counter(r.split+'|'+['Reference','Germline','Somatic'][r.label] for r in records)
        other=Counter(r.chrom for r in records if r.split=='train' and r.chrom not in {'chr'+str(i) for i in range(2,21)})
        item=dict(sample_id=sample.sample_id,input_counts=dict(input_counts),parsed=len(records),
                  omitted_by_model_scope=sum(input_counts.values())-len(records),split_class_counts=dict(counts),
                  non_autosomal_train_counts=dict(other),paths={k:getattr(sample,k) for k in ['vcf','normal_dna','tumor_dna','tumor_rna']})
        result['samples'].append(item)
        print(sample.sample_id+': parser probe completed',flush=True)
        del records
    if any(sha(Path(p))!=h for p,h in hashes.items()):raise ValueError('Probe code changed')
    result.update(status='local_parser_probe_complete_not_training_approved',
        blockers=['Confirm actual execution repo/version with model owner',
                  'Approve the weak-label release and the explicit conflict exclusions',
                  'Resolve split policy: local default trains on every contig except chr1/21/22',
                  'Bind a fresh cache namespace to the label release; do not assume old cache parity',
                  'Inspected loader ignores training-approval flags; launch needs an external approval guard'])
    out.parent.mkdir(parents=True,exist_ok=True)
    with out.open('x') as handle:json.dump(result,handle,indent=2)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--model-root',type=Path,required=True)
    p.add_argument('--bridge',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    print(run(a.model_root,a.bridge,a.out)['status'])
