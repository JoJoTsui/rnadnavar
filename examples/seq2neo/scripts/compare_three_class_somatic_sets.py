#!/usr/bin/env python3
"""Compare every cohort Somatic allele against the previous candidate rerun.

Read-only inputs; explicit fresh report path. Compares full ALT fields, never
splits multiallelic alleles. This is membership parity, not biological approval.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess


def alleles(path):
    result = subprocess.run(['bcftools', 'query', '-i', 'FILTER="Somatic"',
                             '-f', '%CHROM\t%POS\t%REF\t%ALT\n', str(path)],
                            check=True, capture_output=True)
    rows = result.stdout.splitlines(keepends=True)
    if len(rows) != len(set(rows)):
        raise ValueError('Duplicate Somatic allele: ' + str(path))
    return set(rows)


def digest(rows):
    return hashlib.sha256(b''.join(sorted(rows))).hexdigest()


def compare(previous_root, current_root, manifest, report_path):
    if report_path.exists():
        raise ValueError('Use a fresh report path')
    with manifest.open() as handle:
        samples = [r['sample_id'] for r in csv.DictReader(handle, delimiter='\t')]
    if len(samples) != 66 or len(set(samples)) != 66:
        raise ValueError('Require the exact 66-sample cohort')
    report = dict(status='running', training_approved=False, samples=[],
                  previous_root=str(previous_root), current_root=str(current_root),
                  manifest=str(manifest), manifest_sha256=hashlib.sha256(manifest.read_bytes()).hexdigest(),
                  code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  limitations=['Somatic membership parity is not truth accuracy',
                               'Prior VCF declared hashes are recorded, not freshly rehashed here; current VCF hashes are checked by the separate exporter'])
    report_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        for sid in samples:
            states = [json.loads((root/sid/'state.json').read_text()) for root in (previous_root,current_root)]
            old, new = states
            if any(s['status'] != 'candidate_complete_not_training_approved' for s in states):
                raise ValueError('Incomplete sample: '+sid)
            if new['identity']['policy'] != 'separated_three_class_v2':
                raise ValueError('Wrong new candidate policy')
            if old['identity']['policy'] != 'seqc2_refined_v2+seqc2_refined_gate_v1':
                raise ValueError('Wrong previous baseline policy')
            item = dict(sample_id=sid, stages={})
            for stage, name in (('consensus','refined.vcf.gz'), ('rescue','refined.rescue.vcf.gz')):
                a, b = Path(old['output'])/name, Path(new['final_artifacts'][stage])
                if b != Path(new['output'])/f'three_class.{stage}.vcf.gz':
                    raise ValueError('Wrong final candidate path')
                stamps = [(p.stat().st_size,p.stat().st_mtime_ns) for p in (a,b)]
                before, after = alleles(a), alleles(b)
                if stamps != [(p.stat().st_size,p.stat().st_mtime_ns) for p in (a,b)]:
                    raise ValueError('VCF changed during comparison')
                item['stages'][stage] = dict(previous=str(a),current=str(b),
                    previous_declared_sha256=old['outputs'][str(a)],current_declared_sha256=new['outputs'][str(b)],
                    previous_count=len(before),current_count=len(after),
                    added=len(after-before),removed=len(before-after),
                    previous_allele_sha256=digest(before),current_allele_sha256=digest(after))
            report['samples'].append(item)
            if any(s['added'] or s['removed'] for s in item['stages'].values()):
                raise ValueError('Somatic membership changed: '+sid)
            print(f'[{len(report["samples"])}/66] {sid}: DNA and rescue exact Somatic parity', flush=True)
        if (hashlib.sha256(manifest.read_bytes()).hexdigest()!=report['manifest_sha256']
                or hashlib.sha256(Path(__file__).read_bytes()).hexdigest()!=report['code_sha256']):
            raise ValueError('Manifest/comparison code changed during execution')
        report['status']='exact_somatic_membership_parity_pass'
    except Exception as exc:
        report.update(status='failed',error=str(exc))
        raise
    finally:
        with report_path.open('x') as handle:
            json.dump(report,handle,indent=2)
            handle.write('\n')
    return report


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('previous-root','current-root','manifest','report'):
        parser.add_argument('--'+name,type=Path,required=True)
    args=parser.parse_args()
    compare(args.previous_root,args.current_root,args.manifest,args.report)


if __name__=='__main__':
    main()
