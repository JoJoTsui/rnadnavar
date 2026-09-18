#!/usr/bin/env python3
"""Attribute baseline Somatic losses without changing any labels or thresholds.

Full-query som.py replay must match the recorded consensus metrics before its
partitions can attribute losses. Unmatched representations remain unresolved.
This is not a Germline/Reference accuracy estimate or training approval.
"""
import argparse
from collections import Counter
import json
from pathlib import Path
import sqlite3
import subprocess
import sys
import zlib

import pysam

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'bin'))
from apply_refined_rescue import info_dict, transition
from vcf_utils.refined_rescue_policy import biological_veto
from aggregate_benchmark import parse_metrics_json
from validate_refined_native_integration import alleles, digest


def explain(dna_raw, rescue_raw):
    dna = dna_raw.split('\t')
    if dna[6] != 'Somatic':
        raise ValueError('Expected Somatic DNA baseline')
    rescue = rescue_raw.split('\t') if rescue_raw else None
    di = info_dict(dna[7])
    info = info_dict((rescue or dna)[7])
    verification = info.get('DNA_VERIFICATION')
    if di.get('DNA_VERIFICATION') in {'rejected', 'inconclusive'}:
        verification = di['DNA_VERIFICATION']
    veto = biological_veto(info) or biological_veto(di)
    legacy = rescue[6] if rescue else None
    causes = []
    if veto:
        causes.append(veto)
    if verification in {'rejected', 'inconclusive'}:
        causes.append('verification_' + verification)
    if legacy in {'Germline', 'Reference', 'Artifact', 'RNAedit'}:
        causes.append('inherited_label_' + legacy)
    label, reason = transition('Somatic', legacy, False, '', verification,
                               three_class=True, annotation_veto=veto)
    return dict(label=label, reason=reason, causes=causes,
                source_rescue_label=legacy, verification=verification,
                population_af={k:v for k,v in info.items() if k.casefold() == 'gnomad_af'})


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--validation', required=True, type=Path)
    ap.add_argument('--outdir', required=True, type=Path)
    args = ap.parse_args()
    source = json.loads(args.validation.read_text())
    if source['status'] != 'complete_not_training_approved' or not source['sources_unchanged']:
        raise ValueError('Require completed rescue validation with intact sources')
    consensus_path = Path(source['consensus_validation'])
    consensus = json.loads(consensus_path.read_text())
    # Keep truth domains identical: HG008 sensitivity runs require their own
    # corresponding consensus comparison, not a silent historical substitution.
    if source['truth'] != consensus['manifest']['truth']:
        raise ValueError('Rescue/consensus truth differs; explicit matching baseline required')
    args.outdir.mkdir(parents=True, exist_ok=False)
    tracked = {str(p.resolve()):digest(p) for p in (args.validation, consensus_path)}
    # Include the actual truth, reference, intervals and caller sources from
    # the completed baseline run, not just its metadata file.
    tracked.update(consensus['sources'])
    result = dict(scope=__doc__, status='running', training_approved=False,
                  validation=str(args.validation.resolve()), sources=tracked,
                  rounds={}, benchmarks={}, commands=[],
                  code={str(p.resolve()):digest(p) for p in (
                      Path(__file__), ROOT/'bin/apply_refined_rescue.py',
                      ROOT/'bin/vcf_utils/refined_rescue_policy.py')})

    def save():
        (args.outdir/'audit.json').write_text(json.dumps(result, indent=2)+'\n')

    try:
        for job in source['jobs']:
            rd = job['alignment_round']
            dest = Path(job['outdir'])
            dna = Path(job['command'][job['command'].index('--dna-consensus')+1])
            if digest(dna) != consensus['output_sha256']:
                raise ValueError('DNA source no longer matches validation')
            tracked[str(dna)] = consensus['output_sha256']
            adapter = source['rounds'][rd]['adapter']
            tracked.update(adapter['sources'])
            for p,sha in adapter['code'].items():
                if digest(Path(p)) != sha:
                    raise ValueError('Cannot replay changed adapter code: '+p)
            output = dest/'refined.rescue.vcf.gz'
            if digest(output) != adapter['output_sha256']:
                raise ValueError('Rescue output changed')
            tracked[str(output)] = adapter['output_sha256']
            positive = alleles(dna, True)
            observed = {}
            with pysam.VariantFile(str(output)) as reader:
                for r in reader:
                    k = (r.contig, r.pos, r.ref, ','.join(r.alts or []))
                    if k in positive:
                        observed[k] = ';'.join(r.filter)
            counts, losses = Counter(), []
            dbpath = dest/'union.sqlite'
            with sqlite3.connect(dbpath.resolve().as_uri()+'?mode=ro', uri=True) as db:
                for k in sorted(positive):
                    row = db.execute('SELECT dna,rescue FROM variants WHERE chrom=? AND pos=? AND ref=? AND alt=?', k).fetchone()
                    if row is None or row[0] is None:
                        raise ValueError('Missing baseline in staged union: '+str(k))
                    evidence = explain(zlib.decompress(row[0]).decode(), zlib.decompress(row[1]).decode() if row[1] else None)
                    if observed.get(k) != evidence['label']:
                        raise ValueError('Replay/output disagreement: '+str(k))
                    group = '+'.join(evidence['causes']) or 'retained'
                    counts[group] += 1
                    if evidence['label'] != 'Somatic':
                        losses.append(dict(allele=k, **evidence))
            if len(losses) != adapter['counts'].get('baseline_conflict_requires_review', 0):
                raise ValueError('Loss count differs from adapter report')
            result['rounds'][rd] = dict(causes=dict(counts), losses=losses,
                                        all_baseline_records_accounted_for=len(observed)==len(positive))
        save()
        for region in consensus['manifest']['targets']:
            original = next(c for c in consensus['commands'] if c[-1].endswith('/Somatic.'+region))
            cmd = list(original)
            cmd[-1] = str(args.outdir/region)
            scratch = args.outdir/(region+'_scratch')
            cmd += ['--scratch-prefix', str(scratch)]
            result['commands'].append(cmd)
            with (args.outdir/(region+'.log')).open('x') as log:
                subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=True)
            metrics = parse_metrics_json(args.outdir/(region+'.metrics.json'))
            if metrics != consensus['metrics']['Somatic/'+region]['values']:
                raise ValueError('Baseline metric replay mismatch: '+region)
            partitions = {label:alleles(scratch/path) for label,path in {
                'TP':'tpfn/0003.vcf.gz', 'FP':'fp.vcf.gz',
                'ambiguous':'ambi.vcf.gz', 'unscored':'unk.vcf.gz'}.items()}
            result['benchmarks'][region] = dict(metrics=metrics, baseline_parity=True)
            for rd in result['rounds'].values():
                grouped = Counter()
                for row in rd['losses']:
                    labels = [label for label,keys in partitions.items() if tuple(row['allele']) in keys]
                    labels = labels or ['representation_unresolved']
                    row.setdefault('partitions', {})[region] = labels
                    for label in labels:
                        grouped['+'.join(row['causes'])+':'+label] += 1
                rd.setdefault('regional_loss_causes', {})[region] = dict(grouped)
            save()
        result['sources_unchanged'] = all(digest(Path(p)) == sha for p,sha in tracked.items())
        result['code_unchanged'] = all(digest(Path(p)) == sha for p,sha in result['code'].items())
        if not result['sources_unchanged'] or not result['code_unchanged']:
            raise ValueError('Audit integrity failed')
        result['status'] = 'complete_read_only_not_training_approved'
    except Exception as exc:
        result.update(status='failed', error=str(exc))
        raise
    finally:
        save()
    print(json.dumps({k:v['regional_loss_causes'] for k,v in result['rounds'].items()}, indent=2))


if __name__ == '__main__':
    main()
