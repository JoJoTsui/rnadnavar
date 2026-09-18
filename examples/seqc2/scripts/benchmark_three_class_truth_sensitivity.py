#!/usr/bin/env python3
"""Rebenchmark unchanged queries and baselines against an explicit alternate truth.

No consensus rerun or threshold fitting. Intended to distinguish historical
HG008 multi-sample truth from NIST's recommended tumorvariants representation.
"""
import argparse
import json
from pathlib import Path
import subprocess

from aggregate_benchmark import parse_metrics_json
from validate_frozen_hybrid_policy import pass_query
from validate_refined_native_integration import digest


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--validation', type=Path, required=True)
    ap.add_argument('--truth', type=Path, required=True)
    ap.add_argument('--outdir', type=Path, required=True)
    args = ap.parse_args()
    current = json.loads(args.validation.read_text())
    if current['status'] != 'complete_somatic_truth_screen_not_training_approved':
        raise ValueError('Require completed primary screen')
    frozen = json.loads(Path(current['frozen_validation']).read_text())
    config = current['manifest']
    args.outdir.mkdir(parents=True, exist_ok=False)
    queries = {label: args.validation.parent / f'{label}.query.vcf.gz'
               for label in ('Somatic', 'Germline', 'Reference')}
    sources = {str(args.truth): digest(args.truth), str(args.validation):digest(args.validation)}
    for name, biological in (('refined_consensus', True), ('deepsomatic', False)):
        source = Path(frozen['outputs'][name]['vcf'])
        sources[str(source)] = digest(source)
        queries[name] = args.outdir / f'{name}.query.vcf.gz'
        pass_query(source, queries[name], biological)
    sources.update({str(path): digest(path) for path in queries.values()})
    for path in [config['fasta'], config['hc'], *config['targets'].values()]:
        sources[path] = digest(Path(path))
    report = {'scope':__doc__, 'truth':str(args.truth.resolve()), 'sources':sources,
              'commands':[], 'metrics':{}, 'training_approved':False, 'status':'running'}

    def save():
        (args.outdir / 'validation.json').write_text(json.dumps(report, indent=2)+'\n')

    save()
    for name, query in queries.items():
        for domain, target in config['targets'].items():
            prefix = args.outdir / f'{name}.{domain}'
            cmd = ['micromamba','run','-n','happy','som.py',str(args.truth),str(query),
                   '-R',config['hc'],'-T',target,'-r',config['fasta'],'-N','-o',str(prefix)]
            report['commands'].append(cmd)
            save()
            with Path(str(prefix)+'.log').open('x') as log:
                subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=True)
            metrics = parse_metrics_json(Path(str(prefix)+'.metrics.json'))
            report['metrics'][f'{name}/{domain}'] = metrics
            save()
            print(name, domain, metrics['records'], flush=True)
    report['sources_unchanged'] = all(digest(Path(p)) == sha for p,sha in sources.items())
    report['status'] = 'complete_not_training_approved' if report['sources_unchanged'] else 'integrity_failed'
    save()


if __name__ == '__main__':
    main()
