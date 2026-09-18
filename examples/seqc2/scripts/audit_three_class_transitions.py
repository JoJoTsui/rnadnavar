#!/usr/bin/env python3
"""Record exact-allele collisions and lost Somatic labels, not haplotype accuracy."""
import argparse
import json
from pathlib import Path

import pysam
from validate_refined_native_integration import digest


def key(r):
    return r.contig, r.pos, r.ref, r.alts


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--validation', required=True, type=Path)
    ap.add_argument('--out', required=True, type=Path)
    args = ap.parse_args()
    if args.out.exists():
        raise ValueError('Refuse overwrite')
    current = json.loads(args.validation.read_text())
    frozen = json.loads(Path(current['frozen_validation']).read_text())
    baseline_path = frozen['outputs']['refined_consensus']['vcf']
    with pysam.VariantFile(baseline_path) as v:
        baseline = {key(r) for r in v if set(r.filter) == {'Somatic'}}
    with pysam.VariantFile(current['manifest']['truth']) as v:
        truth = {key(r) for r in v}
    output = args.validation.parent / 'three_class.vcf.gz'
    report = {'scope':__doc__, 'validation':str(args.validation.resolve()),
              'output_sha256':digest(output), 'negative_truth_exact_matches':[], 'somatic_label_losses':[],
              'warning':'Exact REF/ALT only, all input loci; use som.py regional metrics for counts. Not a germline/reference truth benchmark.'}
    with pysam.VariantFile(str(output)) as v:
        for r in v:
            k = key(r)
            label = next(iter(r.filter),'')
            row = {'allele':k,'label':label,'rationale':r.info.get('CLASSIFICATION_RATIONALE'),
                   'exact_somatic_truth_match':k in truth}
            if label in {'Reference','Germline'} and k in truth:
                report['negative_truth_exact_matches'].append(row)
            if label != 'Somatic' and k in baseline:
                report['somatic_label_losses'].append(row)
    with args.out.open('x') as f:
        json.dump(report, f, indent=2)
        f.write('\n')
    print({k:len(report[k]) for k in ('negative_truth_exact_matches','somatic_label_losses')})


if __name__ == '__main__':
    main()
