#!/usr/bin/env python3
"""Archive lightweight validation evidence; never copy BAMs, VCFs or Parquet."""
import argparse
import json
from pathlib import Path
import shutil
import subprocess

from validate_refined_native_integration import digest

DATASETS = ('seqc2_wes_ll', 'seqc2_wgs_il', 'hg008_wgs')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--root', type=Path, required=True)
    ap.add_argument('--outdir', type=Path, required=True)
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=False)
    summary = {'status':'not_training_approved', 'cohort_rerun_started':False,
               'scope':'Initial full screens are pre-indel-fallback-fix; corrected CLI replay is a separate regression check, not a full post-fix benchmark.',
               'current_commit':subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
               'heavy_root':str(args.root.resolve()), 'datasets':{}, 'archive_sha256':{}}
    for dataset in DATASETS:
        source = args.root / dataset
        target = args.outdir / 'evidence' / dataset
        target.mkdir(parents=True)
        validation = json.loads((source/'validation.json').read_text())
        if validation['status'] != 'complete_somatic_truth_screen_not_training_approved':
            raise ValueError('Incomplete dataset: '+dataset)
        replay = json.loads((source/'fallback_fix_replay/replay.json').read_text())
        pilot = json.loads((source/'bam_pilot.json').read_text())
        cell = {'pre_fix_class_counts':validation['class_counts'],
                'pre_fix_somatic_metrics':{}, 'somatic_truth_collision_counts':{},
                'bam_pilot':pilot['outcomes'],
                'fallback_replay':{k:replay[k] for k in ('restored_somatic','expected_losses','missing','sources_unchanged')}}
        for key, metric in validation['metrics'].items():
            label, domain = key.split('/')
            if label == 'Somatic':
                cell['pre_fix_somatic_metrics'][domain] = metric['values']
            else:
                cell['somatic_truth_collision_counts'][key] = {kind:m['tp'] for kind,m in metric['values'].items()}
        if dataset == 'hg008_wgs':
            alternate = json.loads((source/'recommended_truth/validation.json').read_text())
            if alternate['status'] != 'complete_not_training_approved':
                raise ValueError('Incomplete alternate truth validation')
            cell['recommended_truth_metrics_pre_fix'] = alternate['metrics']
            cell['orthogonal_normal'] = json.loads((source/'orthogonal_normal.json').read_text())['counts']
        summary['datasets'][dataset] = cell
        paths = ['validation.json','transitions.json','bam_pilot.json','fallback_fix_replay/replay.json']
        if dataset == 'hg008_wgs':
            paths += ['recommended_truth/validation.json','orthogonal_normal.json']
        for relative in paths:
            src, dest = source/relative, target/relative
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(src, dest)
            summary['archive_sha256'][str(dest.relative_to(args.outdir))] = digest(dest)
    (args.outdir/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(args.outdir/'summary.json')


if __name__ == '__main__':
    main()
