#!/usr/bin/env python3
"""Archive lightweight corrected-consensus and negative-evidence checkpoints."""
import argparse
import json
from pathlib import Path
import shutil

from validate_refined_native_integration import digest


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--root',type=Path,required=True)
    ap.add_argument('--outdir',type=Path,required=True)
    args=ap.parse_args()
    args.outdir.mkdir(parents=True,exist_ok=False)
    result={'training_approved':False,'cohort_rerun_started':False,
            'heavy_root':str(args.root.resolve()),'scope':__doc__,'datasets':{},'archived_sha256':{}}
    for dataset in ('seqc2_wes_ll','seqc2_wgs_il','hg008_wgs'):
        source=args.root/dataset
        validation=json.loads((source/'validation.json').read_text())
        if validation['status']!='complete_somatic_truth_screen_not_training_approved':
            raise ValueError('Consensus still incomplete: '+dataset)
        cell={'class_counts':validation['class_counts'],'metrics':validation['metrics'],
              'baseline_somatic_metrics':validation['baseline_somatic_metrics'],
              'sources_unchanged':validation['sources_unchanged'],'code_unchanged':validation['code_unchanged']}
        for name in ('evidence_gate/report.json','recommended_truth/validation.json',
                     'rescue_plan/validation.json','rescue_validation/validation.json'):
            p=source/name
            if p.exists():
                cell[name]=json.loads(p.read_text())
        result['datasets'][dataset]=cell
        names=['validation.json','bam_pilot.json','evidence_gate/report.json','recommended_truth/validation.json',
               'rescue_plan/validation.json','rescue_validation/validation.json']
        for name in names:
            src=source/name
            if not src.exists():continue
            dest=args.outdir/'evidence'/dataset/name
            dest.parent.mkdir(parents=True,exist_ok=True)
            shutil.copyfile(src,dest)
            result['archived_sha256'][str(dest.relative_to(args.outdir))]=digest(dest)
    (args.outdir/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    print(args.outdir/'summary.json')


if __name__=='__main__':main()
