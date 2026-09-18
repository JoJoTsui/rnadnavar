#!/usr/bin/env python3
"""Archive completed candidate validations; never turn them into training approval."""
import argparse
import json
from pathlib import Path
import shutil

from validate_refined_native_integration import digest

ROOT=Path(__file__).resolve().parents[3]


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--root',type=Path,required=True)
    ap.add_argument('--outdir',type=Path,required=True)
    args=ap.parse_args()
    reports={}
    for dataset in ('seqc2_wes_ll','seqc2_wgs_il','hg008_wgs'):
        report=json.loads((args.root/dataset/'validation.json').read_text())
        if (report['status']!='complete_candidate_validation_not_training_approved'
                or not report['sources_unchanged'] or not report['code_unchanged']):
            raise ValueError('Incomplete/failed validation: '+dataset)
        if set(report['stages']) != {'consensus','first','realignment'}:
            raise ValueError('Incomplete stage coverage')
        for stage in report['stages'].values():
            if (not stage['somatic_parity'] or stage['adapter']['somatic_membership_mismatches']
                    or stage['structural']['issues'] or not stage['structural']['sources_unchanged']
                    or stage['negative_evidence']['status']!='complete_not_training_approved'):
                raise ValueError('Candidate validation failed')
        collision=json.loads((args.root/dataset/'negative_collision_check.json').read_text())
        if collision['status']!='pass_targeted_screen_not_accuracy' or collision['supported_known_somatic_collisions']:
            raise ValueError('Known Somatic collision acquired negative evidence support')
        reports[dataset]=report
    args.outdir.mkdir(parents=True,exist_ok=False)
    summary=dict(policy='separated_three_class_v2',candidate_execution_checks_pass=True,
                 biological_training_approved=False,cohort_executed=False,
                 heavy_root=str(args.root.resolve()),datasets={},archived_sha256={},validated_code={})
    for dataset,report in reports.items():
        code={str(Path(p).relative_to(ROOT)):sha for p,sha in report['code'].items()}
        native=json.loads(Path(report['native_validation']).read_text())
        code.update({str(Path(p).relative_to(ROOT)):sha for p,sha in native['code'].items()})
        loss_path=Path(report['native_validation']).parent/('rescue_loss_audit_v3' if dataset=='seqc2_wes_ll' else 'rescue_loss_audit')/'audit.json'
        loss=json.loads(loss_path.read_text())
        if (loss['status']!='complete_read_only_not_training_approved'
                or not loss['sources_unchanged'] or not loss['code_unchanged']):
            raise ValueError('Incomplete loss attribution')
        # The replay independently bound the rescue implementation to the old
        # result. Do not silently stamp the current source as previously tested.
        for p,sha in loss['code'].items():
            if '/bin/' in p:
                code[str(Path(p).relative_to(ROOT))]=sha
        for rel in ('bin/apply_refined_rescue.py', 'examples/seqc2/scripts/audit_refined_label_contract.py'):
            code.setdefault(rel,digest(ROOT/rel))
        if summary['validated_code'] and summary['validated_code']!=code:
            raise ValueError('Datasets used different code')
        if any(digest(ROOT/p)!=sha for p,sha in code.items()):
            raise ValueError('Code changed since validation')
        summary['validated_code']=code
        summary['datasets'][dataset]=dict(truth=report['truth'],stages={})
        for stage,result in report['stages'].items():
            summary['datasets'][dataset]['stages'][stage]=dict(
                class_counts=result['class_counts'],somatic_parity=result['somatic_parity'],
                somatic={domain:result['metrics']['Somatic/'+domain] for domain in ('ukb','medexome')},
                negative_evidence=result['negative_evidence']['counts'])
        names=['validation.json','negative_collision_check.json']
        native_dest=args.outdir/'evidence'/dataset/'native_validation.json'
        native_dest.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(report['native_validation'],native_dest)
        summary['archived_sha256'][str(native_dest.relative_to(args.outdir))]=digest(native_dest)
        loss_dest=native_dest.with_name('loss_attribution.json')
        shutil.copyfile(loss_path,loss_dest)
        summary['archived_sha256'][str(loss_dest.relative_to(args.outdir))]=digest(loss_dest)
        for stage in ('consensus','first','realignment'):
            names += [stage+'/report.json',stage+'/structural.audit.json',stage+'/bam_pilot.json',stage+'/evidence_gate/report.json']
            if (args.root/dataset/stage/'normal_gvcf.json').exists():
                names.append(stage+'/normal_gvcf.json')
        for name in names:
            src=args.root/dataset/name
            dest=args.outdir/'evidence'/dataset/name
            dest.parent.mkdir(parents=True,exist_ok=True)
            shutil.copyfile(src,dest)
            summary['archived_sha256'][str(dest.relative_to(args.outdir))]=digest(dest)
    (args.outdir/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(args.outdir/'summary.json')


if __name__=='__main__':main()
