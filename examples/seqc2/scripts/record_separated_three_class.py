#!/usr/bin/env python3
"""Archive completed candidate validations; never turn them into training approval."""
import argparse
import json
from pathlib import Path
import shutil
import sys

from validate_refined_native_integration import digest

ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'bin'))
from assess_negative_label_evidence import valid_counts


def markdown(summary):
    lines=['# Separated three-class v2: completed candidate validation', '',
           '**Execution checks passed; biological training approval did not occur.**', '',
           'All nine stage outputs preserve their declared established Somatic allele sets exactly.',
           'SNP, indel and aggregate metrics match the corresponding baseline in both regions.',
           'Original inputs and validated source code passed integrity checks. No cohort was executed.', '',
           'See [policy](../../SEPARATED_THREE_CLASS_V2.md) and',
           '[loss diagnosis](../three_class_postfix_20260918/RESCUE_LOSS_DIAGNOSIS.md).', '',
           '## Somatic results', '',
           'These restore the established refined policy; they are not a newly tuned improvement over it.',
           'HG008 uses the recommended tumorvariants truth. Comparisons are paired within each domain,',
           'not across different truths or target regions. HG008 has already been inspected and is not',
           'an untouched holdout. Negative-class queries against Somatic truth are collision screens,',
           'not Germline/Reference precision estimates.', '']
    for domain in ('ukb','medexome'):
        for kind in ('snp','indel','records'):
            lines += [f'### {domain}: {kind}', '',
                      '| Dataset | Stage | TP | FP | FN | Precision | Recall | F1 |',
                      '| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |']
            for dataset,result in summary['datasets'].items():
                for stage,row in result['stages'].items():
                    m=row['somatic'][domain][kind]
                    lines.append(f'| {dataset} | {stage} | {m["tp"]} | {m["fp"]} | {m["fn"]} | '
                                 f'{m["precision"]:.6f} | {m["recall"]:.6f} | {m["f1"]:.6f} |')
            lines.append('')
    lines += ['## Provisional paired-read evidence', '',
              '| Dataset | Stage | Germline supported / usable / selected SNPs | Reference supported / usable / selected SNPs |',
              '| --- | --- | ---: | ---: |']
    for dataset,result in summary['datasets'].items():
        for stage,row in result['stages'].items():
            p=row['pilot_yield']
            lines.append(f'| {dataset} | {stage} | {p["Germline"]["supported"]} / {p["Germline"]["usable"]} / {p["Germline"]["selected"]} | '
                         f'{p["Reference"]["supported"]} / {p["Reference"]["usable"]} / {p["Reference"]["selected"]} |')
    lines += ['', 'These are deterministic, truth-blind SNP pilot retention rates, not class accuracy.',
              'Usable means valid, nonzero paired-DNA read counts, not adequate depth for class support.',
              'Stages reuse many sites and must not be pooled as independent samples. Unassessed',
              'records remain withheld; negative indels require haplotype-aware validation.',
              'No supported negative in the targeted exact-allele challenge overlapped known Somatic truth.',
              'This bounded challenge is not a genome-wide guarantee or a haplotype-equivalence test.', '',
              'Reference uses zero ALT/other alleles and at least 299 observations in each DNA sample',
              '(1% detection limit, 95% confidence per sample under the idealized independence model).',
              'Zero supported Reference yield does not justify relaxing this approved threshold.',
              'HG008 N-P normal gVCF corroboration is different-tissue evidence, not Germline truth',
              'or paired tumor/normal Reference approval. Every output remains TRAINING_ELIGIBLE=NO.', '',
              '## Artifacts and reproduction', '',
              f'Heavy outputs, retained outside Git: `{summary["heavy_root"]}`.',
              'Each dataset contains `consensus`, `first` and `realignment` directories:',
              '`candidates.vcf.gz` is the separated candidate output; `evidence_gate/candidate.evidence.vcf.gz`',
              'adds pilot evidence status without changing FILTER. Neither is an approved training VCF.',
              'Lightweight evidence, source/output hashes, commands and all metrics are archived here.',
              'Earlier Sept18 attempts are failed/interrupted development artifacts, not canonical results.', '',
              'Run each validation in a **fresh** destination; existing outputs are never overwritten:', '', '```bash']
    for dataset in summary['datasets']:
        extra=''
        if dataset=='hg008_wgs': extra=' --truth '+summary['datasets'][dataset]['truth']
        lines.append(f'.venv/bin/python examples/seqc2/scripts/validate_separated_three_class.py '
                     f'--native-validation examples/seqc2/comparison/three_class_postfix_20260918/{dataset}/validation.json '
                     f'--samplesheet examples/seqc2/hybrid/csv/{dataset}_hybrid.csv '
                     f'--outdir /path/to/fresh_validation/{dataset}{extra}')
    lines += ['```', '', 'The validation reports also record exact benchmark and evidence-assessment commands.',
              'The cohort wrapper is a separate, candidate-only preparation step; biological approval',
              'and a reviewed cohort pilot remain necessary before any model-training use.', '']
    return '\n'.join(lines)


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
            pilot=json.loads((args.root/dataset/stage/'bam_pilot.json').read_text())
            yields={label:dict(selected=sum(r['class']==label for r in pilot['results']),
                               usable=sum(r['class']==label and valid_counts(r.get('normal'))
                                          and valid_counts(r.get('tumor')) for r in pilot['results']),
                               supported=sum(n for k,n in result['negative_evidence']['counts'].items()
                                             if k.startswith(label+':SUPPORTED:')))
                    for label in ('Germline','Reference')}
            summary['datasets'][dataset]['stages'][stage]=dict(
                class_counts=result['class_counts'],somatic_parity=result['somatic_parity'],
                somatic={domain:result['metrics']['Somatic/'+domain] for domain in ('ukb','medexome')},
                negative_evidence=result['negative_evidence']['counts'],pilot_yield=yields)
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
    (args.outdir/'README.md').write_text(markdown(summary))
    print(args.outdir/'summary.json')


if __name__=='__main__':main()
