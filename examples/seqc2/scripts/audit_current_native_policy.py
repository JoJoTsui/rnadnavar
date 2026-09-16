#!/usr/bin/env python3
"""Reproduce current native consensus on verified SEQC2 sources, VCF-only.

Uses unrestricted caller candidates; HC and target intervals apply at scoring.
Creates a fresh output directory and never launches Nextflow or reads HG008.
"""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import sys

from aggregate_benchmark import parse_metrics_json
from replay_historical_native_gate import write_query


def digest(path):
    value = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            value.update(block)
    return value.hexdigest()


def alleles(path, accepted):
    result = set()
    with gzip.open(path, 'rt') as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            row = line.rstrip().split('\t')
            if row[6] in accepted:
                for alt in row[4].split(','):
                    result.add((row[0], int(row[1]), row[3], alt))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir', required=True, type=Path)
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[3]
    root = repo / 'examples/seqc2'
    reference = Path('/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/bio_db/references/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta')
    hc = Path('/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/data/giab/data/seqc2/truth/High-Confidence_Regions_v1.2.bed')
    truth = root / 'comparison/common_policy_20260914/wgs_il/ukb/benchmark_truth.vcf.gz'
    targets = {
        'ukb': Path('/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/bio_db/intervals/ukb.pad50.broad.pad50.union.bed'),
        'medexome': root / 'data/SeqCap_EZ_MedExome_hg38_empirical_targets.authoritative.bed',
    }
    sources = [truth, hc, *targets.values()]
    for dataset in ('wes_ll', 'wgs_il'):
        sources.extend((root / 'verified/20260914' / dataset / 'dna_callers').glob('*.vcf.gz'))
    for path in [reference, *sources]:
        if not path.is_file():
            raise FileNotFoundError(path)
    args.outdir = args.outdir.resolve()
    args.outdir.mkdir(parents=True, exist_ok=False)
    report = {
        'status': 'running', 'policy': 'current_native_snv_threshold_indel',
        'git_head': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=repo, text=True).strip(),
        'reference': str(reference), 'reference_size': reference.stat().st_size,
        'source_sha256': {str(p): digest(p) for p in sources},
        'code_sha256': {str(p.relative_to(repo)): digest(p) for p in
                        [repo / 'bin/run_consensus_vcf.py', *sorted((repo / 'bin/vcf_utils').glob('*.py')),
                         *sorted((repo / 'bin/common').glob('*.py'))]},
        'commands': [], 'metrics': {}, 'allele_deltas': {},
        'scope': 'unrestricted candidates; som.py -N with HC -R and target -T; no -P',
    }
    def save():
        (args.outdir / 'audit.json').write_text(json.dumps(report, indent=2) + '\n')
    def run(command, log):
        report['commands'].append(command)
        save()
        with log.open('w') as handle:
            result = subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT)
        if result.returncode:
            report['status'] = 'failed'
            report['failed_log'] = str(log)
            save()
            result.check_returncode()
    save()
    for dataset in ('wes_ll', 'wgs_il'):
        dest = args.outdir / dataset
        dest.mkdir()
        bundle = root / 'verified/20260914' / dataset
        staged = dest / 'inputs'
        staged.mkdir()
        for caller in ('deepsomatic', 'mutect2', 'strelka'):
            source = bundle / 'dna_callers' / (caller + '.vcf.gz')
            for suffix in ('', '.tbi'):
                original = Path(str(source) + suffix)
                if original.is_file():
                    (staged / (dataset + '.' + caller + '.vcf.gz' + suffix)).symlink_to(original.resolve())
        prefix = dest / 'current_native'
        run([sys.executable, '-u', str(repo / 'bin/run_consensus_vcf.py'),
             '--input_dir', str(staged),
             '--expected_callers', 'deepsomatic,mutect2,strelka',
             '--out_prefix', str(prefix), '--snv_thr', '2', '--indel_thr', '2',
             '--min_alt_support', '3', '--native-evidence-snv'], dest / 'consensus.log')
        current = alleles(Path(str(prefix) + '.vcf.gz'), {'Somatic'})
        historical = alleles(bundle / 'native_consensus.vcf.gz', {'Somatic', 'PASS', '.'})
        report['allele_deltas'][dataset] = {
            'added': sorted(current - historical), 'removed': sorted(historical - current),
            'note': 'Exact alleles across differing output scopes, not TP/FP attribution',
        }
        queries = {'current_native': current,
                   'deepsomatic': alleles(bundle / 'dna_callers/deepsomatic.vcf.gz', {'PASS', '.'})}
        for name, selected in queries.items():
            query = write_query(dest / (name + '.query.vcf'), selected)
            for region, target in targets.items():
                cell = dest / region
                cell.mkdir(exist_ok=True)
                output = cell / name
                run(['micromamba', 'run', '-n', 'happy', 'som.py', str(truth), query,
                     '-R', str(hc), '-T', str(target), '-r', str(reference), '-N',
                     '-o', str(output)], cell / (name + '.log'))
                result = parse_metrics_json(Path(str(output) + '.metrics.json'))
                report['metrics'][f'{dataset}/{region}/{name}'] = result
                save()
                print(dataset, region, name, json.dumps(result), flush=True)
    changed = [str(p) for p in sources if digest(p) != report['source_sha256'][str(p)]]
    report['changed_sources'] = changed
    report['status'] = 'failed_source_integrity' if changed else 'complete'
    save()
    if changed:
        raise RuntimeError(f'Source files changed during evaluation: {changed}')


if __name__ == '__main__':
    main()
