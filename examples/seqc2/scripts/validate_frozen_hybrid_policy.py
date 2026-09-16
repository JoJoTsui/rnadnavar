#!/usr/bin/env python3
"""Evaluate frozen opt-in policies using existing full-domain caller VCFs.

No workflow launch, truth-driven selection, default promotion or input writes.
The output remains experimental: historical annotation provenance and label QC
must be reviewed separately before training use.
"""
import argparse
import json
from pathlib import Path
import subprocess
import sys

import pysam

from aggregate_benchmark import parse_metrics_json
from validate_refined_native_integration import digest

CALLERS = ('deepsomatic', 'mutect2', 'strelka')


def load_manifest(path):
    data = json.loads(path.read_text())
    for key in ('dna', 'rna_first', 'rna_realignment'):
        if set(data[key]) != set(CALLERS):
            raise ValueError(f'{key}: require exactly three callers')
        if len(set(data[key].values())) != 3:
            raise ValueError(f'{key}: caller files must be distinct')
    if set(data['rescues']) != {'first', 'realignment'}:
        raise ValueError('Require both rescue rounds')
    if set(data['targets']) != {'ukb', 'medexome'}:
        raise ValueError('Require UKB and MedExome targets')
    sources = [data[k] for k in ('truth', 'hc', 'fasta')]
    for key in ('dna', 'rna_first', 'rna_realignment', 'rescues', 'targets'):
        sources.extend(data[key].values())
    for value in sources:
        p = Path(value)
        if not p.is_absolute() or not p.is_file():
            raise ValueError(f'Require existing absolute input path: {p}')
    return data, sorted(set(sources))


def pass_query(source, dest, biological):
    """Stream an allele-only benchmark copy; never relabel the source VCF."""
    with pysam.VariantFile(str(source)) as reader:
        header = pysam.VariantHeader()
        for name, contig in reader.header.contigs.items():
            if contig.length is None:
                header.contigs.add(name)
            else:
                header.contigs.add(name, length=contig.length)
        with pysam.VariantFile(str(dest), 'wz', header=header) as writer:
            for record in reader:
                filters = set(record.filter)
                selected = filters == {'Somatic'} if biological else filters <= {'PASS', '.'}
                if selected:
                    copied = writer.new_record(contig=record.contig, start=record.start,
                                               alleles=record.alleles)
                    copied.filter.add('PASS')
                    writer.write(copied)
    # Inputs are coordinate-sorted; index only the newly written query.
    pysam.tabix_index(str(dest), preset='vcf', force=False)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--manifest', required=True, type=Path)
    ap.add_argument('--outdir', required=True, type=Path)
    args = ap.parse_args()
    data, sources = load_manifest(args.manifest)
    repo = Path(__file__).resolve().parents[3]
    out = args.outdir.resolve()
    out.mkdir(parents=True, exist_ok=False)
    report = {'status': 'freezing', 'manifest': data, 'commands': [], 'metrics': {},
              'scope': __doc__, 'code': {}, 'sources': {}, 'outputs': {}}

    def save():
        (out / 'validation.json').write_text(json.dumps(report, indent=2) + '\n')

    def run(cmd, log):
        report['commands'].append([str(x) for x in cmd])
        save()
        with (out / log).open('w') as handle:
            subprocess.run(cmd, stdout=handle, stderr=subprocess.STDOUT, check=True)

    def score(name, source, biological):
        query = out / f'{name}.query.vcf.gz'
        pass_query(source, query, biological)
        report['outputs'][name] = {'vcf': str(source), 'sha256': digest(Path(source))}
        for region, bed in data['targets'].items():
            prefix = out / f'{name}.{region}'
            run(['micromamba', 'run', '-n', 'happy', 'som.py', data['truth'], str(query),
                 '-R', data['hc'], '-T', bed, '-r', data['fasta'], '-N', '-o', str(prefix)],
                f'{name}.{region}.log')
            metrics = parse_metrics_json(Path(str(prefix) + '.metrics.json'))
            report['metrics'][f'{name}/{region}'] = metrics
            save()
            print(name, region, metrics['records'], flush=True)

    try:
        code = [repo / 'bin/run_consensus_vcf.py', repo / 'bin/apply_refined_rescue.py',
                Path(__file__).resolve(), Path(__file__).with_name('aggregate_benchmark.py'),
                Path(__file__).with_name('validate_refined_native_integration.py'),
                *sorted((repo / 'bin/vcf_utils').glob('*.py')),
                *sorted((repo / 'bin/common').glob('*.py'))]
        report['code'] = {str(p): digest(p) for p in code}
        # Hash the FASTA as well: slow but unambiguous reference provenance.
        report['sources'] = {p: digest(Path(p)) for p in sources}
        report['manifest_sha256'] = digest(args.manifest)
        report['status'] = 'running'
        save()
        inputs = out / 'dna_inputs'
        inputs.mkdir()
        for caller in CALLERS:
            dest = inputs / f'sample.{caller}.vcf.gz'
            # No HC/target/truth filtering. Reject reference mismatches rather
            # than silently discard them. Match development exact deduplication.
            run(['bcftools', 'norm', '-f', data['fasta'], '-c', 'e', '-d', 'exact',
                 data['dna'][caller], '-Oz', '-o', str(dest)], f'normalize.{caller}.log')
            run(['bcftools', 'index', '-t', str(dest)], f'index.{caller}.log')
        run([sys.executable, str(repo / 'bin/run_consensus_vcf.py'), '--input_dir', str(inputs),
             '--expected_callers', ','.join(CALLERS), '--experimental-refined-native',
             '--out_prefix', str(out / 'refined')], 'consensus.log')
        baseline = out / 'refined.vcf.gz'
        if not baseline.is_file():
            raise RuntimeError('Expected consensus output missing')
        score('refined_consensus', baseline, True)
        score('deepsomatic', data['dna']['deepsomatic'], False)
        for round_name in ('first', 'realignment'):
            dest = out / round_name
            cmd = [sys.executable, str(repo / 'bin/apply_refined_rescue.py'),
                   '--dna-consensus', str(baseline), '--annotated-rescue', data['rescues'][round_name],
                   '--alignment-round', round_name, '--outdir', str(dest)]
            for modality, panel in [('dna', data['dna']), ('rna', data[f'rna_{round_name}'])]:
                for caller in CALLERS:
                    cmd += [f'--{modality}-vcf', f'{caller}={panel[caller]}']
            run(cmd, f'{round_name}.log')
            adapter = json.loads((dest / 'report.json').read_text())
            if not adapter['sources_unchanged']:
                raise RuntimeError('Adapter input integrity failure')
            score(f'refined_{round_name}_rescue', dest / 'refined.rescue.vcf.gz', True)
            score(f'workflow_{round_name}_rescue', data['rescues'][round_name], True)
        report['sources_unchanged'] = all(digest(Path(p)) == h for p, h in report['sources'].items())
        report['code_unchanged'] = all(digest(Path(p)) == h for p, h in report['code'].items())
        report['manifest_unchanged'] = digest(args.manifest) == report['manifest_sha256']
        if not all(report[k] for k in ('sources_unchanged', 'code_unchanged', 'manifest_unchanged')):
            raise RuntimeError('Frozen evaluation integrity failure')
        report['status'] = 'complete_experimental_not_training_ready'
        save()
    except Exception as exc:
        report.update(status='failed', error=f'{type(exc).__name__}: {exc}')
        save()
        raise


if __name__ == '__main__':
    main()
