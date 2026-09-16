#!/usr/bin/env python3
"""Read-only DNA template audit of all frozen WES rescue additions.

Counts reads and unique templates separately; reports conflicting mate bases.
This descriptive assay neither selects variants nor changes the rescue gate.
"""
import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path

import pysam

from audit_current_native_policy import alleles, digest


def summarize_reads(reads, position, ref, alt):
    templates = defaultdict(set)
    observations = []
    counts = Counter()
    for read in reads:
        if read.flag & 0xF04:
            counts['excluded_flags'] += 1
            continue
        qpos = next((q for q, r in read.get_aligned_pairs(matches_only=True) if r == position - 1), None)
        if qpos is None or read.query_sequence is None:
            continue
        base = read.query_sequence[qpos].upper()
        bq = read.query_qualities[qpos] if read.query_qualities is not None else None
        if read.mapping_quality < 20 or bq is None or bq < 20:
            counts['excluded_quality'] += 1
            continue
        counts['passing_reads'] += 1
        counts['alt_reads' if base == alt else 'ref_reads' if base == ref else 'other_reads'] += 1
        identity = (read.get_tag('RG') if read.has_tag('RG') else '', read.query_name)
        templates[identity].add(base)
        if base == alt:
            observations.append({'mapq': read.mapping_quality, 'baseq': bq,
                                 'read_position_fraction': (qpos + 1) / len(read.query_sequence),
                                 'reverse': read.is_reverse, 'proper_pair': read.is_proper_pair,
                                 'edit_distance': read.get_tag('NM') if read.has_tag('NM') else None,
                                 'cigar': read.cigarstring})
    counts['passing_templates'] = len(templates)
    counts['alt_templates'] = sum(bases == {alt} for bases in templates.values())
    counts['ref_templates'] = sum(bases == {ref} for bases in templates.values())
    counts['discordant_templates'] = sum(len(bases) > 1 for bases in templates.values())
    return {'counts': dict(counts), 'alt_observations': observations}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--outdir', required=True, type=Path)
    args = ap.parse_args()
    root = Path(__file__).resolve().parents[3] / 'examples/seqc2'
    native = root / 'verified/20260914/wes_ll/native_consensus.vcf.gz'
    gated = root / 'comparison/rescue_fp_investigation_20260914/gate_tests/wes_ll/nomination_biological/query.vcf.gz'
    evidence = root / 'comparison/rescue_fp_investigation_20260914/evidence.json'
    accepted = {'Somatic', 'PASS', '.'}
    keys = alleles(gated, accepted) - alleles(native, accepted)
    labels = {tuple(x[n] for n in ('chrom', 'pos', 'ref', 'alt')): x['truth_status']
              for x in json.loads(evidence.read_text())['sites'] if x['dataset'] == 'wes_ll'}
    args.outdir.mkdir(parents=True, exist_ok=False)
    report = {'status': 'descriptive_only', 'mapq_min': 20, 'baseq_min': 20, 'excluded_flags': '0xF04',
              'script_sha256': digest(Path(__file__)),
              'sources': {str(p.resolve()): digest(p) for p in (native, gated, evidence)},
              'alignments': {}, 'sites': {}}
    for role in ('dna_tumor', 'dna_normal'):
        bam = root / 'verified/20260914/wes_ll/alignments' / (role + '.bam')
        before = bam.stat()
        report['alignments'][role] = {'path': str(bam.resolve()), 'size': before.st_size, 'mtime_ns': before.st_mtime_ns}
        with pysam.AlignmentFile(str(bam), 'rb') as handle:
            for chrom, pos, ref, alt in sorted(keys):
                if len(ref) != 1 or len(alt) != 1:
                    raise ValueError('Rescue addition is not an SNV')
                key = (chrom, pos, ref, alt)
                name = ':'.join(map(str, key))
                site = report['sites'].setdefault(name, {'key': key, 'truth': labels.get(key, 'unassessed')})
                site[role] = summarize_reads(handle.fetch(chrom, pos - 1, pos), pos, ref, alt)
        after = bam.stat()
        if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
            raise RuntimeError('Alignment changed during read-only audit')
    (args.outdir / 'templates.json').write_text(json.dumps(report, indent=2) + '\n')
    for site in report['sites'].values():
        if site['truth'] in ('TP', 'FP'):
            print(site['truth'], site['key'], site['dna_tumor']['counts'], flush=True)


if __name__ == '__main__':
    main()
