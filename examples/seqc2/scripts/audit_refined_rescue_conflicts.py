#!/usr/bin/env python3
"""Read-only SEQC2 baseline conflict audit; no rule or label changes.

Compare flagged population annotations with exact gnomAD alleles and retained
benchmark partitions. Does not inspect HG008 or treat unscored records as FP.
"""
import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import sqlite3
import subprocess
import sys
import zlib

import pysam

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / 'bin'))
from apply_refined_rescue import info_dict, digest
from vcf_utils.aggregation import resolve_tumor_sample_index, resolve_normal_sample_index, _normal_sample_from_header
from vcf_utils.refined_rescue_policy import biological_veto
from validate_refined_native_integration import alleles


def file_stat(path):
    st = path.stat()
    return {'size': st.st_size, 'mtime_ns': st.st_mtime_ns}


def paired_evidence(reader, caller, key):
    chrom, pos, ref, alt = key
    names = list(reader.header.samples)
    normal = _normal_sample_from_header(str(reader.header))
    indices = {'tumor': resolve_tumor_sample_index(names, caller, normal),
               'normal': resolve_normal_sample_index(names, caller, normal)}
    observations = []
    for record in reader.fetch(chrom, pos-1, pos):
        if record.pos != pos or record.ref != ref or alt not in (record.alts or []):
            continue
        alt_index = record.alts.index(alt) + 1
        observation = {'filter': ';'.join(record.filter) or '.', 'qual': record.qual, 'samples': {}}
        for role, index in indices.items():
            if index is None or not names:
                observation['samples'][role] = None
                continue
            sample = record.samples[names[index]]
            ad = sample.get('AD')
            count = ad[alt_index] if ad and len(ad) > alt_index else None
            if ad is None and caller == 'strelka' and len(ref) == len(alt) == 1:
                value = sample.get(alt + 'U')
                count = value[0] if value else None
            observation['samples'][role] = {'sample': names[index], 'dp': sample.get('DP'), 'alt_count': count}
        observations.append(observation)
    return observations


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--root', type=Path, required=True)
    ap.add_argument('--outdir', type=Path, required=True)
    ap.add_argument('--gnomad', type=Path, default=Path('/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/bio_db/gnomAD/exomes'))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=False)
    baseline = args.root / 'driver_parity_v1/wes/ukb/refined.vcf.gz'
    output_root = args.root / 'rescue_output_wes_ukb_v1/realignment'
    dbpath = (output_root / 'union.sqlite').resolve()
    scratch = args.root / 'wes_indel_attribution/ukb/current_native_scratch'
    partition_paths = {'TP': scratch / 'tpfn/0003.vcf.gz', 'FP': scratch / 'fp.vcf.gz',
                       'ambiguous': scratch / 'ambi.vcf.gz', 'unscored': scratch / 'unk.vcf.gz'}
    partitions = {label: alleles(path) for label, path in partition_paths.items()}
    sources = [baseline, output_root / 'report.json', *partition_paths.values()]
    report = {'status': 'running', 'scope': __doc__, 'sources': {str(p.resolve()): digest(p) for p in sources},
              'sqlite_stat': file_stat(dbpath), 'database_files': {}, 'conflicts': [],
              'script_sha256': digest(Path(__file__))}
    query = subprocess.run(['bcftools', 'query', '-i', 'FILTER="Somatic"', '-f', '%CHROM\t%POS\t%REF\t%ALT\n', str(baseline)],
                           check=True, capture_output=True, text=True)
    handles = {}
    with sqlite3.connect('file:' + str(dbpath) + '?mode=ro', uri=True) as db:
        for line in query.stdout.splitlines():
            chrom, pos, ref, alt = line.split('\t')
            key = chrom, int(pos), ref, alt
            row = db.execute('SELECT rescue FROM variants WHERE chrom=? AND pos=? AND ref=? AND alt=?', key).fetchone()
            if not row or row[0] is None:
                continue
            raw = zlib.decompress(row[0]).decode()
            parts = raw.split('\t'); info = info_dict(parts[7])
            veto = biological_veto(info)
            if not veto:
                continue
            if chrom not in handles:
                path = args.gnomad / f'gnomad.exomes.v4.1.sites.{chrom}.vcf.bgz'
                report['database_files'][str(path.resolve())] = file_stat(path)
                handles[chrom] = pysam.VariantFile(str(path))
            matches = []
            for record in handles[chrom].fetch(chrom, int(pos)-1, int(pos)):
                if record.pos != int(pos) or record.ref != ref:
                    continue
                for i, allele in enumerate(record.alts or []):
                    if allele != alt:
                        continue
                    af = record.info.get('AF')
                    value = af[i] if isinstance(af, tuple) else af
                    matches.append({'af': value, 'record_sha256': hashlib.sha256(str(record).encode()).hexdigest(),
                                    'filters': list(record.filter)})
            source_af = float(info['GNOMAD_AF']) if info.get('GNOMAD_AF') not in (None,'.') else None
            entry = {'allele': key, 'reason': veto, 'source_label': parts[6], 'source_af': source_af,
                     'database_matches': matches,
                     'af_exact_allele_confirmed': any(m['af'] is not None and source_af is not None
                                                   and math.isclose(m['af'],source_af,rel_tol=1e-5,abs_tol=1e-8) for m in matches),
                     'benchmark_partitions': [label for label, keys in partitions.items() if key in keys],
                     'native_filters': info.get('FILTERS_ORIGINAL'), 'source_rationale': info.get('CLASSIFICATION_RATIONALE')}
            report['conflicts'].append(entry)
    for handle in handles.values():
        handle.close()
    for caller in ('deepsomatic','mutect2','strelka'):
        path = REPO / 'examples/seqc2/verified/20260914/wes_ll/dna_callers' / f'{caller}.vcf.gz'
        path = path.resolve(strict=True)
        report['sources'][str(path)] = digest(path)
        with pysam.VariantFile(str(path)) as reader:
            for row in report['conflicts']:
                row.setdefault('paired_caller_evidence', {})[caller] = paired_evidence(reader,caller,row['allele'])
    report['counts'] = dict(Counter(label for row in report['conflicts'] for label in row['benchmark_partitions']))
    report['all_af_confirmed'] = all(row['af_exact_allele_confirmed'] for row in report['conflicts'])
    report['sources_unchanged'] = all(digest(Path(p)) == sha for p,sha in report['sources'].items())
    report['sqlite_unchanged'] = file_stat(dbpath) == report['sqlite_stat']
    report['database_stats_unchanged'] = all(file_stat(Path(p)) == stat for p,stat in report['database_files'].items())
    report['status'] = 'complete_read_only' if all(report[k] for k in ('sources_unchanged','sqlite_unchanged','database_stats_unchanged')) else 'integrity_failure'
    (args.outdir / 'audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ('status','counts','all_af_confirmed')},indent=2))


if __name__ == '__main__':
    main()
