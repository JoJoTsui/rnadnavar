#!/usr/bin/env python3
"""Apply the frozen optimized rescue gate without materializing VCF records.

The evaluator is intentionally VCF-only and writes an allele-only PASS query.
It keeps compact variant keys and performs separate streaming passes over the
rescue, RNA, and DNA VCFs. It never changes workflow VCFs or caller caches.
"""
import argparse
import gzip
import hashlib
import json
import subprocess
from pathlib import Path


def fields(path):
    samples = []
    with gzip.open(path, 'rt') as handle:
        for line in handle:
            if line.startswith('#CHROM'):
                samples = line.rstrip().split('\t')[9:]
                continue
            if line.startswith('#'):
                continue
            row = line.rstrip().split('\t')
            if len(row) < 8:
                continue
            fmt = row[8].split(':') if len(row) > 8 else []
            values = []
            for name, raw in zip(samples, row[9:]):
                values.append((name, dict(zip(fmt, raw.split(':')))))
            for index, alt in enumerate(row[4].split(','), 1):
                yield (row, alt, index, values)


def key(row, alt):
    return row[0], int(row[1]), row[3], alt


def info_map(raw):
    result = {}
    for item in raw.split(';'):
        if not item or item == '.':
            continue
        k, sep, v = item.partition('=')
        result[k] = v if sep else True
    return result


def somatic_filter(value):
    return value in ('Somatic', 'PASS', '.')


def tumor_alt(row, alt_index, alt, caller):
    samples = [(name, value) for name, value in row[3] if name != 'NORMAL' and not name.endswith('_N_1')]
    if len(samples) != 1:
        return False
    sample = samples[0][1]
    if 'AD' in sample:
        try:
            return int(sample['AD'].split(',')[alt_index]) > 0
        except (ValueError, IndexError):
            return False
    if caller == 'strelka':
        for field in (f'{alt}U', 'TIR', 'AO'):
            if field in sample:
                try:
                    return int(sample[field].split(',')[0]) > 0
                except (ValueError, IndexError):
                    return False
    return False


def read_keys(path, filters=None, snv_only=False, exclude=None):
    output = set()
    exclude = exclude or set()
    for row, alt, _, _ in fields(path):
        if filters is not None and row[6] not in filters:
            continue
        if snv_only and (len(row[3]) != 1 or len(alt) != 1):
            continue
        item = key(row, alt)
        if item not in exclude:
            output.add(item)
    return output


def count_rna(path, candidates, counts):
    for row, alt, _, _ in fields(path):
        item = key(row, alt)
        if item in candidates and somatic_filter(row[6]):
            counts[item] = counts.get(item, 0) + 1


def add_dna_nominations(path, caller, candidates, nominated):
    for raw, alt, index, values in fields(path):
        item = key(raw, alt)
        if item not in candidates:
            continue
        if caller == 'deepsomatic' and raw[6] not in ('PASS', '.'):
            continue
        if tumor_alt((raw, alt, index, values), index, alt, caller):
            nominated.add(item)


def biological_veto(path, candidates):
    vetoed = set()
    for row, alt, _, _ in fields(path):
        item = key(row, alt)
        if item not in candidates:
            continue
        info = info_map(row[7])
        af = info.get('GNOMAD_AF') or info.get('gnomAD_AF')
        if af not in (None, '', '.'):
            try:
                if max(float(value) for value in str(af).split(',')) > 0.001:
                    vetoed.add(item)
                    continue
            except ValueError as exc:
                raise ValueError(f'Malformed GNOMAD_AF: {af}') from exc
        if (info.get('REDI_ACCESSION') not in (None, '', '.') and
                info.get('REDI_CANONICAL') == 'YES' and
                info.get('N_DNA_CALLERS_SOMATIC') == '0'):
            vetoed.add(item)
    return vetoed


def write_query(path, keys):
    plain = path.with_suffix('.vcf')
    with plain.open('w') as handle:
        handle.write('##fileformat=VCFv4.2\n')
        for chrom in sorted({item[0] for item in keys}):
            handle.write(f'##contig=<ID={chrom}>\n')
        handle.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
        for chrom, pos, ref, alt in sorted(keys):
            handle.write(f'{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\n')
    subprocess.run(['bcftools', 'view', '-Oz', '-o', str(path), str(plain)], check=True)
    subprocess.run(['bcftools', 'index', '-t', '-f', str(path)], check=True)
    plain.unlink()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--native', type=Path, required=True)
    parser.add_argument('--rescue', type=Path, required=True)
    parser.add_argument('--dna-vcf', action='append', required=True)
    parser.add_argument('--rna-vcf', action='append', required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--report', type=Path)
    args = parser.parse_args()

    native = read_keys(args.native, filters={'Somatic', 'PASS', '.'})
    candidates = read_keys(args.rescue, filters={'Somatic', 'PASS', '.'},
                           snv_only=True, exclude=native)
    rna_counts = {}
    for path in args.rna_vcf:
        count_rna(path, candidates, rna_counts)
    rna_kept = {item for item, count in rna_counts.items() if count >= 2}
    nominated = set()
    for path in args.dna_vcf:
        name = Path(path).name.lower()
        caller = 'deepsomatic' if 'deepsomatic' in name else 'strelka' if 'strelka' in name else 'mutect2'
        add_dna_nominations(path, caller, rna_kept, nominated)
    vetoed = biological_veto(args.rescue, nominated)
    kept = nominated - vetoed
    args.out.parent.mkdir(parents=True, exist_ok=True)
    write_query(args.out, native | kept)
    report = {
        'status': 'exploratory_not_promoted', 'native_records': len(native),
        'candidate_additions': len(candidates), 'rna_supported': len(rna_kept),
        'dna_nominated': len(nominated), 'kept_additions': len(kept),
        'biological_veto': len(vetoed), 'output': str(args.out),
        'sources': {str(Path(p)): hashlib.sha256(Path(p).read_bytes()).hexdigest()
                    for p in [args.native, args.rescue, *args.dna_vcf, *args.rna_vcf]},
        'original_inputs_unchanged': True,
    }
    report_path = args.report or args.out.with_suffix('.report.json')
    report_path.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
