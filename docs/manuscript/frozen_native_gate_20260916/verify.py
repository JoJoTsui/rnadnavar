#!/usr/bin/env python3
"""Verify light archive integrity and frozen code before reproduction (no rerun)."""
import hashlib
import json
import csv
from pathlib import Path


def main():
    archive = Path(__file__).resolve().parent
    repo = archive.parents[2]
    sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    checksums = json.loads((archive/'evidence_checksums.json').read_text())
    for relative, expected in checksums.items():
        if sha(archive/relative) != expected:
            raise ValueError(f'Changed archived evidence: {relative}')
    provenance = json.loads((archive/'provenance.json').read_text())
    with (archive/'metrics.csv').open() as handle:
        rows = list(csv.DictReader(handle))
    keys = {(r['dataset'], r['region'], r['method'], r['variant_type']) for r in rows}
    if len(rows) != 108 or len(keys) != 108:
        raise ValueError('Missing or duplicate metric cells')
    for row in rows:
        report = json.loads((archive/'evidence'/row['dataset']/'validation.json').read_text())
        expected = report['metrics'][row['method']+'/'+row['region']][row['variant_type']]
        if any(float(row[k]) != expected[k] for k in ('tp','fp','fn','precision','recall','f1')):
            raise ValueError(f'Changed metric row: {row}')
    for dataset, entry in provenance['datasets'].items():
        stem = {'seqc2_wes_ll':'seqc2_wes', 'seqc2_wgs_il':'seqc2_wgs', 'hg008_wgs':'hg008'}[dataset]
        manifest = repo/'examples/seqc2/hybrid'/f'frozen_{stem}_validation.json'
        if json.loads(manifest.read_text()) != entry['manifest']:
            raise ValueError(f'Manifest differs from frozen {dataset}')
        code = entry['code_hashes']
        old_repo = next(Path(p).parents[1] for p in code if p.endswith('/bin/run_consensus_vcf.py'))
        for old_path, expected in code.items():
            current = repo/Path(old_path).relative_to(old_repo)
            if sha(current) != expected:
                raise ValueError(f'Code differs from frozen {dataset}: {current}')
    print(f'Verified {len(checksums)} evidence files and all recorded policy code hashes')


if __name__ == '__main__':
    main()
