#!/usr/bin/env python3
"""Read-only shared-path cohort preflight; never approves labels or trains."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

EXCLUDED = {'PRJNA298330_4032', 'PRJNA298376_4081', 'PRJNA298376_4255'}


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def check(root, expected_report_sha256):
    report_path = root / 'report.json'
    if sha(report_path) != expected_report_sha256:
        raise ValueError('Report hash mismatch; require the recorded release identity')
    report = json.loads(report_path.read_text())
    if (report['status'] != 'review_bridge_verified_not_training_approved'
            or report['training_approved'] is not False):
        raise ValueError('Unexpected review status')
    # Export hashes already bind the exact fully roundtripped labels. Rehashing
    # compressed outputs avoids parsing 14 million rows again on every host.
    for path, digest in report['outputs'].items():
        if sha(path) != digest:
            raise ValueError('Output hash mismatch: ' + path)
    # Code paths identify the exporter host, not the destination model code.
    # Verify shared data sources here; historical code identity stays in report.
    for path, digest in report['sources'].items():
        if Path(path).suffix != '.py' and sha(path) != digest:
            raise ValueError('Source data hash mismatch: ' + path)

    def manifest(name):
        path = root / name
        if str(path) not in report['outputs']:
            raise ValueError('Manifest absent from verified report: ' + name)
        data = json.loads(path.read_text())
        if data['schema_version'] != 1 or data['training_approved'] is not False:
            raise ValueError('Unexpected manifest schema/status')
        rows = data['samples']
        if len({r['sample_id'] for r in rows}) != len(rows):
            raise ValueError('Duplicate manifest sample')
        return {r['sample_id']: r for r in rows}

    all_samples = manifest('samples.review.json')
    expected_ids = {key.split('|')[0] for key in report['counts']}
    if (len(all_samples) != 63 or set(all_samples) != expected_ids
            or set(all_samples) & EXCLUDED):
        raise ValueError('Incorrect cohort membership')
    for pool in ('train_pool', 'reserved'):
        subset = manifest(pool + '.review.json')
        expected = {sid: r for sid, r in all_samples.items() if r['sample_pool'] == pool}
        if subset != expected:
            raise ValueError('Pool manifest disagreement: ' + pool)
    pools = Counter(r['sample_pool'] for r in all_samples.values())
    if dict(pools) != report['pool_counts'] or set(pools) != {'train_pool', 'reserved'}:
        raise ValueError('Pool counts mismatch')
    inventory = {(r['sample_id'], r['modality']): r for r in report['alignment_inventory']}
    if len(inventory) != 189 or len(report['alignment_inventory']) != 189:
        raise ValueError('Alignment inventory must have 189 unique entries')
    for sid, row in all_samples.items():
        for suffix in ('', '.tbi'):
            if row['vcf'] + suffix not in report['outputs']:
                raise ValueError('Unverified label/index path')
        for modality, field in [('dn', 'normal_dna'), ('dt', 'tumor_dna'), ('rt', 'tumor_rna')]:
            entry = inventory[(sid, modality)]
            if row[field] != entry['path'] or Path(row[field]).stat().st_size != entry['size']:
                raise ValueError('Alignment path/size changed: ' + sid)
            for path in (row[field], entry['index']):
                with Path(path).open('rb') as handle:
                    if not handle.read(1):
                        raise ValueError('Empty alignment/index: ' + path)
    if sha(report_path) != expected_report_sha256:
        raise ValueError('Report changed during preflight')
    counts = Counter()
    for key, value in report['counts'].items():
        label = key.split('|')[1]
        if label not in {'Somatic', 'Germline', 'Reference'}:
            raise ValueError('Unexpected class')
        counts[label] += value
    return dict(status='shared_data_preflight_passed_not_training_approved',
                training_approved=False, report_sha256=expected_report_sha256,
                sample_count=len(all_samples), pool_counts=dict(pools),
                class_counts=dict(counts), total_records=sum(counts.values()),
                alignment_paths_checked=len(inventory),
                limitations=['Alignment accessibility/size only, not full BAM integrity or identity',
                             'Actual destination loader, reference and feature cache not validated',
                             'Weak-label selection approval remains separate'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bridge', type=Path, required=True)
    parser.add_argument('--expected-report-sha256', required=True)
    args = parser.parse_args()
    print(json.dumps(check(args.bridge, args.expected_report_sha256), indent=2))
