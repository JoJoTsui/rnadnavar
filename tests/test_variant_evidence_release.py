"""Approval gate tests use small fixtures, never mutate the real release."""
import csv
import json

import polars as pl
import pytest

import export_variant_evidence as e


@pytest.fixture
def release(tmp_path, monkeypatch):
    monkeypatch.setattr(e, 'EXPECTED', {'Reference': 61, 'Germline': 1, 'Somatic': 1})
    bridge = tmp_path / 'bridge'
    bridge.mkdir()
    approval_dir = tmp_path / 'cohort63_handoff_tools_20260921'
    approval_dir.mkdir()
    samples = sorted(e.RESERVED) + [f'PRJNATEST_{i:04d}' for i in range(58)]
    manifests = {}
    tsv_rows = []
    for pool, sids in [('reserved', samples[:5]), ('train_pool', samples[5:])]:
        rows = []
        for sid in sids:
            row = {'sample_id': sid, 'sample_pool': pool}
            tsv = {'sample_id': sid}
            for tag, field in e.MODALITIES.items():
                row[field] = f'/registered/{sid}{tag}.bam'
                tsv['bam_' + tag.lower()] = row[field]
            rows.append(row)
            tsv_rows.append(tsv)
        manifests[pool] = rows
        (bridge / f'{pool}.json').write_text(json.dumps({'training_approved': False, 'samples': rows}))
    with (tmp_path / 'manifest.tsv').open('w') as f:
        w = csv.DictWriter(f, fieldnames=list(tsv_rows[0]), delimiter='\t')
        w.writeheader()
        w.writerows(tsv_rows)
    pl.DataFrame({'sample_id': samples, 'CHROM': ['chr1'] * 63, 'POS': [1] * 63,
                  'REF': ['A'] * 63, 'ALT': ['C'] * 63,
                  'FILTER': ['Germline', 'Somatic'] + ['Reference'] * 61,
                  'training_eligible': [False] * 63, 'TRAINING_ELIGIBLE': ['NO'] * 63}).write_parquet(tmp_path / 'variants.parquet')
    (bridge / 'report.json').write_text('{}')
    approval = {'release_id': e.RELEASE_ID, 'training_approved': True, 'class_counts': e.EXPECTED,
                'total_records': 63, 'excluded_samples': ['EXCLUDED'], 'bridge_directory': 'bridge',
                'variant_parquet': 'variants.parquet', 'manifest_tsv': 'manifest.tsv',
                'train_manifest': 'bridge/train_pool.json', 'reserved_manifest': 'bridge/reserved.json'}
    def pin():
        for field in ['variant_parquet', 'manifest_tsv', 'train_manifest', 'reserved_manifest', 'report']:
            p = tmp_path / (approval[field] if field != 'report' else 'bridge/report.json')
            approval[field + '_sha256'] = e.sha256(p)
        (approval_dir / 'RELEASE_APPROVAL.json').write_text(json.dumps(approval))
    pin()
    return tmp_path, approval, pin


def test_sidecar_approves_immutable_stale_flags(release):
    root, approval, _ = release
    before = e.sha256(root / 'variants.parquet')
    _, _, pools, counts = e.load_release(root)
    assert len(pools) == 63 and sum(row['len'] for row in counts) == 63
    assert sum(v['pool'] == 'train_pool' for v in pools.values()) == 58
    assert {s for s, v in pools.items() if v['pool'] == 'reserved'} == e.RESERVED
    assert e.sha256(root / 'variants.parquet') == before


def test_hash_mismatch_aborts(release):
    root, _, _ = release
    with (root / 'manifest.tsv').open('a') as f:
        f.write('\n')
    with pytest.raises(ValueError, match='hash mismatch'):
        e.load_release(root)


def test_unapproved_sidecar_aborts(release):
    root, approval, pin = release
    approval['training_approved'] = False
    pin()
    with pytest.raises(ValueError, match='explicitly approved'):
        e.load_release(root)


def test_duplicate_pool_membership_aborts(release):
    root, _, pin = release
    path = root / 'bridge/train_pool.json'
    d = json.loads(path.read_text())
    d['samples'].append(d['samples'][0])
    path.write_text(json.dumps(d))
    pin()
    with pytest.raises(ValueError, match='Duplicate'):
        e.load_release(root)


@pytest.mark.parametrize('column,value', [('sample_id', 'EXCLUDED'), ('FILTER', 'Artifact')])
def test_source_membership_and_class_mismatch_aborts(release, column, value):
    root, _, pin = release
    path = root / 'variants.parquet'
    df = pl.read_parquet(path)
    df = df.with_columns(pl.when(pl.int_range(pl.len()) == 62).then(pl.lit(value)).otherwise(pl.col(column)).alias(column))
    df.write_parquet(path)
    pin()
    with pytest.raises(ValueError, match='membership/class totals'):
        e.load_release(root)
