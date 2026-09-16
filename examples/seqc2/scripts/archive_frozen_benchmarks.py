#!/usr/bin/env python3
"""Export light evidence and an ignored heavy-artifact index; never move inputs."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--archive', required=True, type=Path)
    ap.add_argument('--heavy-links', type=Path)
    args = ap.parse_args()
    provenance = json.loads((args.archive / 'provenance.json').read_text())
    evidence = args.archive / 'evidence'
    evidence.mkdir(exist_ok=False)
    inventory = []
    for dataset, entry in provenance['datasets'].items():
        report = Path(entry['report'])
        if sha(report) != entry['report_sha256']:
            raise ValueError(f'Changed report: {report}')
        target = evidence / dataset
        target.mkdir()
        paths = [report, *report.parent.glob('*.metrics.json'),
                 *report.parent.glob('*.stats.csv'), *report.parent.glob('*.structural_audit_v1.json')]
        for src in paths:
            shutil.copyfile(src, target / src.name)
        for round_name in ('first', 'realignment'):
            src = report.parent / round_name / 'report.json'
            shutil.copyfile(src, target / f'{round_name}.adapter.json')
        for path in sorted(report.parent.rglob('*')):
            if path.is_file() and (path.name.endswith(('.vcf.gz', '.tbi', '.csi', '.sqlite'))):
                inventory.append({'dataset': dataset, 'path': str(path.resolve()),
                                  'bytes': path.stat().st_size, 'role': 'retained_generated_artifact'})
        if args.heavy_links:
            args.heavy_links.mkdir(parents=True, exist_ok=True)
            (args.heavy_links / dataset).symlink_to(report.parent.resolve(), target_is_directory=True)
    (args.archive / 'heavy_artifacts.json').write_text(json.dumps(inventory, indent=2) + '\n')
    checksums = {str(p.relative_to(args.archive)): sha(p) for p in sorted(evidence.rglob('*')) if p.is_file()}
    (args.archive / 'evidence_checksums.json').write_text(json.dumps(checksums, indent=2) + '\n')
    print(f'Archived {len(checksums)} light files; indexed {len(inventory)} heavy artifacts')


if __name__ == '__main__':
    main()
