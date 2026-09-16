#!/usr/bin/env python3
"""Read-only, bounded-memory structural audit; not biological label-QC approval."""
import argparse
from collections import Counter
import gzip
import json
from pathlib import Path
from urllib.parse import unquote

from validate_refined_native_integration import digest

LABELS = {'Somatic', 'Germline', 'Reference', 'Artifact', 'NoConsensus', 'RNAedit'}


def audit(path):
    before = digest(path)
    counts, issues, examples = Counter(), Counter(), {}
    def flag(reason, key):
        issues[reason] += 1
        if len(examples.setdefault(reason, [])) < 10:
            examples[reason].append(key)
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            row = line.rstrip('\n').split('\t')
            counts['records'] += 1
            if len(row) != 8:
                flag('not_eight_columns', row[:5])
                continue
            key = [row[0], row[1], row[3], row[4]]
            label = row[6]
            counts[label] += 1
            info = {x.partition('=')[0]: x.partition('=')[2] for x in row[7].split(';')}
            if label not in LABELS:
                flag('invalid_label', key)
            if info.get('UNIFIED_FILTER') != label:
                flag('unified_filter_mismatch', key)
            if not info.get('CLASSIFICATION_RATIONALE'):
                flag('missing_rationale', key)
            if 'GATE_POLICY' not in info:
                continue
            if f'class:{label}' not in info.get('CLASSIFICATION_RATIONALE', '').split('|'):
                flag('rationale_class_mismatch', key)
            snapshot = unquote(info.get('GATE_SOURCE_RECORD', '')).split('\t')
            if len(snapshot) != 8 or [snapshot[0], snapshot[1], snapshot[3], snapshot[4]] != key:
                flag('invalid_source_snapshot', key)
            promoted = info.get('RESCUE_PROMOTED') == 'YES'
            baseline = info.get('PASSES_CONSENSUS_DNA') == 'YES'
            if (promoted or baseline) and label != 'Somatic':
                flag('positive_flag_on_negative', key)
            if promoted and baseline:
                flag('baseline_marked_promoted', key)
            if promoted:
                rna = set(info.get('GATE_RNA_ELIGIBLE', '').split('|')) - {'', '.'}
                dna = set(info.get('GATE_DNA_NOMINATORS', '').split('|')) - {'', '.'}
                if len(rna) < 2 or not dna or len(row[3]) != 1 or len(row[4]) != 1:
                    flag('promotion_evidence_contract', key)
    unchanged = digest(path) == before
    return dict(status='structural_pass_not_training_approval' if not issues and unchanged else 'failed',
                source=str(path.resolve()), sha256=before, sources_unchanged=unchanged,
                counts=dict(counts), issues=dict(issues), examples=examples)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--vcf', required=True, type=Path)
    ap.add_argument('--report', required=True, type=Path)
    args = ap.parse_args()
    # Refuse overwriting a previous audit, including the input itself.
    with args.report.open('x') as handle:
        report = audit(args.vcf)
        json.dump(report, handle, indent=2)
        handle.write('\n')
    return 0 if report['status'].startswith('structural_pass') else 1


if __name__ == '__main__':
    raise SystemExit(main())
