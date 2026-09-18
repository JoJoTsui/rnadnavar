#!/usr/bin/env python3
"""Evidence-only gate for native negative candidates, never training approval.

Consumes a hash-bound paired-BAM pilot report. Preserves every FILTER and source
record. Indels/complex contexts abstain pending haplotype-aware adjudication.
No caller votes are summed and missing coverage is never interpreted as zero.
"""
import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path

import pysam

POLICY = "negative_evidence_gate_v1"
BIOLOGICAL_LABELS = {'Somatic','Germline','Reference','Artifact','NoConsensus','RNAedit'}


def digest(path):
    value = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            value.update(chunk)
    return value.hexdigest()


def required_zero_alt_depth(max_af, confidence):
    if (not math.isfinite(max_af) or not math.isfinite(confidence)
            or not 0 < max_af < 1 or not 0 < confidence < 1):
        raise ValueError('Require finite AF/confidence strictly between zero and one')
    return math.ceil(math.log1p(-confidence) / math.log1p(-max_af))


def valid_counts(value):
    if not isinstance(value, dict) or value.get('status'):
        return False
    fields = [value.get(k) for k in ('depth', 'ref', 'alt', 'other')]
    return (all(type(x) is int and x >= 0 for x in fields)
            and fields[0] == sum(fields[1:]) and fields[0] > 0)


def assess(label, ref, alt, normal, tumor, max_af=.01, confidence=.95,
           context_review=None):
    required = required_zero_alt_depth(max_af, confidence)
    if label not in {'Germline', 'Reference'}:
        return 'NOT_APPLICABLE', 'not_a_negative_candidate'
    if len(ref) != 1 or len(alt) != 1 or set(ref + alt) - set('ACGT'):
        return 'WITHHELD', 'requires_haplotype_aware_validation'
    if context_review is not None and context_review != 'clear':
        return 'WITHHELD', 'unresolved_biological_context'
    if not valid_counts(normal) or not valid_counts(tumor):
        return 'WITHHELD', 'missing_or_invalid_paired_read_evidence'
    if label == 'Reference':
        if normal['alt'] or tumor['alt']:
            return 'WITHHELD', 'observed_alt_not_zero'
        if normal['other'] or tumor['other']:
            return 'WITHHELD', 'other_allele_context_requires_review'
        if min(normal['depth'], tumor['depth']) < required:
            return 'WITHHELD', 'insufficient_depth_for_reference_detection_limit'
        return 'SUPPORTED', 'paired_zero_alt_with_stated_detection_limit'
    if normal['depth'] >= 60 and normal['alt'] == 0:
        return 'CONFLICT', 'normal_reference_evidence'
    if normal['other'] or tumor['other']:
        return 'WITHHELD', 'other_allele_context_requires_review'
    if (normal['depth'] < 20 or normal['alt'] < 5
            or normal['alt'] / normal['depth'] < .2
            or tumor['depth'] < 20 or tumor['alt'] < 3):
        return 'WITHHELD', 'insufficient_paired_germline_evidence'
    # Large shifts can reflect purity/CNV/LOH or somatic-on-germline. They do
    # not imply Somatic, but cannot be automatically approved as simple labels.
    if abs(normal['alt'] / normal['depth'] - tumor['alt'] / tumor['depth']) > .3:
        return 'WITHHELD', 'allele_fraction_shift_requires_context_review'
    return 'SUPPORTED', 'paired_germline_read_corroboration_not_truth'


def evidence_index(report, vcf_sha):
    if report.get('vcf_sha256') != vcf_sha:
        raise ValueError('BAM evidence is not bound to this exact candidate VCF')
    thresholds = report.get('thresholds', {})
    if (thresholds.get('MAPQ', 0) < 20 or thresholds.get('BQ', 0) < 20
            or thresholds.get('BAQ') is not True):
        raise ValueError('Require MAPQ/BQ >=20 and BAQ-enabled BAM evidence')
    result = {}
    for row in report['results']:
        key = tuple(row['site'])
        if len(key) != 4 or key in result:
            raise ValueError('Malformed or duplicate evidence allele')
        result[key] = row
    return result


def native_nomination(info):
    trace = info.get('CLASSIFICATION_RATIONALE', '')
    if isinstance(trace, tuple):
        trace = '|'.join(trace)
    return ('three_class_policy:native_three_class_v1' in trace.split('|')
            or info.get('GATE_POLICY') == 'native_three_class_gate_v1')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--vcf', type=Path, required=True)
    ap.add_argument('--bam-evidence', type=Path, required=True)
    ap.add_argument('--outdir', type=Path, required=True)
    ap.add_argument('--max-reference-af', type=float, default=.01)
    ap.add_argument('--confidence', type=float, default=.95)
    args = ap.parse_args()
    depth = required_zero_alt_depth(args.max_reference_af, args.confidence)
    sources = {str(p.resolve()): digest(p) for p in (args.vcf, args.bam_evidence)}
    raw = json.loads(args.bam_evidence.read_text())
    evidence = evidence_index(raw, sources[str(args.vcf.resolve())])
    args.outdir.mkdir(parents=True, exist_ok=False)
    out = args.outdir / 'candidate.evidence.vcf.gz'
    partial = args.outdir / 'candidate.evidence.partial.vcf.gz'
    definitions = {
        'NEGATIVE_EVIDENCE_POLICY': 'Provisional evidence gate identifier; not training approval',
        'NEGATIVE_EVIDENCE_STATUS': 'SUPPORTED, WITHHELD, CONFLICT or NOT_APPLICABLE; read support is not truth',
        'NEGATIVE_EVIDENCE_REASON': 'Reason for evidence status; FILTER is unchanged',
        'TRAINING_ELIGIBLE': 'NO: this tool cannot grant biological training approval',
        'NEGATIVE_NORMAL_COUNTS': 'Quality-filtered ref:alt:other:depth counts; absent evidence is missing',
        'NEGATIVE_TUMOR_COUNTS': 'Quality-filtered ref:alt:other:depth counts; absent evidence is missing',
    }
    summary = {'policy':POLICY, 'training_approved':False, 'sources':sources,
               'parameters':{'max_reference_af':args.max_reference_af,'confidence':args.confidence,
                             'required_zero_alt_depth_each_sample':depth},
               'scope':__doc__, 'counts':{}, 'status':'running'}
    counts = Counter()
    try:
        with pysam.VariantFile(str(args.vcf)) as reader:
            header = reader.header.copy()
            for name, description in definitions.items():
                if name in header.info:
                    raise ValueError('Refuse recursive evidence gating or overwritten eligibility: '+name)
                header.info.add(name, 1, 'String', description)
            with pysam.VariantFile(str(partial), 'wz', header=header) as writer:
                for record in reader:
                    if len(record.filter) != 1 or not set(record.filter) <= BIOLOGICAL_LABELS:
                        raise ValueError('Require a single biological candidate FILTER, not a raw caller VCF')
                    label = next(iter(record.filter), '')
                    key = (record.contig, record.pos, record.ref, ','.join(record.alts or []))
                    row = evidence.get(key)
                    if row is not None and row.get('class') != label:
                        raise ValueError('Evidence candidate class mismatch: '+str(key))
                    if label in {'Germline','Reference'} and not native_nomination(record.info):
                        status, reason = 'WITHHELD', 'no_native_negative_nomination_provenance'
                    else:
                        row = row or {}
                        status, reason = assess(label, key[2], key[3], row.get('normal'), row.get('tumor'),
                                                args.max_reference_af, args.confidence, row.get('context_review'))
                    copied = record.copy()
                    copied.translate(header)
                    for name, value in zip(definitions, (POLICY, status, reason, 'NO')):
                        copied.info[name] = value
                    for role in ('normal', 'tumor'):
                        measured = (row or {}).get(role)
                        if valid_counts(measured):
                            copied.info['NEGATIVE_'+role.upper()+'_COUNTS'] = ':'.join(
                                str(measured[k]) for k in ('ref','alt','other','depth'))
                    writer.write(copied)
                    counts[label+':'+status+':'+reason] += 1
        summary['sources_unchanged'] = all(digest(Path(p)) == sha for p,sha in sources.items())
        if not summary['sources_unchanged']:
            raise ValueError('Sources changed during evidence gating')
        partial.rename(out)
        summary.update(status='complete_not_training_approved', counts=dict(counts),
                       output=str(out.resolve()), output_sha256=digest(out))
    except Exception as error:
        summary.update(status='failed', error=str(error), counts=dict(counts))
        raise
    finally:
        (args.outdir/'report.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary['counts']), flush=True)


if __name__ == '__main__':
    main()
