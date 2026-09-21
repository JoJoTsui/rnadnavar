#!/usr/bin/env python3
"""Bounded candidate audit and paired-read pilot, never training approval."""
import argparse
from collections import Counter
from contextlib import ExitStack
import csv
import hashlib
import heapq
import json
from pathlib import Path
import sys

import pyarrow.parquet as pq
import pysam

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT / 'bin'), str(ROOT / 'examples/seqc2/scripts')]
from assess_negative_label_evidence import assess, digest
from pilot_three_class_bam_evidence import evidence, FLAG_FILTER


def choose(heap, site, limit):
    rank = int(hashlib.sha256(':'.join(map(str, site)).encode()).hexdigest(), 16)
    item = (-rank, site)
    if len(heap) < limit:
        heapq.heappush(heap, item)
    elif item > heap[0]:
        heapq.heapreplace(heap, item)


def fingerprint(path):
    path = Path(path)
    stat = path.stat()
    return dict(path=str(path.resolve()), size=stat.st_size, mtime_ns=stat.st_mtime_ns)


def run(config_path, outdir, limit):
    if limit < 1:
        raise ValueError('Require a positive per-class limit')
    cfg = json.loads(config_path.read_text())
    root = Path(cfg['output_root'])
    summary_path = root / 'exports/summary.json'
    verification_path = root / 'review_20260921/parquet_verification.json'
    manifest_path = root / 'exports/manifest.tsv'
    summary = json.loads(summary_path.read_text())
    verified = json.loads(verification_path.read_text())
    parquet = Path(summary['parquet'])
    if (verified['status'] != 'all_parquet_rows_verified'
            or verified['summary_sha256'] != digest(summary_path)
            or verified['manifest_sha256'] != digest(manifest_path)
            or verified['parquet_sha256'] != digest(parquet)):
        raise ValueError('Export verification binding mismatch')
    with manifest_path.open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    manifest = {r['sample_id']: r for r in rows}
    if len(manifest) != 66 or len(rows) != 66:
        raise ValueError('Require the complete 66-sample manifest')
    pilots = cfg['pilot_samples']
    heaps = {(sid, label): [] for sid in pilots for label in ('Germline', 'Reference')}
    if len(pilots) != 3 or len(set(pilots)) != 3 or not set(pilots) <= set(manifest):
        raise ValueError('Require the three configured pilots')
    outdir.mkdir(parents=True, exist_ok=False)
    sources = {str(p): digest(p) for p in (config_path, summary_path, manifest_path,
               verification_path, Path(__file__), ROOT/'bin/assess_negative_label_evidence.py',
               ROOT/'examples/seqc2/scripts/pilot_three_class_bam_evidence.py')}
    report = dict(status='running', training_approved=False, sources=sources,
                  parquet=str(parquet), parquet_sha256=verified['parquet_sha256'],
                  selection='smallest sha256(CHROM:POS:REF:ALT) per pilot/class; all native conflicts separately',
                  parameters=dict(per_class=limit, MAPQ=20, BQ=20, BAQ=True,
                                  max_depth=8000, flag_filter=FLAG_FILTER,
                                  ignore_overlaps=True, ignore_orphans=True,
                                  reference_max_af=.01, confidence=.95, required_depth_each_dna=299),
                  limitations=['Pilot support is not class accuracy or training approval',
                               'BAMs and reference are stat-bound, not whole-file rehashed',
                               'RNA counts are exploratory using the same quality and pairing filters',
                               'No haplotype-aware indel validation or model ablation is performed',
                               'Biological context is not independently adjudicated'], results=[])
    columns = ['sample_id','CHROM','POS','REF','ALT','FILTER','variant_type',
               'review_reason','THREE_CLASS_NATIVE_FILTER','THREE_CLASS_REVIEW_REASON',
               'training_eligible']
    counts, reasons = Counter(), Counter()
    conflicts = []
    stamps = {}
    try:
        for batch in pq.ParquetFile(parquet).iter_batches(batch_size=32768, columns=columns, use_threads=False):
            for row in batch.to_pylist():
                if row['training_eligible'] is not False:
                    raise ValueError('Unexpected training approval')
                sid, label = row['sample_id'], row['FILTER']
                counts[label] += 1
                reasons[row['review_reason']] += 1
                site = (row['CHROM'], row['POS'], row['REF'], row['ALT'])
                if (sid, label) in heaps and row['variant_type'] == 'SNP':
                    choose(heaps[sid, label], site, limit)
                if 'native_negative_conflict' in (row['THREE_CLASS_REVIEW_REASON'] or ''):
                    conflicts.append((sid, label, row['THREE_CLASS_NATIVE_FILTER'], site, 'targeted_conflict'))
        if dict(counts) != summary['class_counts'] or dict(reasons) != summary['review_reasons']:
            raise ValueError('Candidate audit disagrees with verified export')
        report.update(class_counts=dict(counts), review_reasons=dict(reasons))
        selected = [(sid, label, label, site, 'negative_pilot')
                    for (sid, label), heap in sorted(heaps.items()) for _, site in sorted(heap, reverse=True)]
        selected += conflicts
        for sid in sorted({item[0] for item in selected}):
            row = manifest[sid]
            with ExitStack() as stack:
                fasta = stack.enter_context(pysam.FastaFile(cfg['fasta']))
                bams = {role: stack.enter_context(pysam.AlignmentFile(row['bam_'+suffix],
                        reference_filename=cfg['fasta']))
                        for role, suffix in [('normal','dn'),('tumor','dt'),('rna','rt')]}
                for path in [cfg['fasta'], cfg['fasta']+'.fai'] + [bam.filename.decode() for bam in bams.values()]:
                    stamps[path] = fingerprint(path)
                for bam in bams.values():
                    if not bam.has_index():
                        raise ValueError('Missing existing BAM/CRAM index')
                for sample, label, native, site, group in selected:
                    if sample != sid:
                        continue
                    measured = {role: evidence(bam, fasta, site) for role, bam in bams.items()}
                    status, reason = assess(native, site[2], site[3], measured['normal'], measured['tumor'])
                    report['results'].append(dict(sample_id=sid, site=site, final_class=label,
                        native_class=native, group=group, provisional_status=status, reason=reason, **measured))
                print(f'{sid}: paired DNA + RNA pilot complete', flush=True)
        if any(fingerprint(path) != stamp for path, stamp in stamps.items()):
            raise ValueError('BAM/reference changed during assessment')
        if any(digest(Path(path)) != sha for path, sha in sources.items()) or digest(parquet) != report['parquet_sha256']:
            raise ValueError('Evaluation inputs/code changed')
        report.update(status='pilot_complete_not_training_approved', input_fingerprints=stamps,
                      outcomes=dict(Counter(r['group']+':'+r['native_class']+':'+r['provisional_status']+':'+r['reason']
                                            for r in report['results'])))
    except Exception as exc:
        report.update(status='failed', error=str(exc))
        raise
    finally:
        with (outdir/'report.json').open('x') as handle:
            json.dump(report, handle, indent=2)
            handle.write('\n')
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--outdir', type=Path, required=True)
    parser.add_argument('--per-class', type=int, default=32)
    args = parser.parse_args()
    print(json.dumps(run(args.config, args.outdir, args.per_class)['outcomes'], indent=2))
