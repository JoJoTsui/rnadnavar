#!/usr/bin/env python3
"""Opt-in separated three-class candidates; preserve the declared Somatic baseline.

This tool never grants training approval. No truth VCF is accepted. Native
negative nominations and a stage-specific Somatic baseline are explicit,
hash-bound inputs. SQLite bounds memory for multi-million-record unions.
"""
import argparse
from collections import Counter
import json
from pathlib import Path
import sqlite3
import zlib

import pysam

from apply_refined_rescue import digest, info_dict, LABELS
from vcf_utils.refined_rescue_policy import biological_veto

POLICY = 'separated_three_class_v2'


def decide_label(baseline, native):
    if baseline == 'Somatic':
        return 'Somatic', 'retained_stage_somatic_baseline'
    if native in {'Germline', 'Reference'}:
        return native, 'native_dna_negative_candidate'
    return 'NoConsensus', 'no_supported_native_negative_nomination'


def stage(db, path, column):
    with pysam.VariantFile(str(path)) as reader:
        header = reader.header.copy()
        if header.samples or 'THREE_CLASS_POLICY' in header.info:
            raise ValueError('Require sampleless, non-recursively labelled input')
        for record in reader:
            labels = set(record.filter)
            if len(labels) != 1 or not labels <= LABELS or len(record.alts or []) != 1:
                raise ValueError('Require biallelic biological labels')
            label = next(iter(labels))
            if column == 'native' and label in {'Germline', 'Reference'}:
                trace = record.info.get('CLASSIFICATION_RATIONALE', '')
                if 'three_class_policy:native_three_class_v1' not in str(trace).split('|'):
                    raise ValueError('Negative nomination lacks native policy provenance')
            key = (record.contig, record.pos, record.ref, record.alts[0])
            if record.contig not in header.contigs:
                raise ValueError('Undeclared contig')
            db.execute('INSERT OR IGNORE INTO variants(chrom,pos,ref,alt) VALUES (?,?,?,?)', key)
            if db.execute(f'SELECT {column} FROM variants WHERE chrom=? AND pos=? AND ref=? AND alt=?', key).fetchone()[0] is not None:
                raise ValueError('Duplicate input allele')
            db.execute(f'UPDATE variants SET {column}=? WHERE chrom=? AND pos=? AND ref=? AND alt=?',
                       (zlib.compress(str(record).rstrip('\n').encode()), *key))
            if column == 'baseline' and label == 'Somatic':
                db.execute('UPDATE variants SET baseline_somatic=1 WHERE chrom=? AND pos=? AND ref=? AND alt=?', key)
        db.commit()
    return header


def run(baseline, native, outdir, baseline_sha, native_sha, stage_name):
    sources = {str(p.resolve()):digest(p) for p in (baseline, native)}
    if sources[str(baseline.resolve())] != baseline_sha or sources[str(native.resolve())] != native_sha:
        raise ValueError('Input checksum does not match the declared baseline/candidates')
    outdir.mkdir(parents=True, exist_ok=False)
    report = dict(status='running', policy=POLICY, stage=stage_name, sources=sources,
                  training_approved=False, counts={}, code_sha256=digest(Path(__file__)))
    partial = outdir/'candidates.partial.vcf.gz'
    output = outdir/'candidates.vcf.gz'
    definitions = {
        'THREE_CLASS_POLICY':'Separated candidate policy; not training approval',
        'THREE_CLASS_BASELINE_FILTER':'Original stage-specific baseline class',
        'THREE_CLASS_NATIVE_FILTER':'Native DNA candidate class; missing means unavailable',
        'THREE_CLASS_BASELINE_RATIONALE':'Original baseline rationale; historical provenance',
        'THREE_CLASS_NATIVE_RATIONALE':'Native DNA nomination rationale; not genotype truth',
        'THREE_CLASS_REVIEW_REASON':'Review flags do not silently change Somatic membership',
        'TRAINING_ELIGIBLE':'NO: separate biological approval remains required',
    }
    counts = Counter()
    try:
        with sqlite3.connect(outdir/'union.sqlite') as db:
            db.execute('PRAGMA cache_size=-32768')
            db.execute('CREATE TABLE variants(chrom TEXT,pos INTEGER,ref TEXT,alt TEXT,baseline BLOB,native BLOB,baseline_somatic INTEGER DEFAULT 0,output_somatic INTEGER DEFAULT 0,PRIMARY KEY(chrom,pos,ref,alt))')
            header = stage(db, baseline, 'baseline')
            other = stage(db, native, 'native')
            for name in set(header.contigs) & set(other.contigs):
                a,b = header.contigs[name].length, other.contigs[name].length
                if a and b and a != b:
                    raise ValueError('Reference dictionary mismatch')
            for name in set(header.info) & set(other.info):
                a,b = header.info[name], other.info[name]
                if (a.number,a.type) != (b.number,b.type):
                    raise ValueError('Incompatible INFO: '+name)
            header.merge(other)
            for name,description in definitions.items():
                if name in header.info:
                    raise ValueError('Eligibility/provenance already exists: '+name)
                header.info.add(name,1,'String',description)
            if 'CLASSIFICATION_RATIONALE' not in header.info:
                header.info.add('CLASSIFICATION_RATIONALE',1,'String','Decision trace')
            for label in LABELS:
                if label not in header.filters:
                    header.filters.add(label,None,None,'Biological candidate class')
            with pysam.BGZFile(str(partial),'w') as writer:
                writer.write(str(header).encode())
                for chrom in header.contigs:
                    for pos,ref,alt,b,n in db.execute('SELECT pos,ref,alt,baseline,native FROM variants WHERE chrom=? ORDER BY pos,ref,alt',(chrom,)):
                        b = zlib.decompress(b).decode().split('\t') if b else None
                        n = zlib.decompress(n).decode().split('\t') if n else None
                        bl,nl = b[6] if b else None, n[6] if n else None
                        label,reason = decide_label(bl,nl)
                        bi,ni = info_dict(b[7]) if b else {}, info_dict(n[7]) if n else {}
                        # Negative evidence comes from DNA, not from the old
                        # rescue's failed-Somatic interpretation. Retain other
                        # baseline annotations only when absent from DNA.
                        parts = list(n if label in {'Germline','Reference'} else (b or n))
                        info = dict(ni if label in {'Germline','Reference'} else bi if b else ni)
                        if label in {'Germline','Reference'}:
                            for key,value in bi.items():
                                info.setdefault(key,value)
                        review = []
                        if label == 'Somatic' and nl in {'Germline','Reference'}:
                            review.append('native_negative_conflict')
                        veto = biological_veto(bi) or biological_veto(ni)
                        if veto:
                            review.append(veto)
                        if label == 'Somatic' and bi.get('DNA_VERIFICATION') in {'rejected','inconclusive'}:
                            review.append('dna_verification_'+bi['DNA_VERIFICATION'])
                        info.update(THREE_CLASS_POLICY=POLICY,
                                    THREE_CLASS_BASELINE_FILTER=bl or 'MISSING',
                                    THREE_CLASS_NATIVE_FILTER=nl or 'MISSING',
                                    THREE_CLASS_BASELINE_RATIONALE=bi.get('CLASSIFICATION_RATIONALE','missing'),
                                    THREE_CLASS_NATIVE_RATIONALE=ni.get('CLASSIFICATION_RATIONALE','missing'),
                                    THREE_CLASS_REVIEW_REASON='|'.join(review) or 'none',
                                    TRAINING_ELIGIBLE='NO',
                                    CLASSIFICATION_RATIONALE=f'policy:{POLICY}|decision:{reason}|class:{label}')
                        parts[6]=label
                        parts[7]=';'.join(k if v is None else k+'='+str(v) for k,v in info.items())
                        writer.write(('\t'.join(parts)+'\n').encode())
                        counts[label]+=1
                        for reason in review:
                            counts['review:'+reason]+=1
                        if label=='Somatic':
                            db.execute('UPDATE variants SET output_somatic=1 WHERE chrom=? AND pos=? AND ref=? AND alt=?',(chrom,pos,ref,alt))
            mismatches = db.execute('SELECT count(*) FROM variants WHERE baseline_somatic != output_somatic').fetchone()[0]
            report['somatic_membership_mismatches']=mismatches
            if mismatches:
                raise ValueError('Somatic preservation invariant failed')
            db.commit()
        report['sources_unchanged']=all(digest(Path(p))==h for p,h in sources.items())
        if not report['sources_unchanged'] or digest(Path(__file__))!=report['code_sha256']:
            raise ValueError('Inputs/code changed during execution')
        partial.rename(output)
        pysam.tabix_index(str(output),preset='vcf')
        report.update(status='candidate_complete_not_training_approved',counts=dict(counts),
                      output=str(output.resolve()),output_sha256=digest(output))
    except Exception as exc:
        report.update(status='failed',error=str(exc))
        raise
    finally:
        (outdir/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    return report


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--somatic-baseline',required=True,type=Path)
    ap.add_argument('--native-candidates',required=True,type=Path)
    ap.add_argument('--expected-baseline-sha256',required=True)
    ap.add_argument('--expected-native-sha256',required=True)
    ap.add_argument('--stage',choices=['consensus','first','realignment'],required=True)
    ap.add_argument('--outdir',required=True,type=Path)
    args=ap.parse_args()
    print(json.dumps(run(args.somatic_baseline,args.native_candidates,args.outdir,
                         args.expected_baseline_sha256,args.expected_native_sha256,args.stage)['counts']))


if __name__=='__main__':
    main()
