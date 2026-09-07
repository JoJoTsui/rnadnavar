#!/usr/bin/env python3
"""Build auditable training labels from DeepSomatic DNA PASS and verified RNA nominations."""
import argparse, gzip, json
from pathlib import Path

def op(path, mode='rt'):
    return gzip.open(path, mode) if str(path).endswith('.gz') else open(path, mode)

def main():
    ap=argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument('--deepsomatic-vcf',required=True); ap.add_argument('--rna-nominations'); ap.add_argument('--verification-json'); ap.add_argument('--out',required=True)
    args=ap.parse_args(); verified={}
    if args.verification_json:
        data=json.loads(Path(args.verification_json).read_text())
        verified={(str(r['chrom']),str(r['pos']),r['ref'].upper(),r['alt'].upper()):r for r in data.get('results',[]) if r.get('status')=='confirmed'}
    rows=[]
    with op(args.deepsomatic_vcf) as fh:
        for line in fh:
            if line.startswith('#'): rows.append(line); continue
            f=line.rstrip('\n').split('\t')
            if len(f)<8: continue
            filt=f[6].upper(); info=f[7]
            if filt in {'PASS','.','SOMATIC'} or 'PASS' in {x.upper() for x in filt.split(';')}:
                f[6]='Somatic'; info += ';CLASSIFICATION_RATIONALE=rule:deepsomatic_pass_starting_set|class:Somatic'
                f[7]=info
            rows.append('\t'.join(f)+'\n')
    if args.rna_nominations:
        with op(args.rna_nominations) as fh:
            for line in fh:
                if line.startswith('#'): continue
                f=line.rstrip('\n').split('\t')
                if len(f)<8: continue
                key=(f[0],f[1],f[3].upper(),f[4].upper()); ev=verified.get(key)
                if not ev: continue
                f[6]='Somatic'; f[7]=f[7]+('; ' if f[7] not in ('','.','') else '')+'CLASSIFICATION_RATIONALE=rule:verified_dna_rna_nomination|class:Somatic|tumor_alt:%s|normal_alt:%s' % (ev.get('tumor_alt'),ev.get('normal_alt'))
            rows.append('\t'.join(f)+'\n')
    out=Path(args.out); out.parent.mkdir(parents=True,exist_ok=True)
    if str(out).endswith('.gz'):
        with gzip.open(out,'wt') as fh: fh.writelines(rows)
    else: out.write_text(''.join(rows))
if __name__=='__main__': main()
