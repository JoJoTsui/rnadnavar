"""Read-only exact-allele diagnostic, not normalized benchmark scoring.

Prints JSON; merges BED intervals and uses POS-1 membership.
Unsplit multi-ALT records and indel representations limit exact attribution.
"""
import bisect
import collections as C
import gzip
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / 'examples/seqc2/hybrid/output/seqc2.wes.ll.hybrid.realign.full'
BENCH = ROOT / 'examples/seqc2/hybrid/comparison/comprehensive_realign/WES_LL_T_1_vs_WES_LL_N_1'

def bed(path):
    raw, out = C.defaultdict(list), {}
    with open(path) as f:
        for line in f:
            a = line.split()
            if len(a) >= 3 and a[1].isdigit():
                raw[a[0]].append((int(a[1]), int(a[2])))
    for c, rows in raw.items():
        merged = []
        for s, e in sorted(rows):
            if merged and s <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(e, merged[-1][1]))
            else:
                merged.append((s, e))
        out[c] = merged
    return out

BEDS = [bed(p) for p in ['/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/data/giab/data/seqc2/truth/High-Confidence_Regions_v1.2.bed', '/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/bio_db/intervals/ukb.pad50.broad.pad50.union.bed']]

def inside(c, p):
    for b in BEDS:
        rows = b.get(c, [])
        i = bisect.bisect_right(rows, (p, float('inf'))) - 1
        if i < 0 or not rows[i][0] <= p < rows[i][1]:
            return False
    return True

def read(path):
    rows, duplicates = {}, []
    with gzip.open(path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            a = line.rstrip().split('\t')
            k = (a[0], int(a[1]), a[3], a[4])
            if not inside(k[0], k[1]-1):
                continue
            if k in rows:
                duplicates.append(k)
            rows[k] = (a[6], dict(x.split('=', 1) if '=' in x else (x, True) for x in a[7].split(';')))
    return rows, duplicates

def main():
    patterns = {'dna':'consensus/WES_LL_T_1_vs*/*.consensus.vcf.gz', 'first':'rescue/*/*.filtered.vcf.gz', 'second':'vcf_realignment/rescue/*/*.rescue.filtered.stripped.vep.vcf.gz', 'ds':'variant_calling/deepsomatic/WES_LL_T_1_vs*/*.deepsomatic.vcf.gz', 'mutect':'variant_calling/mutect2/WES_LL_T_1_vs*/*.mutect2.filtered.vcf.gz', 'strelka':'variant_calling/strelka/WES_LL_T_1_vs*/*.variants.vcf.gz'}
    paths = {'truth': BENCH/'high-confidence_sSNV+INDEL_in_HC_regions_v1.2.1.vcf.gz'}
    for name, pattern in patterns.items():
        matches = list(BASE.glob(pattern))
        if len(matches) != 1:
            raise ValueError((name, matches))
        paths[name] = matches[0]
    data, out = {}, {'sources':{k:str(v) for k,v in paths.items()}, 'duplicates':{}, 'method':__doc__}
    for name,p in paths.items():
        data[name], out['duplicates'][name] = read(p)
    truth = set(data['truth'])
    selected = {n:{k for k,(f,_) in d.items() if f == ('PASS' if n in ('ds','mutect','strelka') else 'Somatic')} for n,d in data.items() if n != 'truth'}
    out['truth_distinct'] = len(truth)
    out['transitions'] = {}
    for before,after in [('ds','dna'),('dna','first'),('dna','second'),('first','second')]:
        added,lost = selected[after]-selected[before], selected[before]-selected[after]
        out['transitions'][before+'->'+after] = dict(added_truth=len(added&truth), added_nonmatch=len(added-truth), lost_truth=len(lost&truth), removed_nonmatch=len(lost-truth))
    any_caller = set().union(*(set(data[n]) for n in ('ds','mutect','strelka')))
    any_pass = set().union(*(selected[n] for n in ('ds','mutect','strelka')))
    out['fn'] = {}
    for n in ('dna','first','second'):
        fn = truth-selected[n]
        out['fn'][n] = dict(absent_all_dna_callers=len(fn-any_caller), caller_present_none_pass=len((fn&any_caller)-any_pass), caller_pass=len(fn&any_pass), labels=dict(C.Counter(data[n].get(k,('ABSENT',{}))[0] for k in fn)))
    out['gnomad'] = {}
    for n in ('first','second'):
        counts = C.Counter()
        for k in selected[n]:
            try:
                af = float(data[n][k][1].get('GNOMAD_AF'))
                bucket = '>=.01' if af >= .01 else ('>=.001' if af >= .001 else '<.001')
            except (TypeError,ValueError):
                bucket = 'unknown'
            counts[('truth' if k in truth else 'nonmatch')+':'+bucket] += 1
        out['gnomad'][n] = dict(counts)
    print(json.dumps(out, indent=2, sort_keys=True))

if __name__ == '__main__':
    main()
