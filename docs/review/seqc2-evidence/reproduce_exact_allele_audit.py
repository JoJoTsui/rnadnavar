import gzip,collections,bisect,json,hashlib,csv
from pathlib import Path
root=Path('/t9k/mnt/hdd/work/Vax/pipeline/rnadnavar/examples/seqc2'); ds={}; out={}
def bed(path):
 d=collections.defaultdict(list)
 for l in open(path):
  a=l.split()
  if len(a)>=3 and a[1].isdigit():d[a[0]].append((int(a[1]),int(a[2])))
 return {c:sorted(v) for c,v in d.items()}
beds=[bed('/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/data/giab/data/seqc2/truth/High-Confidence_Regions_v1.2.bed'),bed('/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/bio_db/intervals/ukb.pad50.broad.pad50.union.bed')]
def inside(c,p):
 for b in beds:
  v=b.get(c,[]);i=bisect.bisect_right(v,(p,10**15))-1
  if i<0 or not v[i][0]<=p<v[i][1]:return False
 return True
def read(p,filts):
 d={}
 for l in gzip.open(p,'rt'):
  if l.startswith('#'):continue
  a=l.rstrip().split('\t'); k=(a[0],int(a[1]),a[3],a[4])
  if (filts is None or a[6] in filts) and inside(k[0],k[1]-1):d[k]=a[7]
 return d
truth=read(next((root/'comparison/WES_LL_T_1_vs_WES_LL_N_1').glob('high*.vcf.gz')),None)
base=root/'hybrid/output/seqc2.wes.ll.hybrid.pooling_fix'
for name,glob in [('DNA','consensus/WES_LL_T*/*.vcf.gz'),('RT','consensus/WES_LL_RT*/*.vcf.gz'),('rescue','rescue/*/*.filtered.vcf.gz'),('DS','variant_calling/deepsomatic/WES_LL_T*/*.vcf.gz')]:
 p=next(base.glob(glob));d=read(p,{'Somatic'} if name!='DS' else {'PASS'});ds[name]=d
 out[name]={'source':str(p),'total':len(d),'exactTP':len(d.keys()&truth.keys()),'chr1_exactTP':sum(k[0]=='chr1' for k in d.keys()&truth.keys()),'chr1_total':sum(k[0]=='chr1' for k in d),'bad_pass':dict(collections.Counter('TP' if k in truth else 'nonmatch' for k,v in d.items() if 'PASSES_CONSENSUS=NO;' in v))}
for a,b in [('rescue','DNA'),('DNA','DS'),('RT','DNA')]:
 add=ds[a].keys()-ds[b].keys();drop=ds[b].keys()-ds[a].keys()
 md=collections.Counter((('exactTP' if k in truth else 'nonmatch'),tuple(x for x in ds[a][k].split(';') if x.startswith(('RESCUE_PROMOTED=','PASSES_CONSENSUS_DNA=','PASSES_CONSENSUS_RNA=','CROSS_MODALITY=','CLASSIFICATION_RATIONALE=')))) for k in add)
 out[a+'_vs_'+b]={'add':len(add),'add_exactTP':len(add&truth.keys()),'drop':len(drop),'drop_exactTP':len(drop&truth.keys()),'add_metadata':[[list(k),v] for k,v in md.items()]}
out['truth_chr1']=sum(k[0]=='chr1' for k in truth)
out['method']='Exact CHROM/POS/REF/ALT dictionary sets; BED membership uses POS-1 point within HC and target. Does not normalize haplotypes; differs som.py around overlapping indel boundaries. Not replacement official metrics.'
Path('/tmp/rnadnavar-seqc2-review/docs/review/seqc2-evidence/exact_allele_audit.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2)[:18000])
