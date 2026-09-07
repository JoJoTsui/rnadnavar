import gzip, importlib.util, json, subprocess, sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
def load(name,path):
 s=importlib.util.spec_from_file_location(name,ROOT/path); m=importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
def write_vcf(path, rows):
 text='##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'+''.join(r+'\n' for r in rows)
 if str(path).endswith('.gz'):
  with gzip.open(path,'wt') as f:f.write(text)
 else:path.write_text(text)
def test_label_builder_keeps_all_verified_nominations(tmp_path):
 m=load('labels','bin/build_deepsomatic_labels.py'); ds=tmp_path/'ds.vcf'; rna=tmp_path/'rna.vcf'; out=tmp_path/'out.vcf'
 write_vcf(ds,['1\t10\t.\tA\tG\t.\tPASS\t.'])
 write_vcf(rna,['1\t20\t.\tC\tT\t.\t.\t.','1\t30\t.\tG\tA\t.\t.\t.'])
 ver=tmp_path/'verification.json'; ver.write_text(json.dumps({'results':[{'chrom':'1','pos':20,'ref':'C','alt':'T','status':'confirmed','tumor_alt':5,'normal_alt':0},{'chrom':'1','pos':30,'ref':'G','alt':'A','status':'confirmed','tumor_alt':4,'normal_alt':0}]}))
 subprocess.run([sys.executable,str(ROOT/'bin/build_deepsomatic_labels.py'),'--deepsomatic-vcf',str(ds),'--rna-nominations',str(rna),'--verification-json',str(ver),'--out',str(out)],check=True)
 body=out.read_text(); assert sum(1 for line in body.splitlines() if line and not line.startswith('#') and line.split('\t')[6]=='Somatic')==3

def test_scorer_reports_variant_types_and_transitions(tmp_path):
 truth=tmp_path/'truth.vcf'; base=tmp_path/'base.vcf'; calls=tmp_path/'calls.vcf'; out=tmp_path/'score.json'
 write_vcf(truth,['1\t10\t.\tA\tG\t.\tPASS\t.','1\t20\t.\tA\tAT\t.\tPASS\t.'])
 write_vcf(base,['1\t10\t.\tA\tG\t.\tSomatic\t.'])
 write_vcf(calls,['1\t10\t.\tA\tG\t.\tSomatic\t.','1\t30\t.\tC\tT\t.\tSomatic\t.'])
 subprocess.run([sys.executable,str(ROOT/'examples/seqc2/scripts/score_label_artifact.py'),'--truth',str(truth),'--calls',str(calls),'--baseline',str(base),'--out',str(out)],check=True)
 scored=json.loads(out.read_text()); rows={r['variant_type']:r for r in scored['rows']}; assert rows['SNV']['gained_fp']==1; assert rows['indel']['fn']==1; assert scored['selectors']=={'truth':'all_records','calls':'Somatic','baseline':'PASS'}; assert scored['baseline_sha256']


def test_scorer_accepts_distinct_baseline_selector(tmp_path):
 truth=tmp_path/'truth.vcf'; base=tmp_path/'base.vcf'; calls=tmp_path/'calls.vcf'; out=tmp_path/'score.json'
 write_vcf(truth,['1\t10\t.\tA\tG\t.\tPASS\t.'])
 write_vcf(base,['1\t10\t.\tA\tG\t.\tPASS\t.'])
 write_vcf(calls,['1\t10\t.\tA\tG\t.\tSomatic\t.'])
 subprocess.run([sys.executable,str(ROOT/'examples/seqc2/scripts/score_label_artifact.py'),'--truth',str(truth),'--calls',str(calls),'--baseline',str(base),'--baseline-label','PASS','--out',str(out)],check=True)
 result=json.loads(out.read_text()); assert result['selectors']['baseline']=='PASS'; assert result['baseline_calls_sha256']


def test_label_builder_deduplicates_sorts_and_indexes_output(tmp_path):
    ds=tmp_path/'ds.vcf'; rna=tmp_path/'rna.vcf'; ver=tmp_path/'verification.json'; out=tmp_path/'labels.vcf.gz'
    write_vcf(ds,['1\t30\t.\tG\tA\t.\tPASS\t.','1\t10\t.\tA\tG\t.\tPASS\t.'])
    write_vcf(rna,['1\t30\t.\tG\tA\t.\t.\t.','1\t20\t.\tC\tT\t.\t.\t.'])
    ver.write_text(json.dumps({'results':[{'chrom':'1','pos':20,'ref':'C','alt':'T','status':'confirmed','tumor_alt':5,'normal_alt':0}]}))
    subprocess.run([sys.executable,str(ROOT/'bin/build_deepsomatic_labels.py'),'--deepsomatic-vcf',str(ds),'--rna-nominations',str(rna),'--verification-json',str(ver),'--out',str(out),'--index'],check=True)
    import gzip
    records=[line for line in gzip.open(out,'rt') if line and not line.startswith('#')]
    assert [line.split('\t')[1] for line in records] == ['10','20','30']
    assert sum(line.split('\t')[1]=='30' for line in records)==1
    assert (Path(str(out)+'.tbi')).is_file()


def test_scorer_binds_provenance_manifest(tmp_path):
 truth=tmp_path/'truth.vcf'; calls=tmp_path/'calls.vcf'; provenance=tmp_path/'provenance.json'; out=tmp_path/'score.json'
 write_vcf(truth,['1\t10\t.\tA\tG\t.\tPASS\t.']); write_vcf(calls,['1\t10\t.\tA\tG\t.\tSomatic\t.'])
 import hashlib
 provenance.write_text(json.dumps({'schema':'seqc2-artifact-provenance.v1','artifact':{'sha256':hashlib.sha256(calls.read_bytes()).hexdigest()},'stage':'final_second_rescue'}))
 subprocess.run([sys.executable,str(ROOT/'examples/seqc2/scripts/score_label_artifact.py'),'--truth',str(truth),'--calls',str(calls),'--provenance',str(provenance),'--out',str(out)],check=True)
 result=json.loads(out.read_text()); assert result['provenance_sha256']; assert result['provenance']['stage']=='final_second_rescue'


def test_label_builder_excludes_unknown_deepsomatic_filter(tmp_path):
    ds=tmp_path/'ds.vcf'; out=tmp_path/'labels.vcf'
    write_vcf(ds,['1\t10\t.\tA\tG\t.\t.\t.','1\t15\t.\tA\tT\t.\tLowQual;PASS\t.','1\t20\t.\tC\tT\t.\tLowQual\t.'])
    subprocess.run([sys.executable,str(ROOT/'bin/build_deepsomatic_labels.py'),'--deepsomatic-vcf',str(ds),'--out',str(out)],check=True)
    records=[line for line in out.read_text().splitlines() if line and not line.startswith('#')]
    assert records == []


def test_scorer_rejects_unrelated_provenance(tmp_path):
    truth=tmp_path/'truth.vcf'; calls=tmp_path/'calls.vcf'; provenance=tmp_path/'provenance.json'; out=tmp_path/'score.json'
    write_vcf(truth,['1\t10\t.\tA\tG\t.\tPASS\t.']); write_vcf(calls,['1\t10\t.\tA\tG\t.\tSomatic\t.'])
    provenance.write_text(json.dumps({'schema':'seqc2-artifact-provenance.v1','artifact':{'sha256':'wrong'}}))
    result=subprocess.run([sys.executable,str(ROOT/'examples/seqc2/scripts/score_label_artifact.py'),'--truth',str(truth),'--calls',str(calls),'--provenance',str(provenance),'--out',str(out)],capture_output=True,text=True)
    assert result.returncode != 0 and 'sha256' in result.stderr.lower()


def test_scorer_requires_provenance_for_declared_stage(tmp_path):
    truth=tmp_path/'truth.vcf'; calls=tmp_path/'calls.vcf'; out=tmp_path/'score.json'
    write_vcf(truth,['1\t10\t.\tA\tG\t.\tPASS\t.']); write_vcf(calls,['1\t10\t.\tA\tG\t.\tSomatic\t.'])
    result=subprocess.run([sys.executable,str(ROOT/'examples/seqc2/scripts/score_label_artifact.py'),'--truth',str(truth),'--calls',str(calls),'--stage','final_second_rescue','--out',str(out)],capture_output=True,text=True)
    assert result.returncode != 0 and 'provenance' in result.stderr.lower()
