"""Native class semantics, including DeepSomatic's recoded Germline GT."""
from copy import deepcopy
import gzip
import importlib.util
from pathlib import Path
import sqlite3
import subprocess
import sys

import pysam
import pytest
from vcf_utils.three_class_policy import evaluate, classify
from vcf_utils.classification import compute_unified_classification_consensus
from vcf_utils.aggregation import read_variants_from_vcf

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("rescue_three_test", ROOT/"bin/apply_refined_rescue.py")
rescue = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rescue)

def data(label):
    ad = [20,20] if label == "Germline" else [60,0]
    pl = [40,0,60] if label == "Germline" else [0,40,60]
    return dict(REF="A",ALT="T",is_snv=True,
        callers=["deepsomatic","mutect2"],
        filters_original=["GERMLINE" if label == "Germline" else "RefCall","germline"],
        native_evidence={"deepsomatic":{"tumor_GT":[0,0],"tumor_GQ":[40],
            "tumor_PL":pl,"tumor_DP":[sum(ad)],"tumor_AD":ad}})

@pytest.mark.parametrize("label",["Germline","Reference"])
def test_native_classes_do_not_require_normal_gq_or_generic_gt(label):
    d=data(label)
    assert evaluate(d)[0] == label
    assert classify(d,False)[0] == label
    assert classify(d,True)[0] == "NoConsensus"
    del d["native_evidence"]["deepsomatic"]["tumor_GQ"]
    assert evaluate(d)[0] is None

@pytest.mark.parametrize("filt",["germline","haplotype","panel_of_normals","contamination","possible_numt","LowEVS"])
def test_failed_somatic_is_not_a_negative_class(filt):
    d=data("Germline")
    d["callers"]=["mutect2","strelka"]
    d["filters_original"]=[filt,filt]
    d["native_evidence"].pop("deepsomatic")
    assert classify(d,False)[0] == "NoConsensus"

@pytest.mark.parametrize("quality",[None,[float("nan")],[-1],[29],[30,40]])
def test_bad_native_confidence_abstains(quality):
    d=data("Germline")
    d["native_evidence"]["deepsomatic"]["tumor_GQ"]=quality
    assert evaluate(d)[0] is None

def test_native_filter_must_be_exact_and_pl_consistent():
    d=data("Germline");d["filters_original"][0]="GERMLINE;LowQual"
    assert evaluate(d)[0] is None
    d=data("Germline");d["native_evidence"]["deepsomatic"]["tumor_PL"]=[0,40,60]
    assert evaluate(d)[0] is None

def test_normal_read_counts_corroborate_without_normal_gq():
    d=data("Germline")
    d["native_evidence"]["mutect2"]={"normal_AD":[20,20],"normal_DP":[40]}
    assert "corroborated:mutect2" in evaluate(d)[1]
    d["native_evidence"]["mutect2"]={"normal_AD":[60,0],"normal_DP":[60]}
    assert classify(d,True)[0] == "NoConsensus"
    d=data("Reference")
    d["native_evidence"]["strelka"]={"tumor_AD":[17,3],"tumor_DP":[20]}
    assert classify(d,False)[0] == "NoConsensus"

def test_native_refcall_has_its_own_alt_limit():
    d=data("Reference")
    d["native_evidence"]["deepsomatic"]["tumor_AD"]=[57,3]
    assert classify(d,False)[0]=="NoConsensus"
    d=data("Reference")
    assert d["native_evidence"]["deepsomatic"]["tumor_AD"][1] == 0
    assert classify(d,False)[0] == "Reference"

def test_opt_in_only():
    d=data("Germline")
    d.update(refined_native_enabled=True,refined_native_branch=None,
             refined_native_trace="rule:seqc2_refined_v2",three_class_enabled=True)
    assert compute_unified_classification_consensus(d,2,2) == "Germline"
    d["three_class_enabled"]=False
    assert compute_unified_classification_consensus(d,2,2) == "NoConsensus"

def test_rescue_conflicts_and_inheritance():
    assert rescue.transition("Somatic","Germline",False,"none",None)[0]=="Somatic"
    assert rescue.transition("Somatic","Germline",False,"none",None,True)[0]=="NoConsensus"
    assert rescue.transition("Somatic","Somatic",False,"none",None,True,"common_population_af")[0]=="NoConsensus"
    assert rescue.transition("NoConsensus","Reference",False,"none",None,True)[0]=="NoConsensus"
    assert rescue.transition("Reference","Reference",False,"none",None,True)[0]=="Reference"

def test_strelka_tier1_read_evidence(tmp_path):
    p=tmp_path/"strelka.vcf"
    p.write_text('##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n'
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="depth">\n'
        '##FORMAT=<ID=AU,Number=2,Type=Integer,Description="tier counts">\n'
        '##FORMAT=<ID=TU,Number=2,Type=Integer,Description="tier counts">\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\tTUMOR\n'
        'chr1\t20\t.\tA\tT\t.\tPASS\t.\tDP:AU:TU\t40:20,25:20,30\t40:21,25:19,30\n')
    old=next(iter(read_variants_from_vcf(p,"strelka",refined_native=True).values()))
    new=next(iter(read_variants_from_vcf(p,"strelka",refined_native=True,three_class=True).values()))
    assert "normal_AD" not in old["native_evidence"]
    assert new["native_evidence"]["normal_AD"] == [20,20]
    assert new["native_evidence"]["tumor_AD"] == [21,19]

def test_consensus_cli_three_classes(tmp_path):
    inputs=tmp_path/"inputs";inputs.mkdir()
    header=('##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n'
        '##FILTER=<ID=RefCall,Description="reference">\n##FILTER=<ID=GERMLINE,Description="non somatic">\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="genotype">\n'
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="depth">\n'
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="depth">\n'
        '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="quality">\n')
    for caller in ("mutect2","strelka","deepsomatic"):
        paired=caller!="deepsomatic"
        cols='#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+("NORMAL\tTUMOR" if paired else "TUMOR")+"\n"
        rows=[]
        for pos,ad,dp,gq in [(20,"20,20",40,40),(30,"60,0",60,40),(40,"60,0",60,".")]:
            filt=("GERMLINE" if pos==20 else "RefCall") if not paired else "PASS"
            fmt=f"0/0:{ad}:{dp}:{gq}"
            rows.append(f"chr1\t{pos}\t.\tA\tT\t0\t{filt}\t.\tGT:AD:DP:GQ\t"+(fmt+"\t"+fmt if paired else fmt)+"\n")
        (inputs/f"sample.{caller}.vcf").write_text(header+cols+"".join(rows))
    result=subprocess.run([sys.executable,str(ROOT/"bin/run_consensus_vcf.py"),"--input_dir",str(inputs),
        "--out_prefix",str(tmp_path/"new"),"--experimental-refined-native","--experimental-three-class",
        "--expected_callers","mutect2,strelka,deepsomatic"],capture_output=True,text=True)
    assert result.returncode==0,result.stdout+result.stderr
    with gzip.open(tmp_path/"new.vcf.gz","rt") as f:rows=[r.split("\t") for r in f if not r.startswith("#")]
    assert {int(r[1]):r[6] for r in rows}=={20:"Germline",30:"Reference",40:"NoConsensus"}

def test_rescue_three_class_conflict_output(tmp_path):
    header=('##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n'
        '##INFO=<ID=CLASSIFICATION_RATIONALE,Number=1,Type=String,Description="trace">\n'
        '##FILTER=<ID=Somatic,Description="somatic">\n##FILTER=<ID=Germline,Description="germline">\n')
    cols='#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
    dna,old=tmp_path/"dna.vcf",tmp_path/"old.vcf"
    dna.write_text(header+cols+'chr1\t10\t.\tA\tT\t.\tSomatic\tCLASSIFICATION_RATIONALE=three_class_policy:native_three_class_v1|class:Somatic\n')
    old.write_text(header+cols+'chr1\t10\t.\tA\tT\t.\tGermline\t.\n')
    out=tmp_path/"rescue.vcf.gz"
    with sqlite3.connect(":memory:") as db:
        db.execute("CREATE TABLE variants(chrom TEXT,pos INTEGER,ref TEXT,alt TEXT,dna BLOB,rescue BLOB,dna_votes TEXT DEFAULT '',rna_votes TEXT DEFAULT '',PRIMARY KEY(chrom,pos,ref,alt))")
        h=rescue.stage(db,dna,"dna");h.merge(rescue.stage(db,old,"rescue"))
        result=rescue.write_output(db,h,out,"realignment",three_class=True)
    assert result["NoConsensus"]==1
    with pysam.VariantFile(out) as reader:
        r=next(reader)
        assert r.info["GATE_POLICY"]=="native_three_class_gate_v1"
        assert r.info["RESCUE_PROMOTED"]=="NO"
        assert r.info["PASSES_CONSENSUS_DNA"]=="NO"
        assert r.info["UNIFIED_FILTER"]=="NoConsensus"
