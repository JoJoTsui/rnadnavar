#!/usr/bin/env python3
"""Read-only bounded evidence inventory; prefix sampling is not accuracy validation."""
import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
import csv
from itertools import islice
import json
from pathlib import Path
import sys
import pysam

ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/"bin"))
from vcf_utils.aggregation import resolve_tumor_sample_index, _normal_sample_from_header
from vcf_utils.three_class_policy import native_negative, POLICY

def inspect(item):
    sid,caller,path,limit=item
    with pysam.VariantFile(path) as reader:
        names=list(reader.header.samples)
        if len(names) not in (1,2):
            raise ValueError(f"Unsupported sample count: {path}")
        tumor=resolve_tumor_sample_index(names,caller,_normal_sample_from_header(str(reader.header)))
        normal=names[1-tumor] if len(names)==2 else None
        counts=Counter()
        for r in islice(reader,limit):
            counts["records_scanned"]+=1
            if caller == "deepsomatic":
                label = ";".join(r.filter.keys())
                counts["native_filter:" + label] += 1
                sample = r.samples[names[tumor]]
                evidence = {"tumor_" + key: sample.get(key) for key in ("GQ", "DP", "AD", "PL")}
                nominated, reason = native_negative({
                    "REF": r.ref, "ALT": ",".join(r.alts or []),
                    "callers": [caller], "filters_original": [label],
                    "native_evidence": {caller: evidence},
                })
                counts["native_nomination:" + (nominated or "none")] += 1
                counts["native_reason:" + reason] += 1
                if label == "GERMLINE" and sample.get("GT") == (0, 0):
                    counts["germline_recoded_GT_00"] += 1
            if normal:
                sample=r.samples[normal]
                for key in ("GT","GQ","DP","AD"):
                    val=sample.get(key)
                    present=val is not None and (not isinstance(val,tuple) or all(x is not None for x in val))
                    counts["normal_"+key+"_present"]+=present
        return dict(sample_id=sid,caller=caller,path=path,normal_sample=normal,
                    header_formats=list(reader.header.formats),counts=dict(counts))

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest",type=Path,default=ROOT/"examples/seq2neo/data/processed/sample_manifest.tsv")
    ap.add_argument("--records",type=int,default=200)
    ap.add_argument("--out",type=Path,required=True)
    args=ap.parse_args()
    if args.records<1:raise ValueError("records must be positive")
    with args.manifest.open() as f: rows=list(csv.DictReader(f,delimiter="\t"))
    jobs=[(s["sample_id"],c,s["caller_dna_"+c],args.records)
          for s in rows for c in ("mutect2","strelka","deepsomatic")]
    with ThreadPoolExecutor(max_workers=4) as pool: reports=list(pool.map(inspect,jobs))
    totals={}
    for caller in ("mutect2","strelka","deepsomatic"):
        subset=[r for r in reports if r["caller"]==caller]
        counts=Counter()
        for r in subset:counts.update(r["counts"])
        totals[caller]=dict(files=len(subset),files_with_normal=sum(r["normal_sample"] is not None for r in subset),
            files_with_observed_normal_GQ=sum(r["counts"].get("normal_GQ_present",0)>0 for r in subset),counts=dict(counts))
    args.out.parent.mkdir(parents=True,exist_ok=True)
    with args.out.open("x") as f:
        json.dump(dict(scope="first records per DNA caller file; not random or exhaustive",
                       native_policy=POLICY, native_scope="DeepSomatic nomination only; not cross-caller conflict checks or accuracy",
                       limit=args.records,totals=totals,files=reports),f,indent=2)
        f.write("\n")
    print(json.dumps(totals,indent=2))

if __name__=="__main__": main()
