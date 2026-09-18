#!/usr/bin/env python3
"""Corroborate preselected SNP candidates against an orthogonal normal gVCF.

HG008 N-P is a different normal tissue from the N-D caller input. This is not
a validated germline benchmark, and normal-only evidence cannot approve a
paired Reference label. Missing calls are never treated as reference truth.
"""
import argparse
from collections import Counter
import json
from pathlib import Path

import pysam
from validate_refined_native_integration import digest


def outcome(site, records):
    chrom, pos, ref, alt = site
    evidence = []
    for r in records:
        sample = r.samples[next(iter(r.samples))]
        gt, gq = sample.get("GT"), sample.get("GQ")
        dp = sample.get("MIN_DP")
        if dp is None:
            dp = sample.get("DP")
        entry = dict(pos=r.pos, stop=r.stop, ref=r.ref, alts=r.alts, GT=gt, GQ=gq, DP=dp)
        evidence.append(entry)
        if set(r.filter) - {"PASS", "."} or not gt or any(g is None for g in gt):
            continue
        if gq is None or gq < 30 or dp is None or dp < 20:
            continue
        if r.pos == pos and r.ref == ref and alt in (r.alts or []):
            index = r.alleles.index(alt)
            if index in gt:
                return "normal_alt_corroboration", evidence
        if all(g == 0 for g in gt) and r.start <= pos-1 < r.stop:
            return "normal_reference_evidence", evidence
    return "inconclusive", evidence


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pilot", type=Path, required=True)
    ap.add_argument("--normal-gvcf", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    if args.out.exists():
        raise ValueError("Refuse overwrite")
    pilot = json.loads(args.pilot.read_text())
    result = {"scope": __doc__, "pilot": str(args.pilot.resolve()), "pilot_sha256": digest(args.pilot),
              "normal_gvcf": str(args.normal_gvcf.resolve()), "normal_gvcf_sha256": digest(args.normal_gvcf),
              "training_approved": False, "results": []}
    with pysam.VariantFile(str(args.normal_gvcf)) as vcf:
        result["sample_names"] = list(vcf.header.samples)
        if len(vcf.header.samples) != 1:
            raise ValueError("Require one orthogonal normal sample")
        for row in pilot["results"]:
            chrom,pos,ref,alt = row["site"]
            if len(ref) != 1 or len(alt) != 1:
                raise ValueError("Only SNP pilot supported")
            status, records = outcome(row["site"], vcf.fetch(chrom, pos-1, pos))
            result["results"].append(dict(site=row["site"], label=row["class"], outcome=status, evidence=records))
    result["counts"] = dict(Counter(r["label"] + ":" + r["outcome"] for r in result["results"]))
    with args.out.open("x") as f:
        json.dump(result, f, indent=2)
        f.write("\n")
    print(json.dumps(result["counts"]))


if __name__ == "__main__":
    main()
