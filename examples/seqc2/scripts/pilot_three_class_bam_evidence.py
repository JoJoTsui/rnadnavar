#!/usr/bin/env python3
"""Deterministic, truth-blind SNP candidate pilot; read evidence is not truth.

Indels require allele/haplotype-aware validation and are explicitly not assessed
by this SNP pileup screen. Never changes labels or invokes mapping/calling.
"""
import argparse
from collections import Counter
import csv
import hashlib
import heapq
import json
from pathlib import Path

import pysam
from validate_refined_native_integration import digest

FLAG_FILTER = 0xF04  # unmapped, secondary, QC-fail, duplicate, supplementary


def select_sites(path, limit):
    groups = {label: [] for label in ("Germline", "Reference")}
    populations = Counter()
    with pysam.VariantFile(str(path)) as reader:
        for record in reader:
            label = next(iter(record.filter), "")
            if label not in groups:
                continue
            if len(record.alts or []) != 1 or len(record.ref) != 1 or len(record.alts[0]) != 1:
                populations[label + ":not_SNP"] += 1
                continue
            site = (record.contig, record.pos, record.ref, record.alts[0])
            if set(site[2] + site[3]) - set("ACGT"):
                continue
            populations[label + ":SNP"] += 1
            rank = int(hashlib.sha256((':'.join(map(str, site))).encode()).hexdigest(), 16)
            item = (-rank, site)
            heap = groups[label]
            if len(heap) < limit:
                heapq.heappush(heap, item)
            elif item > heap[0]:
                heapq.heapreplace(heap, item)
    return {label: [site for _,site in sorted(heap, reverse=True)] for label,heap in groups.items()}, dict(populations)


def evidence(bam, fasta, site):
    chrom, pos, ref, alt = site
    if fasta.fetch(chrom, pos-1, pos).upper() != ref:
        raise ValueError("Reference mismatch at " + str(site))
    if bam.get_reference_length(chrom) != fasta.get_reference_length(chrom):
        raise ValueError("BAM/reference contig length mismatch: " + chrom)
    counts = Counter()
    for col in bam.pileup(chrom, pos-1, pos, truncate=True, stepper="samtools",
                         fastafile=fasta, min_mapping_quality=20, min_base_quality=20,
                         ignore_overlaps=True, ignore_orphans=True, max_depth=8000,
                         compute_baq=True, flag_filter=FLAG_FILTER):
        if col.reference_pos != pos-1:
            continue
        if col.nsegments >= 8000:
            return {"status": "depth_cap_inconclusive"}
        for base in col.get_query_sequences():
            if base.upper() in ("A", "C", "G", "T"):
                counts[base.upper()] += 1
    return {"depth": sum(counts.values()), "ref": counts[ref], "alt": counts[alt],
            "other": sum(counts.values()) - counts[ref] - counts[alt]}


def assess(label, normal, tumor):
    if "status" in normal or "status" in tumor:
        return "inconclusive_depth_cap"
    if label == "Germline":
        if normal["depth"] >= 20 and normal["alt"] >= 5 and normal["alt"] / normal["depth"] >= .2:
            return "normal_read_corroborated_not_truth"
        if normal["depth"] >= 60 and normal["alt"] == 0:
            return "normal_reference_conflict"
        return "inconclusive_normal_evidence"
    if any(x["depth"] >= 10 and x["alt"] >= 3 and x["alt"] / x["depth"] > .05 for x in (normal, tumor)):
        return "alt_read_conflict"
    if all(x["depth"] >= 60 and x["alt"] == 0 for x in (normal, tumor)):
        return "paired_zero_alt_corroborated_not_truth"
    return "inconclusive_paired_evidence"


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--vcf", type=Path, required=True)
    ap.add_argument("--samplesheet", type=Path, required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--per-class", type=int, default=32)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    if args.per_class < 1 or args.out.exists():
        raise ValueError("Require positive limit and new output path")
    with args.samplesheet.open() as f:
        rows = list(csv.DictReader(f))
    pair = {}
    for status in ("0", "1"):
        matches = [r for r in rows if r["status"] == status]
        if len(matches) != 1:
            raise ValueError("Require one caller-ready DNA BAM for each status")
        pair[status] = matches[0]
    selected, populations = select_sites(args.vcf, args.per_class)
    report = {"scope": __doc__, "selection": "smallest sha256(CHROM:POS:REF:ALT), separately per class, no truth selection",
              "vcf": str(args.vcf.resolve()), "vcf_sha256": digest(args.vcf),
              "samplesheet": str(args.samplesheet.resolve()), "reference": args.fasta,
              "thresholds": {"MAPQ":20,"BQ":20,"BAQ":True,"max_depth":8000,
                             "flag_filter":FLAG_FILTER,"ignore_overlaps":True,"ignore_orphans":True},
              "bam_paths": {s:r["bam"] for s,r in pair.items()}, "populations": populations,
              "results": [], "training_approved":False}
    with pysam.AlignmentFile(pair["0"]["bam"], index_filename=pair["0"]["bai"]) as normal, \
         pysam.AlignmentFile(pair["1"]["bam"], index_filename=pair["1"]["bai"]) as tumor, \
         pysam.FastaFile(args.fasta) as fasta:
        for label, sites in selected.items():
            for site in sites:
                n, t = evidence(normal, fasta, site), evidence(tumor, fasta, site)
                report["results"].append({"class":label, "site":site, "normal":n, "tumor":t,
                                          "outcome":assess(label,n,t)})
    report["outcomes"] = dict(Counter(r["class"]+":"+r["outcome"] for r in report["results"]))
    with args.out.open("x") as f:
        json.dump(report, f, indent=2)
        f.write("\n")
    print(json.dumps(report["outcomes"]), flush=True)


if __name__ == "__main__":
    main()
