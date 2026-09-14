#!/usr/bin/env python3
"""Read-only DNA pileup audit for the frozen rescue evidence sites.

This is descriptive evidence, not a caller and not a policy threshold tuner.
Unavailable or malformed pileups are recorded as inconclusive.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess


def pileup(samtools, reference, bam, chrom, pos):
    region = f"{chrom}:{pos}-{pos}"
    command = [samtools, "mpileup", "-aa", "-f", str(reference), "-r", region, str(bam)]
    try:
        result = subprocess.run(command, text=True, capture_output=True, timeout=30)
    except subprocess.TimeoutExpired:
        return {"status": "inconclusive", "stderr": "pileup timed out", "command": command}
    if result.returncode:
        return {"status": "inconclusive", "stderr": result.stderr.strip(), "command": command}
    rows = [line.split("\t") for line in result.stdout.splitlines() if line.strip()]
    if not rows or len(rows[0]) < 6:
        return {"status": "inconclusive", "stderr": "empty pileup", "command": command}
    row = rows[0]
    return {"status": "observed", "chrom": row[0], "pos": int(row[1]),
            "ref": row[2], "depth": int(row[3]), "bases": row[4],
            "base_qualities": row[5], "command": command}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--max-sites", type=int, default=0,
                        help="limit sites per dataset (0 audits all scored sites)")
    args = parser.parse_args()
    if args.outdir.exists():
        parser.error(f"output directory already exists: {args.outdir}; choose a fresh namespace")
    repo = Path(__file__).resolve().parents[3]
    bundle = repo / "examples/seqc2/verified/20260914"
    evidence_path = repo / "examples/seqc2/comparison/rescue_fp_investigation_20260914/evidence.json"
    evidence = json.loads(evidence_path.read_text())
    reference = Path("/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/bio_db/references/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta")
    samtools = shutil.which("samtools") or "/t9k/mnt/joey/bio_gizmo/samtools"
    if not reference.exists() or not Path(samtools).exists():
        raise SystemExit("reference or samtools unavailable; no evidence was changed")
    report = {"status": "descriptive_not_validation", "sources": {}, "source_stats": {}, "sites": []}
    report["sources"][str(evidence_path.resolve())] = hashlib.sha256(evidence_path.read_bytes()).hexdigest()
    report["source_stats"][str(reference)] = {"size": reference.stat().st_size, "mtime_ns": reference.stat().st_mtime_ns}
    for dataset in ("wes_ll", "wgs_il"):
        alignments = bundle / dataset / "alignments"
        files = {name: (alignments / filename).resolve(strict=True) for name, filename in
                 (("dna_tumor", "dna_tumor.bam"), ("dna_normal", "dna_normal.bam"),
                  ("rna_tumor", "rna_tumor.cram"), ("rna_realign", "rna_realign.cram"))}
        for path in files.values():
            report["source_stats"][str(path)] = {"size": path.stat().st_size, "mtime_ns": path.stat().st_mtime_ns}
        sites = [row for row in evidence["sites"] if row["dataset"] == dataset and row["truth_status"] in {"TP", "FP"}]
        if args.max_sites:
            sites = sites[:args.max_sites]
        for row in sites:
            result = {"dataset": dataset, "allele": [row[k] for k in ("chrom", "pos", "ref", "alt")],
                      "truth_status": row["truth_status"]}
            for name, path in files.items():
                result[name] = pileup(samtools, reference, path, row["chrom"], row["pos"])
            report["sites"].append(result)
    report["summary"] = {dataset: {status: sum(1 for row in report["sites"] if row["dataset"] == dataset
                                                 and row["dna_tumor"]["status"] == status)
                                   for status in ("observed", "inconclusive")}
                         for dataset in ("wes_ll", "wgs_il")}
    args.outdir.mkdir(parents=True, exist_ok=False)
    (args.outdir / "read_evidence.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report["summary"], indent=2))


if __name__ == "__main__":
    main()
