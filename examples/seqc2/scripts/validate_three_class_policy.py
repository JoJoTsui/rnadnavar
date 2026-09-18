#!/usr/bin/env python3
"""Frozen three-class evaluation from existing normalized callers, never Nextflow.

Somatic truth only measures somatic performance and negative-label collisions;
it cannot establish Germline/Reference precision or biological training approval.
"""
import argparse
from collections import Counter
import json
from pathlib import Path
import subprocess
import sys

import pysam

from aggregate_benchmark import parse_metrics_json
from validate_refined_native_integration import digest

ROOT = Path(__file__).resolve().parents[3]
LABELS = ("Somatic", "Germline", "Reference")


def class_queries(source, folder):
    counts, reasons = Counter(), Counter()
    with pysam.VariantFile(str(source)) as reader:
        header = pysam.VariantHeader()
        for name, contig in reader.header.contigs.items():
            header.contigs.add(name, length=contig.length)
        writers = {label: pysam.VariantFile(str(folder / f"{label}.query.vcf.gz"), "wz", header=header)
                   for label in LABELS}
        try:
            for record in reader:
                labels = set(record.filter)
                label = next(iter(labels)) if len(labels) == 1 else "MULTIPLE"
                counts[label] += 1
                trace = record.info.get("CLASSIFICATION_RATIONALE", "")
                if isinstance(trace, tuple):
                    trace = "|".join(trace)
                for token in trace.split("|"):
                    if token.startswith("decision:"):
                        reasons[token] += 1
                if label in writers:
                    copied = writers[label].new_record(contig=record.contig, start=record.start, alleles=record.alleles)
                    copied.filter.add("PASS")
                    writers[label].write(copied)
        finally:
            for writer in writers.values():
                writer.close()
    for label in LABELS:
        pysam.tabix_index(str(folder / f"{label}.query.vcf.gz"), preset="vcf")
    return dict(counts), dict(reasons)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--frozen-validation", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    args = ap.parse_args()
    frozen = json.loads(args.frozen_validation.read_text())
    config = frozen["manifest"]
    old_command = next(c for c in frozen["commands"] if "--input_dir" in c)
    inputs = Path(old_command[old_command.index("--input_dir") + 1])
    sources = [inputs / f"sample.{c}.vcf.gz" for c in ("deepsomatic", "mutect2", "strelka")]
    sources += [Path(config[k]) for k in ("truth", "hc", "fasta")]
    sources += [Path(p) for p in config["targets"].values()]
    if not all(p.is_file() for p in sources):
        raise ValueError("Missing frozen validation source")
    out = args.outdir.resolve()
    out.mkdir(parents=True, exist_ok=False)
    code = [ROOT / "bin/run_consensus_vcf.py", Path(__file__).resolve(),
            *sorted((ROOT / "bin/vcf_utils").glob("*.py")), *sorted((ROOT / "bin/common").glob("*.py"))]
    report = {"status": "running", "training_approved": False, "scope": __doc__,
              "git_head": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
              "frozen_validation": str(args.frozen_validation.resolve()),
              "frozen_validation_sha256": digest(args.frozen_validation),
              "manifest": config, "sources": {}, "code": {}, "commands": [], "metrics": {},
              "baseline_somatic_metrics": {k:v for k,v in frozen["metrics"].items() if k.startswith("refined_consensus/")}}

    def save():
        (out / "validation.json").write_text(json.dumps(report, indent=2) + "\n")

    def run(command, logfile):
        report["commands"].append(command)
        save()
        with (out / logfile).open("x") as log:
            subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)

    save()
    try:
        report["sources"] = {str(p): digest(p) for p in sources}
        report["code"] = {str(p): digest(p) for p in code}
        run([sys.executable, str(ROOT / "bin/run_consensus_vcf.py"), "--input_dir", str(inputs),
             "--expected_callers", "deepsomatic,mutect2,strelka", "--experimental-refined-native",
             "--experimental-three-class", "--out_prefix", str(out / "three_class")], "consensus.log")
        output = out / "three_class.vcf.gz"
        report["class_counts"], report["decisions"] = class_queries(output, out)
        report["output_sha256"] = digest(output)
        save()
        print(config["dataset"], report["class_counts"], flush=True)
        for label in LABELS:
            for domain, target in config["targets"].items():
                prefix = out / f"{label}.{domain}"
                run(["micromamba", "run", "-n", "happy", "som.py", config["truth"],
                     str(out / f"{label}.query.vcf.gz"), "-R", config["hc"], "-T", target,
                     "-r", config["fasta"], "-N", "-o", str(prefix)], f"{label}.{domain}.log")
                metrics = parse_metrics_json(Path(str(prefix) + ".metrics.json"))
                report["metrics"][f"{label}/{domain}"] = {
                    "interpretation": "somatic_performance" if label == "Somatic" else
                    "TP_is_known_somatic_truth_recovered_by_negative_labels_NOT_negative_class_accuracy",
                    "values": metrics}
                save()
                print(config["dataset"], label, domain, metrics["records"], flush=True)
        report["sources_unchanged"] = all(digest(Path(p)) == sha for p, sha in report["sources"].items())
        report["code_unchanged"] = all(digest(Path(p)) == sha for p, sha in report["code"].items())
        if not report["sources_unchanged"] or not report["code_unchanged"]:
            raise RuntimeError("Input/code changed during evaluation")
        report["status"] = "complete_somatic_truth_screen_not_training_approved"
        save()
    except Exception as error:
        report.update(status="failed", error=str(error))
        save()
        raise


if __name__ == "__main__":
    main()
