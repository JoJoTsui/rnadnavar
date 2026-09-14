#!/usr/bin/env python3
"""Validate effective ingress and re-consensus configuration contracts.

This is a static, read-only check. It does not claim live Nextflow scheduling
or cache reuse; those require an execution trace and are reported separately.
"""
import argparse
import json
from pathlib import Path
import re
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from validation_common import check, overall, write_report


def assignment(text, name):
    match = re.search(rf"^\s*{re.escape(name)}\s*(?:=|:)\s*([^\n]+)", text, re.MULTILINE)
    if not match:
        return None
    value = match.group(1).split("//", 1)[0].split("#", 1)[0].strip()
    return value.strip("\"'")



def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    if args.outdir.exists():
        parser.error(f"output directory already exists: {args.outdir}; choose a fresh namespace")
    repo = Path(__file__).resolve().parents[3]
    checks = []
    configs = {
        "dna_bam": repo / "examples/seqc2/seqc2.shared.config",
        "hybrid_wes": repo / "examples/seqc2/seqc2.hybrid.shared.config",
        "hybrid_wgs": repo / "examples/seqc2/hybrid/seqc2.hybrid.wgs.config",
        "rerun": repo / "examples/seq2neo/config/rerun.yaml",
        "cohort_rerun": repo / "examples/seq2neo/config/rerun_cohort.yaml",
    }
    parsed = {}
    for name, path in configs.items():
        text = path.read_text()
        if path.suffix == ".config":
            parsed[name] = {"step": assignment(text, "step"), "tools": assignment(text, "tools"),
                            "input": assignment(text, "input")}
        else:
            parsed[name] = {key: assignment(text, key) for key in ("step", "tools", "outdir_base", "csv_dir", "state_file", "checksum_dir")}
    checks += [
        check("dna-bam-step", parsed["dna_bam"]["step"] == "variant_calling", str(parsed["dna_bam"])),
        check("dna-bam-input", "bam,bai" in (configs["dna_bam"].read_text()).lower(), str(parsed["dna_bam"])),
        check("wes-hybrid-tools", all(token in (parsed["hybrid_wes"]["tools"] or "") for token in ("consensus", "rescue")), str(parsed["hybrid_wes"])),
        check("wgs-hybrid-tools", all(token in (parsed["hybrid_wgs"]["tools"] or "") for token in ("consensus", "rescue", "realignment")), str(parsed["hybrid_wgs"])),
        check("rerun-step", parsed["rerun"]["step"] == "consensus" and parsed["cohort_rerun"]["step"] == "consensus", "both reruns enter at consensus"),
        check("rerun-no-callers", all(not re.search(r"(?:mutect2|strelka|deepsomatic)", parsed[name]["tools"] or "", re.I) for name in ("rerun", "cohort_rerun")), "caller tools absent"),
        check("cohort-vep-contract", "vep" not in (parsed["cohort_rerun"]["tools"] or "").lower(), str(parsed["cohort_rerun"])),
    ]
    for name in ("rerun", "cohort_rerun"):
        values = [parsed[name][key] for key in ("outdir_base", "csv_dir", "state_file", "checksum_dir")]
        checks.append(check(f"{name}-isolated-paths", all(values) and len(set(values)) == len(values), str(values)))
        checks.append(check(f"{name}-new-output", "output_reconsensus" in (parsed[name]["outdir_base"] or ""), parsed[name]["outdir_base"]))
    report = {"status": overall(checks), "static_only": True, "checks": checks,
              "configs": {name: str(path) for name, path in configs.items()}}
    write_report(args.outdir, "ingress_contract.json", report)
    print(json.dumps({"status": report["status"], "checks": len(checks), "failures": [x for x in checks if x["status"] != "pass"]}, indent=2))
    return 0 if report["status"] == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())
