#!/usr/bin/env python3
"""Static reachability contract for supported input stages.

Checks the DSL conditions that decide whether mapping, caller processing,
consensus, rescue, and second-round realignment are reachable. This does not
execute Nextflow or claim cache reuse.
"""
import argparse
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    if args.outdir.exists():
        parser.error(f"output directory already exists: {args.outdir}; choose a fresh namespace")
    repo = Path(__file__).resolve().parents[3]
    workflow = (repo / "workflows/rnadnavar.nf").read_text()
    bam_align = (repo / "subworkflows/local/bam_align/main.nf").read_text()
    sample_channel = (repo / "subworkflows/local/samplesheet_to_channel/main.nf").read_text()
    consensus = (repo / "subworkflows/local/vcf_consensus_workflow/main.nf").read_text()
    checks = [
        ("bam-align-mapping-gate", "if (params.step == 'mapping')" in bam_align),
        ("bam-input-nonmapping-route", "bam files" in sample_channel and "params.step != 'mapping'" in sample_channel),
        ("fastq-input-mapping-route", "fastq files" in sample_channel and "params.step == 'mapping'" in sample_channel),
        ("consensus-stage-route", "params.step in ['consensus', 'annotate', 'filtering', 'rna_filtering']" in consensus),
        ("rescue-process-route", "VCF_RESCUE_WORKFLOW" in consensus and "params.tools.split(',').contains('rescue')" in consensus),
        ("second-rescue-route", "SECOND_RESCUE_WORKFLOW" in workflow and "second_rescued_vcf" in workflow),
        ("realignment-route", "params.tools.split(',').contains('realignment')" in workflow),
    ]
    report = {"status": "pass" if all(ok for _, ok in checks) else "blocked",
              "static_only": True, "checks": [{"name": n, "status": "pass" if ok else "error"} for n, ok in checks]}
    args.outdir.mkdir(parents=True)
    (args.outdir / "process_reachability.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return 0 if report["status"] == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())
