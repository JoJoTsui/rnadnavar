#!/usr/bin/env python3
"""Complete WGS parity replay using the recovered historical candidate universe.

Run after replay_historical_native_gate.py. The older WGS evidence lookup used
positions in DeepSomatic PASS union original consensus, not all raw VCF sites.
The superseded all-site benchmark is retired; this script needs no such output.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path
import subprocess

from replay_historical_native_gate import records, gated, write_query, historical_native, historical_snapshot
from aggregate_benchmark import parse_metrics_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replay-dir", type=Path, required=True)
    args = parser.parse_args()
    root = args.replay_dir
    dest = root / "wgs_historical_scope"
    dest.mkdir(exist_ok=False)
    provenance = json.loads((root / "provenance.json").read_text())
    ds_path = historical_snapshot('wgs_deepsomatic_pass.vcf.gz', '/tmp/consensus_experiments/wgsil/audit/ds.vcf.gz')
    consensus_path = historical_snapshot('wgs_original_consensus.vcf.gz', '/tmp/consensus_experiments/wgsil/audit/consensus.vcf.gz')
    candidate_bed = historical_snapshot('wgs_candidate_positions.bed', '/tmp/wgsil_sites.bed')
    old_query = historical_snapshot('wgs_native_original.vcf.gz', '/tmp/consensus_experiments/wgsil/policy_relaxed/query.vcf.gz')
    positions = {key[:2] for path in (ds_path, consensus_path) for key in records(path)}
    with candidate_bed.open() as handle:
        recovered = {(parts[0], int(parts[1]) + 1) for line in handle if len(parts := line.split()) >= 3}
    if positions != recovered:
        raise ValueError("Historical candidate BED differs from DS/consensus union")
    all_site_native = historical_native(records(root / "wgs_ds.vcf.gz"),
                                        records(root / "wgs_m2.vcf.gz"))
    native = {key for key in all_site_native if key[:2] in positions}
    if native != set(records(old_query)):
        raise ValueError("Current WGS evidence with historical scope does not reproduce the old query")
    rescue_path = next(Path(path) for path in provenance["sources"] if path.endswith(".rescue.filtered.stripped.vep.vcf.gz"))
    queries = {"historical_native": native,
               "historical_native_gated": gated(native, records(rescue_path))}
    sources = {str(path): hashlib.sha256(path.read_bytes()).hexdigest()
               for path in (ds_path, consensus_path, candidate_bed, old_query, rescue_path)}
    jobs = []
    for name, keys in queries.items():
        query = write_query(dest / f"{name}.vcf", keys)
        for domain in ("ukb", "medexome"):
            (dest / domain).mkdir(exist_ok=True)
            template = next(cmd for cmd in provenance["commands"]
                            if cmd[-1].endswith(f"wgs_il/{domain}/deepsomatic"))
            cmd = list(template)
            cmd[6] = query
            cmd[-1] = str(dest / domain / name)
            jobs.append((domain, name, cmd))
    result = {"historical_candidate_positions": len(positions),
              "candidate_bed_exact_union_match": True, "native_exact_allele_match": True,
              "source_sha256": sources, "commands": [job[-1] for job in jobs]}
    (dest / "provenance.json").write_text(json.dumps(result, indent=2) + "\n")

    def run(job):
        domain, name, cmd = job
        with (dest / domain / f"{name}.log").open("w") as log:
            subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=True)
        metrics = parse_metrics_json(dest / domain / f"{name}.metrics.json")
        print(domain, name, metrics["records"], flush=True)
        return f"{domain}/{name}", metrics

    with ThreadPoolExecutor(max_workers=2) as pool:
        result["metrics"] = dict(pool.map(run, jobs))
    (dest / "provenance.json").write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
