#!/usr/bin/env python3
"""Audit the historical manual policy; never use this as a production caller.

The historical policy includes DeepSomatic-derived indels. This replay is for
comparison parity, not endorsement of that indel rule. Outputs are allele-only
benchmark VCFs, unsuitable for training labels. Original VCFs remain read-only.
WES uses the recovered winning queries. WGS here publishes only the control;
run replay_historical_candidate_scope.py afterward for the exact historical
candidate-scope comparison. Unrestricted WGS metrics are not published.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import gzip
import hashlib
import json
from pathlib import Path
import subprocess

from aggregate_benchmark import parse_metrics_json


def historical_snapshot(name, fallback):
    root = Path(__file__).resolve().parents[3] / 'examples/seqc2/verified/20260914/provenance_inputs'
    durable = root / name
    return durable if durable.exists() else Path(fallback)


def records(path):
    result = {}
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            key = (fields[0], int(fields[1]), fields[3], fields[4])
            if key in result:
                raise ValueError(f"Duplicate allele in {path}: {key}")
            result[key] = fields
    return result


def info(fields):
    return dict(value.split("=", 1) for value in fields[7].split(";") if "=" in value)


def historical_native(ds, m2):
    veto = {"contamination;germline;haplotype;panel_of_normals",
            "contamination;orientation;weak_evidence"}
    selected = set()
    for key, row in ds.items():
        other = m2.get(key)
        if other is not None and other[6] in veto:
            continue
        if row[6] == "PASS":
            selected.add(key)
        evidence = info(other) if other is not None else {}
        try:
            qualifies = (len(key[2]) == len(key[3]) == 1 and float(row[5]) > 0
                         and float(evidence.get("TLOD", 0)) >= 12
                         and float(evidence.get("GERMQ", 0)) >= 60)
        except ValueError:
            qualifies = False
        if qualifies:
            selected.add(key)
    return selected


def gated(native, rescue):
    selected = set(native)
    for key, row in rescue.items():
        if row[6] not in {"Somatic", "PASS"} or not len(key[2]) == len(key[3]) == 1:
            continue
        evidence = info(row)
        # Missing measurements are unavailable evidence, never a passing vote.
        try:
            qualifies = (int(evidence.get("N_DNA_CALLERS_SUPPORT", "0")) >= 1
                         and int(evidence.get("N_RNA_CALLERS_SOMATIC", "0")) >= 2)
        except ValueError:
            qualifies = False
        if qualifies:
            selected.add(key)
    return selected


def write_query(path, keys):
    with path.open("x") as handle:
        handle.write('##fileformat=VCFv4.2\n')
        for chrom in sorted({key[0] for key in keys}):
            handle.write(f'##contig=<ID={chrom}>\n')
        handle.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
        for chrom, pos, ref, alt in sorted(keys):
            handle.write(f'{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\n')
    subprocess.run(["bcftools", "view", "-Oz", "-o", str(path) + ".gz", str(path)], check=True)
    subprocess.run(["bcftools", "index", "-t", str(path) + ".gz"], check=True)
    return str(path) + ".gz"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--shared-repo", type=Path, default=Path(
        "/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/rnadnavar"))
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=False)
    repo = Path(__file__).resolve().parents[3]
    root = repo / "examples/seqc2"
    reference = Path("/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/bio_db/references/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta")
    hc = Path("/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/data/giab/data/seqc2/truth/High-Confidence_Regions_v1.2.bed")
    ukb = Path("/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/bio_db/intervals/ukb.pad50.broad.pad50.union.bed")
    targets = {"ukb": ukb, "medexome": root / "data/SeqCap_EZ_MedExome_hg38_empirical_targets.authoritative.bed"}
    truth = root / "comparison/common_policy_20260914/wgs_il/ukb/benchmark_truth.vcf.gz"
    historical = root / "hybrid/comparison/authoritative_medexome_bed_20260910"
    wgs = args.shared_repo / "examples/seqc2/hybrid/output/seqc2.wgs.il.hybrid"
    pair = "WGS_IL_T_1_vs_WGS_IL_N_1"
    provenance = {"purpose": "historical manual policy reproduction, not policy selection",
                  "sources": {}, "checks": {}, "metrics": {},
                  "domain_note": "Historical native inputs were preselected to HC and UKB. Replay applies the same HC/UKB preselection to every query; MedExome additionally uses -T MedExome.",
                  "indels": "Historical DeepSomatic PASS baseline subject to exact Mutect2 veto combinations; not ordinary consensus indels."}

    def track(path):
        path = Path(path).resolve()
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
        provenance["sources"][str(path)] = digest.hexdigest()
        return path

    def subset(source, dest):
        track(source)
        intermediate = dest.with_suffix(".subset.vcf.gz")
        subprocess.run(["bcftools", "view", "-R", str(hc), "-T", str(ukb),
                        "-Oz", "-o", str(intermediate), str(source)], check=True)
        subprocess.run(["bcftools", "norm", "-f", str(reference), "-m", "-any",
                        "-Oz", "-o", str(dest), str(intermediate)], check=True)
        return records(dest)

    datasets = {}
    wes_native = set(records(track(historical_snapshot('wes_native_original.vcf.gz', '/tmp/consensus_experiments/wesll/policy_relaxed/query.vcf.gz'))))
    wes_rescue = records(track(historical / "realignment_rescue.pass.vcf.gz"))
    wes_gate = gated(wes_native, wes_rescue)
    old_gate = set(records(track(historical / "native_gated_rescue.vcf.gz")))
    if wes_gate != old_gate:
        raise ValueError("Historical WES gate was not exactly reproduced")
    provenance["checks"]["wes_gate_exact_reproduction"] = True
    wes_ds_source = root / "hybrid/output/seqc2.wes.ll.hybrid.realign.policy-default/normalized/deepsomatic/WES_LL_T_1_vs_WES_LL_N_1/WES_LL_T_1_vs_WES_LL_N_1.deepsomatic.variants.dec.norm.vcf.gz"
    wes_ds = subset(wes_ds_source, args.outdir / "wes_ds.vcf.gz")
    datasets["wes_ll"] = {"deepsomatic": {k for k, v in wes_ds.items() if v[6] == "PASS"},
                          "historical_native": wes_native, "historical_native_gated": wes_gate}
    ds = subset(wgs / f"variant_calling/deepsomatic/{pair}/{pair}.deepsomatic.vcf.gz", args.outdir / "wgs_ds.vcf.gz")
    m2 = subset(wgs / f"variant_calling/mutect2/{pair}/{pair}.mutect2.filtered.vcf.gz", args.outdir / "wgs_m2.vcf.gz")
    wgs_native = historical_native(ds, m2)
    old_native = set(records(track(historical_snapshot('wgs_native_original.vcf.gz', '/tmp/consensus_experiments/wgsil/policy_relaxed/query.vcf.gz'))))
    provenance["checks"]["wgs_current_vs_historical_native"] = {
        "added": len(wgs_native - old_native), "removed": len(old_native - wgs_native)}
    rescue_pair = pair + "_rescued_WGS_IL_RT_1_realign_vs_WGS_IL_N_1"
    rescue = records(track(wgs / f"vcf_realignment/rescue/{rescue_pair}/{rescue_pair}.rescue.filtered.stripped.vep.vcf.gz"))
    # Only publish the WGS control here. The second script applies the exact
    # historical candidate universe before publishing native/gated metrics.
    datasets["wgs_il"] = {"deepsomatic": {k for k, v in ds.items() if v[6] == "PASS"}}
    jobs = []
    for dataset, queries in datasets.items():
        folder = args.outdir / dataset
        folder.mkdir()
        for name, keys in queries.items():
            query = write_query(folder / f"{name}.vcf", keys)
            for domain, target in targets.items():
                dest = folder / domain
                dest.mkdir(exist_ok=True)
                cmd = ["micromamba", "run", "-n", "happy", "som.py", str(truth), query,
                       "-R", str(hc), "-T", str(target), "-r", str(reference), "-N",
                       "-o", str(dest / name)]
                jobs.append((dataset, domain, name, dest, cmd))
    for path in (truth, hc, ukb, targets["medexome"]):
        track(path)
    provenance["reference"] = str(reference)
    provenance["commands"] = [job[-1] for job in jobs]
    (args.outdir / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")

    def benchmark(job):
        dataset, domain, name, dest, cmd = job
        with (dest / f"{name}.log").open("w") as log:
            subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=True)
        metrics = parse_metrics_json(dest / f"{name}.metrics.json")
        print(dataset, domain, name, metrics["records"], flush=True)
        return f"{dataset}/{domain}/{name}", metrics

    with ThreadPoolExecutor(max_workers=2) as pool:
        provenance["metrics"] = dict(pool.map(benchmark, jobs))
    (args.outdir / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")


if __name__ == "__main__":
    main()
