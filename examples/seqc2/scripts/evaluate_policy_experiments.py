#!/usr/bin/env python3
"""Evaluate ordinary/native/strict policy metrics against the frozen matrix."""
import argparse
import csv
import json
from pathlib import Path


TYPES = ("snp", "indel", "records")


def metrics(path):
    doc = json.loads(Path(path).read_text())
    data = doc["metrics"][0]["data"]
    by_id = {row["id"]: row["values"] for row in data}
    labels = [str(x).lower() for x in by_id["type"]]
    labels = [{"snvs": "snp", "indels": "indel"}.get(x, x) for x in labels]
    result = {}
    for index, label in enumerate(labels):
        tp, fp, fn = (int(by_id[key][index]) for key in ("tp", "fp", "fn"))
        precision = tp / (tp + fp) if tp + fp else 0.0
        recall = tp / (tp + fn) if tp + fn else 0.0
        f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
        result[label] = {"tp": tp, "fp": fp, "fn": fn, "precision": precision, "recall": recall, "f1": f1}
    return result


def load_csv(path):
    return {(row["query"], row["variant_type"]): row for row in csv.DictReader(path.open())}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--experiments", required=True, type=Path)
    ap.add_argument("--comparisons", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()
    failures = []
    rows = []
    for cell in sorted(args.experiments.glob("*/*")):
        if not cell.is_dir() or cell.name == ".caller_inputs":
            continue
        relative = cell.relative_to(args.experiments)
        comparison = args.comparisons / relative / f"{('WES_LL_T_1_vs_WES_LL_N_1' if relative.parts[0] == 'wes_ll' else 'WGS_IL_T_1_vs_WGS_IL_N_1')}.benchmark_comparison.csv"
        comparator = load_csv(comparison)
        for policy in ("ordinary", "native", "strict"):
            candidate = metrics(cell / policy / "benchmark.metrics.json")
            for variant_type in TYPES:
                c = candidate[variant_type]
                rows.append({"cell": str(relative), "policy": policy, "variant_type": variant_type, **c})
                for caller in ("mutect2", "strelka", "deepsomatic"):
                    base = comparator[(caller, variant_type)]
                    if c["f1"] <= float(base["f1"]):
                        failures.append({"cell": str(relative), "policy": policy, "variant_type": variant_type, "reason": f"f1_not_above_{caller}"})
                ds = comparator[("deepsomatic", variant_type)]
                if c["precision"] < float(ds["precision"]):
                    failures.append({"cell": str(relative), "policy": policy, "variant_type": variant_type, "reason": "precision_below_deepsomatic"})
    status = "qualified" if not failures else "no_qualifying_policy"
    result = {"status": status, "policies": ["ordinary", "native", "strict"], "rows": rows, "failures": failures, "criteria": "F1 strictly above Mutect2/Strelka2/DeepSomatic and precision >= DeepSomatic in every cell and variant class"}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": status, "rows": len(rows), "failures": len(failures)}, indent=2))


if __name__ == "__main__":
    main()
