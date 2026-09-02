#!/usr/bin/env python3
"""Aggregate som.py metrics.json outputs into one benchmark comparison table.

For each query name Q, reads <metrics-dir>/Q.metrics.json (som.py native table
format: each metric has values = [SNV, INDEL, records]) and emits one row per
(query, variant_type) with TP/FP/FN/Precision/Recall/F1.

Stdlib only. Example:
    python3 aggregate_benchmark.py --metrics-dir comparison/WES_LL_T_1_vs_WES_LL_N_1 \
        --queries consensus mutect2 deepsomatic strelka \
        --output comparison/WES_LL_T_1_vs_WES_LL_N_1/WES_LL_T_1_vs_WES_LL_N_1.benchmark_comparison.csv
"""

import argparse
import csv
import json
import sys
from pathlib import Path

# som.py metrics.json rows are ordered [SNV, INDEL, records(=all)]
TYPE_LABELS = ["snv", "indel", "all"]
METRIC_IDS = ("tp", "fp", "fn", "precision", "recall")


def prf(tp, fp, fn):
    precision = tp / (tp + fp) if (tp + fp) else 0.0
    recall = tp / (tp + fn) if (tp + fn) else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) else 0.0
    return precision, recall, f1


def parse_metrics_json(path):
    """Return {variant_type: {tp, fp, fn, precision, recall, f1}}."""
    with open(path) as fh:
        doc = json.load(fh)

    if "metrics" not in doc or not isinstance(doc["metrics"], list):
        raise ValueError(f"{path}: not som.py table-format metrics.json")

    per_metric = {}
    for item in doc["metrics"][0].get("data", []):
        mid, values = item.get("id"), item.get("values", [])
        if mid in METRIC_IDS and values:
            per_metric[mid] = values

    missing = [m for m in ("tp", "fp", "fn") if m not in per_metric]
    if missing:
        raise ValueError(f"{path}: missing metric ids: {missing}")

    n_rows = len(per_metric["tp"])
    # Prefer explicit column labels when som.py provides them; else positional.
    columns = doc["metrics"][0].get("columns")
    if columns and len(columns) == n_rows:
        def _norm(c):
            c = str(c).lower()
            if c in ("records", "total", "all"):
                return "all"
            return "snv" if c in ("snp", "snv") else c
        labels = [_norm(c) for c in columns]
    else:
        labels = TYPE_LABELS if n_rows == 3 else (
            TYPE_LABELS[: n_rows - 1] + ["all"] if n_rows > 1 else ["all"])

    rows = {}
    for i, label in enumerate(labels):
        tp = int(per_metric["tp"][i])
        fp = int(per_metric["fp"][i])
        fn = int(per_metric["fn"][i])
        # Prefer som.py's own precision/recall when present; recompute otherwise.
        precision = float(per_metric["precision"][i]) if "precision" in per_metric else None
        recall = float(per_metric["recall"][i]) if "recall" in per_metric else None
        if precision is None or recall is None:
            precision, recall, _ = prf(tp, fp, fn)
        _, _, f1 = prf(tp, fp, fn)
        rows[label] = dict(tp=tp, fp=fp, fn=fn,
                           precision=precision, recall=recall, f1=f1)
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--metrics-dir", required=True, type=Path,
                    help="Directory containing <query>.metrics.json files")
    ap.add_argument("--queries", required=True, nargs="+",
                    help="Query names, e.g. consensus mutect2 deepsomatic strelka")
    ap.add_argument("--output", required=True, type=Path, help="Output CSV path")
    args = ap.parse_args(argv)

    fieldnames = ["query", "variant_type", "tp", "fp", "fn",
                  "precision", "recall", "f1"]
    table = []
    for query in args.queries:
        path = args.metrics_dir / f"{query}.metrics.json"
        if not path.exists():
            print(f"WARNING: skipping {query}: {path} not found", file=sys.stderr)
            continue
        for vtype, m in parse_metrics_json(path).items():
            table.append({"query": query, "variant_type": vtype, **m})

    if not table:
        sys.exit("ERROR: no metrics parsed, nothing to write")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for row in table:
            writer.writerow({k: (f"{v:.4f}" if isinstance(v, float) else v)
                             for k, v in row.items()})

    # Human-readable echo
    header = ("query", "type", "tp", "fp", "fn", "prec", "rec", "f1")
    print(("{:<12} {:<6} {:>7} {:>6} {:>6} {:>8} {:>8} {:>8}").format(*header))
    for r in table:
        print("{:<12} {:<6} {:>7} {:>6} {:>6} {:>8.4f} {:>8.4f} {:>8.4f}".format(
            r["query"], r["variant_type"], r["tp"], r["fp"], r["fn"],
            r["precision"], r["recall"], r["f1"]))
    print(f"\nWrote {args.output}")


if __name__ == "__main__":
    main()
