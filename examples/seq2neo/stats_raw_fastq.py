#!/usr/bin/env python3
"""
Raw FASTQ statistics for all running seq2neo samples using seqkit.

Generates two output files:
  1. FASTQ-level TSV (--fastq-out) — per-file stats with set/sample/modality context
  2. Sample-level TSV (--sample-out) — per sample with standalone r1/r2 columns
     for primary pairs, extra pairs aggregated together.

Progress is printed to stderr.
"""

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
MERGED_JSON = HERE / "data" / "processed" / "merged.json"
REQUIRED_MODALITIES = ["DN", "DT", "RT"]

FASTQ_COLS = [
    "set", "patient_id", "disease", "status",
    "modality", "pair_type", "pair_idx", "r_label",
    "file",
    "format", "type", "num_seqs", "sum_len",
    "min_len", "avg_len", "max_len", "Q1", "Q2", "Q3",
    "sum_gap", "N50", "N50_num", "Q20(%)", "Q30(%)",
    "AvgQual", "GC(%)", "sum_n",
]
# seqkit output columns (in order): file, format, type, num_seqs, sum_len,
#   min_len, avg_len, max_len, Q1, Q2, Q3, sum_gap, N50, N50_num,
#   Q20(%), Q30(%), AvgQual, GC(%), sum_n
# FASTQ_COLS[9:] aligns with seqkit columns starting from "format"
# (skipping seqkit's "file" since we write our own via row.append)

SAMPLE_STAT_KEYS = ["reads", "bases", "GC", "Q20", "Q30", "AvgQual"]


def build_sample_header():
    """TSV header for sample-level output."""
    h = ["set", "patient_id", "disease", "status", "n_extra_pairs"]
    for mod in REQUIRED_MODALITIES:
        for r in ("r1", "r2"):
            for k in SAMPLE_STAT_KEYS:
                h.append(f"{mod}_{r}_{k}")
        for k in SAMPLE_STAT_KEYS:
            h.append(f"{mod}_extra_{k}")
    h.extend(["total_reads", "total_bases"])
    return h


def build_sample_rows(entries, stats_by_file):
    """
    Aggregate FASTQ-level stats to sample level.

    Returns a list of dicts, one per sample, sorted by set then patient_id.
    Primary pairs: r1 and r2 are standalone columns.
    Extra pairs: all extra pairs in a modality are aggregated together.
    """
    # Intermediate accumulator: sample_key -> mod -> (r1 stats, r2 stats, extra stats list)
    samples = {}

    for e in entries:
        sample_key = (e["set"], e["patient"], e["disease"], e["status"])
        mod = e["modality"]
        pt = e["pair_type"]
        rl = e["r_label"]

        stats = stats_by_file.get(e["file"])
        if stats is None:
            continue

        num_seqs = int(stats.get("num_seqs", 0))
        sum_len = int(stats.get("sum_len", 0))
        gc = float(stats.get("GC(%)", 0))
        q20 = float(stats.get("Q20(%)", 0))
        q30 = float(stats.get("Q30(%)", 0))
        avgqual = float(stats.get("AvgQual", 0))

        if sample_key not in samples:
            samples[sample_key] = {"_n_extra": 0}

        smp = samples[sample_key]

        if pt == "extra":
            smp["_n_extra"] = max(smp["_n_extra"], e["pair_idx"] - 1)
            prefix = f"{mod}_extra"
        else:
            prefix = f"{mod}_{rl}"

        # Accumulate reads/bases (sum) and quality metrics (weighted by num_seqs)
        for key, val in [("reads", num_seqs), ("bases", sum_len),
                          ("w", num_seqs)]:
            smp.setdefault(f"{prefix}_{key}", 0)
            smp[f"{prefix}_{key}"] += val

        for key, val in [("gc_w", gc * num_seqs), ("q20_w", q20 * num_seqs),
                          ("q30_w", q30 * num_seqs), ("avgqual_w", avgqual * num_seqs)]:
            smp.setdefault(f"{prefix}_{key}", 0.0)
            smp[f"{prefix}_{key}"] += val

    # Build output rows
    out = []
    for (pset, patient, disease, status), smp in sorted(
        samples.items(), key=lambda x: ((x[0][0] or 99), x[0][1])
    ):
        row = {
            "set": pset,
            "patient_id": patient,
            "disease": disease,
            "status": status,
            "n_extra_pairs": smp["_n_extra"],
        }
        total_reads = 0
        total_bases = 0

        for mod in REQUIRED_MODALITIES:
            for suffix in ("r1", "r2", "extra"):
                prefix = f"{mod}_{suffix}"
                w = smp.get(f"{prefix}_w", 0)
                reads = smp.get(f"{prefix}_reads", 0)
                bases = smp.get(f"{prefix}_bases", 0)

                if w > 0:
                    row[f"{prefix}_reads"] = reads
                    row[f"{prefix}_bases"] = bases
                    row[f"{prefix}_GC"] = round(smp[f"{prefix}_gc_w"] / w, 2)
                    row[f"{prefix}_Q20"] = round(smp[f"{prefix}_q20_w"] / w, 2)
                    row[f"{prefix}_Q30"] = round(smp[f"{prefix}_q30_w"] / w, 2)
                    row[f"{prefix}_AvgQual"] = round(smp[f"{prefix}_avgqual_w"] / w, 2)

                total_reads += reads
                total_bases += bases

        row["total_reads"] = total_reads
        row["total_bases"] = total_bases
        out.append(row)

    return out


def collect_entries(data):
    """Parse merged.json into entries list and set of file paths."""
    entries = []
    all_paths = set()

    for proj in data["projects"]:
        project_id = proj["project_id"]
        for s in proj["samples"]:
            if s["status"] not in ("standard", "extra"):
                continue

            pset = s.get("partition_set")
            patient = f"{project_id}_{s['patient_id']}"
            disease = s["disease_normalized"]
            status = s["status"]

            for mod in REQUIRED_MODALITIES:
                mod_data = s["modalities"].get(mod)
                if not mod_data or not mod_data["pairs"]:
                    continue

                for i, pair in enumerate(mod_data["pairs"]):
                    pair_type = "primary" if i == 0 else "extra"
                    pair_idx = i + 1
                    for label in ("r1", "r2"):
                        fpath = pair[label]
                        all_paths.add(fpath)
                        entries.append({
                            "set": pset,
                            "patient": patient,
                            "disease": disease,
                            "status": status,
                            "modality": mod,
                            "pair_type": pair_type,
                            "pair_idx": pair_idx,
                            "r_label": label,
                            "file": fpath,
                        })
    return entries, all_paths


def main():
    parser = argparse.ArgumentParser(
        description="Raw FASTQ statistics for seq2neo samples using seqkit"
    )
    parser.add_argument("--fastq-out", type=Path, required=True,
                        help="Write FASTQ-level per-file stats to this file")
    parser.add_argument("--sample-out", type=Path, default=None,
                        help="Write sample-level aggregated stats to this file")
    parser.add_argument("--max-samples", type=int, default=None,
                        help="Limit to first N eligible samples (for testing)")
    args = parser.parse_args()

    with open(MERGED_JSON) as fh:
        data = json.load(fh)

    entries, all_paths = collect_entries(data)

    if args.max_samples is not None:
        # Keep only the first N unique patients (preserves 1 std + 1 extra ordering)
        seen = set()
        keep_patients = set()
        for e in entries:
            if e["patient"] not in seen:
                seen.add(e["patient"])
                keep_patients.add(e["patient"])
                if len(keep_patients) >= args.max_samples:
                    break
        entries = [e for e in entries if e["patient"] in keep_patients]
        all_paths = {e["file"] for e in entries}
        print(f"Testing with {len(keep_patients)} sample(s): "
              f"{len(entries)} FASTQ files", file=sys.stderr)

    # Pre-check file existence
    print(f"Checking {len(all_paths)} FASTQ files...", file=sys.stderr)
    existing = []
    missing = []
    for p in sorted(all_paths):
        if os.path.exists(p):
            existing.append(p)
        else:
            missing.append(p)
    if missing:
        for p in missing:
            print(f"WARNING: file not found, skipping: {p}", file=sys.stderr)
        print(f"WARNING: {len(missing)} of {len(all_paths)} files missing",
              file=sys.stderr)

    if not existing:
        print("ERROR: no FASTQ files found", file=sys.stderr)
        sys.exit(1)

    # Run seqkit on existing files
    print(f"Running seqkit stats on {len(existing)} files "
          f"({len(entries)} entries, {len({e['patient'] for e in entries})} samples)...",
          file=sys.stderr)
    path_list = "\n".join(existing)
    result = subprocess.run(
        ["seqkit", "stats", "--all", "--tabular", "--skip-err",
         "--skip-file-check", "--threads", "32", "--infile-list", "-"],
        input=path_list, capture_output=True, text=True,
    )

    if result.returncode != 0:
        print(result.stderr, file=sys.stderr)
        sys.exit(result.returncode)

    lines = result.stdout.strip().split("\n")
    if not lines:
        print("No output from seqkit", file=sys.stderr)
        sys.exit(1)

    seqkit_header = lines[0].split("\t")
    stats_by_file = {}
    for line in lines[1:]:
        if not line.strip():
            continue
        cols = line.split("\t")
        stats_by_file[cols[0]] = dict(zip(seqkit_header, cols))

    # --- FASTQ-level output ---
    print("Writing FASTQ-level stats...", file=sys.stderr)
    mod_order = {"DN": 0, "DT": 1, "RT": 2}
    type_order = {"primary": 0, "extra": 1}
    entries.sort(key=lambda e: (
        e["set"] or 99,
        e["patient"],
        mod_order.get(e["modality"], 9),
        type_order.get(e["pair_type"], 9),
        e["pair_idx"],
        e["r_label"],
    ))

    fastq_missing = 0
    with open(args.fastq_out, "w") as fh:
        fh.write("\t".join(FASTQ_COLS) + "\n")
        for e in entries:
            stats = stats_by_file.get(e["file"])
            if stats is None:
                fastq_missing += 1
                continue
            row = [str(e[k]) for k in ("set", "patient", "disease", "status",
                                         "modality", "pair_type", "pair_idx", "r_label")]
            row.append(e["file"])
            row.extend(stats.get(c, "") for c in FASTQ_COLS[9:])
            fh.write("\t".join(row) + "\n")

    print(f"FASTQ-level stats written to {args.fastq_out} "
          f"({len(entries) - fastq_missing} rows)", file=sys.stderr)
    if fastq_missing:
        print(f"Warning: {fastq_missing} files missing from seqkit output",
              file=sys.stderr)

    # --- Sample-level output ---
    if args.sample_out:
        print("Writing sample-level stats...", file=sys.stderr)
        sample_rows = build_sample_rows(entries, stats_by_file)
        header = build_sample_header()
        with open(args.sample_out, "w") as fh:
            fh.write("\t".join(header) + "\n")
            for row in sample_rows:
                fh.write("\t".join(str(row.get(c, "")) for c in header) + "\n")
        print(f"Sample-level stats written to {args.sample_out} "
              f"({len(sample_rows)} rows)", file=sys.stderr)


if __name__ == "__main__":
    main()
