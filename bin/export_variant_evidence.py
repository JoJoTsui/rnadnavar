#!/usr/bin/env python3
"""Content-bound, read-only BAM evidence exporter for the approved seq2neo release.

See examples/seq2neo/docs/VARIANT_EVIDENCE.md for measurement semantics.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
import csv
import fcntl
import hashlib
import json
import os
from pathlib import Path
import resource
import subprocess
import sys
import time
from bisect import bisect_left

import polars as pl
import pyarrow as pa
import pyarrow.parquet as pq
import pysam

VERSION = "seq2neo-info-v1"
RELEASE_ID = "seq2neo_separated_three_class_v2_cohort63_20260921"
EXPECTED = {"Reference": 12791828, "Germline": 1522528, "Somatic": 28241}
RESERVED = {"PRJNA298310_3812", "PRJNA298376_4166", "PRJNA298376_4214",
            "PRJNA298376_4231", "PRJNA298376_4242"}
KEY = ["sample_id", "CHROM", "POS", "REF", "ALT"]
LABEL = KEY + ["FILTER"]
MODALITIES = {"DN": "normal_dna", "DT": "tumor_dna", "RT": "tumor_rna"}
COUNTS = ["depth", "ref_count", "alt_count", "other_count", "vaf_denominator",
          "ref_forward", "ref_reverse", "alt_forward", "alt_reverse",
          "ref_f1r2", "ref_f2r1", "alt_f1r2", "alt_f2r1",
          "ref_orientation_unknown", "alt_orientation_unknown", "mq_count"]
FLOATS = ["vaf", "mean_bq", "mean_mq"]
SETTINGS = {"min_bq": 20, "min_mq": 20, "max_depth": 10000,
            "window_bp": 100000, "exclude_flags": 0x4 | 0x100 | 0x200 | 0x400 | 0x800,
            "overlapping_mates": "count_both_reads", "baq": False,
            "star_mq255": "accept_only_NH1_effective_MQ60",
            "indels": "exact_CIGAR_with_right_flank_no_normalization",
            "split": "reserved_first_chr1_test_chr21_chr22_val_other_train"}


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def stat_identity(path):
    s = Path(path).stat()
    return {"realpath": str(Path(path).resolve()), "size": s.st_size,
            "mtime_ns": s.st_mtime_ns, "ctime_ns": s.st_ctime_ns,
            "device": s.st_dev, "inode": s.st_ino}


def identity(path):
    before = stat_identity(path)
    result = {"path": str(path), **before, "sha256": sha256(path)}
    if before != stat_identity(path):
        raise RuntimeError(f"Input changed while hashing: {path}")
    return result


def atomic_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + f".tmp.{os.getpid()}")
    with tmp.open("w") as f:
        json.dump(value, f, indent=2, sort_keys=True)
        f.write("\n")
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp, path)


def split(pool, chrom):
    if pool == "reserved":
        return "reserved"
    if pool != "train_pool":
        raise ValueError(f"Unknown pool: {pool}")
    return "test" if chrom == "chr1" else "val" if chrom in {"chr21", "chr22"} else "train"


def unique(df):
    if df.select(KEY).null_count().sum_horizontal().item():
        raise ValueError("Null allele key")
    if df.select(KEY).is_duplicated().any():
        raise ValueError("Duplicate full allele key")


def reconcile(expected, actual, pool):
    unique(expected)
    unique(actual)
    if not expected.select(LABEL).sort(KEY).equals(actual.select(LABEL).sort(KEY)):
        raise ValueError("Missing/extra allele keys or changed FILTER")
    if actual["pool"].to_list() != [pool] * len(actual):
        raise ValueError("Sample pool isolation failure")
    if actual["split"].to_list() != [split(pool, c) for c in actual["CHROM"]]:
        raise ValueError("Incorrect chromosome split")


def load_release(root):
    root = Path(root).resolve()
    ap = root / "cohort63_handoff_tools_20260921/RELEASE_APPROVAL.json"
    a = json.loads(ap.read_text())
    if a["release_id"] != RELEASE_ID or a["training_approved"] is not True:
        raise ValueError("This exporter requires the explicitly approved release")
    if a["class_counts"] != EXPECTED or a["total_records"] != sum(EXPECTED.values()):
        raise ValueError("Unexpected release totals")
    hashes = {"approval": identity(ap)}
    for field in ["variant_parquet", "manifest_tsv", "train_manifest", "reserved_manifest", "report"]:
        rel = a[field] if field != "report" else a["bridge_directory"] + "/report.json"
        hashes[field] = identity(root / rel)
        if hashes[field]["sha256"] != a[field + "_sha256"]:
            raise ValueError(f"Approval hash mismatch: {field}")
    pools = {}
    for pool, field in [("train_pool", "train_manifest"), ("reserved", "reserved_manifest")]:
        for sample in json.loads((root / a[field]).read_text())["samples"]:
            sid = sample["sample_id"]
            if sid in pools or sample.get("sample_pool") != pool:
                raise ValueError("Duplicate/conflicting pool membership")
            pools[sid] = {**sample, "pool": pool}
    if len(pools) != 63 or {s for s, v in pools.items() if v["pool"] == "reserved"} != RESERVED:
        raise ValueError("Incorrect approved sample pools")
    if set(a["excluded_samples"]) & pools.keys():
        raise ValueError("Excluded sample restored")
    with (root / a["manifest_tsv"]).open() as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    if len(rows) != 63 or {r["sample_id"] for r in rows} != pools.keys():
        raise ValueError("TSV / pool membership mismatch")
    for row in rows:
        for tag, field in MODALITIES.items():
            if row["bam_" + tag.lower()] != pools[row["sample_id"]][field]:
                raise ValueError("TSV / registered BAM mismatch")
    scan = pl.scan_parquet(root / a["variant_parquet"])
    counts = scan.group_by(["sample_id", "FILTER"]).len().collect()
    totals = counts.group_by("FILTER").agg(pl.col("len").sum())
    if dict(totals.iter_rows()) != EXPECTED or set(counts["sample_id"]) != pools.keys():
        raise ValueError("Parquet membership/class totals do not match approval")
    return a, hashes, pools, counts.to_dicts()


def allele_kind(ref, alt):
    if not ref or not alt or set(ref + alt) - set("ACGT") or ref == alt:
        return None
    if len(ref) == len(alt) == 1:
        return "snv"
    if len(ref) == 1 and len(alt) > 1 and alt.startswith(ref):
        return "insertion"
    if len(alt) == 1 and len(ref) > 1 and ref.startswith(alt):
        return "deletion"
    return None


def observation(read, pos0, ref, alt, settings, star=False):
    """Return (ref|alt|other, min event BQ, effective MQ, reverse, orientation).

    Read CIGAR projects the complete reference interval and indel right flank.
    No anchor-only indel inference and no N (splice) treated as D (deletion).
    """
    if read.flag & settings["exclude_flags"]:
        return None
    mq = read.mapping_quality
    if mq == 255:
        if not star or not read.has_tag("NH") or read.get_tag("NH") != 1:
            return None
        mq = 60
    if mq < settings["min_mq"] or read.query_sequence is None or read.query_qualities is None:
        return None
    kind = allele_kind(ref, alt)
    if kind is None:
        return None
    # Indels require a sequenced right flank, including REF observations.
    end = pos0 + len(ref) + (kind != "snv")
    r, q = read.reference_start, 0
    bases, qualities = [], []
    anchor = flank = False
    for op, n in read.cigartuples or []:
        if op in (0, 7, 8):
            lo, hi = max(pos0, r), min(end, r + n)
            if lo < hi:
                qs, qe = q + lo - r, q + hi - r
                bases.append(read.query_sequence[qs:qe])
                qualities.extend(read.query_qualities[qs:qe])
                anchor |= lo == pos0
                flank |= hi == end
            r += n
            q += n
        elif op == 1:
            if pos0 < r < end:
                bases.append(read.query_sequence[q:q+n])
                qualities.extend(read.query_qualities[q:q+n])
            q += n
        elif op in (2, 3):
            if op == 3 and r < end and r + n > pos0:
                return None
            r += n
        elif op == 4:
            q += n
        elif op not in (5, 6):
            return None
        if r >= end:
            break
    if not anchor or not flank or not qualities or min(qualities) < settings["min_bq"]:
        return None
    seq = "".join(bases).upper()
    if set(seq) - set("ACGT"):
        return None
    if kind != "snv":
        seq = seq[:-1]  # flank proves event is spanned; not part of allele
    state = "ref" if seq == ref else "alt" if seq == alt else "other"
    orientation = "unknown"
    if read.is_paired and read.is_read1 != read.is_read2:
        orientation = "f1r2" if read.is_read1 != read.is_reverse else "f2r1"
    return state, min(qualities), mq, read.is_reverse, orientation


def unavailable(reason, detail=None):
    return {**{k: None for k in COUNTS + FLOATS}, "status": reason,
            "missing_reason": detail or reason}


def finish(acc, settings):
    n = acc["depth"]
    if n > settings["max_depth"]:
        return unavailable("depth_cap_exceeded")
    out = {k: acc[k] for k in COUNTS}
    out.update(vaf_denominator=n, vaf=acc["alt_count"] / n if n else None,
               mean_bq=acc["bq_sum"] / n if n else None,
               mean_mq=acc["mq_sum"] / n if n else None,
               status="ok" if n else "zero_usable_depth",
               missing_reason=None if n else "no_read_spans_allele_after_filters")
    return out


def add(acc, obs, settings):
    if acc["depth"] > settings["max_depth"]:
        return
    state, bq, mq, reverse, orientation = obs
    acc["depth"] += 1
    acc[state + "_count"] += 1
    acc["bq_sum"] += bq
    acc["mq_sum"] += mq
    acc["mq_count"] += 1
    if state != "other":
        acc[state + ("_reverse" if reverse else "_forward")] += 1
        acc[state + ("_orientation_unknown" if orientation == "unknown" else "_" + orientation)] += 1


def measure(rows, bam, reference, settings, star=False, input_error=None):
    result = [None] * len(rows)
    windows = defaultdict(list)
    for i, row in enumerate(rows):
        chrom, pos, ref, alt = (row[k] for k in ["CHROM", "POS", "REF", "ALT"])
        if input_error:
            result[i] = unavailable(*input_error)
        elif allele_kind(ref, alt) is None:
            result[i] = unavailable("unsupported_allele")
        elif chrom not in reference.references:
            result[i] = unavailable("reference_contig_absent")
        elif pos < 1 or pos + len(ref) - 1 > reference.get_reference_length(chrom):
            result[i] = unavailable("invalid_coordinate")
        elif reference.fetch(chrom, pos - 1, pos - 1 + len(ref)).upper() != ref:
            result[i] = unavailable("reference_mismatch")
        elif chrom not in bam.references:
            result[i] = unavailable("bam_contig_absent")
        elif bam.get_reference_length(chrom) != reference.get_reference_length(chrom):
            result[i] = unavailable("reference_length_mismatch")
        else:
            windows[(chrom, (pos - 1) // settings["window_bp"])].append(i)
    for (chrom, _), indices in sorted(windows.items()):
        indices.sort(key=lambda i: rows[i]["POS"])
        positions = [rows[i]["POS"] - 1 for i in indices]
        acc = {i: Counter() for i in indices}
        try:
            for read in bam.fetch(chrom, positions[0], positions[-1] + 1):
                if read.flag & settings["exclude_flags"] or read.reference_end is None:
                    continue
                lo = bisect_left(positions, read.reference_start)
                hi = bisect_left(positions, read.reference_end)
                for j in range(lo, hi):
                    i = indices[j]
                    obs = observation(read, positions[j], rows[i]["REF"], rows[i]["ALT"], settings, star)
                    if obs is not None:
                        add(acc[i], obs, settings)
            for i in indices:
                result[i] = finish(acc[i], settings)
        except (OSError, ValueError) as e:
            for i in indices:
                result[i] = unavailable("bam_query_error", str(e))
    return result


def schema():
    fields = [pa.field(k, pa.int64() if k == "POS" else pa.string()) for k in LABEL]
    fields += [pa.field(k, pa.string()) for k in ["pool", "split", "expression_status"]]
    for tag in MODALITIES:
        fields += [pa.field(f"bam_{tag}_{k}", pa.int64()) for k in COUNTS]
        fields += [pa.field(f"bam_{tag}_{k}", pa.float64()) for k in FLOATS]
        fields += [pa.field(f"bam_{tag}_{k}", pa.string()) for k in ["status", "missing_reason"]]
    return pa.schema(fields, metadata={b"schema_version": VERSION.encode()})


def sample_labels(parquet, sid, pilot_rows=0):
    df = pl.scan_parquet(parquet).filter(pl.col("sample_id") == sid).select(LABEL).collect().sort(KEY)
    unique(df)
    if pilot_rows:
        # Deterministic allele-stratified pilot; no label-driven selection.
        indel = (pl.col("REF").str.len_chars() != 1) | (pl.col("ALT").str.len_chars() != 1)
        n = min(pilot_rows // 2, df.filter(indel).height)
        df = pl.concat([df.filter(indel).head(n), df.filter(~indel).head(pilot_rows - n)]).sort(KEY)
    return df


def open_bam(path, sid, tag):
    info = {"registered_path": path}
    if not Path(path).is_file():
        return None, info, ("bam_absent", path)
    print(f"{sid} {tag}: hashing BAM", flush=True)
    try:
        info["bam"] = identity(path)
        candidates = [Path(path + ".bai"), Path(path).with_suffix(".bai"), Path(path + ".csi")]
        index = next((p for p in candidates if p.is_file()), None)
        if index is None:
            return None, info, ("index_absent", path)
        info["index"] = identity(index)
        bam = pysam.AlignmentFile(path, "rb", index_filename=str(index), require_index=True)
        info["header"] = bam.header.to_dict()
        info["index_older_than_bam"] = index.stat().st_mtime_ns < Path(path).stat().st_mtime_ns
        samples = {r.get("SM") for r in info["header"].get("RG", [])}
        if samples != {sid + tag}:
            bam.close()
            return None, info, ("bam_sample_identity_mismatch", repr(sorted(str(s) for s in samples)))
        info["star"] = any(p.get("PN") == "STAR" for p in info["header"].get("PG", []))
        return bam, info, None
    except (OSError, ValueError) as e:
        return None, info, ("bam_unreadable", str(e))


def cache_valid(path, receipt, expected_identity):
    if not path.exists() or not receipt.exists():
        return False
    data = json.loads(receipt.read_text())
    if data["identity"] != expected_identity:
        raise ValueError(f"Existing output has a different input/extractor identity: {path}")
    if sha256(path) != data["sha256"]:
        raise ValueError(f"Output checksum mismatch: {path}")
    return True


def verify_measurements(df):
    """Check numerical invariants without using labels to judge evidence."""
    for tag in MODALITIES:
        prefix = f"bam_{tag}_"
        measured = df.filter(pl.col(prefix + "status").is_in(["ok", "zero_usable_depth"]))
        missing = df.filter(~pl.col(prefix + "status").is_in(["ok", "zero_usable_depth"]))
        if any(missing[prefix + k].null_count() != len(missing) for k in COUNTS + FLOATS):
            raise ValueError("Unavailable measurement contains fabricated values")
        if any(measured[prefix + k].null_count() for k in COUNTS):
            raise ValueError("Measured count is null")
        if len(measured):
            depth = measured[prefix + "depth"]
            if not (depth == measured[prefix + "ref_count"] + measured[prefix + "alt_count"] + measured[prefix + "other_count"]).all():
                raise ValueError("Depth does not equal REF+ALT+other")
            if not (depth == measured[prefix + "vaf_denominator"]).all():
                raise ValueError("VAF denominator mismatch")
            for state in ["ref", "alt"]:
                count = measured[prefix + state + "_count"]
                if not (count == measured[prefix + state + "_forward"] + measured[prefix + state + "_reverse"]).all():
                    raise ValueError("Strand count mismatch")
                if not (count == measured[prefix + state + "_f1r2"] + measured[prefix + state + "_f2r1"] + measured[prefix + state + "_orientation_unknown"]).all():
                    raise ValueError("Orientation count mismatch")
            nonzero = measured.filter(pl.col(prefix + "depth") > 0)
            if len(nonzero) and not ((nonzero[prefix + "vaf"] - nonzero[prefix + "alt_count"] / nonzero[prefix + "depth"]).abs() < 1e-12).all():
                raise ValueError("VAF mismatch")


def process_sample(task):
    started = time.monotonic()
    root, sample, context, args = task
    root = Path(root)
    sid, pool = sample["sample_id"], sample["pool"]
    if args["deadline"] and time.time() >= args["deadline"]:
        return {"sample_id": sid, "pool": pool, "rows": 0, "expected_rows": None,
                "classes": {}, "complete": False, "runtime_seconds": 0,
                "stop_reason": "configured_wall_time_budget_exhausted_before_sample"}
    outdir = root / ("reserved" if pool == "reserved" else "development") / sid
    outdir.mkdir(parents=True, exist_ok=True)
    labels = sample_labels(context["parquet"], sid, args["pilot_rows"])
    handles, bam_info, errors = {}, {}, {}
    for tag, field in MODALITIES.items():
        handles[tag], bam_info[tag], errors[tag] = open_bam(sample[field], sid, tag)
    base_id = digest({"context": context, "bam_info": bam_info, "errors": errors,
                      "sample_id": sid, "pool": pool, "pilot_rows": args["pilot_rows"]})
    inputs_path = outdir / "inputs.json"
    if inputs_path.exists():
        if not args["resume"] or json.loads(inputs_path.read_text())["identity"] != base_id:
            for bam in handles.values():
                if bam is not None:
                    bam.close()
            raise ValueError("Existing sample input identity differs, or --resume missing; provenance preserved")
    else:
        atomic_json(inputs_path, {"identity": base_id, "bams": bam_info,
                                 "input_errors": errors, "expected_rows": len(labels)})
    qc = {"sample_id": sid, "pool": pool, "rows": 0, "expected_rows": len(labels),
          "classes": Counter(), "splits": Counter(), "modalities": {t: Counter() for t in MODALITIES},
          "covered": Counter(), "alt_positive": Counter(), "parts": [], "complete": False}
    reference = pysam.FastaFile(context["reference"]["path"])
    try:
        for batch_no, batch in enumerate(labels.iter_slices(args["batch_rows"])):
            if args["deadline"] and time.time() >= args["deadline"]:
                qc["stop_reason"] = "configured_wall_time_budget_exhausted"
                break
            path = outdir / f"part-{batch_no:05d}.parquet"
            receipt = path.with_suffix(".json")
            rows = batch.to_dicts()
            part_id = digest({"sample": base_id, "keys_labels": rows})
            if args["resume"] and cache_valid(path, receipt, part_id):
                actual = pl.read_parquet(path)
            else:
                if receipt.exists() or (path.exists() and not args["resume"]):
                    raise ValueError(f"Refusing to overwrite existing output; use --resume: {path}")
                for row in rows:
                    row.update(pool=pool, split=split(pool, row["CHROM"]),
                               expression_status="unavailable_no_versioned_quantification")
                for tag in MODALITIES:
                    values = measure(rows, handles[tag], reference, SETTINGS,
                                     bam_info[tag].get("star", False), errors[tag])
                    for row, measurements in zip(rows, values):
                        row.update({f"bam_{tag}_{k}": v for k, v in measurements.items()})
                table = pa.Table.from_pylist(rows, schema=schema())
                tmp = path.with_name(path.name + f".tmp.{os.getpid()}")
                pq.write_table(table, tmp, compression="zstd")
                actual = pl.read_parquet(tmp)
                reconcile(batch, actual, pool)
                verify_measurements(actual)
                with tmp.open("rb") as f:
                    os.fsync(f.fileno())
                os.replace(tmp, path)
                atomic_json(receipt, {"identity": part_id, "sha256": sha256(path), "rows": len(batch)})
            reconcile(batch, actual, pool)
            verify_measurements(actual)
            qc["rows"] += len(actual)
            qc["classes"].update(actual["FILTER"].to_list())
            qc["splits"].update(actual["split"].to_list())
            for tag in MODALITIES:
                qc["modalities"][tag].update(actual[f"bam_{tag}_status"].to_list())
                qc["covered"][tag] += int((actual[f"bam_{tag}_depth"] > 0).sum())
                qc["alt_positive"][tag] += int((actual[f"bam_{tag}_alt_count"] > 0).sum())
            qc["parts"].append(str(path.relative_to(root)))
            qc["runtime_seconds"] = time.monotonic() - started
            atomic_json(outdir / "qc.json", qc)
            print(f"{sid}: {qc['rows']}/{len(labels)} rows; {qc['runtime_seconds']:.1f}s", flush=True)
        # Detect source modification during extraction (hash was taken before).
        for info in bam_info.values():
            for field in ["bam", "index"]:
                if field in info:
                    recorded = info[field]
                    if stat_identity(recorded["path"]) != {k: recorded[k] for k in stat_identity(recorded["path"])}:
                        raise RuntimeError("BAM/index changed during extraction")
        qc["complete"] = qc["rows"] == len(labels)
        expected_parts = {f"part-{i:05d}.parquet" for i in range((len(labels) + args["batch_rows"] - 1) // args["batch_rows"])}
        if {p.name for p in outdir.glob("part-*.parquet")} - expected_parts:
            raise ValueError("Unexpected extra output parts")
        qc["runtime_seconds"] = time.monotonic() - started
        usage = resource.getrusage(resource.RUSAGE_SELF)
        qc["resource_use"] = {"max_rss_kib": usage.ru_maxrss, "user_cpu_seconds": usage.ru_utime,
                               "system_cpu_seconds": usage.ru_stime, "block_inputs": usage.ru_inblock}
        atomic_json(outdir / "qc.json", qc)
        return qc
    finally:
        reference.close()
        for bam in handles.values():
            if bam is not None:
                bam.close()


def code_identity():
    repo = Path(__file__).resolve().parents[1]
    def git(*args):
        return subprocess.check_output(["git", "-C", str(repo), *args], text=True).strip()
    paths = [Path(__file__).resolve(), repo / "tests/test_export_variant_evidence.py",
             repo / "examples/seq2neo/docs/VARIANT_EVIDENCE.md",
             repo / "examples/seq2neo/docs/TRAINING_DATA_GUIDE.md", repo / "uv.lock"]
    return {"revision": git("rev-parse", "HEAD"), "git_status": git("status", "--porcelain"),
            "relevant_file_hashes": {str(p.relative_to(repo)): sha256(p) for p in paths if p.exists()},
            "python": sys.version, "pysam": pysam.__version__, "htslib": pysam.__samtools_version__,
            "polars": pl.__version__, "pyarrow": pa.__version__}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--release-root", type=Path, required=True)
    ap.add_argument("--reference", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--sample", action="append")
    ap.add_argument("--pilot-rows", type=int, default=0)
    ap.add_argument("--batch-rows", type=int, default=10000)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--max-runtime-seconds", type=int, default=0,
                    help="Optional resumable scheduling budget, checked between batches")
    ap.add_argument("--resume", action="store_true")
    args = ap.parse_args()
    if args.batch_rows < 1 or args.workers < 1 or args.pilot_rows < 0:
        ap.error("Invalid batch size, concurrency or pilot size")
    args.output = args.output.resolve()
    release_root = args.release_root.resolve()
    if args.output == release_root or release_root in args.output.parents:
        ap.error("Output must be separate from the immutable approved release")
    historical = Path(__file__).resolve().parents[1] / "examples/seq2neo/stats"
    if args.output == historical or historical in args.output.parents:
        ap.error("Output must be separate from historical statistics")
    args.output.mkdir(parents=True, exist_ok=True)
    lock = (args.output / ".lock").open("w")
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    start = time.time()
    approval, hashes, pools, counts = load_release(release_root)
    selected = sorted(args.sample or pools)
    if len(selected) != len(set(selected)) or set(selected) - pools.keys():
        ap.error("Samples must be unique full approved sample IDs")
    if args.pilot_rows and (len(selected) != 1 or pools[selected[0]]["pool"] != "train_pool"):
        ap.error("A bounded pilot must select one training-pool sample")
    print("Release hashes, 63 sample memberships, and class totals verified", flush=True)
    context = {"schema_version": VERSION, "release_id": RELEASE_ID, "release_hashes": hashes,
               "reference": identity(args.reference.resolve()),
               "reference_index": identity(Path(str(args.reference.resolve()) + ".fai")),
               "parquet": str(release_root / approval["variant_parquet"]),
               "settings": SETTINGS, "code": code_identity(), "batch_rows": args.batch_rows,
               "pilot_rows": args.pilot_rows}
    # Git status is recorded but unrelated dirty files must not invalidate evidence.
    cache_context = json.loads(json.dumps(context))
    del cache_context["code"]["git_status"]
    manifest = args.output / "manifest.json"
    run_id = digest(cache_context)
    if manifest.exists():
        if not args.resume or json.loads(manifest.read_text())["identity"] != run_id:
            raise ValueError("Output identity differs, or --resume missing; choose a new version directory")
    else:
        atomic_json(manifest, {"identity": run_id, **context, "approval_policy":
                    "Hash-bound sidecar acknowledged; immutable eligibility flags intentionally not used",
                    "expression": "unavailable: no versioned quantification source registered",
                    "historical_measurements_reused": False, "approved_counts": counts})
    atomic_json(args.output / "schema.json", {"version": VERSION,
                "fields": [{"name": f.name, "type": str(f.type), "nullable": f.nullable} for f in schema()],
                "model_feature_allowlist": [f"bam_{t}_{k}" for t in MODALITIES for k in COUNTS + FLOATS],
                "provenance_columns": LABEL + ["pool", "split", "expression_status"],
                "missingness_columns": [f"bam_{t}_{k}" for t in MODALITIES for k in ["status", "missing_reason"]]})
    run_args = {"pilot_rows": args.pilot_rows, "batch_rows": args.batch_rows, "resume": args.resume,
                "deadline": start + args.max_runtime_seconds if args.max_runtime_seconds else 0}
    tasks = [(str(args.output), pools[s], cache_context, run_args) for s in selected]
    if args.workers == 1:
        results = [process_sample(t) for t in tasks]
    else:
        # spawn avoids forking a live Polars/Rayon thread pool.
        import multiprocessing
        with ProcessPoolExecutor(max_workers=args.workers, mp_context=multiprocessing.get_context("spawn")) as ex:
            results = list(ex.map(process_sample, tasks))
    total_classes = Counter()
    for q in results:
        total_classes.update(q["classes"])
    full = set(selected) == pools.keys() and not args.pilot_rows and all(q["complete"] for q in results)
    if full and dict(total_classes) != EXPECTED:
        raise ValueError("Final class reconciliation failed")
    # Recheck all approval-bound inputs after generation.
    load_release(release_root)
    if stat_identity(args.reference) != {k: context["reference"][k] for k in stat_identity(args.reference)}:
        raise RuntimeError("Reference changed during extraction")
    report = {"full_release_complete": full, "selected_samples_complete": all(q["complete"] for q in results),
              "pilot_rows": args.pilot_rows, "rows": sum(q["rows"] for q in results),
              "classes": total_classes, "samples": results, "runtime_seconds": time.time() - start,
              "command": [sys.executable, *sys.argv], "workers": args.workers,
              "reconciliation": "exact full-key/FILTER equality checked for every written or resumed part",
              "remaining_release_rows": sum(EXPECTED.values()) - sum(q["rows"] for q in results)}
    atomic_json(args.output / "qc.json", report)
    atomic_json(args.output / f"run-{time.time_ns()}.json", report)
    print(json.dumps({k: v for k, v in report.items() if k != "samples"}, indent=2), flush=True)


if __name__ == "__main__":
    main()
