#!/usr/bin/env python3
"""Audit an evidence dataset against the approved source, with optional BAM SNV checks."""
import argparse
from collections import Counter
import json
from pathlib import Path
import time

import polars as pl
import pysam

import export_variant_evidence as e


def independent_indel_counts(row, bam, star):
    """Use htslib-expanded aligned pairs instead of the exporter's CIGAR walk."""
    start, end = row["POS"] - 1, row["POS"] + len(row["REF"])
    counts = Counter(depth=0, ref_count=0, alt_count=0)
    for read in bam.fetch(row["CHROM"], start, start + 1):
        if read.flag & e.SETTINGS["exclude_flags"] or read.query_sequence is None or read.query_qualities is None:
            continue
        mq = read.mapping_quality
        if mq == 255:
            if not star or not read.has_tag("NH") or read.get_tag("NH") != 1:
                continue
            mq = 60
        if mq < e.SETTINGS["min_mq"]:
            continue
        pairs = read.get_aligned_pairs(with_cigar=True)
        anchor = flank = False
        spliced = False
        query_positions = []
        previous_reference = None
        for q, r, op in pairs:
            if r is not None:
                previous_reference = r
                if start <= r < end:
                    spliced |= op == 3
                    if q is not None:
                        query_positions.append(q)
                        anchor |= r == start
                        flank |= r == end - 1
            elif op == 1 and previous_reference is not None and start <= previous_reference < end - 1:
                query_positions.append(q)
        if spliced or not anchor or not flank or not query_positions:
            continue
        if min(read.query_qualities[q] for q in query_positions) < e.SETTINGS["min_bq"]:
            continue
        sequence = ''.join(read.query_sequence[q] for q in query_positions).upper()
        if set(sequence) - set('ACGT'):
            continue
        counts['depth'] += 1
        sequence = sequence[:-1]
        if sequence == row['REF']:
            counts['ref_count'] += 1
        elif sequence == row['ALT']:
            counts['alt_count'] += 1
    return dict(counts)


def independent_snv_check(rows, bam, star, count):
    """Use htslib count_coverage, independent of the exporter's CIGAR walker."""
    def keep(read):
        if read.flag & e.SETTINGS["exclude_flags"]:
            return False
        mq = read.mapping_quality
        if mq == 255:
            return star and read.has_tag("NH") and read.get_tag("NH") == 1
        return mq >= e.SETTINGS["min_mq"]
    tested = 0
    for row in rows:
        if len(row["REF"]) != 1 or len(row["ALT"]) != 1 or set(row["REF"] + row["ALT"]) - set("ACGT"):
            continue
        counts = bam.count_coverage(row["CHROM"], row["POS"] - 1, row["POS"],
                                    quality_threshold=e.SETTINGS["min_bq"], read_callback=keep)
        counts = {base: int(c[0]) for base, c in zip("ACGT", counts)}
        yield row, {"depth": sum(counts.values()), "ref_count": counts[row["REF"]],
                    "alt_count": counts[row["ALT"]]}
        tested += 1
        if tested >= count:
            return


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--allow-partial", action="store_true")
    ap.add_argument("--check-bam-snvs", type=int, default=0,
                    help="Independent htslib count_coverage checks per sample and modality (development only)")
    ap.add_argument("--check-bam-indels", type=int, default=0,
                    help="Independent htslib aligned-pair indel checks per development sample/modality")
    args = ap.parse_args()
    start = time.time()
    root = args.output.resolve()
    manifest = json.loads((root / "manifest.json").read_text())
    parquet = Path(manifest["parquet"])
    approval, hashes, pools, approved_counts = e.load_release(parquet.parent.parent)
    if hashes != manifest["release_hashes"]:
        raise ValueError("Release identity changed")
    rows_total, classes, splits = 0, Counter(), Counter()
    records = []
    found_samples = set()
    observed_paths = set()
    for pooldir, pool in [("development", "train_pool"), ("reserved", "reserved")]:
        for directory in sorted((root / pooldir).glob("*")):
            if not directory.is_dir():
                raise ValueError(f"Unexpected dataset entry: {directory}")
            sid = directory.name
            if sid not in pools or pools[sid]["pool"] != pool or sid in found_samples:
                raise ValueError("Unapproved, duplicate or cross-pool sample")
            found_samples.add(sid)
            expected = e.sample_labels(parquet, sid, manifest["pilot_rows"])
            info = json.loads((directory / "inputs.json").read_text())
            parts = sorted(directory.glob("part-*.parquet"))
            frames = []
            for part in parts:
                observed_paths.add(part)
                receipt = json.loads(part.with_suffix(".json").read_text())
                if e.sha256(part) != receipt["sha256"]:
                    raise ValueError(f"Output hash mismatch: {part}")
                index = int(part.stem.split("-")[1])
                expected_part = expected.slice(index * manifest["batch_rows"], manifest["batch_rows"])
                expected_id = e.digest({"sample": info["identity"], "keys_labels": expected_part.to_dicts()})
                if expected_id != receipt["identity"]:
                    raise ValueError("Part identity does not match original approved keys/labels")
                actual = pl.read_parquet(part)
                e.reconcile(expected_part, actual, pool)
                e.verify_measurements(actual)
                frames.append(actual)
            if not frames:
                continue
            actual = pl.concat(frames)
            e.unique(actual)  # catches duplication across part files
            complete = len(actual) == len(expected)
            if complete:
                e.reconcile(expected, actual, pool)
            elif not args.allow_partial:
                raise ValueError(f"Incomplete sample: {sid}")
            sample = {"sample_id": sid, "pool": pool, "rows": len(actual),
                      "expected_selected_rows": len(expected), "complete_selected_sample": complete,
                      "classes": dict(Counter(actual["FILTER"])), "modalities": {}}
            for tag in e.MODALITIES:
                prefix = f"bam_{tag}_"
                sample["modalities"][tag] = {
                    "status_counts": dict(Counter(actual[prefix + "status"])),
                    "covered_rows": int((actual[prefix + "depth"] > 0).sum()),
                    "alt_positive_rows": int((actual[prefix + "alt_count"] > 0).sum()),
                    "independent_snv_checks": 0, "independent_indel_checks": 0}
                if (args.check_bam_snvs or args.check_bam_indels) and pool == "train_pool" and info["input_errors"][tag] is None:
                    bam_info = info["bams"][tag]
                    for field in ["bam", "index"]:
                        recorded = bam_info[field]
                        now = e.stat_identity(recorded["path"])
                        if now != {k: recorded[k] for k in now}:
                            raise ValueError("BAM/index identity changed since extraction")
                    measurable = actual.filter(pl.col(prefix + "status").is_in(["ok", "zero_usable_depth"]))
                    with pysam.AlignmentFile(bam_info["bam"]["path"], "rb", index_filename=bam_info["index"]["path"]) as bam:
                        for row, counts in independent_snv_check(measurable.to_dicts(), bam, bam_info["star"], args.check_bam_snvs) if args.check_bam_snvs else []:
                            if any(row[prefix + k] != v for k, v in counts.items()):
                                raise ValueError(f"Independent SNV disagreement: {sid}/{tag}/{row['CHROM']}:{row['POS']}")
                            sample["modalities"][tag]["independent_snv_checks"] += 1
                        if args.check_bam_indels:
                            indels = measurable.filter((pl.col('REF').str.len_chars() != 1) | (pl.col('ALT').str.len_chars() != 1))
                            # Alternate insertion/deletion, positive/zero support; never use FILTER.
                            selected = []
                            for is_ins in [True, False]:
                                for positive in [True, False]:
                                    selected.extend(indels.filter(
                                        (pl.col('ALT').str.len_chars() > pl.col('REF').str.len_chars()) == is_ins
                                    ).filter((pl.col(prefix + 'alt_count') > 0) == positive).head(max(1, args.check_bam_indels // 4)).to_dicts())
                            for row in selected[:args.check_bam_indels]:
                                counts = independent_indel_counts(row, bam, bam_info['star'])
                                if any(row[prefix + k] != v for k, v in counts.items()):
                                    raise ValueError(f"Independent indel disagreement: {sid}/{tag}/{row['CHROM']}:{row['POS']}/{row['REF']}/{row['ALT']}: {counts}")
                                sample['modalities'][tag]['independent_indel_checks'] += 1
            rows_total += len(actual)
            classes.update(actual["FILTER"])
            splits.update(actual["split"])
            records.append(sample)
    if set(root.glob("**/*.parquet")) != observed_paths:
        raise ValueError("Unregistered parquet outside sample partitions")
    full = rows_total == approval["total_records"] and found_samples == pools.keys() and not manifest["pilot_rows"]
    if full and dict(classes) != e.EXPECTED:
        raise ValueError("Class totals mismatch")
    if not full and not args.allow_partial:
        raise ValueError("Not a complete approved-release evidence dataset")
    result = {"full_release_complete": full, "output_rows": rows_total,
              "approved_release_rows": approval["total_records"],
              "remaining_release_rows": approval["total_records"] - rows_total,
              "classes": classes, "splits": splits, "samples": records,
              "runtime_seconds": time.time() - start, "verifier_sha256": e.sha256(__file__),
              "verification": "PASS: full-key/FILTER, pool, split, output hashes, numerical invariants",
              "independent_snv_checks": sum(v["independent_snv_checks"] for s in records for v in s["modalities"].values()),
              "independent_indel_checks": sum(v["independent_indel_checks"] for s in records for v in s["modalities"].values())}
    e.atomic_json(root / "verification.json", result)
    lines = ["# Variant evidence verification", "",
             f"Full approved release complete: **{full}**. Output rows: **{rows_total:,}** / {approval['total_records']:,}.",
             "", result["verification"], "", f"Independent htslib SNV checks: {result['independent_snv_checks']}; indel checks: {result['independent_indel_checks']}.", "",
             "| Sample | Pool | Rows | Modality | Covered | ALT positive | Status counts |",
             "|---|---|---:|---|---:|---:|---|"]
    for s in records:
        for tag, v in s["modalities"].items():
            lines.append(f"| {s['sample_id']} | {s['pool']} | {s['rows']:,} | {tag} | {v['covered_rows']:,} | {v['alt_positive_rows']:,} | {json.dumps(v['status_counts'], sort_keys=True)} |")
    (root / "verification.md").write_text("\n".join(lines) + "\n")
    print(json.dumps({k: v for k, v in result.items() if k != "samples"}, indent=2))


if __name__ == "__main__":
    main()
