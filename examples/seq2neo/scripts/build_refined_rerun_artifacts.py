#!/usr/bin/env python3
"""Build lightweight provenance and a Somatic-candidate variant Parquet.

Original cohort inputs remain read-only. The output Parquet is intentionally
written outside Git; this script commits only TSV/JSON path and summary metadata.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess
from collections import Counter
import tempfile

import polars as pl
import pysam

# Keep the committed/ignored variant table compact and training-label focused.
# Full FILTER distributions (including Germline/Reference) remain in summary.json.
FILTERS = {"Somatic"}


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def info_value(record, key):
    value = record.info.get(key)
    if isinstance(value, tuple):
        return ",".join(str(x) for x in value)
    return "" if value in (None, ".") else str(value)


def variant_type(ref, alt):
    if len(ref) == 1 and len(alt) == 1:
        return "SNP"
    if len(ref) == len(alt):
        return "MNP"
    return "INDEL"


def scan_vcf(path, collect_rows=False, sample_id="", row_writer=None):
    counts, rows = {}, []
    with pysam.VariantFile(str(path)) as reader:
        for record in reader:
            filters = list(record.filter)
            label = filters[0] if filters else "PASS"
            counts[label] = counts.get(label, 0) + 1
            if collect_rows and label in FILTERS:
                alt = str(record.alts[0]) if record.alts else ""
                row = {"sample_id": sample_id, "CHROM": record.contig,
                    "POS": int(record.pos), "REF": str(record.ref), "ALT": alt,
                    "FILTER": label, "variant_type": variant_type(str(record.ref), alt),
                    "RESCUED": info_value(record, "RESCUED"),
                    "RESCUE_PROMOTED": info_value(record, "RESCUE_PROMOTED"),
                    "N_DNA_CALLERS_SUPPORT": info_value(record, "N_DNA_CALLERS_SUPPORT"),
                    "N_RNA_CALLERS_SUPPORT": info_value(record, "N_RNA_CALLERS_SUPPORT"),
                    "GNOMAD_AF": info_value(record, "GNOMAD_AF"),
                    "candidate_status": "candidate_not_training_approved"}
                if row_writer is None:
                    rows.append(row)
                else:
                    row_writer.writerow(row)
    return counts, rows


def count_filters_fast(path):
    proc = subprocess.Popen(["bcftools", "query", "-f", "%FILTER\\n", str(path)], stdout=subprocess.PIPE, text=True)
    counts = Counter()
    assert proc.stdout is not None
    for line in proc.stdout:
        counts[line.strip() or "PASS"] += 1
    if proc.wait() != 0:
        raise RuntimeError(f"bcftools query failed: {path}")
    return dict(counts)


def stream_somatic_rows(path, sample_id, row_writer):
    proc = subprocess.Popen(
        ["bcftools", "query", "-i", 'FILTER="Somatic"', "-f", "%CHROM\\t%POS\\t%REF\\t%FIRST_ALT\\t%FILTER\\n", str(path)],
        stdout=subprocess.PIPE, text=True,
    )
    assert proc.stdout is not None
    count = 0
    for line in proc.stdout:
        chrom, pos, ref, alt, filt = line.rstrip("\n").split("\t")
        row_writer.writerow({"sample_id": sample_id, "CHROM": chrom, "POS": int(pos),
            "REF": ref, "ALT": alt, "FILTER": filt, "variant_type": variant_type(ref, alt),
            "RESCUED": "", "RESCUE_PROMOTED": "", "N_DNA_CALLERS_SUPPORT": "",
            "N_RNA_CALLERS_SUPPORT": "", "GNOMAD_AF": "", "candidate_status": "candidate_not_training_approved"})
        count += 1
    if proc.wait() != 0:
        raise RuntimeError(f"bcftools Somatic query failed: {path}")
    return count


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source-manifest", type=Path, default=Path("examples/seq2neo/data/processed/sample_manifest.tsv"))
    ap.add_argument("--output-root", type=Path, default=Path("/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/rnadnavar/examples/seq2neo/output_refined_native_v2_20260916"))
    ap.add_argument("--outdir", type=Path, default=Path("examples/seq2neo/data/processed/refined_native_v2_20260917"))
    ap.add_argument("--parquet", type=Path, default=Path("/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/rnadnavar/examples/seq2neo/output_refined_native_v2_20260916/variant_parquet/refined_candidates.parquet"))
    args = ap.parse_args(); args.outdir.mkdir(parents=True, exist_ok=True)
    with args.source_manifest.open() as handle: source_rows = list(csv.DictReader(handle, delimiter="\t"))
    summaries = []
    schema = {"sample_id": pl.Utf8, "CHROM": pl.Utf8, "POS": pl.Int64, "REF": pl.Utf8, "ALT": pl.Utf8, "FILTER": pl.Utf8, "variant_type": pl.Utf8, "RESCUED": pl.Utf8, "RESCUE_PROMOTED": pl.Utf8, "N_DNA_CALLERS_SUPPORT": pl.Utf8, "N_RNA_CALLERS_SUPPORT": pl.Utf8, "GNOMAD_AF": pl.Utf8, "candidate_status": pl.Utf8}
    row_fields = list(schema)
    candidate_rows = 0
    temp = tempfile.NamedTemporaryFile("w", suffix=".tsv", prefix="refined_candidates_", delete=False, newline="")
    temp_path = Path(temp.name)
    row_writer = csv.DictWriter(temp, fieldnames=row_fields, delimiter="\t")
    row_writer.writeheader()
    for source in source_rows:
        sid = source["sample_id"]; state_path = args.output_root / sid / "state.json"
        if not state_path.is_file(): raise SystemExit(f"missing completed state: {sid}")
        state = json.loads(state_path.read_text())
        if state.get("status") != "candidate_complete_not_training_approved": raise SystemExit(f"sample not complete: {sid}")
        new_consensus = Path(state["output"]) / "refined.vcf.gz"; new_rescue = Path(state["output"]) / "refined.rescue.vcf.gz"; original = Path(source["rescue_vcf_path"])
        if not original.is_file(): raise SystemExit(f"missing original rescue: {sid}: {original}")
        old_counts = count_filters_fast(original); consensus_counts = count_filters_fast(new_consensus); new_counts = count_filters_fast(new_rescue)
        candidate_rows += stream_somatic_rows(new_rescue, sid, row_writer)
        source_hashes = state.get("source_hashes", {})
        output_hashes = state.get("outputs", {})
        original_hash = source_hashes.get(str(original)) or sha(original)
        consensus_hash = output_hashes.get(str(new_consensus)) or sha(new_consensus)
        rescue_hash = output_hashes.get(str(new_rescue)) or sha(new_rescue)
        summaries.append({"sample_id": sid, "status": state["status"], "original_rescue": str(original), "new_consensus": str(new_consensus), "new_rescue": str(new_rescue), "original_sha256": original_hash, "new_consensus_sha256": consensus_hash, "new_rescue_sha256": rescue_hash, "original_counts": old_counts, "consensus_counts": consensus_counts, "new_rescue_counts": new_counts, "original_records": sum(old_counts.values()), "new_rescue_records": sum(new_counts.values()), "comparison_note": "Exact allele overlap omitted from lightweight scan; use VCFs for allele-level comparison"})
    temp.close()
    args.parquet.parent.mkdir(parents=True, exist_ok=True)
    pl.scan_csv(temp_path, separator="\t", schema_overrides=schema).sink_parquet(args.parquet, compression="zstd")
    temp_path.unlink(missing_ok=True)
    manifest = args.outdir / "manifest.tsv"
    fields = ["sample_id", "original_rescue", "new_consensus", "new_rescue", "new_rescue_sha256", "candidate_status"]
    with manifest.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t"); writer.writeheader()
        for item in summaries: writer.writerow({**{k: item[k] for k in fields if k != "candidate_status"}, "candidate_status": "candidate_not_training_approved"})
    summary = {"schema_version": 1, "samples": len(summaries), "candidate_rows": candidate_rows, "parquet": str(args.parquet), "manifest": str(manifest), "variant_filter": "Somatic", "status": "candidate_not_training_approved", "source_manifest": str(args.source_manifest), "samples_detail": summaries}
    (args.outdir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({"samples": len(summaries), "candidate_rows": candidate_rows, "parquet": str(args.parquet), "manifest": str(manifest)}, indent=2))


if __name__ == "__main__": main()
