#!/usr/bin/env python3
"""Export all three candidate label classes, without asserting training approval.

Streams batches into Parquet. Preserves source INFO rather than filling evidence
columns with blanks. Completed audit counts are reused only after hash checks.
"""
import argparse
from collections import Counter
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

import pyarrow as pa
import pyarrow.parquet as pq

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / "bin"))
from vcf_utils.refined_rescue_policy import biological_veto

LABELS = ("Somatic", "Germline", "Reference")
STATUS = "candidate_not_training_approved"
SHARED = Path("/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/rnadnavar/examples/seq2neo")
INFO_FIELDS = (
    "RESCUED", "RESCUE_PROMOTED", "N_DNA_CALLERS_SUPPORT", "N_RNA_CALLERS_SUPPORT",
    "GNOMAD_AF", "GT_BY_CALLER", "DP_BY_CALLER", "AD_BY_CALLER", "VAF_BY_CALLER",
    "NORMAL_GT_BY_CALLER", "NORMAL_DP_BY_CALLER", "NORMAL_AD_BY_CALLER",
    "NORMAL_VAF_BY_CALLER", "CLASSIFICATION_RATIONALE", "GATE_DNA_NOMINATORS",
    "GATE_RNA_ELIGIBLE", "GATE_ALIGNMENT_ROUND", "UNIFIED_FILTER_DNA",
    "PASSES_CONSENSUS_DNA", "DNA_VERIFICATION", "REDI_ACCESSION",
    "REDI_CANONICAL", "N_DNA_CALLERS_SOMATIC",
    "NEGATIVE_EVIDENCE_POLICY", "NEGATIVE_EVIDENCE_STATUS", "NEGATIVE_EVIDENCE_REASON",
    "NEGATIVE_NORMAL_COUNTS", "NEGATIVE_TUMOR_COUNTS",
    "THREE_CLASS_POLICY", "THREE_CLASS_BASELINE_FILTER", "THREE_CLASS_NATIVE_FILTER",
    "THREE_CLASS_BASELINE_RATIONALE", "THREE_CLASS_NATIVE_RATIONALE",
    "THREE_CLASS_REVIEW_REASON", "TRAINING_ELIGIBLE",
)
NUMERIC_FIELDS = ("DP_DNA_MEAN", "DP_RNA_MEAN", "VAF_DNA_MEAN", "VAF_RNA_MEAN")
SCHEMA = pa.schema(
    [(k, pa.string()) for k in ("sample_id", "CHROM", "REF", "ALT", "FILTER", "variant_type",
                               "candidate_status", "label_confidence", "review_reason")]
    + [("POS", pa.int64()), ("training_eligible", pa.bool_())]
    + [(k, pa.string()) for k in INFO_FIELDS]
    + [(k, pa.float64()) for k in NUMERIC_FIELDS]
)

def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()

def variant_type(ref, alt):
    if "," in alt or set(ref + alt) - set("ACGTN"):
        return "OTHER"
    if len(ref) == len(alt) == 1:
        return "SNP"
    return "MNP" if len(ref) == len(alt) else "INDEL"

def review_reason(label, info):
    veto = biological_veto(info)
    if label == "Somatic" and veto:
        return "annotation_conflict:" + veto
    if label == "Somatic" and info.get("DNA_VERIFICATION") in {"rejected", "inconclusive"}:
        return "dna_verification_conflict"
    if label == "Somatic" and info.get("THREE_CLASS_REVIEW_REASON") not in (None,"none"):
        return "three_class_review:" + info["THREE_CLASS_REVIEW_REASON"]
    if label in {"Germline", "Reference"}:
        if info.get("NEGATIVE_EVIDENCE_STATUS") == "SUPPORTED":
            return "negative_read_supported_biological_approval_pending"
        if info.get("NEGATIVE_EVIDENCE_STATUS") in {"WITHHELD", "CONFLICT"}:
            return "negative_evidence_withheld:" + (info.get("NEGATIVE_EVIDENCE_REASON") or "unspecified")
        if info.get("THREE_CLASS_POLICY") == "separated_three_class_v2":
            return "native_negative_requires_paired_validation"
        return "inherited_or_legacy_negative_requires_paired_validation"
    return "somatic_candidate_requires_label_qc"

def stream_rows(path, sid, writer, batch_size=50000):
    """One record per VCF row; retain full ALT, including multi-allelic values."""
    fields = INFO_FIELDS + NUMERIC_FIELDS
    fmt = "%CHROM\\t%POS\\t%REF\\t%ALT\\t%FILTER" + "".join(
        "\\t%INFO/" + name for name in fields
    ) + "\\n"
    expr = 'FILTER="Somatic" || FILTER="Germline" || FILTER="Reference"'
    proc = subprocess.Popen(["bcftools", "query", "-u", "-i", expr, "-f", fmt, str(path)],
                            stdout=subprocess.PIPE, text=True)
    counts, types, reasons = Counter(), Counter(), Counter()
    batch = []
    try:
        for line in proc.stdout:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 5 + len(fields):
                raise ValueError(f"Unexpected query columns: {sid}")
            chrom, pos, ref, alt, label = parts[:5]
            if label not in LABELS:
                raise ValueError(f"Ambiguous FILTER: {label}")
            info = {k: None if v == "." else v for k, v in zip(fields, parts[5:])}
            reason = review_reason(label, info)
            row = dict(sample_id=sid, CHROM=chrom, POS=int(pos), REF=ref, ALT=alt,
                       FILTER=label, variant_type=variant_type(ref, alt),
                       candidate_status=STATUS, label_confidence="unvalidated",
                       training_eligible=False, review_reason=reason,
                       **{k: info[k] for k in INFO_FIELDS})
            for key in NUMERIC_FIELDS:
                value = info[key]
                row[key] = float(value) if value is not None else None
            batch.append(row)
            counts[label] += 1
            types[label + ":" + row["variant_type"]] += 1
            reasons[reason] += 1
            if len(batch) >= batch_size:
                writer.write_table(pa.Table.from_pylist(batch, schema=SCHEMA))
                batch.clear()
        if proc.wait() != 0:
            raise RuntimeError(f"bcftools query failed: {path}")
        if batch:
            writer.write_table(pa.Table.from_pylist(batch, schema=SCHEMA))
    finally:
        proc.stdout.close()
        if proc.poll() is None:
            proc.terminate()
            proc.wait()
    return dict(counts), dict(types), dict(reasons)

def checked_audit(state, path, name):
    expected = state["outputs"][str(path)]
    if sha(path) != expected:
        raise ValueError(f"Changed VCF: {path}")
    audit_path = path.parent / (name + ".audit.json")
    if sha(audit_path) != state["outputs"][str(audit_path)]:
        raise ValueError(f"Changed audit: {audit_path}")
    audit = json.loads(audit_path.read_text())
    if audit["issues"] or not audit["sources_unchanged"] or audit["sha256"] != expected:
        raise ValueError(f"Failed audit: {audit_path}")
    return {k: v for k, v in audit["counts"].items() if k != "records"}


def candidate_paths(state):
    folder=Path(state['output'])
    paths=state.get('final_artifacts')
    if state.get('identity',{}).get('policy')=='separated_three_class_v2' and not paths:
        raise ValueError('Missing separated-policy final artifacts; refuse baseline fallback')
    if paths is None:
        return folder/'refined.rescue.vcf.gz',folder/'refined.vcf.gz'
    if set(paths)!={'consensus','rescue'}:
        raise ValueError('Incomplete final candidate paths')
    selected=[Path(paths[name]) for name in ('rescue','consensus')]
    if any(p.resolve().parent!=folder.resolve() for p in selected):
        raise ValueError('Final artifact outside recorded output directory')
    if state.get('identity',{}).get('policy')=='separated_three_class_v2' and any(
            Path(paths[k])!=folder/f'three_class.{k}.vcf.gz' for k in paths):
        raise ValueError('Separated policy cannot export an old baseline as final')
    return tuple(selected)


def sample_summary(old,state,consensus,rescue,cc,rc,counts,types,reasons):
    """Keep original provenance, but never pair new counts with old result paths."""
    return dict(old,new_consensus=str(consensus),new_rescue=str(rescue),
                new_consensus_sha256=state['outputs'][str(consensus)],
                new_rescue_sha256=state['outputs'][str(rescue)],new_rescue_records=sum(rc.values()),
                consensus_counts=cc,new_rescue_counts=rc,exported_counts=counts,
                exported_types=types,review_reasons=reasons,status=STATUS,
                candidate_policy=state.get('identity',{}).get('policy','unrecorded'),
                prior_comparison_note=old.get('comparison_note'),
                comparison_note='Original counts are a prior snapshot; current candidates are not approved truth')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--source-manifest", type=Path, default=REPO/"examples/seq2neo/data/processed/sample_manifest.tsv")
    ap.add_argument("--prior-summary", type=Path, default=REPO/"examples/seq2neo/data/processed/refined_native_v2_20260917/summary.json")
    ap.add_argument("--output-root", type=Path, default=SHARED/"output_refined_native_v2_20260916")
    ap.add_argument("--outdir", type=Path, default=REPO/"examples/seq2neo/data/processed/refined_three_class_20260918")
    ap.add_argument("--parquet", type=Path, default=SHARED/"output_refined_native_v2_20260916/variant_parquet/refined_three_class_candidates.parquet")
    args = ap.parse_args()
    if args.parquet.exists() or args.outdir.exists():
        raise ValueError("Use fresh destinations; existing exports are never overwritten")
    with args.source_manifest.open() as handle:
        sources = list(csv.DictReader(handle, delimiter="\t"))
    if len({s["sample_id"] for s in sources}) != len(sources):
        raise ValueError("Duplicate sample IDs")
    prior = {x["sample_id"]: x for x in json.loads(args.prior_summary.read_text())["samples_detail"]}
    args.parquet.parent.mkdir(parents=True, exist_ok=True)
    tmp = args.parquet.with_name(args.parquet.name + f".{os.getpid()}.partial")
    summaries, manifests = [], []
    totals, type_totals, reason_totals = Counter(), Counter(), Counter()
    try:
        with pq.ParquetWriter(tmp, SCHEMA, compression="zstd") as writer:
            for i, src in enumerate(sources, 1):
                sid = src["sample_id"]
                state = json.loads((args.output_root/sid/"state.json").read_text())
                if state["status"] != "candidate_complete_not_training_approved":
                    raise ValueError(f"Incomplete sample: {sid}")
                folder = Path(state["output"])
                rescue, consensus = candidate_paths(state)
                rc = checked_audit(state, rescue, "rescue")
                cc = checked_audit(state, consensus, "consensus")
                counts, types, reasons = stream_rows(rescue, sid, writer)
                if any(counts.get(k, 0) != rc.get(k, 0) for k in LABELS):
                    raise ValueError(f"Export/audit class-count mismatch: {sid}")
                if sha(rescue) != state["outputs"][str(rescue)]:
                    raise ValueError(f"Rescue changed during export: {sid}")
                old = prior[sid]
                if old["original_rescue"] != src["rescue_vcf_path"]:
                    raise ValueError(f"Original path changed: {sid}")
                item = sample_summary(old,state,consensus,rescue,cc,rc,counts,types,reasons)
                summaries.append(item)
                manifests.append(dict(src, original_rescue_vcf_path=src["rescue_vcf_path"],
                    rescue_vcf_path=str(rescue), consensus_vcf_path=str(consensus),
                    variant_parquet_path=str(args.parquet), candidate_status=STATUS,
                    candidate_policy=item['candidate_policy'],
                    label_qc_verdict="NOT_APPROVED", training_label_vcf="",
                    new_rescue_sha256=state["outputs"][str(rescue)]))
                totals.update(counts); type_totals.update(types); reason_totals.update(reasons)
                print(f"[{i}/{len(sources)}] {sid}: {counts}", flush=True)
        if pq.ParquetFile(tmp).metadata.num_rows != sum(totals.values()):
            raise ValueError("Parquet row count mismatch")
        args.outdir.mkdir(parents=True)
        manifest = args.outdir/"manifest.tsv"
        with manifest.open("w", newline="") as handle:
            w = csv.DictWriter(handle, fieldnames=list(manifests[0]), delimiter="\t", lineterminator="\n")
            w.writeheader(); w.writerows(manifests)
        # Small sample Parquet is kept next to the heavy variant Parquet, outside Git.
        pq.write_table(pa.Table.from_pylist(manifests), args.parquet.parent/"refined_three_class_manifest.parquet", compression="zstd")
        tmp.rename(args.parquet)
        summary = dict(schema_version=3, samples=len(sources), candidate_rows=sum(totals.values()),
                       class_counts=dict(totals), class_type_counts=dict(type_totals),
                       review_reasons=dict(reason_totals), variant_filters=list(LABELS),
                       parquet=str(args.parquet), parquet_sha256=sha(args.parquet),
                       status=STATUS, prior_summary=str(args.prior_summary),
                       original_counts_provenance="Recorded prior summary; not freshly rescanned",
                       samples_detail=summaries)
        (args.outdir/"summary.json").write_text(json.dumps(summary, indent=2)+"\n")
        print(json.dumps({k: v for k, v in summary.items() if k != "samples_detail"}, indent=2))
    finally:
        tmp.unlink(missing_ok=True)

if __name__ == "__main__":
    main()
