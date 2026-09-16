#!/usr/bin/env python3
"""Experimental post-annotation rescue adapter; never enabled by default.

Uses existing VCFs only. Stores the union on disk, preserving negative records
and source provenance. Outputs require label QC before any training use.
"""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import sqlite3
from urllib.parse import quote
import zlib

import pysam

from vcf_utils.aggregation import resolve_tumor_sample_index, _normal_sample_from_header
from vcf_utils.refined_rescue_policy import CALLERS, biological_veto, decide, nominates, rna_supports

LABELS = {"Somatic", "Germline", "Reference", "Artifact", "NoConsensus", "RNAedit"}


def digest(path):
    value = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


def panel(values):
    result = {}
    for value in values:
        caller, sep, path = value.partition("=")
        if not sep or caller not in CALLERS or caller in result:
            raise ValueError("Require unique caller=path entries for the three callers")
        result[caller] = Path(path).resolve(strict=True)
    if set(result) != CALLERS or len(set(result.values())) != len(CALLERS):
        raise ValueError("Require three distinct caller VCFs")
    return result


def info_dict(raw):
    return {part.partition("=")[0]: part.partition("=")[2] if "=" in part else None
            for part in raw.split(";") if part not in ("", ".")}


def transition(dna_label, rescue_label, allowed, reason, verification):
    if dna_label == "Somatic":
        return "Somatic", "retained_dna_baseline"
    if rescue_label and rescue_label != "Somatic":
        return rescue_label, "retained_rescue_negative"
    if dna_label in {"Artifact", "Germline", "Reference", "RNAedit"}:
        return dna_label, "protected_dna_negative"
    if allowed and verification not in {"rejected", "inconclusive"}:
        return "Somatic", reason
    return dna_label or "NoConsensus", "verification_withheld" if allowed else reason


def stage(db, path, kind):
    with pysam.VariantFile(str(path)) as reader:
        if reader.header.samples:
            raise ValueError(f"Require sampleless consensus/rescue VCF: {path}")
        header = reader.header.copy()
        if 'GATE_POLICY' in header.info:
            raise ValueError('Refuse recursively gating a previous adapter output')
        for record in reader:
            if record.contig not in header.contigs:
                raise ValueError(f'Undeclared contig: {record.contig}')
            raw = str(record).rstrip("\n")
            parts = raw.split("\t")
            if len(parts) != 8 or parts[6] not in LABELS:
                raise ValueError(f"Require biological label VCF: {path}: {parts[:7]}")
            key = (parts[0], int(parts[1]), parts[3], parts[4])
            db.execute("INSERT OR IGNORE INTO variants(chrom,pos,ref,alt) VALUES (?,?,?,?)", key)
            if db.execute(f"SELECT {kind} FROM variants WHERE chrom=? AND pos=? AND ref=? AND alt=?", key).fetchone()[0] is not None:
                raise ValueError(f"Duplicate {kind} allele: {key}")
            db.execute(f"UPDATE variants SET {kind}=? WHERE chrom=? AND pos=? AND ref=? AND alt=?",
                       (zlib.compress(raw.encode()), *key))
        db.commit()
    return header


def scan_panel(db, inputs, modality):
    column = "dna_votes" if modality == "DNA" else "rna_votes"
    for caller, path in inputs.items():
        with pysam.VariantFile(str(path)) as reader:
            names = list(reader.header.samples)
            if not 1 <= len(names) <= 2:
                raise ValueError(f"Require paired or tumor-only caller: {path}")
            idx = resolve_tumor_sample_index(names, caller, _normal_sample_from_header(str(reader.header)))
            for record in reader:
                if len(record.ref) != 1:
                    continue
                sample = record.samples[names[idx]]
                filt = ";".join(record.filter) or "."
                for i, alt in enumerate(record.alts or [], 1):
                    if len(alt) != 1 or alt not in "ACGT":
                        continue
                    ad = sample.get("AD")
                    count = ad[i] if ad and i < len(ad) else None
                    if ad is None and caller == "strelka":
                        values = sample.get(alt + "U")
                        count = values[0] if values else None
                    allowed = nominates(caller, filt, count) if modality == "DNA" else rna_supports(caller, filt, count)
                    if allowed:
                        db.execute(f"UPDATE variants SET {column}={column} || ? WHERE chrom=? AND pos=? AND ref=? AND alt=? AND rescue IS NOT NULL",
                                   (caller + "|", record.contig, record.pos, record.ref, alt))
        db.commit()


def write_output(db, header, out, alignment_round):
    definitions = {
        "GATE_POLICY": "Experimental policy identifier",
        "GATE_ALIGNMENT_ROUND": "RNA round supplied for this evaluation; not a second independent vote",
        "GATE_DNA_NOMINATORS": "Unique DNA callers with positive native nomination evidence",
        "GATE_RNA_ELIGIBLE": "Unique RNA callers with native PASS and at least three tumor alternate reads",
        "GATE_SOURCE_RECORD": "Percent-encoded original chosen source record before relabeling",
        "GATE_DNA_RECORD": "Percent-encoded refined DNA consensus record when rescue is the chosen source",
        "UNIFIED_FILTER": "Current biological label",
        "UNIFIED_FILTER_DNA": "Refined DNA consensus biological label",
        "CLASSIFICATION_RATIONALE": "Current label decision trace",
        "PASSES_CONSENSUS_DNA": "Whether refined DNA consensus is Somatic",
        "RESCUED": "Passed cross-modality or newly admitted rescue",
        "RESCUE_PROMOTED": "New Somatic admission relative to refined DNA consensus",
    }
    for name, description in definitions.items():
        if name not in header.info:
            header.add_meta("INFO", items=[("ID", name), ("Number", 1), ("Type", "String"), ("Description", description)])
        elif header.info[name].number != 1 or header.info[name].type != "String":
            raise ValueError(f"Incompatible INFO definition: {name}")
    for label in LABELS:
        if label not in header.filters:
            header.add_meta("FILTER", items=[("ID", label), ("Description", label)])
    header.add_meta("refined_rescue_status", value="experimental_not_training_ready")
    counts = Counter()
    with pysam.BGZFile(str(out), "w") as output:
        output.write(str(header).encode())
        for chrom in header.contigs:
            cursor = db.execute("SELECT ref,alt,dna,rescue,dna_votes,rna_votes FROM variants WHERE chrom=? ORDER BY pos,ref,alt", (chrom,))
            for ref, alt, dna, rescue, dna_votes, rna_votes in cursor:
                dna_raw = zlib.decompress(dna).decode() if dna else None
                rescue_raw = zlib.decompress(rescue).decode() if rescue else None
                raw = rescue_raw or dna_raw
                parts = raw.split("\t")
                info = info_dict(parts[7])
                dna_label = dna_raw.split("\t")[6] if dna_raw else None
                rescue_label = parts[6] if rescue_raw else None
                dna_set, rna_set = set(dna_votes.split("|")) - {""}, set(rna_votes.split("|")) - {""}
                allowed, reason = decide(ref, alt, rescue_label, dna_set, rna_set, info)
                label, reason = transition(dna_label, rescue_label, allowed, reason, info.get("DNA_VERIFICATION"))
                promoted = label == "Somatic" and dna_label != "Somatic"
                if dna_label == "Somatic" and biological_veto(info):
                    counts["baseline_annotation_conflicts"] += 1
                if dna_label == "Somatic" and not rescue_raw:
                    counts["baseline_somatic_without_rescue_record"] += 1
                if dna_label == "Somatic" and rescue_label and rescue_label != "Somatic":
                    counts["baseline_vs_rescue_negative_conflicts"] += 1
                info.update(GATE_POLICY="seqc2_refined_gate_v1", GATE_ALIGNMENT_ROUND=alignment_round,
                            GATE_DNA_NOMINATORS="|".join(sorted(dna_set)) or ".",
                            GATE_RNA_ELIGIBLE="|".join(sorted(rna_set)) or ".",
                            GATE_SOURCE_RECORD=quote(raw, safe=""), UNIFIED_FILTER=label,
                            UNIFIED_FILTER_DNA=dna_label or "NoConsensus",
                            PASSES_CONSENSUS_DNA="YES" if dna_label == "Somatic" else "NO",
                            CLASSIFICATION_RATIONALE=f"rule:seqc2_refined_gate_v1|branch:{reason}|class:{label}",
                            RESCUE_PROMOTED="YES" if promoted else "NO",
                            RESCUED="YES" if label == "Somatic" and (promoted or info.get("PASSES_CONSENSUS_RNA") == "YES") else "NO")
                if dna_raw and rescue_raw:
                    info["GATE_DNA_RECORD"] = quote(dna_raw, safe="")
                parts[6] = label
                parts[7] = ";".join(k if v is None else f"{k}={v}" for k, v in info.items())
                output.write(("\t".join(parts) + "\n").encode())
                counts[label] += 1
                counts[reason] += 1
                counts['records_written'] += 1
    if counts['records_written'] != db.execute('SELECT COUNT(*) FROM variants').fetchone()[0]:
        raise ValueError('Output did not cover the entire input union')
    pysam.tabix_index(str(out), preset="vcf", force=False)
    return dict(counts)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dna-consensus", type=Path, required=True)
    ap.add_argument("--annotated-rescue", type=Path, required=True)
    ap.add_argument("--dna-vcf", action="append", required=True, help="caller=path")
    ap.add_argument("--rna-vcf", action="append", required=True, help="caller=path")
    ap.add_argument("--alignment-round", choices=["first", "realignment"], required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    args = ap.parse_args()
    dna, rna = panel(args.dna_vcf), panel(args.rna_vcf)
    sources = [args.dna_consensus.resolve(strict=True), args.annotated_rescue.resolve(strict=True), *dna.values(), *rna.values()]
    args.outdir.mkdir(parents=True, exist_ok=False)
    report = {"status": "running", "alignment_round": args.alignment_round,
              "sources": {str(p): digest(p) for p in sources}, "script_sha256": digest(Path(__file__))}
    report['code'] = {str(p): digest(p) for p in (Path(__file__),
                       Path(__file__).parent / 'vcf_utils/refined_rescue_policy.py',
                       Path(__file__).parent / 'vcf_utils/aggregation.py')}
    try:
        with sqlite3.connect(args.outdir / "union.sqlite") as db:
            db.execute("PRAGMA cache_size=-32768")
            db.execute("CREATE TABLE variants(chrom TEXT,pos INTEGER,ref TEXT,alt TEXT,dna BLOB,rescue BLOB,dna_votes TEXT DEFAULT '',rna_votes TEXT DEFAULT '',PRIMARY KEY(chrom,pos,ref,alt))")
            header = stage(db, sources[0], "dna")
            rescue_header = stage(db, sources[1], "rescue")
            for name in set(header.contigs) & set(rescue_header.contigs):
                lengths = header.contigs[name].length, rescue_header.contigs[name].length
                if all(lengths) and lengths[0] != lengths[1]:
                    raise ValueError(f'Conflicting source contig length: {name}')
            for name in set(header.info) & set(rescue_header.info):
                if (header.info[name].number, header.info[name].type) != (rescue_header.info[name].number, rescue_header.info[name].type):
                    raise ValueError(f"Conflicting source INFO definition: {name}")
            header.merge(rescue_header)
            scan_panel(db, dna, "DNA")
            scan_panel(db, rna, "RNA")
            report["counts"] = write_output(db, header, args.outdir / "refined.rescue.vcf.gz", args.alignment_round)
        report["sources_unchanged"] = all(digest(Path(p)) == sha for p, sha in report["sources"].items())
        if not report["sources_unchanged"]:
            raise ValueError("Source integrity changed")
        report["output_sha256"] = digest(args.outdir / "refined.rescue.vcf.gz")
        report["status"] = "experimental_labels_not_training_ready"
    except Exception as exc:
        report.update(status="failed", error=str(exc))
        raise
    finally:
        (args.outdir / "report.json").write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
