#!/usr/bin/env python3
"""Compare frozen metadata rescue with fresh caller-native evidence, no truth.

No VCF labels are rewritten, no workflow is launched, and no HG008 is read.
The historical comparator intentionally retains its recorded support gates.
"""
import argparse
from collections import Counter
import json
from pathlib import Path
import sys

import pysam

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / "bin"))
from vcf_utils.aggregation import resolve_tumor_sample_index, _normal_sample_from_header
from vcf_utils.refined_rescue_policy import biological_veto, decide, integer, nominates, rna_supports
from audit_current_native_policy import digest
from apply_optimized_rescue_gate import read_keys, fields, info_map

CALLERS = ("deepsomatic", "mutect2", "strelka")


def caller_evidence(path, caller, candidates, dna):
    observed = set()
    eligible = set()
    with pysam.VariantFile(str(path)) as reader:
        samples = list(reader.header.samples)
        if not 1 <= len(samples) <= 2:
            raise ValueError(f"Expected tumor-only or paired caller VCF: {path}")
        index = resolve_tumor_sample_index(samples, caller, _normal_sample_from_header(str(reader.header)))
        for record in reader:
            filt = ";".join(record.filter) or "."
            for allele_index, alt in enumerate(record.alts or [], 1):
                key = record.contig, record.pos, record.ref, alt
                if key not in candidates:
                    continue
                sample = record.samples[samples[index]]
                ad = sample.get("AD")
                count = ad[allele_index] if ad and len(ad) > allele_index else None
                if ad is None and caller == "strelka":
                    field = alt + "U" if len(alt) == len(record.ref) == 1 else "TIR"
                    values = sample.get(field)
                    count = values[0] if values else None
                if not dna:
                    if filt in {"PASS", ".", "Somatic"}:
                        observed.add(key)
                    if rna_supports(caller, filt, count):
                        eligible.add(key)
                    continue
                if nominates(caller, filt, count):
                    observed.add(key)
    return observed if dna else (observed, eligible)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--outdir", type=Path, required=True)
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=False)
    root = REPO / "examples/seqc2"
    report = {"status": "running", "sources": {}, "datasets": {}, "code": {}}
    for path in (Path(__file__), REPO / "bin/vcf_utils/refined_rescue_policy.py",
                 REPO / "bin/vcf_utils/aggregation.py"):
        report["code"][str(path)] = digest(path)

    def track(path):
        path = path.resolve(strict=True)
        report["sources"][str(path)] = digest(path)
        return path

    for dataset in ("wes_ll", "wgs_il"):
        bundle = root / "verified/20260914" / dataset
        native = read_keys(track(bundle / "native_consensus.vcf.gz"), {"Somatic", "PASS", "."})
        expected = read_keys(track(root / "comparison/rescue_fp_investigation_20260914/gate_tests"
                                   / dataset / "nomination_biological/query.vcf.gz")) - native
        pattern = "*/*.filtered.vcf.gz" if dataset == "wes_ll" else "*/*.rescue.filtered.stripped.vep.vcf.gz"
        paths = list((bundle / "original_workflow_rescue").glob(pattern))
        if len(paths) != 1:
            raise ValueError(f"Ambiguous rescue source: {paths}")
        candidates = {}
        for row, alt, _, _ in fields(track(paths[0])):
            key = row[0], int(row[1]), row[3], alt
            if len(row[3]) != 1 or len(alt) != 1 or row[6] not in {"Somatic", "PASS", "."}:
                continue
            if key in candidates:
                raise ValueError(f"Duplicate rescue candidate: {key}")
            candidates[key] = {"filter": row[6], "info": info_map(row[7])}
        dna, rna, eligible_rna = {}, {}, {}
        for caller in CALLERS:
            dna[caller] = caller_evidence(track(bundle / f"dna_callers/{caller}.vcf.gz"), caller, candidates, True)
            rna[caller], eligible_rna[caller] = caller_evidence(track(bundle / f"rna_realign_callers/{caller}.vcf.gz"), caller, candidates, False)
        recorded, recomputed, raw_with_recorded_dna = set(), set(), set()
        eligible_recomputed = set()
        details, reasons = [], Counter()
        for key, row in candidates.items():
            info = row["info"]
            nominated = {c for c in CALLERS if key in dna[c]}
            passed = {c for c in CALLERS if key in rna[c]}
            eligible = {c for c in CALLERS if key in eligible_rna[c]}
            old_dna = integer(info.get("N_DNA_CALLERS_SUPPORT"))
            old_rna = integer(info.get("N_RNA_CALLERS_SOMATIC"))
            old_gate = (old_dna is not None and old_dna >= 1 and old_rna is not None
                        and old_rna >= 2 and nominated and not biological_veto(info))
            if old_gate:
                recorded.add(key)
            allowed, reason = decide(key[2], key[3], row["filter"], nominated, passed, info)
            eligible_allowed, eligible_reason = decide(key[2], key[3], row["filter"], nominated, eligible, info)
            if eligible_allowed:
                eligible_recomputed.add(key)
            reasons[reason] += 1
            if allowed:
                recomputed.add(key)
                if old_dna is not None and old_dna >= 1:
                    raw_with_recorded_dna.add(key)
            if bool(old_gate) != allowed or key in expected:
                details.append({"allele": key, "recorded_dna_support": old_dna,
                                "recorded_rna_somatic": old_rna, "dna_nominators": sorted(nominated),
                                "raw_rna_pass_callers": sorted(passed), "decision": reason,
                                "eligible_rna_callers": sorted(eligible), "eligible_decision": eligible_reason,
                                "historical_expected": key in expected, "recomputed": allowed})
        # Retain all accepted candidates for transfer onto the NEW baseline.
        # Excluding the historical baseline too early can hide re-admissions
        # of SNPs deliberately rejected by the refined consensus.
        all_accepted = sorted(recomputed)
        all_eligible = sorted(eligible_recomputed)
        eligible_recomputed -= native
        recorded -= native
        recomputed -= native
        raw_with_recorded_dna -= native
        cell = {"candidate_count": len(candidates), "expected_additions": len(expected),
                "recomputed_all_accepted": all_accepted,
                "eligible_all_accepted": all_eligible,
                "eligible_count": len(eligible_recomputed),
                "eligible_extra": sorted(eligible_recomputed - expected),
                "eligible_missing": sorted(expected - eligible_recomputed),
                "recorded_gate_count": len(recorded), "recorded_extra": sorted(recorded - expected),
                "recorded_missing": sorted(expected - recorded), "recomputed_count": len(recomputed),
                "recomputed_extra": sorted(recomputed - expected), "recomputed_missing": sorted(expected - recomputed),
                "raw_with_recorded_dna_count": len(raw_with_recorded_dna),
                "reasons": dict(reasons), "details": details}
        report["datasets"][dataset] = cell
        (args.outdir / "audit.json").write_text(json.dumps(report, indent=2) + "\n")
        print(dataset, {k: v for k, v in cell.items() if k.endswith("count") or k == "expected_additions"}, flush=True)
    report["sources_unchanged"] = all(digest(Path(p)) == sha for p, sha in report["sources"].items())
    report["code_unchanged"] = all(digest(Path(p)) == sha for p, sha in report["code"].items())
    report["status"] = "complete_not_promoted" if report["sources_unchanged"] and report["code_unchanged"] else "integrity_failure"
    (args.outdir / "audit.json").write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
