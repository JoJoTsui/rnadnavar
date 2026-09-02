#!/usr/bin/env python3
"""
build_rerun_manifest.py
───────────────────────
Regenerate data/processed/sample_manifest_rerun.tsv from the first-run
sample manifest (data/processed/sample_manifest.tsv).

What this manifest is: per-sample pointers to the re-consensus RERUN label
VCFs (consensus+rescue re-run from the FIRST run's caller VCFs — the bam_*
and caller_* columns deliberately still reference first-run outputs, which
are the rerun's actual inputs).

Transformations vs the source manifest:
  * rescue_vcf_path → rerun label VCF (rescue filtered stripped), prefixed
    with the rsynced tree root (RSYNC_PREFIX). Local existence is verified
    against LOCAL_PREFIX; the build fails if a label VCF is missing locally.
  * status/status_reason → 'useless' + reason for EXCLUDED_SAMPLES (kept in
    the table, marked not usable for training).
  * appended columns:
      label_qc_verdict     — cohort label_qc gate verdict (PASS/WARN/FAIL),
                             read from runs/label_qc/cohort66_apply/samples_qc.tsv
      label_qc_cleaned_vcf — for WARN samples, the label_qc --apply cleaned
                             VCF (bgzipped) under the rsynced tree; empty else
      known_limitations    — cohort-wide caveats of the rerun labels
                             (chrM records dropped by the gnomAD scatter-gather
                             merge, fixed in d336c7f AFTER the cohort; Rule-2
                             common-AF/prior-artifact vetoes also post-date the
                             cohort — measured zero Rule-2 firings cohort-wide)

Usage:
  python3 scripts/build_rerun_manifest.py            # write the TSV
  python3 scripts/build_rerun_manifest.py --check    # verify vs existing TSV
"""

import argparse
import csv
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
SEQ2NEO_ROOT = SCRIPT_DIR.parent
SOURCE_MANIFEST = SEQ2NEO_ROOT / "data" / "processed" / "sample_manifest.tsv"
DEFAULT_OUT = SEQ2NEO_ROOT / "data" / "processed" / "sample_manifest_rerun.tsv"
LABEL_QC_TSV = SEQ2NEO_ROOT / "runs" / "label_qc" / "cohort66_apply" / "samples_qc.tsv"

LOCAL_PREFIX = SEQ2NEO_ROOT  # where outputs live on this machine
RSYNC_PREFIX = Path(
    "/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/rnadnavar/examples/seq2neo"
)

# Samples excluded from training by the rerun-vs-clean anomaly analysis
# (docs/RERUN_LABEL_ANOMALY_ANALYSIS.md). Kept in the manifest, marked useless.
EXCLUDED_SAMPLES = {
    "PRJNA298330_4032": (
        "excluded: TruthQC suspect sample; rerun label VCF is raw (un-QC'd) and "
        "resurrects 4288 QC-dropped RNA-only common-AF sites "
        "(docs/RERUN_LABEL_ANOMALY_ANALYSIS.md)"
    ),
    "PRJNA298376_4081": (
        "excluded: TruthQC suspect sample (RefCall-heavy inputs, 93% Somatic "
        "drop rate); clean baseline predates 2026-03-18 germline rule "
        "(docs/RERUN_LABEL_ANOMALY_ANALYSIS.md)"
    ),
    "PRJNA298376_4255": (
        "excluded: TruthQC suspect sample (RefCall-heavy inputs, 89% Somatic "
        "drop rate); clean baseline predates 2026-03-18 germline rule "
        "(docs/RERUN_LABEL_ANOMALY_ANALYSIS.md)"
    ),
}

KNOWN_LIMITATIONS = "nuclear_only_chrM_dropped;pre_rule2_veto_zero_firings"

ADDED_COLUMNS = ["label_qc_verdict", "label_qc_cleaned_vcf", "known_limitations"]


def rescue_label_relpath(row: dict) -> Path:
    """Rerun label VCF path relative to the output root."""
    p = row["vcf_prefix"]
    pair = f"{p}DT_vs_{p}DN_rescued_{p}RT_realign_vs_{p}DN"
    return (Path("output_reconsensus") / row["sample_id"] / "rescue" / pair
            / f"{pair}.filtered.vcf.stripped.vcf.gz")


def load_label_qc():
    """sample -> (verdict, cleaned_vcf_relpath_or_empty)."""
    qc = {}
    if not LABEL_QC_TSV.exists():
        print(f"WARN: label_qc verdicts not found at {LABEL_QC_TSV}; "
              f"label_qc_verdict will be empty", file=sys.stderr)
        return qc
    with open(LABEL_QC_TSV, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            sid = row["sample"]
            cleaned = ""
            if row["verdict"] == "WARN":
                rel = (Path("runs/label_qc/cohort66_apply/cleaned_vcf")
                       / f"{sid}.bgz.vcf.gz")
                if not (LOCAL_PREFIX / rel).exists():
                    print(f"WARN: cleaned VCF missing locally: {rel}",
                          file=sys.stderr)
                else:
                    cleaned = str(rel)
            qc[sid] = (row["verdict"], cleaned)
    return qc


def build_rows():
    with open(SOURCE_MANIFEST, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    qc = load_label_qc()
    out_rows = []
    for row in rows:
        sid = row["sample_id"]
        rel = rescue_label_relpath(row)
        if not (LOCAL_PREFIX / rel).exists():
            sys.exit(f"ERROR: rerun label VCF missing locally: {LOCAL_PREFIX / rel}")
        row["rescue_vcf_path"] = str(RSYNC_PREFIX / rel)
        if sid in EXCLUDED_SAMPLES:
            row["status"] = "useless"
            row["status_reason"] = EXCLUDED_SAMPLES[sid]
        verdict, cleaned = qc.get(sid, ("", ""))
        row["label_qc_verdict"] = verdict
        row["label_qc_cleaned_vcf"] = (
            str(RSYNC_PREFIX / cleaned) if cleaned else "")
        row["known_limitations"] = KNOWN_LIMITATIONS
        out_rows.append(row)
    return out_rows


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=str(DEFAULT_OUT))
    ap.add_argument("--check", action="store_true",
                    help="Do not write; fail if the existing TSV differs")
    args = ap.parse_args()

    rows = build_rows()
    header = list(rows[0].keys())  # source columns, then ADDED_COLUMNS
    assert header[-len(ADDED_COLUMNS):] == ADDED_COLUMNS

    text = "\t".join(header) + "\n" + "".join(
        "\t".join(row.get(h, "") for h in header) + "\n" for row in rows)

    out = Path(args.out)
    if args.check:
        if not out.exists():
            sys.exit(f"missing: {out}")
        if out.read_text() != text:
            sys.exit(f"STALE: {out} differs from regenerated content — "
                     f"re-run without --check to rewrite")
        print(f"OK: {out} is current ({len(rows)} samples)")
        return

    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(text)
    n_useless = sum(1 for r in rows if r["status"] == "useless")
    n_warn = sum(1 for r in rows if r["label_qc_verdict"] == "WARN")
    print(f"wrote {out}  ({len(rows)} samples; {n_useless} useless, "
          f"{n_warn} WARN with cleaned VCF)")


if __name__ == "__main__":
    main()
