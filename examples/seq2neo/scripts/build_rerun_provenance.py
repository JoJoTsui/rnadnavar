#!/usr/bin/env python3
"""
build_rerun_provenance.py
─────────────────────────
Build the committed provenance table for the consensus+rescue rerun:
which EXACT input files (path + md5 + size) each sample's rerun consumed.

Source: runs/rerun_checksums/<sample>.input_checksums.json — recorded by
run_reconsensus_rerun.py immediately before each sample's Nextflow run and
verified unchanged after it. These checksum manifests are gitignored runtime
artifacts; this script distills them into a small committable TSV.

Output (default: data/processed/rerun_input_provenance.tsv), long format:
  sample_id  input  path  md5  size_bytes  recorded
where input ∈ {dna,rna}_{mutect2,strelka,deepsomatic} (.tbi/.csi excluded —
they are derived files; the VCF checksum is the identity proof).

Usage:
  python3 scripts/build_rerun_provenance.py            # write the TSV
  python3 scripts/build_rerun_provenance.py --check    # verify vs existing TSV
"""

import argparse
import csv
import json
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
SEQ2NEO_ROOT = SCRIPT_DIR.parent
CHECKSUM_DIR = SEQ2NEO_ROOT / "runs" / "rerun_checksums"
DEFAULT_OUT = SEQ2NEO_ROOT / "data" / "processed" / "rerun_input_provenance.tsv"

CALLER_SUFFIX = {
    ".mutect2.filtered.vcf.gz": "mutect2",
    ".strelka.variants.vcf.gz": "strelka",
    ".deepsomatic.vcf.gz": "deepsomatic",
}

HEADER = ["sample_id", "input", "path", "md5", "size_bytes", "recorded"]


def classify_input(path_str: str):
    """Map an input VCF path to its '<modality>_<caller>' input name."""
    if path_str.endswith((".tbi", ".csi")):
        return None
    caller = next((c for sfx, c in CALLER_SUFFIX.items() if path_str.endswith(sfx)), None)
    if caller is None:
        return None
    modality = "rna" if "/vcf_realignment/" in path_str else "dna"
    return f"{modality}_{caller}"


def build_rows():
    rows = []
    manifests = sorted(CHECKSUM_DIR.glob("*.input_checksums.json"))
    if not manifests:
        sys.exit(f"no checksum manifests found in {CHECKSUM_DIR}")
    for mf in manifests:
        sid = mf.name.replace(".input_checksums.json", "")
        data = json.loads(mf.read_text())
        recorded = data.get("recorded", "")
        for path_str, entry in sorted(data["files"].items()):
            name = classify_input(path_str)
            if name is None:
                continue
            rows.append({
                "sample_id": sid,
                "input": name,
                "path": path_str,
                "md5": entry["md5"],
                "size_bytes": entry["size"],
                "recorded": recorded,
            })
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=str(DEFAULT_OUT))
    ap.add_argument("--check", action="store_true",
                    help="Do not write; fail if the existing TSV differs")
    args = ap.parse_args()

    rows = build_rows()

    # Sanity: every sample must contribute exactly the 6 caller VCFs.
    per_sample = {}
    for r in rows:
        per_sample.setdefault(r["sample_id"], set()).add(r["input"])
    expected = {f"{m}_{c}" for m in ("dna", "rna")
                for c in ("mutect2", "strelka", "deepsomatic")}
    problems = {s: sorted(expected - got) for s, got in per_sample.items()
                if got != expected}
    if problems:
        for s, missing in sorted(problems.items()):
            print(f"WARN: {s}: incomplete inputs: missing {missing}",
                  file=sys.stderr)

    out = Path(args.out)
    text = "\t".join(HEADER) + "\n" + "".join(
        "\t".join(str(r[h]) for h in HEADER) + "\n" for r in rows)

    if args.check:
        if not out.exists():
            sys.exit(f"missing: {out}")
        if out.read_text() != text:
            sys.exit(f"STALE: {out} differs from runs/rerun_checksums/ — "
                     f"re-run without --check to regenerate")
        print(f"OK: {out} matches runs/rerun_checksums/ "
              f"({len(per_sample)} samples, {len(rows)} inputs)")
        return

    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(text)
    print(f"wrote {out}  ({len(per_sample)} samples, {len(rows)} input rows)")


if __name__ == "__main__":
    main()
