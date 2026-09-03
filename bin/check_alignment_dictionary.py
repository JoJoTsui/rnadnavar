#!/usr/bin/env python3
"""Audit an alignment sequence dictionary against the selected reference.

The checker deliberately does not rewrite alignments.  It establishes whether
an external BAM can be normalized safely with Picard ReorderSam, or whether an
external BAM/CRAM must be rejected because a shared contig has a conflicting
length (an assembly mismatch).  A JSON audit is always written, including for
rejected inputs.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--header", required=True, type=Path)
    parser.add_argument("--idxstats", required=True, type=Path)
    parser.add_argument("--reference-dict", required=True, type=Path)
    parser.add_argument("--alignment", required=True, type=Path)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--stage", choices=("input", "normalized"), required=True)
    parser.add_argument("--policy", choices=("normalize", "strict"), required=True)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_sq(path: Path) -> tuple[list[tuple[str, int]], list[str]]:
    sequence: list[tuple[str, int]] = []
    malformed: list[str] = []
    seen: set[str] = set()
    for line_number, raw in enumerate(path.read_text().splitlines(), start=1):
        if not raw.startswith("@SQ\t"):
            continue
        fields = {}
        for field in raw.split("\t")[1:]:
            if ":" in field:
                key, value = field.split(":", 1)
                fields[key] = value
        name = fields.get("SN")
        length = fields.get("LN")
        if not name or not length:
            malformed.append(f"line {line_number}: missing SN or LN")
            continue
        try:
            parsed_length = int(length)
        except ValueError:
            malformed.append(f"line {line_number}: invalid LN={length!r}")
            continue
        if name in seen:
            malformed.append(f"line {line_number}: duplicate SN={name!r}")
            continue
        seen.add(name)
        sequence.append((name, parsed_length))
    if not sequence:
        malformed.append("no @SQ records found")
    return sequence, malformed


def parse_idxstats(path: Path) -> tuple[dict[str, dict[str, int]], list[str]]:
    stats: dict[str, dict[str, int]] = {}
    malformed: list[str] = []
    for line_number, raw in enumerate(path.read_text().splitlines(), start=1):
        fields = raw.split("\t")
        if len(fields) != 4:
            malformed.append(f"line {line_number}: expected four columns")
            continue
        name, length, mapped, unmapped = fields
        try:
            stats[name] = {
                "length": int(length),
                "mapped": int(mapped),
                "unmapped": int(unmapped),
            }
        except ValueError:
            malformed.append(f"line {line_number}: non-integer count or length")
    return stats, malformed


def build_audit(args: argparse.Namespace) -> tuple[dict, bool, str | None]:
    alignment_sq, alignment_errors = parse_sq(args.header)
    reference_sq, reference_errors = parse_sq(args.reference_dict)
    idxstats, idxstats_errors = parse_idxstats(args.idxstats)

    alignment_map = dict(alignment_sq)
    reference_map = dict(reference_sq)
    extra_names = [name for name, _ in alignment_sq if name not in reference_map]
    missing_names = [name for name, _ in reference_sq if name not in alignment_map]
    shared_length_mismatches = [
        {
            "contig": name,
            "alignment_length": alignment_map[name],
            "reference_length": reference_map[name],
        }
        for name, _ in alignment_sq
        if name in reference_map and alignment_map[name] != reference_map[name]
    ]
    shared_alignment_order = [name for name, _ in alignment_sq if name in reference_map]
    shared_reference_order = [name for name, _ in reference_sq if name in alignment_map]
    order_mismatch = shared_alignment_order != shared_reference_order

    extra_contigs = []
    for name in extra_names:
        counts = idxstats.get(name, {"mapped": 0, "unmapped": 0})
        extra_contigs.append(
            {
                "name": name,
                "length": alignment_map[name],
                "mapped_reads": counts["mapped"],
                "unmapped_reads": counts["unmapped"],
            }
        )

    malformed = alignment_errors + reference_errors + idxstats_errors
    exact = (
        not malformed
        and not extra_names
        and not missing_names
        and not shared_length_mismatches
        and not order_mismatch
        and alignment_sq == reference_sq
    )
    unsafe_reason = None
    if malformed:
        unsafe_reason = "malformed sequence dictionary or idxstats"
    elif shared_length_mismatches:
        unsafe_reason = "shared contig lengths differ"

    normalization_required = not exact and unsafe_reason is None
    if unsafe_reason:
        decision = "rejected"
    elif normalization_required and args.policy == "strict":
        decision = "rejected"
        unsafe_reason = "dictionary differs and policy is strict"
    elif normalization_required:
        decision = "normalize"
    else:
        decision = "compatible"

    total_mapped = sum(v["mapped"] for name, v in idxstats.items() if name != "*")
    total_unmapped = sum(v["unmapped"] for v in idxstats.values())
    audit = {
        "schema_version": 1,
        "sample_id": args.sample_id,
        "stage": args.stage,
        "policy": args.policy,
        "decision": decision,
        "reason": unsafe_reason,
        "alignment": {
            "path": args.alignment.name,
            "size_bytes": args.alignment.stat().st_size,
            "sha256": sha256(args.alignment),
            "contig_count": len(alignment_sq),
            "mapped_reads": total_mapped,
            "unmapped_reads": total_unmapped,
        },
        "reference": {
            "path": args.reference_dict.name,
            "size_bytes": args.reference_dict.stat().st_size,
            "sha256": sha256(args.reference_dict),
            "contig_count": len(reference_sq),
        },
        "differences": {
            "extra_contigs": extra_contigs,
            "missing_reference_contigs": [
                {"name": name, "length": reference_map[name]} for name in missing_names
            ],
            "shared_length_mismatches": shared_length_mismatches,
            "shared_order_mismatch": order_mismatch,
            "malformed_records": malformed,
        },
        "normalization": {
            "required": normalization_required,
            "candidate_reads_affected": sum(
                c["mapped_reads"] + c["unmapped_reads"] for c in extra_contigs
            ),
            "all_reference_contigs_retained": not missing_names,
        },
    }
    return audit, normalization_required, unsafe_reason


def main() -> int:
    args = parse_args()
    audit, normalization_required, reason = build_audit(args)
    args.output.write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n")
    if audit["decision"] == "rejected":
        print(
            f"ERROR: alignment dictionary rejected for {args.sample_id}: {reason}; "
            f"see {args.output}",
            file=sys.stderr,
        )
        return 2
    print("true" if normalization_required else "false")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
