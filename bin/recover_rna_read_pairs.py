#!/usr/bin/env python3
"""Recover selected reads from one original RNA library's paired FASTQs.

The caller must supply read IDs and FASTQs scoped to the same library. Paired
outputs contain synchronized mates only; orphan reads have separate outputs.
"""
import argparse
import gzip
import json
from pathlib import Path


def read_id(name):
    """Remove only an explicit FASTQ mate suffix, preserving the template ID."""
    name = name.split()[0]
    return name[:-2] if name.endswith(("/1", "/2")) else name


def read_fastq(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            sequence = handle.readline()
            separator = handle.readline()
            quality = handle.readline()
            if not (header.startswith("@") and header[1:].strip()
                    and sequence and separator.startswith("+") and quality):
                raise ValueError(
                    f"Malformed or truncated FASTQ record in {path}: {header.strip()!r}"
                )
            if len(sequence.rstrip("\r\n")) != len(quality.rstrip("\r\n")):
                raise ValueError(
                    f"FASTQ sequence/quality length mismatch in {path}: {header.strip()!r}"
                )
            yield read_id(header[1:]), header + sequence + separator + quality


def selected_reads(path, requested_ids, library):
    reads = {}
    for name, record in read_fastq(path):
        if name not in requested_ids:
            continue
        if name in reads:
            raise ValueError(
                f"Duplicate selected read name {name!r} in library {library!r}: {path}"
            )
        reads[name] = record
    return reads


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--read-ids", required=True, help="Library-scoped template IDs")
    parser.add_argument("--fastq-1", required=True)
    parser.add_argument("--fastq-2", required=True)
    parser.add_argument("--library", required=True)
    parser.add_argument("--out-prefix", required=True)
    parser.add_argument("--stats", required=True)
    args = parser.parse_args()

    try:
        requested_ids = {
            line.split()[0]
            for line in Path(args.read_ids).read_text().splitlines() if line.strip()
        }
        mate1 = selected_reads(args.fastq_1, requested_ids, args.library)
        mate2 = selected_reads(args.fastq_2, requested_ids, args.library)
    except (OSError, ValueError, EOFError) as error:
        parser.error(str(error))

    pairs = sorted(mate1.keys() & mate2.keys())
    singletons1 = sorted(mate1.keys() - mate2.keys())
    singletons2 = sorted(mate2.keys() - mate1.keys())
    for suffix, names, records in (
        ("_R1.fastq", pairs, mate1),
        ("_R2.fastq", pairs, mate2),
        ("_singleton_R1.fastq", singletons1, mate1),
        ("_singleton_R2.fastq", singletons2, mate2),
    ):
        with open(args.out_prefix + suffix, "w") as output:
            for name in names:
                output.write(records[name])
    stats = {
        "library": args.library,
        "requested": len(requested_ids),
        "recovered_r1": len(mate1),
        "recovered_r2": len(mate2),
        "paired": len(pairs),
        "singleton_r1": len(singletons1),
        "singleton_r2": len(singletons2),
        "missing": len(requested_ids - mate1.keys() - mate2.keys()),
    }
    Path(args.stats).write_text(json.dumps(stats, indent=2) + "\n")


if __name__ == "__main__":
    main()
