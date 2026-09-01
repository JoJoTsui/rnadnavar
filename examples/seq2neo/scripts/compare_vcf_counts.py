#!/usr/bin/env python3
"""compare_vcf_counts.py — fast per-sample variant-count comparison between two VCF sets.

Counts records (non-header lines) in gzipped or plain VCFs, matched across two
directories by a sample ID extracted with per-side regexes. Parallelized with
multiprocessing; counting uses bulk newline counting (headers are contiguous at
the top of a VCF, so only the header block is scanned line-by-line).

Stdlib-only.

Example (re-consensus rerun vs QC-cleaned labels):

    python3 compare_vcf_counts.py \
      --left-dir  /path/to/cleaned_vcf --left-glob '*.vcf.gz' \
      --left-id-regex '^(?:(PRJNA[0-9]+)_)?([0-9]+?)(?:_rnadnavar)?\\.vcf\\.gz$' \
      --right-dir output_reconsensus \
      --right-glob '*/rescue/*_rescued_*/*.filtered.vcf.stripped.vcf.gz' \
      --right-id-regex '(PRJNA[0-9]+)_([0-9]+)' \
      --default-project PRJNA298376 \
      --jobs 16 --out rerun_vs_clean_counts.tsv
"""

import argparse
import gzip
import re
import sys
from multiprocessing import Pool
from pathlib import Path

CHUNK = 8 * 1024 * 1024  # 8 MiB read chunks


def count_vcf_records(path: str) -> int:
    """Count non-header lines in a (possibly gzipped) VCF, fast path.

    Scans the header line-by-line until the first data line, then counts
    newlines in the remaining chunks in bulk.
    """
    opener = gzip.open if path.endswith(".gz") else open
    count = 0
    last_byte = b"\n"
    with opener(path, "rb") as fh:
        carry = b""
        in_header = True
        while True:
            chunk = fh.read(CHUNK)
            if not chunk:
                break
            last_byte = chunk[-1:]
            if in_header:
                data = carry + chunk
                lines = data.split(b"\n")
                # last element may be a partial line; keep it in carry
                carry = lines.pop()
                for i, line in enumerate(lines):
                    if not line.startswith(b"#"):
                        # first data line found: count it + the rest in bulk
                        in_header = False
                        count += 1
                        count += len(lines) - 1 - i
                        break
            else:
                count += chunk.count(b"\n")
        if not in_header and last_byte != b"\n":
            count += 1  # final line without trailing newline
    return count


def count_vcf_records_by_filter(path: str) -> dict:
    """Count records grouped by FILTER value (column 7). Slower, full parse."""
    opener = gzip.open if path.endswith(".gz") else open
    counts: dict[str, int] = {}
    with opener(path, "rb") as fh:
        for line in fh:
            if line.startswith(b"#"):
                continue
            filt = line.split(b"\t", 7)[6].decode()
            counts[filt] = counts.get(filt, 0) + 1
    counts["__total__"] = sum(counts.values())
    return counts


def extract_id(path: Path, pattern: re.Pattern, root: Path, default_project: str) -> str | None:
    """Extract sample key. The regex may carry one group (sample id) or two
    (project, sample id); a missing project group falls back to
    default_project. The regex is searched against the path relative to the
    scan root, so it can also match directory components."""
    m = pattern.search(str(path.relative_to(root)))
    if not m:
        return None
    groups = m.groups()
    if len(groups) >= 2 and groups[-2]:
        return f"{groups[-2]}_{groups[-1]}"
    return f"{default_project}_{groups[-1]}" if default_project else groups[-1]


def collect(directory: str, glob_pat: str, id_regex: str, default_project: str) -> dict[str, str]:
    pat = re.compile(id_regex)
    root = Path(directory)
    out = {}
    for p in sorted(root.glob(glob_pat)):
        sid = extract_id(p, pat, root, default_project)
        if sid is None:
            print(f"WARNING: no ID match: {p}", file=sys.stderr)
            continue
        if sid in out:
            print(f"WARNING: duplicate ID {sid}: {p} (kept {out[sid]})", file=sys.stderr)
            continue
        out[sid] = str(p)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--left-dir", required=True)
    ap.add_argument("--left-glob", required=True)
    ap.add_argument("--left-id-regex", required=True, help="regex with one capture group = sample ID")
    ap.add_argument("--right-dir", required=True)
    ap.add_argument("--right-glob", required=True)
    ap.add_argument("--right-id-regex", required=True)
    ap.add_argument("--default-project", default="", help="project prefix used when the regex has no project group")
    ap.add_argument("--left-name", default="left")
    ap.add_argument("--right-name", default="right")
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--by-filter", action="store_true", help="also tally counts per FILTER value (slower)")
    ap.add_argument("--out", default=None, help="output TSV (default: stdout)")
    args = ap.parse_args()

    left = collect(args.left_dir, args.left_glob, args.left_id_regex, args.default_project)
    right = collect(args.right_dir, args.right_glob, args.right_id_regex, args.default_project)
    samples = sorted(set(left) | set(right))

    paths = sorted({*left.values(), *right.values()})
    counter = count_vcf_records_by_filter if args.by_filter else count_vcf_records
    with Pool(args.jobs) as pool:
        results = dict(zip(paths, pool.map(counter, paths)))

    rows = ["sample\t{}_count\t{}_count\tdiff({}-{})\tpct_diff".format(
        args.left_name, args.right_name, args.right_name, args.left_name)]
    n_eq = n_diff = 0
    for sid in samples:
        lp, rp = left.get(sid), right.get(sid)
        if lp is None:
            rows.append(f"{sid}\tMISSING_{args.left_name.upper()}\t{results[rp]['__total__'] if args.by_filter else results[rp]}\t-\t-")
            continue
        if rp is None:
            rows.append(f"{sid}\t{results[lp]['__total__'] if args.by_filter else results[lp]}\tMISSING_{args.right_name.upper()}\t-\t-")
            continue
        lc, rc = results[lp], results[rp]
        if args.by_filter:
            lc, rc = lc["__total__"], rc["__total__"]
        diff = rc - lc
        pct = (100.0 * diff / lc) if lc else 0.0
        rows.append(f"{sid}\t{lc}\t{rc}\t{diff}\t{pct:+.4f}%")
        if diff == 0:
            n_eq += 1
        else:
            n_diff += 1

    if args.by_filter:
        rows.append("")
        rows.append("# per-FILTER breakdown (sample, side, FILTER, count)")
        for sid in samples:
            for side, m in ((args.left_name, left), (args.right_name, right)):
                p = m.get(sid)
                if p is None:
                    continue
                for filt, n in sorted(results[p].items()):
                    if filt != "__total__":
                        rows.append(f"{sid}\t{side}\t{filt}\t{n}")

    text = "\n".join(rows) + "\n"
    if args.out:
        Path(args.out).write_text(text)
    else:
        sys.stdout.write(text)

    paired = n_eq + n_diff
    print(f"compared {paired} paired samples: {n_eq} identical, {n_diff} differ; "
          f"{len(samples) - paired} unpaired (see MISSING rows)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
