#!/usr/bin/env python3
"""
build_split_manifest.py
───────────────────────
Build the seq2neo split manifest from the reconsensus rerun truth-label VCFs.

Reads data/processed/sample_manifest_rerun.tsv, streams each working-cohort
sample's training_label_vcf (gzipped, 8-column VCF, no sample column) at its
LITERAL path (the /t9k/mnt/WorkSpace/... rsync tree is mounted on this
machine — no prefix remapping, unlike build_rerun_manifest.py), selects the
whitelist variants, assigns each a split, computes sub-pool tags from VCF
INFO, and writes three artifacts into --outdir:

  sample_split.tsv            sample-level split assignment (63 rows)
  selected_variants.parquet   variant-level split + tags (zstd)
  split_manifest.provenance.txt  generation provenance + totals + caveats

No tensors are created here; that stays downstream in neo_var.

Working cohort: manifest rows with non-empty training_label_vcf (63 samples:
56 PASS + 7 WARN). The build FAILS if any such VCF is missing/unreadable, if
any sample_id is duplicated, or if the reserved-sample validation fails
(every RESERVED_SAMPLE_IDS id must be in the cohort with verdict PASS and the
expected disease_normalized).

Selection: FILTER in {Somatic, Germline, Reference}; all other records
(Artifact, NoConsensus, RNAedit) are dropped. This absorbs the whitelist
that neo_var previously applied at tensor-extraction time into split-build
time, so downstream extraction never sees a non-whitelist row.

Split routing per record:
  * sample in RESERVED_SAMPLE_IDS            → "reserved" (all chromosomes)
  * else sample label_qc_verdict == "WARN"   → "train"    (all chromosomes;
    cleaned-label noise does the least harm in train — WARN samples never
    enter val/test/reserved evaluation pools)
  * else deepsomatic chromosome map: chr1 → test, chr21–22 → val,
    chr2–20 → train, anything else (chrX/chrY/chrM/unrecognized) → train.

Port notes (explicit):
  * DEEPSOMATIC_SPLIT and _normalize_chrom are ported verbatim from
    neo_var/src/neo_var/data/split_dataset.py (accepts "1", "chr1", "Chr1",
    "CHR1", "chromosome1", X/Y/M aliases; chromosomes not in the map → train).
  * RESERVED_SAMPLE_IDS, _EXPECTED_RESERVED_DISEASES and the sub-pool tag
    predicates are ported verbatim in semantics from
    neo_var/src/neo_var/data/reserve_downstream.py (compute_sub_pool_tags),
    with one substitution: neo_var's is_low_dp uses BAM_DT_DP (BAM pileup
    depth); here DP_DNA_MEAN (caller-mean DNA depth from the label VCF INFO)
    stands in for it, including the null semantics (missing DP_DNA_MEAN →
    is_low_dp is False, not True).
  * Sub-pool tags are computed for EVERY selected row (not just reserved);
    is_zero_shot is True iff split == "reserved".

Usage:
  python3 scripts/build_split_manifest.py                 # write artifacts
  python3 scripts/build_split_manifest.py --workers 16    # parallel parse
  python3 scripts/build_split_manifest.py --check         # verify vs existing
"""

import argparse
import csv
import gzip
import hashlib
import os
import subprocess
import sys
import time
from datetime import datetime, timezone
from multiprocessing import Pool
from pathlib import Path

import polars as pl

SCRIPT_DIR = Path(__file__).resolve().parent
SEQ2NEO_ROOT = SCRIPT_DIR.parent
DEFAULT_MANIFEST = SEQ2NEO_ROOT / "data" / "processed" / "sample_manifest_rerun.tsv"
DEFAULT_OUTDIR = SEQ2NEO_ROOT / "data" / "processed"
FASTQ_MD5SUMS = SEQ2NEO_ROOT / "data" / "raw" / "md5sums.txt"

SAMPLE_TSV_NAME = "sample_split.tsv"
PARQUET_NAME = "selected_variants.parquet"
PROVENANCE_NAME = "split_manifest.provenance.txt"

# Extraction-time FILTER whitelist, absorbed into split-build time.
WHITELIST_FILTERS = frozenset({"Somatic", "Germline", "Reference"})

KNOWN_LIMITATIONS = "nuclear_only_chrM_dropped;pre_rule2_veto_zero_firings"

# ── Ported verbatim from neo_var/src/neo_var/data/reserve_downstream.py ─────
RESERVED_SAMPLE_IDS = {
    "PRJNA298310_3812",  # melanoma (skin cancer)
    "PRJNA298376_4166",  # tumor dna sample from lung of a human female participant
    "PRJNA298376_4214",  # lung cancer
    "PRJNA298376_4231",  # ampullary cancer (rare, non-colorectal)
    "PRJNA298376_4242",  # gastric cancer (stomach, non-colorectal)
}

_EXPECTED_RESERVED_DISEASES = {
    "PRJNA298310_3812": "melanoma",
    "PRJNA298376_4166": "tumor dna sample from lung of a human female participant",
    "PRJNA298376_4214": "lung cancer",
    "PRJNA298376_4231": "ampullary cancer",
    "PRJNA298376_4242": "gastric cancer",
}

# ── Ported verbatim from neo_var/src/neo_var/data/split_dataset.py ──────────
# DeepSomatic-style: train on chr2-20, val on chr21-22, test on chr1
DEEPSOMATIC_SPLIT = {
    "train": [f"chr{i}" for i in range(2, 21)],  # chr2–chr20
    "val": ["chr21", "chr22"],
    "test": ["chr1"],
}

_ALT_NAMES: dict = {}
for _base, _names in [
    ("chr1", ["1", "chr1", "Chr1", "CHR1", "chromosome1"]),
    ("chr2", ["2", "chr2", "Chr2", "CHR2", "chromosome2"]),
    ("chr3", ["3", "chr3", "Chr3", "CHR3", "chromosome3"]),
    ("chr4", ["4", "chr4", "Chr4", "CHR4", "chromosome4"]),
    ("chr5", ["5", "chr5", "Chr5", "CHR5", "chromosome5"]),
    ("chr6", ["6", "chr6", "Chr6", "CHR6", "chromosome6"]),
    ("chr7", ["7", "chr7", "Chr7", "CHR7", "chromosome7"]),
    ("chr8", ["8", "chr8", "Chr8", "CHR8", "chromosome8"]),
    ("chr9", ["9", "chr9", "Chr9", "CHR9", "chromosome9"]),
    ("chr10", ["10", "chr10", "Chr10", "CHR10", "chromosome10"]),
    ("chr11", ["11", "chr11", "Chr11", "CHR11", "chromosome11"]),
    ("chr12", ["12", "chr12", "Chr12", "CHR12", "chromosome12"]),
    ("chr13", ["13", "chr13", "Chr13", "CHR13", "chromosome13"]),
    ("chr14", ["14", "chr14", "Chr14", "CHR14", "chromosome14"]),
    ("chr15", ["15", "chr15", "Chr15", "CHR15", "chromosome15"]),
    ("chr16", ["16", "chr16", "Chr16", "CHR16", "chromosome16"]),
    ("chr17", ["17", "chr17", "Chr17", "CHR17", "chromosome17"]),
    ("chr18", ["18", "chr18", "Chr18", "CHR18", "chromosome18"]),
    ("chr19", ["19", "chr19", "Chr19", "CHR19", "chromosome19"]),
    ("chr20", ["20", "chr20", "Chr20", "CHR20", "chromosome20"]),
    ("chr21", ["21", "chr21", "Chr21", "CHR21", "chromosome21"]),
    ("chr22", ["22", "chr22", "Chr22", "CHR22", "chromosome22"]),
    ("chrX", ["X", "chrX", "ChrX", "CHRX", "chromosomeX", "23"]),
    ("chrY", ["Y", "chrY", "ChrY", "CHRY", "chromosomeY", "24"]),
    ("chrM", ["M", "MT", "chrM", "chrMT", "ChrM", "CHRM", "chromosomeM", "mitochondrial"]),
]:
    _ALT_NAMES.update({name: _base for name in _names})


def _normalize_chrom(name: str) -> str:
    """Normalise a chromosome name to the ``chrN`` canonical form."""
    if not name:
        return ""
    name = str(name).strip()
    if name in _ALT_NAMES:
        return _ALT_NAMES[name]
    if not name.lower().startswith("chr"):
        cand = f"chr{name}"
        if cand in _ALT_NAMES:
            return _ALT_NAMES[cand]
    return name


# Canonical chromosome → split; anything not in this map routes to train.
_CHROM_TO_SPLIT = {c: s for s, chroms in DEEPSOMATIC_SPLIT.items() for c in chroms}

# Natural-sort rank: chr1–22 → 1–22, then X, Y, M, then anything else.
_CHROM_RANK = {"X": 23, "Y": 24, "M": 25}


def _chrom_rank(normalized: str) -> int:
    if normalized.startswith("chr"):
        rest = normalized[3:]
        if rest.isdigit():
            n = int(rest)
            return n if 1 <= n <= 22 else 26
        return _CHROM_RANK.get(rest, 26)
    return 26


SAMPLE_TSV_COLUMNS = [
    "sample_id", "project_id", "patient_id", "set_number", "disease_normalized",
    "label_qc_verdict", "sample_pool", "training_label_vcf", "n_variants_selected",
]

PARQUET_COLUMNS = [
    "sample_id", "CHROM", "POS", "REF", "ALT", "FILTER", "split", "label_verdict",
    "variant_type", "is_zero_shot", "is_low_vaf_a", "is_low_vaf_b", "is_low_dp",
    "is_rescued", "is_non_rescued", "is_indel",
    "VAF_DNA_MEAN", "VAF_RNA_MEAN", "DP_DNA_MEAN", "DP_RNA_MEAN", "RESCUED",
]

TAG_COLUMNS = [
    "is_zero_shot", "is_low_vaf_a", "is_low_vaf_b", "is_low_dp",
    "is_rescued", "is_non_rescued", "is_indel",
]

SPLIT_ORDER = ["train", "val", "test", "reserved"]
FILTER_ORDER = ["Somatic", "Germline", "Reference"]

# Worker frame schema (pre-tag; tags are added vectorized in the parent).
_WORKER_SCHEMA = {
    "sample_id": pl.Utf8, "CHROM": pl.Utf8, "POS": pl.Int64,
    "REF": pl.Utf8, "ALT": pl.Utf8, "FILTER": pl.Utf8,
    "split": pl.Utf8, "label_verdict": pl.Utf8, "variant_type": pl.Utf8,
    "VAF_DNA_MEAN": pl.Float64, "VAF_RNA_MEAN": pl.Float64,
    "DP_DNA_MEAN": pl.Float64, "DP_RNA_MEAN": pl.Float64,
    "RESCUED": pl.Utf8, "_chrom_rank": pl.Int64,
}


def _fail(msg: str):
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(1)


def _md5(path: Path) -> str:
    return hashlib.md5(Path(path).read_bytes()).hexdigest()


def _git_head() -> str:
    try:
        out = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=SCRIPT_DIR,
            capture_output=True, text=True, check=True,
        )
        return out.stdout.strip()
    except Exception:
        return "unknown"


def _load_working_cohort(manifest_path: Path) -> list:
    """Manifest rows with non-empty training_label_vcf, with build-time validation."""
    if not manifest_path.exists():
        _fail(f"manifest not found: {manifest_path}")
    with open(manifest_path, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    cohort = [r for r in rows if (r.get("training_label_vcf") or "").strip()]
    if not cohort:
        _fail(f"no working-cohort rows (non-empty training_label_vcf) in {manifest_path}")

    seen, dups = set(), set()
    for r in cohort:
        sid = r["sample_id"]
        if sid in seen:
            dups.add(sid)
        seen.add(sid)
    if dups:
        _fail(f"duplicate sample_id in working cohort: {sorted(dups)}")

    missing = [r["training_label_vcf"] for r in cohort
               if not (os.path.isfile(r["training_label_vcf"])
                       and os.access(r["training_label_vcf"], os.R_OK))]
    if missing:
        _fail("training label VCFs missing/unreadable:\n  " + "\n  ".join(missing))

    by_id = {r["sample_id"]: r for r in cohort}
    problems = []
    for sid in sorted(RESERVED_SAMPLE_IDS):
        r = by_id.get(sid)
        if r is None:
            problems.append(f"{sid}: not present in working cohort")
            continue
        if r["label_qc_verdict"] != "PASS":
            problems.append(
                f"{sid}: label_qc_verdict={r['label_qc_verdict']!r}, expected 'PASS'")
        expected = _EXPECTED_RESERVED_DISEASES[sid]
        if r["disease_normalized"] != expected:
            problems.append(
                f"{sid}: disease_normalized={r['disease_normalized']!r}, "
                f"expected {expected!r}")
    if problems:
        _fail("reserved-sample validation failed:\n  " + "\n  ".join(problems))

    return cohort


def _variant_type(ref: str, alt: str) -> str:
    if len(ref) == 1 and len(alt) == 1:
        return "SNV"
    if len(alt) > len(ref):
        return "INS"
    if len(ref) > len(alt):
        return "DEL"
    return "OTHER"


def _process_sample(task: dict) -> dict:
    """Worker: parse one label VCF → polars DataFrame of selected rows."""
    sid = task["sample_id"]
    verdict = task["label_verdict"]
    is_reserved = sid in RESERVED_SAMPLE_IDS

    cols = {name: [] for name in _WORKER_SCHEMA}
    n_records = 0
    n_off_chrom = 0
    off_chroms = set()

    with gzip.open(task["training_label_vcf"], "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            n_records += 1
            parts = line.split("\t", 8)
            filt = parts[6]
            if filt not in WHITELIST_FILTERS:
                continue
            chrom, ref, alt = parts[0], parts[3], parts[4]
            info = parts[7].rstrip("\n")

            vaf_dna = vaf_rna = dp_dna = dp_rna = None
            rescued = None
            for field in info.split(";"):
                if field.startswith("VAF_DNA_MEAN="):
                    vaf_dna = float(field[13:])
                elif field.startswith("VAF_RNA_MEAN="):
                    vaf_rna = float(field[13:])
                elif field.startswith("DP_DNA_MEAN="):
                    dp_dna = float(field[12:])
                elif field.startswith("DP_RNA_MEAN="):
                    dp_rna = float(field[12:])
                elif field.startswith("RESCUED="):
                    v = field[8:]
                    rescued = None if v in ("", ".") else v

            norm = _normalize_chrom(chrom)
            rank = _chrom_rank(norm)
            if rank >= 23:
                n_off_chrom += 1
                off_chroms.add(norm or chrom)

            if is_reserved:
                split = "reserved"
            elif verdict == "WARN":
                split = "train"
            else:
                split = _CHROM_TO_SPLIT.get(norm, "train")

            cols["sample_id"].append(sid)
            cols["CHROM"].append(chrom)
            cols["POS"].append(int(parts[1]))
            cols["REF"].append(ref)
            cols["ALT"].append(alt)
            cols["FILTER"].append(filt)
            cols["split"].append(split)
            cols["label_verdict"].append(verdict)
            cols["variant_type"].append(_variant_type(ref, alt))
            cols["VAF_DNA_MEAN"].append(vaf_dna)
            cols["VAF_RNA_MEAN"].append(vaf_rna)
            cols["DP_DNA_MEAN"].append(dp_dna)
            cols["DP_RNA_MEAN"].append(dp_rna)
            cols["RESCUED"].append(rescued)
            cols["_chrom_rank"].append(rank)

    return {
        "sample_id": sid,
        "df": pl.DataFrame(cols, schema=_WORKER_SCHEMA),
        "n_records": n_records,
        "n_selected": len(cols["POS"]),
        "n_off_chrom": n_off_chrom,
        "off_chroms": sorted(off_chroms),
    }


def _add_sub_pool_tags(df: pl.DataFrame) -> pl.DataFrame:
    """Port of neo_var reserve_downstream.compute_sub_pool_tags semantics.

    DP_DNA_MEAN substitutes for neo_var's BAM_DT_DP (caller-mean depth, not
    BAM pileup depth); null DP_DNA_MEAN → is_low_dp False, as in neo_var.
    Tags are computed for every row; is_zero_shot is True iff reserved.
    """
    rna_nonzero = (pl.col("VAF_RNA_MEAN").fill_null(0.0) > 0.0) & (
        pl.col("DP_RNA_MEAN").fill_null(0.0) > 0.0)
    vaf_dna = pl.col("VAF_DNA_MEAN")
    return df.with_columns(
        (pl.col("split") == "reserved").alias("is_zero_shot"),
        ((vaf_dna < 0.15).fill_null(False) & rna_nonzero).alias("is_low_vaf_a"),
        ((vaf_dna >= 0.15).fill_null(False)
         & (vaf_dna < 0.30).fill_null(False)
         & rna_nonzero).alias("is_low_vaf_b"),
        ((pl.col("DP_DNA_MEAN") < 20).fill_null(False) & rna_nonzero).alias("is_low_dp"),
        (pl.col("RESCUED") == "YES").fill_null(False).alias("is_rescued"),
        (pl.col("RESCUED") == "NO").fill_null(False).alias("is_non_rescued"),
        pl.col("variant_type").is_in(["INS", "DEL"]).alias("is_indel"),
    )


def _run_workers(cohort: list, workers: int) -> list:
    tasks = [
        {"sample_id": r["sample_id"],
         "training_label_vcf": r["training_label_vcf"],
         "label_verdict": r["label_qc_verdict"]}
        for r in cohort
    ]
    n = len(tasks)
    results = []

    def _report(i, res):
        off = (f"; non-1-22 chrom records: {res['n_off_chrom']:,} "
               f"({', '.join(res['off_chroms'])})" if res["n_off_chrom"] else "")
        print(f"[{i}/{n}] {res['sample_id']}: {res['n_selected']:,} selected / "
              f"{res['n_records']:,} records{off}", flush=True)

    if workers > 1 and n > 1:
        with Pool(workers) as pool:
            for i, res in enumerate(pool.imap_unordered(_process_sample, tasks), 1):
                _report(i, res)
                results.append(res)
    else:
        for i, task in enumerate(tasks, 1):
            res = _process_sample(task)
            _report(i, res)
            results.append(res)

    # Deterministic concat order regardless of pool completion order.
    results.sort(key=lambda r: r["sample_id"])
    return results


def _assemble(cohort: list, results: list):
    """→ (sample_split.tsv text, sorted tagged variant frame)."""
    by_id = {r["sample_id"]: r for r in results}
    df = pl.concat([r["df"] for r in results], how="vertical")
    df = _add_sub_pool_tags(df)
    df = df.sort(["sample_id", "_chrom_rank", "CHROM", "POS"], maintain_order=True)
    df = df.select(PARQUET_COLUMNS)

    lines = ["\t".join(SAMPLE_TSV_COLUMNS)]
    for row in sorted(cohort, key=lambda r: r["sample_id"]):
        sid = row["sample_id"]
        pool = "reserved" if sid in RESERVED_SAMPLE_IDS else "train_pool"
        lines.append("\t".join([
            sid, row["project_id"], row["patient_id"], row["set_number"],
            row["disease_normalized"], row["label_qc_verdict"], pool,
            row["training_label_vcf"], str(by_id[sid]["n_selected"]),
        ]))
    return "\n".join(lines) + "\n", df


def _split_filter_totals(df: pl.DataFrame) -> dict:
    totals = {s: {f: 0 for f in FILTER_ORDER} for s in SPLIT_ORDER}
    for row in df.group_by("split", "FILTER").len().iter_rows(named=True):
        totals[row["split"]][row["FILTER"]] = row["len"]
    return totals


def _reserved_somatic_tag_counts(df: pl.DataFrame) -> dict:
    rs = df.filter((pl.col("split") == "reserved") & (pl.col("FILTER") == "Somatic"))
    return {t: int(rs[t].sum()) for t in TAG_COLUMNS}


def _totals_lines(totals: dict) -> list:
    lines = [f"  {'split':<10}{'Somatic':>12}{'Germline':>12}{'Reference':>12}{'total':>12}"]
    for s in SPLIT_ORDER:
        t = totals[s]
        lines.append(f"  {s:<10}{t['Somatic']:>12,}{t['Germline']:>12,}"
                     f"{t['Reference']:>12,}{sum(t.values()):>12,}")
    return lines


def _provenance_text(manifest_path: Path, totals: dict, tag_counts: dict,
                     known_limitations: str) -> str:
    lines = [
        "seq2neo split manifest — provenance",
        f"generated_at: {datetime.now(timezone.utc).isoformat(timespec='seconds')}",
        f"git_commit: {_git_head()}",
        f"input_manifest: {manifest_path}",
        f"input_manifest_md5: {_md5(manifest_path)}",
    ]
    if FASTQ_MD5SUMS.exists():
        lines.append(f"fastq_md5sums: {FASTQ_MD5SUMS}")
        lines.append(f"fastq_md5sums_md5: {_md5(FASTQ_MD5SUMS)}")
    else:
        lines.append(f"fastq_md5sums: absent (expected at {FASTQ_MD5SUMS})")
    lines.append("")
    lines.append("per-split per-FILTER totals (selected_variants.parquet):")
    lines.extend(_totals_lines(totals))
    lines.append("")
    lines.append("reserved-pool Somatic-restricted sub-pool tag counts:")
    for t in TAG_COLUMNS:
        lines.append(f"  {t}: {tag_counts[t]:,}")
    lines.append("")
    lines.append(
        "caveat: the low_vaf_a / low_vaf_b / low_dp / rescued sub-pools are "
        "thin per-sample (low_dp especially: is_low_dp requires a present "
        "DP_DNA_MEAN, so it only fires on variants also detected by DNA "
        "callers); they must be evaluated pooled-only across the 5 reserved "
        "samples, never per-sample.")
    lines.append(f"known_limitations: {known_limitations}")
    return "\n".join(lines) + "\n"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", default=str(DEFAULT_MANIFEST),
                    help="Input rerun sample manifest TSV")
    ap.add_argument("--outdir", default=str(DEFAULT_OUTDIR),
                    help="Output directory for the three artifacts")
    ap.add_argument("--workers", type=int, default=8,
                    help="Parallel VCF-parse workers (default: 8; 1 = sequential)")
    ap.add_argument("--check", action="store_true",
                    help="Do not write; fail (STALE, exit 1) if existing "
                         "sample_split.tsv / selected_variants.parquet differ "
                         "from regenerated content (provenance is volatile and "
                         "excluded)")
    args = ap.parse_args(argv)

    manifest_path = Path(args.manifest)
    outdir = Path(args.outdir)
    tsv_path = outdir / SAMPLE_TSV_NAME
    pq_path = outdir / PARQUET_NAME

    t0 = time.time()
    cohort = _load_working_cohort(manifest_path)
    print(f"working cohort: {len(cohort)} samples from {manifest_path}", flush=True)

    results = _run_workers(cohort, args.workers)
    tsv_text, df = _assemble(cohort, results)
    totals = _split_filter_totals(df)
    tag_counts = _reserved_somatic_tag_counts(df)
    elapsed = time.time() - t0

    n_selected = df.height
    n_records = sum(r["n_records"] for r in results)

    if args.check:
        from polars.testing import assert_frame_equal
        stale = []
        if not tsv_path.exists():
            stale.append(f"STALE: missing {tsv_path}")
        elif tsv_path.read_bytes() != tsv_text.encode("utf-8"):
            stale.append(f"STALE: {tsv_path} differs from regenerated content")
        if not pq_path.exists():
            stale.append(f"STALE: missing {pq_path}")
        else:
            try:
                assert_frame_equal(pl.read_parquet(pq_path), df)
            except AssertionError as exc:
                stale.append(f"STALE: {pq_path} differs from regenerated content "
                             f"({exc})")
        if stale:
            for s in stale:
                print(s, file=sys.stderr)
            print("re-run without --check to rewrite", file=sys.stderr)
            sys.exit(1)
        print(f"OK: {tsv_path} and {pq_path} are current "
              f"({len(cohort)} samples, {n_selected:,} variants)")
        return 0

    outdir.mkdir(parents=True, exist_ok=True)
    tsv_path.write_text(tsv_text)
    df.write_parquet(pq_path, compression="zstd")

    known = sorted({r.get("known_limitations", "") for r in cohort} - {""})
    known_limitations = ";".join(known) if known else KNOWN_LIMITATIONS
    (outdir / PROVENANCE_NAME).write_text(
        _provenance_text(manifest_path, totals, tag_counts, known_limitations))

    print(f"\nwrote {tsv_path}  ({len(cohort)} samples)")
    print(f"wrote {pq_path}  ({n_selected:,} variants)")
    print(f"wrote {outdir / PROVENANCE_NAME}")
    print(f"\nscanned {n_records:,} records, selected {n_selected:,} "
      f"({n_selected / max(n_records, 1):.1%}) in {elapsed:.1f}s "
      f"({args.workers} workers)")
    print("per-split per-FILTER totals:")
    for line in _totals_lines(totals):
        print(line)
    print("reserved-pool Somatic-restricted sub-pool tag counts:")
    for t in TAG_COLUMNS:
        print(f"  {t}: {tag_counts[t]:,}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
