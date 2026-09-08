#!/usr/bin/env python3
"""Validate SEQC2 hybrid-realignment ingress and final run contracts."""

from __future__ import annotations

import argparse
import csv
import gzip
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path


# Require the actual process in its own workflow scope, not a substring
# supplied by a first-pass caller or an enclosing subworkflow name.
SECOND_PASS_PROCESSES = {
    "candidate extraction": ("PREPARE_REALIGNMENT_VCF", "VCF2BED"),
    "paired-read validation": ("PREPARE_REALIGNMENT_VCF", "VALIDATE_READ_IDS"),
    "HISAT2": ("PREPARE_REALIGNMENT_VCF", "HISAT2_ALIGN"),
    "Mutect2": ("RNA_REALIGNMENT_WORKFLOW", "MUTECT2_PAIRED"),
    "Strelka2": ("RNA_REALIGNMENT_WORKFLOW", "STRELKA_SOMATIC"),
    "DeepSomatic": ("RNA_REALIGNMENT_WORKFLOW", "DEEPSOMATIC"),
    "RNA consensus": ("RNA_REALIGNMENT_WORKFLOW", "VCF_CONSENSUS"),
    "second rescue": ("SECOND_RESCUE_WORKFLOW", "VCF_RESCUE"),
    "RNA editing": ("SECOND_RESCUE_WORKFLOW", "RNA_EDITING_ANNOTATION"),
    "COSMIC/gnomAD": ("SECOND_RESCUE_WORKFLOW", "COSMIC_GNOMAD_ANNOTATION"),
    "VEP": ("SECOND_RESCUE_WORKFLOW", "ENSEMBLVEP_VEP"),
}

SECOND_PASS_ARTIFACTS = (
    "vcf_realignment/**/**.deepsomatic.vcf.gz",
    "vcf_realignment/**/**.mutect2.filtered.vcf.gz",
    "vcf_realignment/**/**.strelka.variants.vcf.gz",
    "vcf_realignment/consensus/**/*.consensus.vcf.gz",
    "vcf_realignment/rescue/**/*.rescued.vcf.gz",
)
FINAL_VCF = (
    "vcf_realignment/rescue/*/*.rescue.filtered.stripped.vep.vcf.gz"
)


class ZeroCandidates(ValueError):
    """The first-pass candidate set was validly empty, not missing or corrupt."""


def fail(message: str) -> None:
    raise ValueError(message)


def read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        fail("samplesheet is empty")
    return rows


def validate_manifest(rows: list[dict[str, str]]) -> None:
    by_patient: dict[str, list[dict[str, str]]] = defaultdict(list)
    for number, row in enumerate(rows, 2):
        patient = row.get("patient", "").strip()
        sample = row.get("sample", "").strip()
        stage = row.get("input_stage", "").strip()
        if not patient or not sample:
            fail(f"row {number}: patient and sample are required")
        if stage not in {"raw_reads", "raw_alignment", "caller_ready"}:
            fail(f"row {number}: unsupported input_stage {stage!r}")
        has_fastq = bool(row.get("fastq_1", "").strip() or row.get("fastq_2", "").strip())
        has_bam = bool(row.get("bam", "").strip() or row.get("bai", "").strip())
        if has_fastq == has_bam:
            fail(f"row {number}: exactly one FASTQ or BAM/BAI payload is required")
        if has_fastq:
            if stage != "raw_reads" or not row.get("fastq_1", "").strip() or not row.get("fastq_2", "").strip():
                fail(f"row {number}: raw_reads requires both FASTQs")
        elif stage not in {"caller_ready", "raw_alignment"} or not row.get("bam", "").strip() or not row.get("bai", "").strip():
            fail(f"row {number}: raw_alignment/caller_ready requires BAM and BAI")
        for field in ("fastq_1", "fastq_2", "bam", "bai"):
            value = row.get(field, "").strip()
            if value and not Path(value).is_file():
                fail(f"row {number}: {field} does not exist: {value}")
        by_patient[patient].append(row)

    for patient, patient_rows in by_patient.items():
        roles = defaultdict(list)
        for row in patient_rows:
            try:
                roles[int(row.get("status", ""))].append(row)
            except ValueError:
                fail(f"patient {patient}: status must be 0, 1, or 2")
        for status, label in ((0, "DN"), (1, "DT")):
            if len(roles[status]) != 1:
                fail(f"patient {patient}: expected exactly one {label} row")
        rt_rows = roles[2]
        if not rt_rows:
            fail(f"patient {patient}: expected at least one RT row")
        if len({row["sample"].strip() for row in rt_rows}) != 1:
            fail(f"patient {patient}: RT libraries must share one logical sample")
        libraries = [row.get("library", "").strip() for row in rt_rows]
        if any(not library for library in libraries) or len(libraries) != len(set(libraries)):
            fail(f"patient {patient}: RT libraries must be unique and non-empty")


def fasta_contigs(fasta: Path) -> dict[str, int]:
    fai = Path(str(fasta) + ".fai")
    if not fasta.is_file() or not fai.is_file():
        fail(f"reference FASTA and index are required: {fasta}")
    contigs = {}
    for line in fai.read_text().splitlines():
        fields = line.split("\t")
        if len(fields) < 2:
            fail(f"malformed FASTA index record: {line}")
        contigs[fields[0]] = int(fields[1])
    if not contigs:
        fail("reference FASTA index has no contigs")
    return contigs


def validate_bam_dictionary(rows: list[dict[str, str]], reference: dict[str, int]) -> None:
    samtools = shutil.which("samtools")
    if not samtools:
        fail("samtools is required for caller-ready BAM preflight")
    for row in rows:
        if row.get("input_stage", "").strip() != "caller_ready":
            continue
        bam = row["bam"].strip()
        quickcheck = subprocess.run([samtools, "quickcheck", "-v", bam], text=True, capture_output=True)
        if quickcheck.returncode:
            fail(f"caller-ready BAM failed samtools quickcheck: {bam}: {quickcheck.stderr.strip()}")
        header = subprocess.run([samtools, "view", "-H", bam], text=True, capture_output=True)
        if header.returncode:
            fail(f"cannot read caller-ready BAM header: {bam}: {header.stderr.strip()}")
        observed = []
        for line in header.stdout.splitlines():
            if not line.startswith("@SQ"):
                continue
            fields = dict(field.split(":", 1) for field in line.split("\t")[1:] if ":" in field)
            name, length = fields.get("SN"), fields.get("LN")
            observed.append((name, int(length)))
        observed_by_name = dict(observed)
        shared_mismatch = [name for name, length in observed_by_name.items()
                           if name in reference and reference[name] != length]
        if shared_mismatch:
            fail(f"caller-ready BAM has incompatible shared contigs {shared_mismatch}: {bam}")
        # The pipeline's configured normalize policy safely drops BAM-only
        # decoys and supplies reference-only unused contigs. Record this
        # explicit, auditable case rather than treating file presence as proof.
        if observed != list(reference.items()):
            print(f"normalizable caller-ready dictionary difference: {bam}")


def validate_hisat2_resources(directory: Path, splicesites: Path) -> None:
    small = sorted(directory.glob("*.1.ht2"))
    large = sorted(directory.glob("*.1.ht2l"))
    if small and large:
        fail(f"HISAT2 index mixes .ht2 and .ht2l formats: {directory}")
    candidates = small or large
    basenames = {path.name.rsplit(".1.", 1)[0] for path in candidates}
    if len(basenames) != 1:
        fail(f"HISAT2 directory must contain one index basename, found {sorted(basenames)}")
    prefix = next(iter(basenames))
    suffixes = [directory / f"{prefix}.{number}.ht2" for number in range(1, 9)]
    suffixes_l = [directory / f"{prefix}.{number}.ht2l" for number in range(1, 9)]
    if not all(path.is_file() and path.stat().st_size > 0 for path in suffixes) and not all(path.is_file() and path.stat().st_size > 0 for path in suffixes_l):
        fail(f"HISAT2 directory must contain one complete non-empty .ht2 or .ht2l set: {directory}")
    if not splicesites.is_file() or splicesites.stat().st_size == 0:
        fail(f"HISAT2 splice-site file is missing or empty: {splicesites}")


def preflight(samplesheet: Path, fasta: Path, hisat2_dir: Path, splicesites: Path) -> None:
    rows = read_manifest(samplesheet)
    validate_manifest(rows)
    validate_hisat2_resources(hisat2_dir, splicesites)
    validate_bam_dictionary(rows, fasta_contigs(fasta))


def latest_trace(outdir: Path) -> Path:
    traces = sorted((outdir / "pipeline_info").glob("execution_trace*.txt"), key=lambda path: path.stat().st_mtime)
    if not traces:
        fail("execution trace is missing")
    return traces[-1]


def validate_trace(path: Path) -> None:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        fail("execution trace is empty")
    failed = [row for row in rows if row.get("status", "").upper() == "FAILED"]
    if failed:
        fail("execution trace contains failed tasks")
    names = "\n".join(" ".join((row.get("process", ""), row.get("name", ""))).upper() for row in rows)
    if "VCF2BED" not in names:
        fail("candidate extraction is missing; absence is not evidence of zero candidates")
    candidate_rows = [row for row in rows if "VCF2BED" in
                      (row.get("process", "") + row.get("name", "")).upper()]
    if any(row.get("status", "").upper() not in {"COMPLETED", "CACHED"} for row in candidate_rows):
        fail("candidate extraction did not complete successfully")
    candidate_dir = path.parent.parent / "vcf_realignment" / "vcf2bed"
    candidate_beds = list(candidate_dir.glob("*/*.bed"))
    if (len(candidate_beds) == len(candidate_rows)
            and all(bed.stat().st_size == 0 for bed in candidate_beds)):
        if "HISAT2_ALIGN" in names:
            fail("empty candidate BED contradicts executed HISAT2 alignment")
        raise ZeroCandidates("RT consensus was converted successfully to an empty candidate BED")
    missing = []
    for label, (scope, process) in SECOND_PASS_PROCESSES.items():
        matches = []
        for row in rows:
            qualified = (row.get("process") or row.get("name", "")).upper()
            parts = qualified.split(" (", 1)[0].split(":")
            if scope in parts[:-1] and parts[-1] == process:
                matches.append(row)
        if not matches or any(
            row.get("status", "").upper() not in {"COMPLETED", "CACHED"}
            or row.get("exit", "").strip() not in {"", "0"}
            for row in matches
        ):
            missing.append(label)
    if missing:
        fail("second-pass process groups missing or unsuccessful: " + ", ".join(missing))


def validate_final_vcf(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        fail(f"final VCF is missing or empty: {path}")
    index = Path(str(path) + ".tbi")
    if not index.is_file() or index.stat().st_size == 0:
        fail(f"final VCF index is missing or empty: {index}")
    try:
        with gzip.open(path, "rt") as handle:
            header = [line.rstrip("\n") for line in handle if line.startswith("#")]
    except OSError as error:
        fail(f"final VCF is not readable gzip: {error}")
    chrom = next((line for line in header if line.startswith("#CHROM")), "")
    if len(chrom.split("\t")) != 8:
        fail("final VCF must retain the stripped eight-column contract")
    tabix = shutil.which("tabix")
    if not tabix:
        fail("tabix is required to validate the final VCF index")
    result = subprocess.run([tabix, "-l", str(path)], text=True, capture_output=True)
    if result.returncode != 0:
        fail(f"final VCF index is unreadable: {result.stderr.strip()}")


def complete(outdir: Path) -> None:
    trace = latest_trace(outdir)
    validate_trace(trace)
    missing = [pattern for pattern in SECOND_PASS_ARTIFACTS if not list(outdir.glob(pattern))]
    if missing:
        fail("second-pass artifacts missing: " + ", ".join(missing))
    final_vcfs = list(outdir.glob(FINAL_VCF))
    if len(final_vcfs) != 1:
        fail(f"expected exactly one annotated rescue VCF, found {len(final_vcfs)}")
    trace_text = trace.read_text()
    if final_vcfs[0].parent.name not in trace_text:
        fail("final rescue identity is not present in the current execution trace")
    validate_final_vcf(final_vcfs[0])


def main() -> int:
    parser = argparse.ArgumentParser()
    commands = parser.add_subparsers(dest="command", required=True)
    preflight_parser = commands.add_parser("preflight")
    preflight_parser.add_argument("--samplesheet", type=Path, required=True)
    preflight_parser.add_argument("--fasta", type=Path, required=True)
    preflight_parser.add_argument("--hisat2-dir", type=Path, required=True)
    preflight_parser.add_argument("--splicesites", type=Path, required=True)
    complete_parser = commands.add_parser("complete")
    complete_parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    try:
        if args.command == "preflight":
            preflight(args.samplesheet, args.fasta, args.hisat2_dir, args.splicesites)
        else:
            complete(args.outdir)
    except ZeroCandidates as error:
        print(f"ZERO_CANDIDATES: {error}", file=sys.stderr)
        return 3
    except ValueError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1
    print(f"hybrid realignment {args.command} validation PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
