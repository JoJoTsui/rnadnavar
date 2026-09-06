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


SECOND_PASS_PROCESSES = {
    "candidate extraction": ("VCF2BED",),
    "HISAT2": ("FASTQ_ALIGN_HISAT2",),
    "Mutect2": ("MUTECT2",),
    "Strelka2": ("STRELKA",),
    "DeepSomatic": ("DEEPSOMATIC",),
    "RNA consensus": ("VCF_CONSENSUS",),
    "second rescue": ("SECOND_RESCUE_WORKFLOW",),
    "RNA editing": ("RNA_EDITING_ANNOTATION",),
    "COSMIC/gnomAD": ("COSMIC_GNOMAD_ANNOTATION",),
    "VEP": ("ENSEMBLVEP_VEP",),
}
SECOND_PASS_ARTIFACTS = (
    "vcf_realignment/**/**.deepsomatic.vcf.gz",
    "vcf_realignment/**/**.mutect2.filtered.vcf.gz",
    "vcf_realignment/**/**.strelka.variants.vcf.gz",
    "vcf_realignment/consensus/**/*.consensus.vcf.gz",
    "vcf_realignment/rescue/**/*.rescued.vcf.gz",
)
FINAL_VCF = (
    "vcf_realignment/rescue/WES_LL_RT_1_vs_WES_LL_N_1/"
    "WES_LL_RT_1_vs_WES_LL_N_1.rescue.filtered.stripped.vep.vcf.gz"
)


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
    indexes = sorted(directory.glob("*.ht2"))
    if len(indexes) != 8 or any(index.stat().st_size == 0 for index in indexes):
        fail(f"HISAT2 directory must contain eight non-empty .ht2 indexes: {directory}")
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
    missing = [label for label, patterns in SECOND_PASS_PROCESSES.items()
               if not any(pattern in names for pattern in patterns)]
    if missing:
        fail("second-pass process groups missing: " + ", ".join(missing))
    non_realign = [label for label, patterns in SECOND_PASS_PROCESSES.items()
                   if label not in {"candidate extraction", "HISAT2"}
                   and not any(f"{pattern}_REALIGN" in names for pattern in patterns)]
    if non_realign:
        fail("second-pass realigned process groups missing: " + ", ".join(non_realign))


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
    validate_trace(latest_trace(outdir))
    missing = [pattern for pattern in SECOND_PASS_ARTIFACTS if not list(outdir.glob(pattern))]
    if missing:
        fail("second-pass artifacts missing: " + ", ".join(missing))
    validate_final_vcf(outdir / FINAL_VCF)


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
    except ValueError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1
    print(f"hybrid realignment {args.command} validation PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
