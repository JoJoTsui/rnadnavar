"""Focused public checks for SEQC2 hybrid-realignment validation."""

from __future__ import annotations

import importlib.util
import gzip
from types import SimpleNamespace
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    "hybrid_realign_validation",
    ROOT / "examples/seqc2/hybrid/scripts/validate_realign_run.py",
)
validator = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(validator)


def test_second_pass_artifact_globs_are_pathlib_compatible(tmp_path):
    """Completion checks must not reject valid output with an invalid glob."""
    files = {
        "vcf_realignment/normalized/deepsomatic/sample/sample.deepsomatic.vcf.gz",
        "vcf_realignment/normalized/mutect2/sample/sample.mutect2.filtered.vcf.gz",
        "vcf_realignment/normalized/strelka/sample/sample.strelka.variants.vcf.gz",
        "vcf_realignment/consensus/sample/sample.consensus.vcf.gz",
        "vcf_realignment/rescue/sample/sample.rescued.vcf.gz",
    }
    for relative in files:
        path = tmp_path / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.touch()

    matches = [list(tmp_path.glob(pattern)) for pattern in validator.SECOND_PASS_ARTIFACTS]
    assert all(matches)


def row(tmp_path: Path, status: int, sample: str, library: str) -> dict[str, str]:
    fastq_1 = tmp_path / "r1.fastq.gz"
    fastq_2 = tmp_path / "r2.fastq.gz"
    fastq_1.touch(exist_ok=True)
    fastq_2.touch(exist_ok=True)
    return {
        "patient": "P",
        "status": str(status),
        "sample": sample,
        "input_stage": "raw_reads",
        "library": library,
        "fastq_1": str(fastq_1),
        "fastq_2": str(fastq_2),
        "bam": "",
        "bai": "",
    }


def test_manifest_accepts_one_logical_pooled_rt_sample(tmp_path):
    rows = [row(tmp_path, 0, "DN", "dn"), row(tmp_path, 1, "DT", "dt")]
    rows.extend(row(tmp_path, 2, "RT", library) for library in ("A", "B", "C"))

    validator.validate_manifest(rows)


def test_manifest_rejects_cross_sample_pooled_rt(tmp_path):
    rows = [row(tmp_path, 0, "DN", "dn"), row(tmp_path, 1, "DT", "dt"), row(tmp_path, 2, "RT_A", "A"), row(tmp_path, 2, "RT_B", "B")]

    with pytest.raises(ValueError, match="one logical sample"):
        validator.validate_manifest(rows)


def test_manifest_rejects_duplicate_rna_library(tmp_path):
    rows = [row(tmp_path, 0, "DN", "dn"), row(tmp_path, 1, "DT", "dt"), row(tmp_path, 2, "RT", "A"), row(tmp_path, 2, "RT", "A")]

    with pytest.raises(ValueError, match="unique and non-empty"):
        validator.validate_manifest(rows)


def trace(tmp_path, name='VCF2BED', status='COMPLETED'):
    path = tmp_path / 'pipeline_info/trace.txt'
    path.parent.mkdir()
    path.write_text(f'name\tstatus\n{name}\t{status}\n')
    return path


def test_missing_candidate_process_is_not_zero_candidates(tmp_path):
    with pytest.raises(ValueError, match='absence is not evidence'):
        validator.validate_trace(trace(tmp_path, 'VCF_CONSENSUS'))


def test_verified_empty_bed_is_zero_candidates(tmp_path):
    bed = tmp_path / 'vcf_realignment/vcf2bed/RT/RT.bed'
    bed.parent.mkdir(parents=True)
    bed.touch()
    with pytest.raises(validator.ZeroCandidates, match='empty candidate BED'):
        validator.validate_trace(trace(tmp_path))


def test_nonempty_candidates_without_alignment_fail(tmp_path):
    bed = tmp_path / 'vcf_realignment/vcf2bed/RT/RT.bed'
    bed.parent.mkdir(parents=True)
    bed.write_text('chr7\t9\t10\n')
    with pytest.raises(ValueError, match='second-pass process groups missing'):
        validator.validate_trace(trace(tmp_path))


def test_hisat2_resource_validator_accepts_complete_large_index_set(tmp_path):
    for suffix in range(1, 9):
        path = tmp_path / f"assembly.{suffix}.ht2l"
        path.write_bytes(b"index")
    splice = tmp_path / "splice.txt"
    splice.write_text("chr1\t1\t5\t+\n")

    validator.validate_hisat2_resources(tmp_path, splice)


def test_hisat2_resource_validator_rejects_mixed_or_incomplete_sets(tmp_path):
    for suffix in range(1, 8):
        (tmp_path / f"assembly.{suffix}.ht2").write_bytes(b"index")
    (tmp_path / "other.8.ht2").write_bytes(b"index")
    splice = tmp_path / "splice.txt"
    splice.write_text("chr1\t1\t5\t+\n")

    with pytest.raises(ValueError, match="one complete"):
        validator.validate_hisat2_resources(tmp_path, splice)


def test_hisat2_resource_validator_rejects_mixed_index_formats(tmp_path):
    for suffix in range(1, 9):
        (tmp_path / f"assembly.{suffix}.ht2").write_bytes(b"index")
        (tmp_path / f"other.{suffix}.ht2l").write_bytes(b"index")
    splice = tmp_path / "splice.txt"
    splice.write_text("chr1\t1\t5\t+\n")

    with pytest.raises(ValueError, match="mixes"):
        validator.validate_hisat2_resources(tmp_path, splice)


def complete_trace_rows():
    return [
        ('PREPARE_REALIGNMENT_VCF:VCF2BED', 'COMPLETED'),
        ('PREPARE_REALIGNMENT_VCF:VALIDATE_READ_IDS', 'CACHED'),
        ('PREPARE_REALIGNMENT_VCF:FASTQ_ALIGN_HISAT2:HISAT2_ALIGN', 'COMPLETED'),
        ('RNA_REALIGNMENT_WORKFLOW:CALLING:MUTECT2_PAIRED', 'COMPLETED'),
        ('RNA_REALIGNMENT_WORKFLOW:CALLING:STRELKA_SOMATIC', 'COMPLETED'),
        ('RNA_REALIGNMENT_WORKFLOW:CALLING:DEEPSOMATIC', 'COMPLETED'),
        ('RNA_REALIGNMENT_WORKFLOW:VCF_CONSENSUS_WORKFLOW:VCF_CONSENSUS', 'COMPLETED'),
        ('SECOND_RESCUE_WORKFLOW:VCF_RESCUE_WORKFLOW:VCF_RESCUE', 'COMPLETED'),
        ('SECOND_RESCUE_WORKFLOW:POST:RNA_EDITING_ANNOTATION', 'COMPLETED'),
        ('SECOND_RESCUE_WORKFLOW:POST:COSMIC_GNOMAD_ANNOTATION', 'COMPLETED'),
        ('SECOND_RESCUE_WORKFLOW:VCF_ANNOTATE:VCF_ANNOTATE_ENSEMBLVEP:ENSEMBLVEP_VEP', 'COMPLETED'),
    ]


def write_trace(tmp_path, rows):
    path = tmp_path / 'pipeline_info/trace.txt'
    path.parent.mkdir(exist_ok=True)
    path.write_text('name\tstatus\n' + ''.join(f'ROOT:{name} (sample)\t{status}\n' for name, status in rows))
    return path


def test_complete_scoped_trace_accepts_cached_tasks(tmp_path):
    validator.validate_trace(write_trace(tmp_path, complete_trace_rows()))


@pytest.mark.parametrize('leaf', ['DEEPSOMATIC', 'MUTECT2_PAIRED', 'STRELKA_SOMATIC',
                                  'VCF_CONSENSUS', 'VCF_RESCUE', 'RNA_EDITING_ANNOTATION',
                                  'COSMIC_GNOMAD_ANNOTATION', 'ENSEMBLVEP_VEP'])
def test_first_pass_cannot_replace_second_pass_task(tmp_path, leaf):
    rows = [(('BAM_PROCESSING:' + name.split(':')[-1]) if name.endswith(':' + leaf) else name, status)
            for name, status in complete_trace_rows()]
    with pytest.raises(ValueError, match='second-pass process groups missing'):
        validator.validate_trace(write_trace(tmp_path, rows))


@pytest.mark.parametrize('status', ['ABORTED', 'RUNNING', 'SUBMITTED'])
def test_incomplete_vep_is_rejected_even_with_a_successful_task(tmp_path, status):
    rows = complete_trace_rows()
    rows.append((rows[-1][0], status))
    with pytest.raises(ValueError, match='did not complete successfully'):
        validator.validate_trace(write_trace(tmp_path, rows))


def final_vcf(tmp_path, monkeypatch, csq_header=True, info='CSQ=C|missense_variant', record=True):
    path = tmp_path / 'sample.rescue.filtered.stripped.vep.vcf.gz'
    with gzip.open(path, 'wt') as handle:
        handle.write('##fileformat=VCFv4.2\n')
        if csq_header:
            handle.write('##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence">\n')
        handle.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
        if record:
            handle.write(f'chr1\t1\t.\tA\tC\t.\tSomatic\t{info}\n')
    Path(str(path) + '.tbi').write_bytes(b'fixture-index')
    monkeypatch.setattr(validator.shutil, 'which', lambda name: '/fixture/tabix')
    monkeypatch.setattr(validator.subprocess, 'run', lambda *a, **kw: SimpleNamespace(returncode=0))
    return path


@pytest.mark.parametrize('record', [True, False])
def test_annotated_vcf_including_empty_result_is_valid(tmp_path, monkeypatch, record):
    validator.validate_final_vcf(final_vcf(tmp_path, monkeypatch, record=record))


def test_unannotated_vcf_is_rejected(tmp_path, monkeypatch):
    path = final_vcf(tmp_path, monkeypatch, csq_header=False, info='.')
    with pytest.raises(ValueError, match='CSQ'):
        validator.validate_final_vcf(path)


@pytest.mark.parametrize('info', ['.', 'CSQ=', 'CSQ=.'])
def test_unannotated_record_is_rejected(tmp_path, monkeypatch, info):
    path = final_vcf(tmp_path, monkeypatch, info=info)
    with pytest.raises(ValueError, match='CSQ'):
        validator.validate_final_vcf(path)


@pytest.mark.parametrize('leaf,replacement', [
    ('HISAT2_ALIGN', 'SAMTOOLS_SORT'),
    ('MUTECT2_PAIRED', 'FILTERMUTECTCALLS_REALIGN'),
    ('STRELKA_SOMATIC', 'MERGE_STRELKA'),
])
def test_postprocessing_is_not_evidence_of_caller_execution(tmp_path, leaf, replacement):
    rows = [(name.removesuffix(leaf) + replacement if name.endswith(':' + leaf) else name, status)
            for name, status in complete_trace_rows()]
    with pytest.raises(ValueError, match='second-pass process groups missing'):
        validator.validate_trace(write_trace(tmp_path, rows))


@pytest.mark.parametrize('identity', ['sample', 'stale'])
def test_completion_binds_output_identity_to_vep_task(tmp_path, monkeypatch, identity):
    trace_path = write_trace(tmp_path, complete_trace_rows())
    trace_path.rename(trace_path.with_name('execution_trace.txt'))
    path = final_vcf(tmp_path, monkeypatch)
    dest = tmp_path / 'vcf_realignment/rescue' / identity / path.name
    dest.parent.mkdir(parents=True)
    path.rename(dest)
    Path(str(path) + '.tbi').rename(Path(str(dest) + '.tbi'))
    # Other artifacts have their own existence checks; isolate VEP identity.
    monkeypatch.setattr(validator, 'SECOND_PASS_ARTIFACTS', ())
    if identity == 'sample':
        validator.complete(tmp_path)
    else:
        trace_path = trace_path.with_name('execution_trace.txt')
        with trace_path.open('a') as handle:
            handle.write('ROOT:SECOND_RESCUE_WORKFLOW:VCF_RESCUE (stale)\tCOMPLETED\n')
        with pytest.raises(ValueError, match='successful second-rescue VEP task'):
            validator.complete(tmp_path)
