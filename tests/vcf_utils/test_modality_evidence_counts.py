"""Observed records must not masquerade as eligible RNA-editing evidence."""

import ast
import sys
from pathlib import Path

import pysam
from cyvcf2 import VCF

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "bin"))

from vcf_utils.aggregation import aggregate_genotypes
from vcf_utils.io_utils import write_union_vcf
from vcf_utils.rna_editing_core import classify_rna_editing_evidence


def _extract_for_editing(record):
    # The CLI imports annotation dependencies requiring newer Python syntax.
    # Execute its actual extraction method independently of those optional tools.
    path = Path(__file__).resolve().parents[2] / "bin" / "annotate_rna_editing.py"
    tree = ast.parse(path.read_text())
    cls = next(n for n in tree.body if isinstance(n, ast.ClassDef) and n.name == "RNAEditingAnnotator")
    method = next(n for n in cls.body if isinstance(n, ast.FunctionDef) and n.name == "_extract_variant_data_pysam")
    namespace = {"Dict": dict}
    exec(compile(ast.Module(body=[method], type_ignores=[]), str(path), "exec"), namespace)
    return namespace[method.name](None, record)


def test_writer_distinguishes_observed_eligible_and_somatic(tmp_path):
    template = tmp_path / "template.vcf"
    template.write_text(
        '##fileformat=VCFv4.2\n##contig=<ID=chr1,length=10000>\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
    )
    callers = ["DNA_mutect2", "DNA_strelka", "DNA_deepsomatic", "RNA_mutect2", "RNA_consensus"]
    labels = ["Artifact", "Germline", "Somatic", "Somatic", "Somatic"]
    data = dict(
        CHROM="chr1", POS=100, REF="A", ALT="G", is_snv=True,
        callers=callers, filters_original=labels, filters_normalized=labels,
        filters_category=labels, qualities=[], genotypes={}, ids=[],
        support_callers={"DNA_strelka", "DNA_deepsomatic", "RNA_mutect2"},
        passes_consensus=True, gt_aggregated=aggregate_genotypes({}, callers),
    )
    out = tmp_path / "out.vcf"
    reader = VCF(str(template))
    write_union_vcf(
        {"chr1:100:A:G": data}, reader, "sample", str(out), "vcf", callers[:-1],
        modality_map={c: c.split("_")[0] for c in callers},
    )
    reader.close()
    with pysam.VariantFile(out) as vcf:
        rec = next(vcf)
        assert rec.info["N_DNA_CALLERS_SUPPORT"] == 3  # compatibility
        assert rec.info["DNA_SUPPORT"] == 3
        assert rec.info["N_DNA_CALLERS_OBSERVED"] == 3
        assert rec.info["N_DNA_CALLERS_ELIGIBLE"] == 2
        assert rec.info["N_DNA_CALLERS_SOMATIC"] == 1
        assert rec.info["N_RNA_CALLERS_OBSERVED"] == 1  # consensus excluded
        assert rec.info["N_RNA_CALLERS_ELIGIBLE"] == 1
        assert rec.info["N_RNA_CALLERS_SOMATIC"] == 1
        assert data["final_classification"] == rec.info["UNIFIED_FILTER"]


def _editing_data(eligible=None):
    header = pysam.VariantHeader()
    header.contigs.add("chr1")
    for modality in ("DNA", "RNA"):
        header.info.add(f"N_{modality}_CALLERS_SUPPORT", 1, "Integer", "Legacy observed")
        if eligible is not None:
            header.info.add(f"N_{modality}_CALLERS_ELIGIBLE", 1, "Integer", "Eligible")
    for modality in ("DNA", "RNA"):
        header.info.add(f"VAF_{modality}_MEAN", 1, "Float", "VAF")
    rec = header.new_record(contig="chr1", start=99, alleles=("A", "G"))
    rec.info["N_DNA_CALLERS_SUPPORT"] = 3
    rec.info["N_RNA_CALLERS_SUPPORT"] = 2
    rec.info["VAF_DNA_MEAN"] = 0.0
    rec.info["VAF_RNA_MEAN"] = 0.4
    if eligible is not None:
        rec.info["N_DNA_CALLERS_ELIGIBLE"] = eligible
        rec.info["N_RNA_CALLERS_ELIGIBLE"] = 2
    return _extract_for_editing(rec)


def test_rejected_dna_presence_does_not_suppress_editing():
    data = _editing_data(eligible=0)
    assert data["N_DNA_CALLERS_SUPPORT"] == 0
    assert classify_rna_editing_evidence(data, True, 2) == "VERY_HIGH"


def test_eligible_dna_still_protects_from_editing_relabel():
    assert classify_rna_editing_evidence(_editing_data(eligible=1), True, 2) == "MEDIUM"


def test_legacy_vcf_retains_previous_behavior():
    data = _editing_data()
    assert data["N_DNA_CALLERS_SUPPORT"] == 3
    assert classify_rna_editing_evidence(data, True, 2) == "MEDIUM"
