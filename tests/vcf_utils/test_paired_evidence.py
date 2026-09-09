from vcf_utils.classification import compute_unified_classification_rescue
from vcf_utils.aggregation import (
    aggregate_genotypes,
    aggregate_variants,
    resolve_normal_sample_index,
    resolve_tumor_sample_index,
)

import pytest
import pysam
from cyvcf2 import VCF
from vcf_utils.aggregation import read_variants_from_vcf, extract_genotype_info
from vcf_utils.io_utils import write_union_vcf


@pytest.mark.parametrize("modality", [None, "DNA"])
def test_tumor_only_normal_evidence_round_trips_as_missing(tmp_path, modality):
    source = tmp_path / "sample.deepsomatic.vcf"
    source.write_text(
        '##fileformat=VCFv4.2\n'
        '##contig=<ID=chr1,length=1000>\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">\n'
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depth">\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tWES_LL_T_1\n'
        'chr1\t10\t.\tA\tG\t50\tPASS\t.\tGT:DP:AD\t0/1:20:15,5\n'
    )
    caller = "DNA_deepsomatic" if modality else "deepsomatic"
    variants = read_variants_from_vcf(str(source), caller)
    data = aggregate_variants([(caller, variants, modality)], 1, 1)
    output = tmp_path / "union.vcf"
    template = VCF(str(source))
    write_union_vcf(data, template, "unused", str(output), "vcf", [caller],
                    modality_map={caller: modality} if modality else None,
                    snv_threshold=1, indel_threshold=1)
    template.close()
    with pysam.VariantFile(output) as reader:
        record = next(reader)
        for field in ("GT", "DP", "AD", "VAF", "VAF_SOURCE"):
            value = record.info[f"NORMAL_{field}_BY_CALLER"]
            assert (value[0] if isinstance(value, tuple) else value) == f"{caller}:."
        value = record.info["VAF_BY_CALLER"]
        assert (value[0] if isinstance(value, tuple) else value) == f"{caller}:0.2500"
        assert not reader.header.samples


@pytest.mark.parametrize("indel", [False, True])
def test_strelka_derived_vaf_has_provenance(tmp_path, indel):
    fields = ["TAR", "TIR"] if indel else ["AU", "CU", "GU", "TU"]
    values = "15,0:5,0" if indel else "15,0:0,0:5,0:0,0"
    source = tmp_path / "sample.strelka.vcf"
    source.write_text(
        '##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n'
        + ''.join(f'##FORMAT=<ID={f},Number=2,Type=Integer,Description="Counts">\n' for f in fields)
        + '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\tTUMOR\n'
        + f'chr1\t10\t.\tA\t{"AT" if indel else "G"}\t50\tPASS\t.\t{":".join(fields)}\t{values}\t{values}\n'
    )
    reader = VCF(str(source))
    variant = next(reader)
    for sample_idx in (0, 1):
        evidence = extract_genotype_info(variant, "strelka", sample_idx=sample_idx)
        assert evidence["VAF"] == 0.25
        assert evidence["VAF_SOURCE"] == "derived"
    reader.close()


def test_paired_sample_resolution_is_explicit_and_complementary():
    samples = ["WES_LL_N_1", "WES_LL_T_1"]
    assert resolve_tumor_sample_index(samples, "mutect2") == 1
    assert resolve_normal_sample_index(samples, "mutect2") == 0


def test_tumor_only_vcf_has_unavailable_normal():
    assert resolve_normal_sample_index(["tumor"], "deepsomatic") is None


def test_ambiguous_multi_sample_vcf_does_not_guess_normal():
    assert resolve_normal_sample_index(["A", "B", "C"], "unknown") is None


def test_aggregate_genotypes_keeps_caller_af_distinct_from_ad():
    result = aggregate_genotypes(
        {
            "mutect2": {
                "GT": "0/1",
                "DP": 20,
                "AD": "18,2",
                "VAF": 0.071,
                "GQ": 40,
            }
        },
        ["mutect2"],
    )
    assert result["dp_by_caller"] == [20]
    assert result["vaf_by_caller"] == [0.071]
    assert result["alt_count_by_caller"] == [2]


def test_multiallelic_alt_support_uses_all_alternates():
    result = aggregate_genotypes(
        {"caller": {"GT": "1/2", "ALT_INDICES": [1, 2], "DP": 30, "AD": "10,3,7", "VAF": None, "GQ": 50}},
        ["caller"],
    )
    assert result["alt_count_by_caller"] == [7]


def test_unverified_rna_only_nomination_is_not_somatic():
    data = {
        "callers": ["RNA_consensus"],
        "filters_normalized": ["Somatic"],
        "caller_modality_map": {"RNA_consensus": "RNA"},
        "is_snv": True,
        "support_callers": {"RNA_consensus"},
        "dna_verification_status": "inconclusive",
    }
    assert compute_unified_classification_rescue(
        data, {"RNA_consensus": "RNA"}, snv_threshold=2, indel_threshold=2
    ) == "NoConsensus"


def test_source_evidence_survives_aggregation():
    variants = {
        "chr1:10:A:G": {
            "CHROM": "chr1", "POS": 10, "REF": "A", "ALT": "G",
            "is_snv": True, "caller": "DNA_consensus",
            "filter_original": "Somatic", "filter_normalized": "Somatic",
            "filter_category": "Somatic", "quality": 50.0,
            "genotype": {"GT": "0/1", "DP": 10, "AD": "8,2", "VAF": 0.2},
            "normal_genotype": None, "source_evidence": {"GT_BY_CALLER": "DNA_mutect2:0/1"},
            "id": None, "classification": "Somatic",
        }
    }
    out = aggregate_variants([("DNA_consensus", variants, "DNA")], 1, 1, min_alt_support=0)
    assert out["chr1:10:A:G"]["source_evidence"]["GT_BY_CALLER"] == "DNA_mutect2:0/1"
