"""Regression tests for the rescue contract redesign (audit M1/M2/M3, ticket 06).

Covers:
- M1 promotion: when neither modality has a consensus label, DNA and RNA
  individual callers agreeing on Somatic rescue the site as Somatic and tag
  it (RESCUE_PROMOTED / RESCUED) instead of collapsing to NoConsensus.
  A consensus label of "NoConsensus" means "no consensus" and must not block
  promotion.
- M2 veto: a DNA Artifact consensus label outranks RNA non-Artifact evidence;
  RNA can no longer flip a DNA-flagged artifact to Somatic. Direction is
  configurable (rescue_veto_direction, default "dna").
- M3 truthful flags: RESCUED / CROSS_MODALITY / PASSES_CONSENSUS_DNA /
  PASSES_CONSENSUS_RNA are computed from records that PASSED as Somatic,
  never from mere presence in the union file.
- Config: rescue_promotion_enabled / rescue_promotion_min_dna_callers /
  rescue_promotion_min_rna_callers / rescue_veto_direction knobs and their
  CLI counterparts on run_rescue_vcf.py.

Run with: .venv/bin/python -m pytest tests/vcf_utils/test_rescue_contract.py -v
"""

import subprocess
import sys
from pathlib import Path

import pytest
from cyvcf2 import VCF

# Ensure bin/ is on the Python path
_BIN_DIR = Path(__file__).resolve().parent.parent.parent / "bin"
if str(_BIN_DIR) not in sys.path:
    sys.path.insert(0, str(_BIN_DIR))

from vcf_utils.aggregation import aggregate_genotypes
from vcf_utils.io_utils import write_union_vcf
from vcf_utils.tagging import mark_rescued_variants
from vcf_utils.variant_classifier_unified import UnifiedVariantClassifier

BIN_DIR = Path(__file__).resolve().parents[2] / "bin"

# ---------------------------------------------------------------------------
# Synthetic VCF fixtures
# ---------------------------------------------------------------------------

# Consensus-style VCFs (FILTER holds the biological class). The "failed"
# variant failed within-modality consensus; the "somatic" variant passed.
DNA_CONSENSUS_VCF = """\
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=Somatic,Description="Somatic variant">
##FILTER=<ID=Artifact,Description="Artifact">
##FILTER=<ID=NoConsensus,Description="No consensus">
##contig=<ID=chr1,length=1000000>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t1000\t.\tA\tG\t.\tNoConsensus\t.
chr1\t2000\t.\tC\tT\t.\tSomatic\t.
chr1\t3000\t.\tG\tA\t.\tArtifact\t.
"""

RNA_CONSENSUS_VCF = """\
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=Somatic,Description="Somatic variant">
##FILTER=<ID=Artifact,Description="Artifact">
##FILTER=<ID=NoConsensus,Description="No consensus">
##contig=<ID=chr1,length=1000000>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t1000\t.\tA\tG\t.\tNoConsensus\t.
chr1\t2000\t.\tC\tT\t.\tSomatic\t.
chr1\t3000\t.\tG\tA\t.\tNoConsensus\t.
"""

# Mutect2 (normal CRAM first -> sample 0 = normal), PASS with tumor alt 40.
MUTECT2_PASS_VCF = """\
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=chr1,length=1000000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth">
##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele fractions of alternate alleles">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tCOO8801DN\tCOO8801DT
chr1\t1000\t.\tA\tG\t.\tPASS\t.\tGT:AD:DP:AF\t0/0:80,0:80:0.0\t0/1:60,40:100:0.4
"""

# Strelka (fixed NORMAL,TUMOR order), PASS with tumor alt (G tier1) = 40.
STRELKA_PASS_VCF = """\
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=chr1,length=1000000>
##INFO=<ID=NT,Number=1,Type=String,Description="Genotype of the normal sample">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=AU,Number=2,Type=Integer,Description="A allele counts (tier1,tier2)">
##FORMAT=<ID=CU,Number=2,Type=Integer,Description="C allele counts (tier1,tier2)">
##FORMAT=<ID=GU,Number=2,Type=Integer,Description="G allele counts (tier1,tier2)">
##FORMAT=<ID=TU,Number=2,Type=Integer,Description="T allele counts (tier1,tier2)">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\tTUMOR
chr1\t1000\t.\tA\tG\t.\tPASS\tNT=ref\tGT:DP:AU:CU:GU:TU\t0/0:80:80,0:0,0:0,0:0,0\t0/1:100:60,0:0,0:40,0:0,0
"""

# Consensus VCF pair for the veto scenario: DNA Artifact vs RNA Somatic.
DNA_CONSENSUS_ARTIFACT_VCF = """\
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=Somatic,Description="Somatic variant">
##FILTER=<ID=Artifact,Description="Artifact">
##contig=<ID=chr1,length=1000000>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t1000\t.\tA\tG\t.\tArtifact\t.
"""

RNA_CONSENSUS_SOMATIC_VCF = DNA_CONSENSUS_ARTIFACT_VCF.replace("Artifact\t", "Somatic\t")


def _write_vcf(directory, name, content):
    path = directory / name
    path.write_text(content)
    return str(path)


def _run_rescue(tmp_path, dna_consensus, rna_consensus, extra_args=None, name="rescued"):
    dna = _write_vcf(tmp_path, f"{name}.dna.consensus.vcf", dna_consensus)
    rna = _write_vcf(tmp_path, f"{name}.rna.consensus.vcf", rna_consensus)
    out_prefix = tmp_path / name
    cmd = [
        sys.executable,
        str(BIN_DIR / "run_rescue_vcf.py"),
        "--dna_consensus",
        dna,
        "--rna_consensus",
        rna,
        "--out_prefix",
        str(out_prefix),
        "--output_format",
        "vcf",
    ] + (extra_args or [])
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    return str(out_prefix) + ".vcf"


def _info_scalar(variant, key):
    val = variant.INFO.get(key)
    if isinstance(val, (tuple, list)):
        val = val[0] if val else None
    return val


def _classify(callers, filters, config=None, is_snv=True):
    """Unit seam: run the rescue classifier on a synthetic variant_data dict."""
    modality_map = {}
    for caller in callers:
        if caller.startswith("DNA_"):
            modality_map[caller] = "DNA"
        elif caller.startswith("RNA_"):
            modality_map[caller] = "RNA"
    data = {
        "callers": list(callers),
        "filters_normalized": list(filters),
        "is_snv": is_snv,
    }
    result = UnifiedVariantClassifier(config).classify_rescue_variant(
        data, modality_map
    )
    return result, data


# ---------------------------------------------------------------------------
# M1: cross-modality promotion
# ---------------------------------------------------------------------------


class TestCrossModalityPromotion:
    def test_agreeing_somatic_callers_promote(self):
        """1 DNA + 1 RNA individual caller agreeing on Somatic, no consensus
        labels -> Somatic. Pre-fix this returned NoConsensus (audit M1)."""
        result, _ = _classify(
            ["DNA_mutect2", "RNA_strelka"], ["Somatic", "Somatic"],
            config={"rescue_promotion_min_rna_callers": 1},
        )
        assert result == "Somatic"

    def test_promotion_marks_record_rescued(self):
        """A promoted record is tagged as cross-modality rescued."""
        _, data = _classify(
            ["DNA_mutect2", "RNA_strelka"], ["Somatic", "Somatic"],
            config={"rescue_promotion_min_rna_callers": 1},
        )
        assert data.get("rescued") is True
        assert data.get("rescue_promoted") is True

    def test_noconsensus_consensus_labels_do_not_block_promotion(self):
        """A 'NoConsensus' consensus label means 'no consensus': it must be
        treated as absent so individual callers can still promote the site.
        Pre-fix the both-labels-agree branch returned NoConsensus."""
        result, _ = _classify(
            ["DNA_consensus", "RNA_consensus", "DNA_mutect2", "RNA_strelka"],
            ["NoConsensus", "NoConsensus", "Somatic", "Somatic"],
            config={"rescue_promotion_min_rna_callers": 1},
        )
        assert result == "Somatic"

    def test_indel_is_not_promoted_by_gated_rule(self):
        result, data = _classify(
            ["DNA_mutect2", "RNA_strelka"], ["Somatic", "Somatic"],
            config={"rescue_promotion_min_rna_callers": 1}, is_snv=False,
        )
        assert result == "NoConsensus"
        assert data.get("rescue_promoted") is not True

    def test_germline_agreement_not_promoted(self):
        """Only Somatic agreement promotes; Germline agreement keeps the
        legacy NoConsensus outcome."""
        result, _ = _classify(
            ["DNA_mutect2", "RNA_strelka"], ["Germline", "Germline"]
        )
        assert result == "NoConsensus"

    def test_cross_modality_disagreement_stays_artifact(self):
        """DNA Somatic vs RNA Germline disagreement -> Artifact (unchanged)."""
        result, _ = _classify(
            ["DNA_mutect2", "RNA_strelka"], ["Somatic", "Germline"]
        )
        assert result == "Artifact"

    def test_promotion_requires_both_modalities(self):
        result, _ = _classify(["DNA_mutect2", "DNA_strelka"], ["Somatic", "Somatic"])
        assert result == "NoConsensus"

    def test_promotion_min_caller_thresholds_configurable(self):
        """rescue_promotion_min_dna_callers / rescue_promotion_min_rna_callers
        gate the promotion."""
        # 1 DNA caller is below a min of 2 -> no promotion
        result, _ = _classify(
            ["DNA_mutect2", "RNA_strelka"],
            ["Somatic", "Somatic"],
            config={"rescue_promotion_min_dna_callers": 2},
        )
        assert result == "NoConsensus"
        # 2 DNA callers meet the raised floor
        result, _ = _classify(
            ["DNA_mutect2", "DNA_deepsomatic", "RNA_strelka"],
            ["Somatic", "Somatic", "Somatic"],
            config={"rescue_promotion_min_dna_callers": 2, "rescue_promotion_min_rna_callers": 1},
        )
        assert result == "Somatic"

    def test_promotion_param_gated_legacy_behavior(self):
        """rescue_promotion_enabled=False restores the pre-fix NoConsensus."""
        result, data = _classify(
            ["DNA_mutect2", "RNA_strelka"],
            ["Somatic", "Somatic"],
            config={"rescue_promotion_enabled": False},
        )
        assert result == "NoConsensus"
        assert data.get("rescue_promoted") is not True


def test_default_gated_policy_requires_two_rna_callers():
    result, _ = _classify(["DNA_mutect2", "RNA_strelka"], ["Somatic", "Somatic"])
    assert result == "NoConsensus"


# ---------------------------------------------------------------------------
# M2: DNA-Artifact veto
# ---------------------------------------------------------------------------


class TestDnaArtifactVeto:
    VETO_CALLERS = ["DNA_consensus", "RNA_consensus", "RNA_mutect2", "RNA_strelka"]
    VETO_FILTERS = ["Artifact", "Somatic", "Somatic", "Somatic"]

    def test_dna_artifact_vetoes_rna_somatic(self):
        """DNA Artifact consensus + RNA Somatic consensus with >=2 RNA
        callers -> stays Artifact. Pre-fix the RNA-first branch flipped it
        to Somatic (audit M2)."""
        result, _ = _classify(self.VETO_CALLERS, self.VETO_FILTERS)
        assert result == "Artifact"

    def test_veto_direction_none_restores_legacy(self):
        """rescue_veto_direction='none' restores the legacy RNA override."""
        result, _ = _classify(
            self.VETO_CALLERS,
            self.VETO_FILTERS,
            config={"rescue_veto_direction": "none"},
        )
        assert result == "Somatic"

    def test_veto_direction_rna(self):
        """rescue_veto_direction='rna' mirrors the veto: RNA Artifact outranks
        DNA Somatic even with >=2 DNA callers."""
        result, _ = _classify(
            ["DNA_consensus", "RNA_consensus", "DNA_mutect2", "DNA_strelka"],
            ["Somatic", "Artifact", "Somatic", "Somatic"],
            config={"rescue_veto_direction": "rna"},
        )
        assert result == "Artifact"

    def test_rna_artifact_dna_somatic_keeps_legacy(self):
        """RNA Artifact vs DNA Somatic keeps current behavior: the DNA
        non-Artifact label still wins with >=2 DNA callers under the default
        'dna' veto direction."""
        result, _ = _classify(
            ["DNA_consensus", "RNA_consensus", "DNA_mutect2", "DNA_strelka"],
            ["Somatic", "Artifact", "Somatic", "Somatic"],
        )
        assert result == "Somatic"


# ---------------------------------------------------------------------------
# M3: truthful RESCUED / CROSS_MODALITY / PASSES_CONSENSUS_* flags
# ---------------------------------------------------------------------------


def _consensus_record(filter_value):
    """Minimal per-record dict as produced by read_variants_from_vcf for a
    consensus VCF record (filter_normalized holds the FILTER class)."""
    return {
        "CHROM": "chr1",
        "POS": 1000,
        "REF": "A",
        "ALT": "G",
        "is_snv": True,
        "caller": "DNA_consensus",
        "filter_original": filter_value,
        "filter_normalized": filter_value,
        "filter_category": filter_value,
        "classification": filter_value,
        "quality": None,
        "genotype": {"GT": None, "DP": None, "AD": None, "VAF": None, "GQ": None},
        "id": None,
    }


class TestTruthfulRescueFlags:
    def test_mark_rescued_requires_somatic_pass(self):
        """A site present in both union files but Artifact in DNA and
        NoConsensus in RNA is NOT rescued. Pre-fix mere presence gave
        rescued=True (audit M3)."""
        vkey = "1:1000:A:G"
        dna = {vkey: _consensus_record("Artifact")}
        rna = {vkey: _consensus_record("NoConsensus")}
        aggregated = {vkey: {}}
        marked = mark_rescued_variants(aggregated, dna, rna)
        assert marked[vkey]["rescued"] is False
        assert marked[vkey]["cross_modality"] is False

    def test_mark_rescued_somatic_in_both(self):
        vkey = "1:1000:A:G"
        dna = {vkey: _consensus_record("Somatic")}
        rna = {vkey: _consensus_record("Somatic")}
        aggregated = {vkey: {}}
        marked = mark_rescued_variants(aggregated, dna, rna)
        assert marked[vkey]["rescued"] is True
        assert marked[vkey]["cross_modality"] is True

    def test_mark_rescued_legacy_key_sets_still_work(self):
        """Plain key sets (no record evidence) keep legacy presence semantics."""
        vkey = "1:1000:A:G"
        marked = mark_rescued_variants({vkey: {}}, {vkey}, {vkey})
        assert marked[vkey]["rescued"] is True


def _write_record(tmp_path, callers, filters, rescued=False, name="out.vcf"):
    """Write a single aggregated record via write_union_vcf (rescue mode)."""
    data = {
        "CHROM": "chr1",
        "POS": 1000,
        "REF": "A",
        "ALT": "G",
        "is_snv": True,
        "callers": list(callers),
        "filters_original": list(filters),
        "filters_normalized": list(filters),
        "filters_category": list(filters),
        "qualities": [],
        "genotypes": {},
        "ids": [],
        "support_callers": set(callers),
        "passes_consensus": True,
        "rescued": rescued,
        "gt_aggregated": aggregate_genotypes({}, callers),
    }
    modality_map = {}
    for c in callers:
        if c.startswith("DNA_"):
            modality_map[c] = "DNA"
        elif c.startswith("RNA_"):
            modality_map[c] = "RNA"
    template_path = _write_vcf(tmp_path, "template.vcf", DNA_CONSENSUS_ARTIFACT_VCF)
    out = tmp_path / name
    write_union_vcf(
        {"1:1000:A:G": data},
        VCF(template_path),
        "SAMPLE",
        str(out),
        "vcf",
        list(callers),
        modality_map=modality_map,
    )
    (rec,) = list(VCF(str(out)))
    return rec


class TestTruthfulInfoFlags:
    def test_flags_no_for_failed_consensus_records(self, tmp_path):
        """Record present in both modalities' union files but NoConsensus in
        DNA and Artifact in RNA -> all passed-flags NO. Pre-fix these were
        YES from mere presence (audit M3)."""
        rec = _write_record(
            tmp_path,
            callers=["DNA_consensus", "RNA_consensus"],
            filters=["NoConsensus", "Artifact"],
        )
        assert _info_scalar(rec, "PASSES_CONSENSUS_DNA") == "NO"
        assert _info_scalar(rec, "PASSES_CONSENSUS_RNA") == "NO"
        assert _info_scalar(rec, "CROSS_MODALITY") == "NO"
        assert _info_scalar(rec, "RESCUED") == "NO"

    def test_flags_yes_for_somatic_in_both(self, tmp_path):
        """Genuinely Somatic in both modalities' consensus -> all flags YES."""
        rec = _write_record(
            tmp_path,
            callers=["DNA_consensus", "RNA_consensus"],
            filters=["Somatic", "Somatic"],
            rescued=True,
        )
        assert _info_scalar(rec, "PASSES_CONSENSUS_DNA") == "YES"
        assert _info_scalar(rec, "PASSES_CONSENSUS_RNA") == "YES"
        assert _info_scalar(rec, "CROSS_MODALITY") == "YES"
        assert _info_scalar(rec, "RESCUED") == "YES"


# ---------------------------------------------------------------------------
# End-to-end CLI tests (run_rescue_vcf.py)
# ---------------------------------------------------------------------------


class TestRescueCliContract:
    def test_e2e_promotion_rescues_failed_consensus_site(self, tmp_path):
        """Site failed within-modality consensus (NoConsensus in both) but
        1 DNA + 1 RNA caller call it Somatic -> rescued as Somatic and tagged.
        Pre-fix the output FILTER was NoConsensus."""
        mutect2 = _write_vcf(tmp_path, "sample.mutect2.variants.vcf", MUTECT2_PASS_VCF)
        strelka = _write_vcf(tmp_path, "sample.strelka.variants.vcf", STRELKA_PASS_VCF)
        out = _run_rescue(
            tmp_path,
            DNA_CONSENSUS_VCF,
            RNA_CONSENSUS_VCF,
            extra_args=["--dna_vcf", mutect2, "--rna_vcf", strelka,
                        "--rescue_min_rna_callers", "1"],
        )
        records = {r.POS: r for r in VCF(out)}
        rec = records[1000]
        assert rec.FILTER == "Somatic"
        assert _info_scalar(rec, "RESCUE_PROMOTED") == "YES"
        assert _info_scalar(rec, "RESCUED") == "YES"
        assert _info_scalar(rec, "CROSS_MODALITY") == "YES"

    def test_e2e_promotion_disabled_by_cli(self, tmp_path):
        """--disable_rescue_promotion restores the legacy NoConsensus."""
        mutect2 = _write_vcf(tmp_path, "sample.mutect2.variants.vcf", MUTECT2_PASS_VCF)
        strelka = _write_vcf(tmp_path, "sample.strelka.variants.vcf", STRELKA_PASS_VCF)
        out = _run_rescue(
            tmp_path,
            DNA_CONSENSUS_VCF,
            RNA_CONSENSUS_VCF,
            extra_args=[
                "--dna_vcf",
                mutect2,
                "--rna_vcf",
                strelka,
                "--disable_rescue_promotion",
            ],
        )
        records = {r.POS: r for r in VCF(out)}
        assert records[1000].FILTER == "NoConsensus"

    def test_e2e_dna_artifact_veto(self, tmp_path):
        """DNA Artifact consensus + RNA Somatic consensus with 2 RNA callers
        -> stays Artifact. Pre-fix flipped to Somatic."""
        mutect2 = _write_vcf(tmp_path, "sample.mutect2.variants.vcf", MUTECT2_PASS_VCF)
        strelka = _write_vcf(tmp_path, "sample.strelka.variants.vcf", STRELKA_PASS_VCF)
        out = _run_rescue(
            tmp_path,
            DNA_CONSENSUS_ARTIFACT_VCF,
            RNA_CONSENSUS_SOMATIC_VCF,
            extra_args=["--rna_vcf", mutect2, "--rna_vcf", strelka],
        )
        (rec,) = list(VCF(out))
        assert rec.FILTER == "Artifact"

    def test_e2e_veto_none_restores_legacy(self, tmp_path):
        """--rescue_veto none restores the legacy RNA override."""
        mutect2 = _write_vcf(tmp_path, "sample.mutect2.variants.vcf", MUTECT2_PASS_VCF)
        strelka = _write_vcf(tmp_path, "sample.strelka.variants.vcf", STRELKA_PASS_VCF)
        out = _run_rescue(
            tmp_path,
            DNA_CONSENSUS_ARTIFACT_VCF,
            RNA_CONSENSUS_SOMATIC_VCF,
            extra_args=[
                "--rna_vcf",
                mutect2,
                "--rna_vcf",
                strelka,
                "--rescue_veto",
                "none",
            ],
        )
        (rec,) = list(VCF(out))
        assert rec.FILTER == "Somatic"

    def test_e2e_flags_never_fired_by_mere_presence(self, tmp_path):
        """chr1:3000 is Artifact in DNA consensus and NoConsensus in RNA
        consensus -> RESCUED / PASSES_CONSENSUS_* / CROSS_MODALITY all NO.
        chr1:2000 is Somatic in both -> all YES."""
        out = _run_rescue(tmp_path, DNA_CONSENSUS_VCF, RNA_CONSENSUS_VCF)
        records = {r.POS: r for r in VCF(out)}

        failed = records[3000]
        assert _info_scalar(failed, "RESCUED") == "NO"
        assert _info_scalar(failed, "PASSES_CONSENSUS_DNA") == "NO"
        assert _info_scalar(failed, "PASSES_CONSENSUS_RNA") == "NO"
        assert _info_scalar(failed, "CROSS_MODALITY") == "NO"

        somatic = records[2000]
        assert _info_scalar(somatic, "RESCUED") == "YES"
        assert _info_scalar(somatic, "PASSES_CONSENSUS_DNA") == "YES"
        assert _info_scalar(somatic, "PASSES_CONSENSUS_RNA") == "YES"
        assert _info_scalar(somatic, "CROSS_MODALITY") == "YES"
