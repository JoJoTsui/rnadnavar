"""Regression tests for annotation-stage germline rule ordering (audit M6, ticket 13).

Covers the narrowed M6 fix in ``bin/vcf_utils/variant_classifier.py``:
- Common population frequency (gnomAD AF > ``annotation_germline_freq_threshold``,
  default 0.001) vetoes a Rule-1 Somatic (re)classification, so >=2 DNA-caller
  "Somatic" agreement can no longer override common-AF evidence.
- Germline demotion no longer requires a Germline-labeled caller in BOTH DNA and
  RNA: DNA-side germline evidence at common AF suffices (RNA Mutect2 is sparse,
  so the old both-modalities requirement was nearly unreachable).

Run with: .venv/bin/python -m pytest tests/vcf_utils/test_annotation_germline_ordering.py -v
"""

import sys
from pathlib import Path

# Ensure bin/ is on the Python path
_BIN_DIR = Path(__file__).resolve().parent.parent.parent / "bin"
if str(_BIN_DIR) not in sys.path:
    sys.path.insert(0, str(_BIN_DIR))

from vcf_utils.variant_classifier import create_variant_classifier


def _classifier():
    """Classifier with default thresholds (germline AF threshold 0.001)."""
    return create_variant_classifier()


def _info(filters_normalized, gnomad_af, n_dna, n_rna, original_filter):
    info = {
        "FILTER": original_filter,
        "CHROM": "chr1",
        "POS": 1000,
        "REF": "A",
        "ALT": "G",
        "N_DNA_CALLERS": n_dna,
        "N_RNA_CALLERS": n_rna,
        "FILTERS_NORMALIZED": filters_normalized,
    }
    if gnomad_af is not None:
        info["GNOMAD_AF"] = gnomad_af
    return info


# ---------------------------------------------------------------------------
# M6(a): common-AF evidence vetoes Rule-1 Somatic reclassification
# ---------------------------------------------------------------------------


class TestCommonAfVetoesRule1:
    def test_two_dna_somatic_labels_at_common_af_demote_to_germline(self):
        """Two DNA callers label Somatic a variant at gnomAD AF 0.3 (with a
        third DNA caller labeling Germline): demoted to Germline.
        Pre-fix Rule 1 fired first (2/3 > 50% Somatic) -> stayed Somatic."""
        result = _classifier().classify_variant_from_info(
            _info(
                "DNA_mutect2:Somatic|DNA_strelka:Somatic|DNA_deepsomatic:Germline",
                gnomad_af=0.3,
                n_dna=3,
                n_rna=0,
                original_filter="Somatic",
            )
        )
        assert result.classification == "Germline"

    def test_common_af_without_germline_labels_blocks_rule1(self):
        """Common AF with no Germline-labeled caller at all: Rule 1 must not
        fire, so the original FILTER is preserved (no Somatic rescue).
        Pre-fix Rule 1 reclassified to Somatic."""
        result = _classifier().classify_variant_from_info(
            _info(
                "DNA_mutect2:Somatic|DNA_strelka:Somatic",
                gnomad_af=0.3,
                n_dna=2,
                n_rna=0,
                original_filter="NoConsensus",
            )
        )
        assert result.classification == "NoConsensus"

    def test_subthreshold_af_keeps_rule1_somatic(self):
        """AF below the germline threshold (0.0005 < 0.001): Rule 1 still
        classifies Somatic on majority DNA Somatic-labeled support."""
        result = _classifier().classify_variant_from_info(
            _info(
                "DNA_mutect2:Somatic|DNA_strelka:Somatic|DNA_deepsomatic:Germline",
                gnomad_af=0.0005,
                n_dna=3,
                n_rna=0,
                original_filter="Somatic",
            )
        )
        assert result.classification == "Somatic"

    def test_no_af_annotation_keeps_rule1_somatic(self):
        """No gnomAD annotation at all: Rule 1 behavior is unchanged."""
        result = _classifier().classify_variant_from_info(
            _info(
                "DNA_mutect2:Somatic|DNA_strelka:Somatic",
                gnomad_af=None,
                n_dna=2,
                n_rna=0,
                original_filter="NoConsensus",
            )
        )
        assert result.classification == "Somatic"


# ---------------------------------------------------------------------------
# M6(b): germline demotion must not require RNA germline evidence
# ---------------------------------------------------------------------------


class TestDnaOnlyGermlineDemotion:
    def test_dna_only_germline_evidence_at_common_af_demotes(self):
        """DNA-side Germline-labeled caller at common AF, no RNA caller labels
        at all (sparse RNA Mutect2): demoted to Germline.
        Pre-fix the rule required RNA germline evidence too -> unreachable."""
        result = _classifier().classify_variant_from_info(
            _info(
                "DNA_mutect2:Germline",
                gnomad_af=0.5,
                n_dna=1,
                n_rna=0,
                original_filter="NoConsensus",
            )
        )
        assert result.classification == "Germline"

    def test_common_af_without_any_germline_label_preserves_original(self):
        """Common AF alone, without any Germline-labeled DNA caller, is not
        enough to demote on its own (it only vetoes Rule 1)."""
        result = _classifier().classify_variant_from_info(
            _info(
                "DNA_mutect2:Reference",
                gnomad_af=0.5,
                n_dna=1,
                n_rna=0,
                original_filter="Reference",
            )
        )
        assert result.classification == "Reference"

    def test_artifact_protection_still_blocks_germline_demotion(self):
        """Artifact-labeled callers in BOTH modalities still protect the
        original FILTER despite common AF (existing protection preserved)."""
        result = _classifier().classify_variant_from_info(
            _info(
                "DNA_mutect2:Artifact|RNA_mutect2:Artifact",
                gnomad_af=0.5,
                n_dna=1,
                n_rna=1,
                original_filter="Artifact",
            )
        )
        assert result.classification == "Artifact"
