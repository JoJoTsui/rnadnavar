"""Regression tests for RNA-editing over-masking (audit M7, ticket 13).

Covers the M7 fix in ``bin/vcf_utils/rna_editing_core.py`` +
``bin/vcf_utils/filter_updater.py`` (+ ``evidence_tiering.py`` consistency):

- FILTER is changed to ``RNAedit`` only for tiers with NO DNA support
  (VERY_HIGH/HIGH in the four-tier matrix of ``rna_editing_core.py``; HIGH in
  ``evidence_tiering.py``'s HIGH/MEDIUM/LOW scheme, where HIGH == RNA-only).
- MEDIUM (canonical transition, REDIportal match, RNA consensus, but HAS DNA
  presence) and LOW (non-canonical transition at a REDIportal site) keep their
  original FILTER; the tier is still emitted as the ``REDI_EVIDENCE`` INFO
  annotation by ``annotate_rna_editing.py`` (existing REDI_* INFO convention),
  so the editing evidence remains visible without masking DNA-supported
  somatic calls.

Run with: .venv/bin/python -m pytest tests/vcf_utils/test_rna_editing_overmasking.py -v
"""

import sys
from pathlib import Path

# Ensure bin/ is on the Python path
_BIN_DIR = Path(__file__).resolve().parent.parent.parent / "bin"
if str(_BIN_DIR) not in sys.path:
    sys.path.insert(0, str(_BIN_DIR))

from vcf_utils.evidence_tiering import create_evidence_tiering_processor
from vcf_utils.filter_updater import FilterUpdateValidator, create_filter_updater
from vcf_utils.rna_editing_core import (
    classify_rna_editing_biological_category,
    classify_rna_editing_evidence,
)

# A DNA-supported canonical A>G at a REDIportal site: tumor DNA VAF 0.3,
# two RNA callers support it, two DNA callers support it.
DNA_SUPPORTED_VARIANT = {
    "CHROM": "chr1",
    "POS": 1000,
    "REF": "A",
    "ALT": "G",
    "N_RNA_CALLERS_SUPPORT": 2,
    "N_DNA_CALLERS_SUPPORT": 2,
    "VAF_RNA_MEAN": 0.4,
    "VAF_DNA_MEAN": 0.3,
}

# A canonical A>G at a REDIportal site with NO DNA support (RNA-only).
RNA_ONLY_VARIANT = {
    "CHROM": "chr1",
    "POS": 2000,
    "REF": "A",
    "ALT": "G",
    "N_RNA_CALLERS_SUPPORT": 3,
    "N_DNA_CALLERS_SUPPORT": 0,
    "VAF_RNA_MEAN": 0.15,
    "VAF_DNA_MEAN": 0.0,
}

# A non-canonical transition (G>A) at a REDIportal site with RNA consensus.
NON_CANONICAL_VARIANT = {
    "CHROM": "chr1",
    "POS": 3000,
    "REF": "G",
    "ALT": "A",
    "N_RNA_CALLERS_SUPPORT": 2,
    "N_DNA_CALLERS_SUPPORT": 1,
    "VAF_RNA_MEAN": 0.2,
    "VAF_DNA_MEAN": 0.1,
}


# ---------------------------------------------------------------------------
# Tier semantics (four-tier matrix in rna_editing_core.py)
# ---------------------------------------------------------------------------


class TestEvidenceTierSemantics:
    def test_dna_supported_canonical_site_is_medium_tier(self):
        """Canonical + REDIportal + RNA consensus + DNA presence -> MEDIUM."""
        tier = classify_rna_editing_evidence(DNA_SUPPORTED_VARIANT, True, 2)
        assert tier == "MEDIUM"

    def test_rna_only_canonical_site_is_no_dna_support_tier(self):
        """Canonical + REDIportal + RNA consensus, RNA-only, high RNA VAF,
        low DNA VAF -> VERY_HIGH (no DNA support)."""
        tier = classify_rna_editing_evidence(RNA_ONLY_VARIANT, True, 2)
        assert tier == "VERY_HIGH"

    def test_non_canonical_site_is_low_tier(self):
        """Non-canonical transition at a REDIportal site -> LOW."""
        tier = classify_rna_editing_evidence(NON_CANONICAL_VARIANT, True, 2)
        assert tier == "LOW"


# ---------------------------------------------------------------------------
# M7: FILTER change to RNAedit only for tiers with no DNA support
# ---------------------------------------------------------------------------


class TestBiologicalCategory:
    def test_dna_supported_medium_tier_keeps_somatic_filter(self):
        """DNA-supported A>G at a REDIportal site (tumor DNA VAF 0.3) keeps
        FILTER=Somatic. Pre-fix it was relabeled RNAedit and removed from the
        somatic set."""
        tier = classify_rna_editing_evidence(DNA_SUPPORTED_VARIANT, True, 2)
        assert tier == "MEDIUM"
        category = classify_rna_editing_biological_category(tier, "Somatic")
        assert category == "Somatic"

    def test_low_tier_keeps_original_filter(self):
        """Non-canonical transition at a REDIportal site keeps its FILTER."""
        tier = classify_rna_editing_evidence(NON_CANONICAL_VARIANT, True, 2)
        assert tier == "LOW"
        category = classify_rna_editing_biological_category(tier, "Somatic")
        assert category == "Somatic"

    def test_no_dna_support_tiers_still_become_rnaedit(self):
        """VERY_HIGH/HIGH (no DNA support) are still relabeled RNAedit."""
        for tier in ("VERY_HIGH", "HIGH"):
            assert classify_rna_editing_biological_category(tier, "Somatic") == "RNAedit"

    def test_none_tier_preserves_original(self):
        assert classify_rna_editing_biological_category("NONE", "Somatic") == "Somatic"


class TestFilterUpdater:
    def test_medium_tier_preserves_somatic_filter(self):
        """End-to-end FILTER decision for the DNA-supported case:
        MEDIUM tier + REDIportal match + RNA consensus must NOT update FILTER.
        Pre-fix this returned RNAedit."""
        updater = create_filter_updater()
        new_filter, was_updated = updater.update_variant_filter(
            original_filter="Somatic",
            evidence_tier="MEDIUM",
            has_rediportal_match=True,
            rna_support=2,
            min_rna_support=2,
        )
        assert new_filter == "Somatic"
        assert was_updated is False

    def test_low_tier_preserves_filter(self):
        updater = create_filter_updater()
        new_filter, was_updated = updater.update_variant_filter(
            original_filter="Somatic",
            evidence_tier="LOW",
            has_rediportal_match=True,
            rna_support=2,
            min_rna_support=2,
        )
        assert new_filter == "Somatic"
        assert was_updated is False

    def test_high_tier_still_updates_to_rnaedit(self):
        """HIGH tier (RNA-only, no DNA support) still becomes RNAedit."""
        updater = create_filter_updater()
        new_filter, was_updated = updater.update_variant_filter(
            original_filter="Somatic",
            evidence_tier="HIGH",
            has_rediportal_match=True,
            rna_support=3,
            min_rna_support=2,
        )
        assert new_filter == "RNAedit"
        assert was_updated is True

    def test_very_high_tier_still_updates_to_rnaedit(self):
        updater = create_filter_updater()
        new_filter, was_updated = updater.update_variant_filter(
            original_filter="Somatic",
            evidence_tier="VERY_HIGH",
            has_rediportal_match=True,
            rna_support=3,
            min_rna_support=2,
        )
        assert new_filter == "RNAedit"
        assert was_updated is True

    def test_validator_matches_new_rules(self):
        """FilterUpdateValidator encodes the same no-DNA-support-only rule."""
        assert FilterUpdateValidator.validate_filter_update_logic(
            "Somatic", "Somatic", "MEDIUM", True, 2, 2
        )
        assert FilterUpdateValidator.validate_filter_update_logic(
            "Somatic", "RNAedit", "HIGH", True, 3, 2
        )
        assert not FilterUpdateValidator.validate_filter_update_logic(
            "Somatic", "RNAedit", "MEDIUM", True, 2, 2
        )


class TestEvidenceTieringConsistency:
    def test_rna_only_tier_marks_filter_update(self):
        """evidence_tiering.py HIGH (RNA-only) still flags a FILTER update."""
        processor = create_evidence_tiering_processor(min_rna_support=2)
        result = processor.process_variant(
            {"N_RNA_CALLERS_SUPPORT": 3, "N_DNA_CALLERS_SUPPORT": 0}, True
        )
        assert result["evidence_tier"] == "HIGH"
        assert result["update_filter"] is True

    def test_dna_supported_tier_does_not_mark_filter_update(self):
        """evidence_tiering.py MEDIUM (has DNA support) no longer flags a
        FILTER update; the tier is still reported for the REDI_EVIDENCE INFO
        annotation written by annotate_rna_editing.py."""
        processor = create_evidence_tiering_processor(min_rna_support=2)
        result = processor.process_variant(
            {"N_RNA_CALLERS_SUPPORT": 2, "N_DNA_CALLERS_SUPPORT": 2}, True
        )
        assert result["evidence_tier"] == "MEDIUM"
        assert result["update_filter"] is False
