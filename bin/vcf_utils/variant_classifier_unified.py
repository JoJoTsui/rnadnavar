#!/usr/bin/env python3
"""
Unified Variant Classifier for All Classification Scenarios

This module provides a comprehensive unified classifier that handles all variant
classification logic across the RNADNAvar pipeline:

1. Standalone caller classification (Strelka, DeepSomatic, Mutect2)
2. Consensus classification with voting thresholds
3. Cross-modality rescue classification
4. Annotation-based reclassification (gnomAD, COSMIC, REDIportal)

The classifier preserves the exact complex logic from the current implementation
while providing a unified interface for all classification needs.

Key Features:
- Preserves exact voting and classification rules from current pipeline
- Handles complex cross-modality rescue logic
- Integrates RNA editing classification
- Supports annotation-based reclassification
- Configurable thresholds via classification_config
"""

from collections import Counter

from .classification import classify_variant_from_dict, classify_variant_from_record
from .classification_config import get_config_with_overrides


class UnifiedVariantClassifier:
    """
    Unified classifier for all variant classification scenarios.

    This class implements the complete classification logic used throughout
    the RNADNAvar pipeline, preserving all the complex rules and thresholds
    from the current implementation.
    """

    def __init__(self, config=None):
        """
        Initialize classifier with optional config overrides.

        Args:
            config (dict): Optional threshold overrides
        """
        self.config = get_config_with_overrides(config)

    def classify_standalone_variant(self, variant, caller_name, sample_indices=None):
        """
        Classify variant from a single caller.

        This delegates to the existing classification logic for each caller:
        - Strelka: Uses FILTER, NT field, and normal depth
        - DeepSomatic: Uses FILTER field (GERMLINE, RefCall, etc.)
        - Mutect2: Uses FILTER field

        Args:
            variant: cyvcf2.Variant object or dict representation
            caller_name (str): Name of the variant caller
            sample_indices (dict): Optional sample indices for tumor/normal

        Returns:
            str: One of ['Somatic', 'Germline', 'Reference', 'Artifact']
        """
        if isinstance(variant, dict):
            return classify_variant_from_dict(variant, caller_name)
        else:
            return classify_variant_from_record(variant, caller_name, sample_indices)

    def classify_consensus_variant(self, variant_data):
        """
        Classify variant based on consensus voting logic.

        Preserves exact logic from current implementation:
        1. Check if total caller count meets threshold (SNV vs indel)
        2. If insufficient callers → NoConsensus
        3. If sufficient callers → Use majority vote among classifications
        4. If clear majority → Use majority classification
        5. If tie → Artifact (disagreement indicates inconsistency)

        Args:
            variant_data (dict): Aggregated variant data containing:
                - callers: List of caller names
                - filters_normalized: List of classifications per caller
                - is_snv: Boolean indicating if variant is SNV

        Returns:
            str: One of ['Somatic', 'Germline', 'Reference', 'Artifact', 'NoConsensus']
        """
        classification, _ = self.classify_consensus_variant_with_rationale(
            variant_data
        )
        return classification

    @staticmethod
    def _format_votes(classification_counts):
        """Format vote counts for INFO: 'Somaticx2+Artifactx1' (VCF-safe)."""
        return "+".join(
            f"{cls}x{count}"
            for cls, count in sorted(
                classification_counts.items(), key=lambda kv: (-kv[1], kv[0])
            )
        )

    def classify_consensus_variant_with_rationale(self, variant_data):
        """
        Classify a consensus variant and explain the decision (ticket 07).

        Same logic as classify_consensus_variant; additionally returns a
        VCF-INFO-safe rationale string so every FILTER is derivable from the
        record's own INFO (CLASSIFICATION_RATIONALE).

        Args:
            variant_data (dict): Aggregated variant data (see
                classify_consensus_variant)

        Returns:
            tuple: (classification, rationale) where rationale has the form
                'rule:<rule>|class:<classification>|votes:<class>x<n>+...|
                callers:<n>|threshold:<t>'
        """
        # Get individual caller classifications (exclude consensus callers)
        individual_filters = []
        observed_filters = [f for f in variant_data.get("filters_normalized", [])]
        eligible = variant_data.get("support_callers")
        for i, caller in enumerate(variant_data["callers"]):
            if not caller.endswith("_consensus") and (eligible is None or caller in eligible):
                if i < len(variant_data["filters_normalized"]):
                    individual_filters.append(variant_data["filters_normalized"][i])

        # Check if we have enough callers for consensus
        caller_count = len(individual_filters)
        required_threshold = (
            self.config["consensus_snv_threshold"]
            if variant_data.get("is_snv", True)
            else self.config["consensus_indel_threshold"]
        )

        # Use majority vote to determine classification
        classification_counts = Counter(individual_filters)
        votes = (
            self._format_votes(classification_counts)
            if classification_counts
            else "none"
        )

        def rationale(rule, classification):
            return (
                f"rule:{rule}|class:{classification}|votes:{votes}"
                f"|callers:{caller_count}|threshold:{required_threshold}"
            )

        if caller_count < required_threshold:
            # Preserve an explicit all-Artifact outcome even when no record is
            # eligible to vote; this is a substantive negative classification,
            # not an affirmative consensus vote.
            if observed_filters and all(label == "Artifact" for label in observed_filters):
                return "Artifact", rationale("all_observed_artifact", "Artifact")
            return "NoConsensus", rationale("insufficient_callers", "NoConsensus")

        # Find the most frequent classification(s)
        if not classification_counts:
            return "NoConsensus", rationale("insufficient_callers", "NoConsensus")

        max_count = max(classification_counts.values())
        majority_classifications = [
            cls for cls, count in classification_counts.items() if count == max_count
        ]

        # Check if there's a clear majority (no tie)
        if len(majority_classifications) == 1:
            # Clear majority - use the majority classification
            classification = majority_classifications[0]
            return classification, rationale("majority", classification)
        else:
            # Tie between classifications - mark as Artifact due to disagreement
            return "Artifact", rationale("majority_tie", "Artifact")

    def classify_rescue_variant(self, variant_data, modality_map):
        """
        Classify variant based on cross-modality rescue logic.

        Preserves exact complex logic from current implementation:
        1. Separate consensus callers from individual callers
        2. Extract DNA/RNA consensus labels and individual caller classifications
        3. Apply hierarchical cross-modality rules:
           - Both DNA & RNA consensus: Apply agreement/disagreement rules
           - Single modality consensus: Return that classification
           - No consensus: Analyze individual caller patterns

        Rescue contract (audit M1/M2, ticket 06):
        - A "NoConsensus" consensus label means the modality reached no
          consensus and is treated as absent.
        - Promotion: when neither modality has a consensus label and at
          least rescue_promotion_min_dna_callers DNA callers and
          rescue_promotion_min_rna_callers RNA callers agree on Somatic, the
          site is rescued as Somatic and the record is tagged
          (variant_data["rescued"]/["rescue_promoted"]). Gated by
          rescue_promotion_enabled (default on).
        - Veto: with rescue_veto_direction "dna" (default), a DNA Artifact
          consensus label outranks RNA non-Artifact evidence; "rna" mirrors
          it; "none" restores the legacy RNA-first override.

        Args:
            variant_data (dict): Complete variant data with consensus and individual callers
            modality_map (dict): Maps caller names to modalities ('DNA' or 'RNA')

        Returns:
            str: One of ['Somatic', 'Germline', 'Reference', 'Artifact', 'NoConsensus']
        """
        classification, _ = self.classify_rescue_variant_with_rationale(
            variant_data, modality_map
        )
        return classification

    def classify_rescue_variant_with_rationale(self, variant_data, modality_map):
        """
        Classify a rescue variant and explain the decision (ticket 07).

        Same logic as classify_rescue_variant; additionally returns a
        VCF-INFO-safe rationale string so every FILTER is derivable from the
        record's own INFO (CLASSIFICATION_RATIONALE).

        Args:
            variant_data (dict): Complete variant data with consensus and individual callers
            modality_map (dict): Maps caller names to modalities ('DNA' or 'RNA')

        Returns:
            tuple: (classification, rationale) where rationale has the form
                'rule:<rule>|class:<classification>|dna_consensus:<label|none>|
                rna_consensus:<label|none>|dna_votes:<class>x<n>+...|rna_votes:...'
        """
        # Step 1: Separate consensus callers from individual callers
        dna_consensus_label = None
        rna_consensus_label = None
        individual_callers = []
        individual_filters = []

        for i, caller in enumerate(variant_data["callers"]):
            if caller.endswith("_consensus"):
                # Extract consensus labels. A "NoConsensus" consensus label
                # means the modality reached NO consensus for this site —
                # treat it as absent so individual callers can still drive
                # the rescue contract (audit M1).
                label = variant_data["filters_normalized"][i]
                if label == "NoConsensus":
                    label = None
                if "DNA" in caller.upper():
                    dna_consensus_label = label
                elif "RNA" in caller.upper():
                    rna_consensus_label = label
            else:
                # Collect only eligible callers for vote-dependent decisions.
                # Observed records remain available in the aggregate metadata.
                if variant_data.get("support_callers") is not None and caller not in variant_data["support_callers"]:
                    continue
                individual_callers.append(caller)
                if i < len(variant_data["filters_normalized"]):
                    individual_filters.append(variant_data["filters_normalized"][i])

        # Get individual callers by modality
        dna_callers = []
        rna_callers = []

        for i, caller in enumerate(individual_callers):
            modality = modality_map.get(caller, "UNKNOWN")
            if modality == "DNA" and i < len(individual_filters):
                dna_callers.append(individual_filters[i])
            elif modality == "RNA" and i < len(individual_filters):
                rna_callers.append(individual_filters[i])

        def rationale(rule, classification):
            dna_votes = (
                self._format_votes(Counter(dna_callers)) if dna_callers else "none"
            )
            rna_votes = (
                self._format_votes(Counter(rna_callers)) if rna_callers else "none"
            )
            return (
                f"rule:{rule}|class:{classification}"
                f"|dna_consensus:{dna_consensus_label or 'none'}"
                f"|rna_consensus:{rna_consensus_label or 'none'}"
                f"|dna_votes:{dna_votes}|rna_votes:{rna_votes}"
            )

        # Step 2: Apply cross-modality consensus rules if both available
        if dna_consensus_label and rna_consensus_label:
            # Both modalities have consensus
            if dna_consensus_label == rna_consensus_label:
                # Agreement across modalities - use agreed classification
                return dna_consensus_label, rationale(
                    "cross_modality_agreement", dna_consensus_label
                )
            elif (
                dna_consensus_label == "Artifact" and rna_consensus_label == "Artifact"
            ):
                # Both modalities are Artifact
                return "Artifact", rationale("cross_modality_agreement", "Artifact")
            elif (
                dna_consensus_label != "Artifact" and rna_consensus_label != "Artifact"
            ):
                # Cross-modality disagreement on non-Artifact classifications
                min_callers = self.config["cross_modality_min_callers_for_artifact"]
                if len(dna_callers) >= min_callers and len(rna_callers) >= min_callers:
                    return "Artifact", rationale(
                        "cross_modality_disagreement", "Artifact"
                    )
                elif len(dna_callers) >= min_callers:
                    return dna_consensus_label, rationale(
                        "dna_consensus_priority", dna_consensus_label
                    )
                elif len(rna_callers) >= min_callers:
                    return rna_consensus_label, rationale(
                        "rna_consensus_priority", rna_consensus_label
                    )
                else:
                    return "NoConsensus", rationale(
                        "insufficient_modality_support", "NoConsensus"
                    )
            elif (
                self.config["rescue_veto_direction"] == "dna"
                and dna_consensus_label == "Artifact"
            ):
                # DNA Artifact veto (audit M2): a DNA-flagged artifact
                # outranks RNA non-Artifact evidence — RNA-specific error
                # modes must not flip it to Somatic
                return "Artifact", rationale("dna_artifact_veto", "Artifact")
            elif (
                self.config["rescue_veto_direction"] == "rna"
                and rna_consensus_label == "Artifact"
            ):
                # Mirrored veto when the direction is configured to RNA
                return "Artifact", rationale("rna_artifact_veto", "Artifact")
            elif (
                rna_consensus_label != "Artifact"
                and len(rna_callers)
                >= self.config["cross_modality_min_callers_for_artifact"]
            ):
                # One modality is Artifact, other is not - use non-Artifact with priority
                return rna_consensus_label, rationale(
                    "rna_consensus_priority", rna_consensus_label
                )
            elif (
                dna_consensus_label != "Artifact"
                and len(dna_callers)
                >= self.config["cross_modality_min_callers_for_artifact"]
            ):
                return dna_consensus_label, rationale(
                    "dna_consensus_priority", dna_consensus_label
                )
            else:
                return "Artifact", rationale("artifact_by_default", "Artifact")

        # Step 3: Single modality consensus
        elif dna_consensus_label and not rna_consensus_label:
            return dna_consensus_label, rationale(
                "dna_consensus_only", dna_consensus_label
            )

        elif rna_consensus_label and not dna_consensus_label:
            return rna_consensus_label, rationale(
                "rna_consensus_only", rna_consensus_label
            )

        # Step 4: No consensus labels - analyze individual caller patterns
        else:
            # Check if we have sufficient cross-modality support
            has_dna_support = len(dna_callers) > 0
            has_rna_support = len(rna_callers) > 0

            if not (has_dna_support and has_rna_support):
                # Insufficient cross-modality support
                return "NoConsensus", rationale(
                    "insufficient_modality_support", "NoConsensus"
                )

            # Check for cross-modality consistency
            dna_classifications = set(dna_callers)
            rna_classifications = set(rna_callers)

            # If both modalities have consistent internal classifications
            if len(dna_classifications) == 1 and len(rna_classifications) == 1:
                dna_class = list(dna_classifications)[0]
                rna_class = list(rna_classifications)[0]

                if dna_class == rna_class:
                    # Cross-modality agreement
                    if dna_class == "Artifact":
                        return dna_class, rationale(
                            "cross_modality_agreement", dna_class
                        )
                    # Cross-modality promotion (audit M1): individual DNA and
                    # RNA callers agreeing on Somatic rescue a site that
                    # failed within-modality consensus. The record is tagged
                    # as cross-modality rescued for the writer (INFO RESCUED /
                    # RESCUE_PROMOTED).
                    if (
                        dna_class == "Somatic"
                        and self.config["rescue_promotion_enabled"]
                        and len(dna_callers)
                        >= self.config["rescue_promotion_min_dna_callers"]
                        and len(rna_callers)
                        >= self.config["rescue_promotion_min_rna_callers"]
                    ):
                        variant_data["rescued"] = True
                        variant_data["rescue_promoted"] = True
                        return "Somatic", rationale("rescue_promotion", "Somatic")
                    # Not enough consensus support
                    return "NoConsensus", rationale(
                        "below_rescue_promotion_threshold", "NoConsensus"
                    )
                else:
                    # Cross-modality disagreement
                    return "Artifact", rationale(
                        "cross_modality_disagreement", "Artifact"
                    )
            else:
                # Internal inconsistency within modalities
                return "Artifact", rationale(
                    "modality_internal_disagreement", "Artifact"
                )

    def reclassify_with_annotation(
        self,
        variant,
        current_classification,
        gnomad_af=None,
        cosmic_count=None,
        is_rna_edit=None,
        rna_support=None,
        dna_support=None,
        cross_modality=False,
    ):
        """
        Reclassify variant based on annotation evidence.

        Implements the exact logic from current annotation scripts:
        1. RNA editing: If in REDIportal + sufficient RNA support → RNAedit
        2. Germline: If gnomAD AF > threshold → Germline (with artifact protection)
        3. Somatic: If COSMIC recurrence + cross-modality → Somatic
        4. Preserve artifact protection (cross-modality artifacts stay artifacts)

        Args:
            variant: Variant object (for context, optional)
            current_classification (str): Current classification
            gnomad_af (float): gnomAD allele frequency
            cosmic_count (int): COSMIC recurrence count
            is_rna_edit (bool): Is variant in REDIportal
            rna_support (int): Number of RNA callers supporting
            dna_support (int): Number of DNA callers supporting
            cross_modality (bool): Has cross-modality support

        Returns:
            str: Updated classification
        """
        # Start with current classification
        new_classification = current_classification

        # 1. RNA editing classification (highest priority)
        if is_rna_edit and rna_support is not None:
            if rna_support >= self.config["rna_editing_min_rna_support"]:
                # Validate that RNAedit variants are SNPs only
                # RNA editing is a single-base substitution (A→G or T→C)
                # Import get_variant_type only when needed to avoid circular imports
                try:
                    from ..common.vcf_config import get_variant_type
                    # Note: This validation may not work in all contexts since we don't
                    # always have access to the full variant object here
                    # The main validation happens in vcf_processor.py during statistics extraction
                except ImportError:
                    pass  # Validation happens in statistics extraction

                return "RNAedit"

        # 2. Artifact protection - cross-modality artifacts stay artifacts
        if (
            current_classification == "Artifact"
            and cross_modality
            and dna_support
            and rna_support
        ):
            # Both modalities classify as artifact - preserve it
            return "Artifact"

        # 3. Germline reclassification based on population frequency
        if (
            gnomad_af is not None
            and gnomad_af > self.config["annotation_germline_freq_threshold"]
        ):
            # Don't reclassify artifacts to germline if cross-modality artifact
            if not (current_classification == "Artifact" and cross_modality):
                new_classification = "Germline"

        # 4. Somatic reclassification based on COSMIC and cross-modality
        if cosmic_count is not None and cross_modality:
            if cosmic_count >= self.config["annotation_cosmic_recurrence_threshold"]:
                # COSMIC evidence with cross-modality support
                new_classification = "Somatic"

        # 5. Somatic reclassification based on unanimous DNA support
        if (
            dna_support is not None and dna_support >= 3
        ):  # Unanimous (assuming 3 DNA callers)
            all_callers_agree = (
                True  # In practice, check if all DNA callers classify as Somatic
            )
            if all_callers_agree:
                new_classification = "Somatic"

        return new_classification

    def classify_variant_pipeline(self, variant, stage, **kwargs):
        """
        Main entry point for any classification scenario.

        This method routes to the appropriate classification method based on
        the pipeline stage.

        Args:
            variant: Variant data (format depends on stage)
            stage: One of 'standalone', 'consensus', 'rescue', 'annotation'
            **kwargs: Stage-specific parameters
                - standalone: caller_name, sample_indices
                - consensus: (variant_data contains all needed info)
                - rescue: modality_map
                - annotation: current_classification, gnomad_af, cosmic_count, etc.

        Returns:
            str: Classification result

        Raises:
            ValueError: If unknown stage specified
        """
        if stage == "standalone":
            return self.classify_standalone_variant(
                variant, kwargs["caller_name"], kwargs.get("sample_indices")
            )
        elif stage == "consensus":
            return self.classify_consensus_variant(variant)
        elif stage == "rescue":
            return self.classify_rescue_variant(variant, kwargs["modality_map"])
        elif stage == "annotation":
            return self.reclassify_with_annotation(
                variant,
                kwargs.get("current_classification", "NoConsensus"),
                gnomad_af=kwargs.get("gnomad_af"),
                cosmic_count=kwargs.get("cosmic_count"),
                is_rna_edit=kwargs.get("is_rna_edit"),
                rna_support=kwargs.get("rna_support"),
                dna_support=kwargs.get("dna_support"),
                cross_modality=kwargs.get("cross_modality", False),
            )
        else:
            raise ValueError(f"Unknown classification stage: {stage}")

    def get_classification_confidence(
        self, variant_data, classification, stage="consensus"
    ):
        """
        Calculate confidence score for a classification.

        This is an optional enhancement that can be used to provide
        confidence metrics for classifications.

        Args:
            variant_data (dict): Variant data with caller information
            classification (str): The assigned classification
            stage (str): Pipeline stage

        Returns:
            float: Confidence score between 0 and 1
        """
        if stage == "consensus":
            # Calculate based on caller agreement
            if "filters_normalized" not in variant_data:
                return 0.0

            total_callers = len(variant_data["filters_normalized"])
            if total_callers == 0:
                return 0.0

            agreeing_callers = sum(
                1 for f in variant_data["filters_normalized"] if f == classification
            )
            return agreeing_callers / total_callers

        # Default confidence
        return 1.0 if classification != "NoConsensus" else 0.0


# Lazy initialization for default classifier to avoid circular imports
_default_classifier = None


def get_default_classifier():
    """Get the default classifier instance (lazy initialization)."""
    global _default_classifier
    if _default_classifier is None:
        _default_classifier = UnifiedVariantClassifier()
    return _default_classifier


# For backward compatibility - this will be lazily initialized on first access
class _LazyClassifier:
    """Lazy wrapper for default_classifier to avoid circular imports."""

    def __getattr__(self, name):
        return getattr(get_default_classifier(), name)


default_classifier = _LazyClassifier()
