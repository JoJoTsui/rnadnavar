#!/usr/bin/env python3
"""
Annotation Configuration Module

Centralized configuration for COSMIC, gnomAD, and REDIportal annotation pipelines.
Provides a single source of truth for all annotation thresholds and parameters.

This module ensures consistency across all annotation modules and simplifies
parameter management for the entire annotation infrastructure.

USAGE:
    from annotation_config import COSMIC_THRESHOLDS, GNOMAD_THRESHOLDS, REDIPORTAL_THRESHOLDS
    
    min_recurrence = COSMIC_THRESHOLDS['recurrence_minimum']
    germline_freq = GNOMAD_THRESHOLDS['germline_frequency']
    min_rna_support = REDIPORTAL_THRESHOLDS['min_rna_support']

Author: COSMIC/gnomAD/REDIportal Enhancement Pipeline
Date: 2025-12-22
Version: 1.0
"""

from typing import Dict, Any, List, Tuple
from dataclasses import dataclass
from enum import Enum

# ============================================================================
# COSMIC ANNOTATION CONFIGURATION
# ============================================================================

COSMIC_THRESHOLDS = {
    # Minimum mutation recurrence in COSMIC database for rescue consideration
    # Rationale: Filters singleton mutations, focuses on recurrent cancer drivers
    # Range: 1-20 (1 = all matches, 10+ = well-established drivers)
    # Recommendation: 5 for balanced sensitivity/specificity
    "recurrence_minimum": 5,
    
    # Minimum number of samples with mutation in COSMIC screens
    # Rationale: Filters artifacts, requires multiple sample validation
    # Range: 1-100 (1 = minimum, 2+ = cross-sample validation)
    # Recommendation: 2 for basic validation
    "sample_count_minimum": 2,
    
    # COSMIC annotation is database lookup - no hard thresholds
    # All matches returned, interpretation at variant classification stage
}

COSMIC_DESCRIPTION = """
COSMIC Annotation Configuration

Purpose: Annotate somatic mutations with cancer mutation database evidence

Annotation Strategy:
  - VCF-to-VCF exact matching (CHROM, POS, REF, ALT)
  - Extract mutation recurrence (COSMIC_CNT)
  - Extract sample frequency (COSMIC_SAMPLE_COUNT)
  - Extract COSMIC identifiers (COSMIC_ID)

Thresholds:
  - recurrence_minimum: Minimum COSMIC occurrences to consider significant
  - sample_count_minimum: Minimum samples supporting mutation

Interpretation:
  - COSMIC_CNT ≥ recurrence_minimum: Known cancer driver candidate
  - COSMIC_CNT < recurrence_minimum: Rare or novel somatic mutation
  - COSMIC_CNT = 0: Not in database (novel variant)

Database Information:
  - Source: COSMIC v103 or later
  - Format: VCF (bgzipped and indexed)
  - Coverage: ~22 million somatic mutations from cancer screens
  - Update frequency: Monthly

Clinical Use:
  - Somatic variant prioritization
  - Cancer driver identification
  - Recurrence validation
"""

# ============================================================================
# GNOMAD ANNOTATION CONFIGURATION
# ============================================================================

GNOMAD_THRESHOLDS = {
    # Population frequency threshold for germline classification
    # Rationale: Clinical standard (ACMG) for benign variant cutoff
    # Range: 0.001-0.05 (0.1%-5%)
    # Recommendation: 0.001 (ACMG standard for rare disease)
    "germline_frequency": 0.001,
    
    # Allele count minimum to avoid singleton errors
    # Rationale: Single observations in large databases can be artifacts
    # Range: 1-10 (1 = any count, 2+ = cross-chromosome validation)
    # Recommendation: 2 for basic validation
    "allele_count_minimum": 2,
    
    # Maximum allele count for rare variant definition
    # Rationale: Optional threshold to define "common" variants
    # Range: 1000-10000 (higher = more permissive)
    # Recommendation: Not enforced (use AF threshold instead)
    "allele_count_maximum": None,
    
    # Filter allele frequency (95% CI upper bound)
    # Rationale: Conservative estimate of AF for filtering
    # Range: 0.001-0.05
    # Recommendation: 0.001 (same as germline_frequency)
    "filter_allele_frequency_95ci": 0.001,
}

GNOMAD_DESCRIPTION = """
gnomAD Annotation Configuration

Purpose: Annotate variants with population frequency data for germline classification

Annotation Strategy:
  - Scatter-gather parallelization by chromosome
  - VCF-to-VCF annotation with chromosome-split gnomAD files
  - Extract population allele frequencies
  - Auto-optimize parallel workers

Database Information:
  - Source: gnomAD v4.0 or later
  - Coverage: ~141,000 genomes + ~7.6 million exomes
  - Populations: 7 major populations + all combined
  - Format: VCF per chromosome (bgzipped and indexed)

Thresholds:
  - germline_frequency: Frequency cutoff for benign variants (ACMG standard)
  - allele_count_minimum: Minimum observation count
  - filter_allele_frequency_95ci: Conservative upper bound estimate

Interpretation:
  - AF > germline_frequency: Likely benign/common variant (germline)
  - AF ≤ germline_frequency: Rare or novel variant (potential disease-causing)
  - AC < allele_count_minimum: Potential singleton error (exclude)

Clinical Use:
  - Germline variant classification
  - Population frequency annotation
  - Benign variant filtering (ACMG standard)
  - Hardy-Weinberg equilibrium assessment

Note on 0.001 Threshold:
  - ACMG recommendation for rare disease genetic testing
  - Filters out variants with >0.1% population frequency
  - Conservative cutoff that prioritizes sensitivity
  - Can be adjusted for specific disease prevalence
"""

# ============================================================================
# REDIPORTAL ANNOTATION CONFIGURATION
# ============================================================================

REDIPORTAL_THRESHOLDS = {
    # Minimum number of RNA callers supporting variant for classification
    # Rationale: Requires cross-caller validation to reduce artifacts
    # Range: 1-3 (1 = any caller, 2 = dual validation, 3 = high confidence)
    # Recommendation: 2 for balanced sensitivity/specificity
    "min_rna_support": 2,
    
    # Canonical RNA editing transitions (A-to-I editing signature)
    # Rationale: A>G and T>C are hallmark RNA editing signatures
    # Recommendation: Keep canonical transitions only (avoid spurious calls)
    "canonical_transitions": [
        ("A", "G"),  # A-to-G (A-to-I)
        ("T", "C"),  # T-to-C (reverse complement of A-to-G)
    ],
    
    # Minimum RNA VAF for VERY_HIGH/HIGH evidence classification
    # Rationale: Clear signal above sequencing noise
    # Range: 0.05-0.25 (5%-25%)
    # Recommendation: 0.10 (10% - standard VAF threshold)
    "rna_vaf_minimum": 0.10,
    
    # Maximum DNA VAF to distinguish RNA-specific editing from DNA mutation
    # Rationale: RNA editing should be RNA-specific, not present in DNA
    # Range: 0.01-0.10 (1%-10%)
    # Recommendation: 0.05 (5% - allows minor DNA contamination)
    "dna_vaf_maximum": 0.05,
    
    # REDIportal exact match requirement (CHROM, POS)
    # Rationale: Database-validated sites only
    # Recommendation: Keep True (exact matches)
    "require_exact_match": True,
}

REDIPORTAL_EVIDENCE_LEVELS = {
    # Evidence classification matrix (4-tier system + NONE)
    # Higher level = stronger evidence for RNA editing site
    
    "VERY_HIGH": {
        "criteria": "Canonical + RNA-only + High RNA VAF + Low DNA VAF + REDIportal + RNA consensus",
        "description": "Strongest evidence for RNA editing site",
        "action": "Reclassify to RNAedit",
        "quality_score": 100,
    },
    
    "HIGH": {
        "criteria": "Canonical + RNA-only + High RNA VAF + REDIportal + RNA consensus",
        "description": "Strong evidence for RNA editing site",
        "action": "Reclassify to RNAedit",
        "quality_score": 80,
    },
    
    "MEDIUM": {
        "criteria": "Canonical + REDIportal + RNA consensus + DNA presence",
        "description": "Moderate evidence for RNA editing (DNA present - possible contamination)",
        "action": "Reclassify to RNAedit",
        "quality_score": 60,
    },
    
    "LOW": {
        "criteria": "Non-canonical + REDIportal + RNA consensus",
        "description": "Weak evidence for RNA editing (non-canonical transition)",
        "action": "Reclassify to RNAedit",
        "quality_score": 40,
    },
    
    "NONE": {
        "criteria": "No exact REDIportal match OR insufficient RNA support",
        "description": "Insufficient evidence for RNA editing",
        "action": "Preserve original classification",
        "quality_score": 0,
    },
}

REDIPORTAL_DESCRIPTION = """
REDIportal Annotation Configuration

Purpose: Identify and classify RNA editing sites with evidence tiering

Annotation Strategy:
  - REDIportal database lookup (exact coordinate matching)
  - Canonical transition detection (A>G, T>C)
  - Evidence-based classification (VERY_HIGH → NONE)
  - Dynamic FILTER update based on evidence level

Database Information:
  - Source: REDIportal v3.0 or later
  - Coverage: ~4.6 million RNA editing sites in humans
  - Format: VCF (bgzipped and indexed)
  - Coverage: 55 datasets across different tissues/conditions

Thresholds:
  - min_rna_support: Minimum RNA callers (cross-validation)
  - canonical_transitions: A>G and T>C only
  - rna_vaf_minimum: Clear signal above noise (10%)
  - dna_vaf_maximum: RNA-specific sites (≤5% DNA)

Evidence Classification (4-tier):
  1. VERY_HIGH: Maximum confidence RNA editing
  2. HIGH: Strong RNA editing evidence
  3. MEDIUM: Moderate evidence with DNA presence
  4. LOW: Weak evidence (non-canonical)
  5. NONE: Insufficient evidence (preserve original)

Interpretation:
  - VERY_HIGH/HIGH: Reclassify to RNAedit
  - MEDIUM: Reclassify to RNAedit (likely editing)
  - LOW: Reclassify to RNAedit (weak support)
  - NONE: Keep original classification (genetic mutation)

Clinical Use:
  - RNA editing site validation
  - RNA mutation vs DNA mutation distinction
  - Tissue-specific variant interpretation
  - Functional consequence prediction

Canonical Transitions Rationale:
  - A>G (A-to-I) editing: Most common in humans
  - T>C (reverse complement): Same biological process
  - Non-canonical transitions: Ignored (too rare in humans)
"""

# ============================================================================
# MULTI-MODAL VARIANT CLASSIFICATION CONFIGURATION
# ============================================================================

CLASSIFICATION_THRESHOLDS = {
    # Minimum DNA caller support for somatic variant confidence
    # Rationale: Cross-caller validation reduces artifacts
    # Range: 1-3 (1 = any caller, 2 = dual validation, 3 = high confidence)
    # Recommendation: 2 for balanced somatic classification
    "somatic_consensus_threshold": 2,
    
    # Whether to require dual DNA+RNA support for rescue variants
    # Rationale: Cross-modality validation increases confidence
    # Range: True/False
    # Recommendation: True (stronger evidence)
    "require_dna_rna_dual_support": True,
    
    # Minimum DNA support when RNA support exists
    # Rationale: RNA alone may indicate RNA-specific variant
    # Range: 0-3
    # Recommendation: 1 (at least one DNA caller)
    "min_dna_support_with_rna": 1,
    
    # Minimum RNA support when DNA support exists
    # Rationale: DNA alone may indicate germline/common variant
    # Range: 0-3
    # Recommendation: 1 (RNA-specific validation)
    "min_rna_support_with_dna": 1,
    
    # Whether germline variants must have database support
    # Rationale: Germline in population database increases confidence
    # Range: True/False
    # Recommendation: True (clinical standard)
    "require_germline_database_support": True,
}

CLASSIFICATION_DESCRIPTION = """
Multi-Modal Variant Classification Configuration

Purpose: Define thresholds for comprehensive variant quality assessment

Classification Strategy:
  - DNA caller consensus (Somatic, Germline classification)
  - RNA caller support (RNA-specific validation)
  - Database evidence (gnomAD, COSMIC, REDIportal)
  - Cross-modality validation (DNA + RNA)

Thresholds:
  - somatic_consensus_threshold: Minimum DNA callers for somatic
  - require_dna_rna_dual_support: Cross-modality validation
  - min_dna_support_with_rna: DNA requirement when RNA present
  - min_rna_support_with_dna: RNA requirement when DNA present
  - require_germline_database_support: Germline validation in gnomAD

Classification Rules:
  1. Somatic: ≥2 DNA callers + (0 or ≥1 RNA caller) + COSMIC support (optional)
  2. Germline: ≥1 DNA caller + gnomAD AF > 0.001 (ACMG standard)
  3. RNA_Edit: REDIportal match + HIGH/VERY_HIGH evidence + RNA consensus
  4. Artifact: <2 DNA callers + no database support
  5. NoConsensus: DNA/RNA callers disagree on classification

Quality Scoring:
  - Somatic + COSMIC + RNA support: Highest confidence
  - Somatic + RNA support: High confidence
  - Somatic alone: Moderate confidence
  - Germline + gnomAD: High confidence
  - RNA_Edit + VERY_HIGH evidence: Highest for RNA editing
  - Artifact or NoConsensus: Lowest confidence (review required)
"""

# ============================================================================
# DATABASE PATH CONFIGURATION
# ============================================================================

DATABASE_PATHS = {
    # COSMIC database path (template, override in config)
    "cosmic_vcf": "/t9k/mnt/joey/bio_db/COSMIC/ucsc_chr_cosmic_GenomeScreensMutant_Normal_v103_GRCh38.vcf.gz",
    "cosmic_tbi": "/t9k/mnt/joey/bio_db/COSMIC/ucsc_chr_cosmic_GenomeScreensMutant_Normal_v103_GRCh38.vcf.gz.tbi",
    
    # gnomAD database path (template, override in config)
    "gnomad_dir": "/t9k/mnt/joey/bio_db/gnomAD/exomes",
    "gnomad_version": "4.0",
    
    # REDIportal database path (template, override in config)
    "rediportal_vcf": "/t9k/mnt/joey/bio_db/rna_editing/REDIportal/REDIportal_hg38_v3.vcf.gz",
    "rediportal_tbi": "/t9k/mnt/joey/bio_db/rna_editing/REDIportal/REDIportal_hg38_v3.vcf.gz.tbi",
}

# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

def get_threshold(category: str, threshold_name: str) -> Any:
    """
    Get threshold value by category and name.
    
    Args:
        category: "COSMIC", "GNOMAD", "REDIPORTAL", "CLASSIFICATION"
        threshold_name: Name of specific threshold
    
    Returns:
        Threshold value
    
    Raises:
        KeyError: If category or threshold not found
    """
    config_map = {
        "COSMIC": COSMIC_THRESHOLDS,
        "GNOMAD": GNOMAD_THRESHOLDS,
        "REDIPORTAL": REDIPORTAL_THRESHOLDS,
        "CLASSIFICATION": CLASSIFICATION_THRESHOLDS,
    }
    
    if category not in config_map:
        raise KeyError(f"Unknown configuration category: {category}")
    
    if threshold_name not in config_map[category]:
        raise KeyError(f"Unknown threshold '{threshold_name}' in {category}")
    
    return config_map[category][threshold_name]


def get_all_thresholds() -> Dict[str, Dict[str, Any]]:
    """
    Get all configured thresholds as dictionary.
    
    Returns:
        Dictionary with all threshold categories and values
    """
    return {
        "COSMIC": COSMIC_THRESHOLDS,
        "GNOMAD": GNOMAD_THRESHOLDS,
        "REDIPORTAL": REDIPORTAL_THRESHOLDS,
        "CLASSIFICATION": CLASSIFICATION_THRESHOLDS,
    }


def print_configuration_summary() -> None:
    """Print comprehensive configuration summary for debugging."""
    print("=" * 80)
    print("ANNOTATION CONFIGURATION SUMMARY")
    print("=" * 80)
    
    print("\n📊 COSMIC THRESHOLDS:")
    for key, value in COSMIC_THRESHOLDS.items():
        print(f"  {key}: {value}")
    
    print("\n📊 GNOMAD THRESHOLDS:")
    for key, value in GNOMAD_THRESHOLDS.items():
        print(f"  {key}: {value}")
    
    print("\n📊 REDIPORTAL THRESHOLDS:")
    for key, value in REDIPORTAL_THRESHOLDS.items():
        if key != "canonical_transitions":
            print(f"  {key}: {value}")
        else:
            print(f"  {key}: {value}")
    
    print("\n📊 CLASSIFICATION THRESHOLDS:")
    for key, value in CLASSIFICATION_THRESHOLDS.items():
        print(f"  {key}: {value}")
    
    print("\n📊 DATABASE PATHS:")
    for key, value in DATABASE_PATHS.items():
        print(f"  {key}: {value}")
    
    print("\n" + "=" * 80)


def validate_thresholds() -> Dict[str, Any]:
    """
    Validate all configured thresholds for consistency.
    
    Returns:
        Dictionary with validation results
    """
    issues = []
    warnings = []
    
    # Validate COSMIC
    if COSMIC_THRESHOLDS["recurrence_minimum"] < 1:
        issues.append("COSMIC recurrence_minimum must be ≥ 1")
    if COSMIC_THRESHOLDS["sample_count_minimum"] < 1:
        issues.append("COSMIC sample_count_minimum must be ≥ 1")
    
    # Validate gnomAD
    if GNOMAD_THRESHOLDS["germline_frequency"] < 0 or GNOMAD_THRESHOLDS["germline_frequency"] > 0.5:
        issues.append("gnomAD germline_frequency must be 0-0.5")
    if GNOMAD_THRESHOLDS["allele_count_minimum"] < 1:
        issues.append("gnomAD allele_count_minimum must be ≥ 1")
    
    # Validate REDIportal
    if REDIPORTAL_THRESHOLDS["min_rna_support"] < 1:
        issues.append("REDIportal min_rna_support must be ≥ 1")
    if REDIPORTAL_THRESHOLDS["rna_vaf_minimum"] < 0 or REDIPORTAL_THRESHOLDS["rna_vaf_minimum"] > 1:
        issues.append("REDIportal rna_vaf_minimum must be 0-1")
    if REDIPORTAL_THRESHOLDS["dna_vaf_maximum"] < 0 or REDIPORTAL_THRESHOLDS["dna_vaf_maximum"] > 1:
        issues.append("REDIportal dna_vaf_maximum must be 0-1")
    
    # Warn about threshold values
    if GNOMAD_THRESHOLDS["germline_frequency"] != 0.001:
        warnings.append(f"gnomAD germline_frequency is {GNOMAD_THRESHOLDS['germline_frequency']} (ACMG standard is 0.001)")
    
    if REDIPORTAL_THRESHOLDS["min_rna_support"] == 1:
        warnings.append("REDIportal min_rna_support = 1 allows single-caller artifacts")
    
    return {
        "valid": len(issues) == 0,
        "issues": issues,
        "warnings": warnings,
    }


if __name__ == "__main__":
    # Print configuration when run directly
    print_configuration_summary()
    
    # Validate thresholds
    validation = validate_thresholds()
    if validation["issues"]:
        print("\n⚠️  VALIDATION ISSUES:")
        for issue in validation["issues"]:
            print(f"  ✗ {issue}")
    else:
        print("\n✓ All thresholds validated successfully")
    
    if validation["warnings"]:
        print("\n⚠️  WARNINGS:")
        for warning in validation["warnings"]:
            print(f"  ! {warning}")
