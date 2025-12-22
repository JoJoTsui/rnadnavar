#!/usr/bin/env python3
"""
Annotation-Specific Classification Functions (Notebook Analysis Only)

** ANALYSIS-SPECIFIC MODULE **
This module contains functions for detecting special annotations that are only
needed for Jupyter notebook analysis, NOT for the main Nextflow pipeline.

Specifically, this module provides:
- detect_rna_edit_variant(): Detection of RNA editing annotations
- detect_cosmic_gnomad_annotation(): Detection of COSMIC/gnomAD annotations
- detect_no_consensus_variant(): Detection of non-consensus variants
- classify_annotated_variant(): Combined annotation classification

For BASE variant classification (Somatic/Germline/Reference/Artifact), use:
- bin/vcf_utils/classification.py (canonical - used by production pipeline)

This module wraps those canonical functions for notebook-specific analysis.
"""

from typing import Dict, Optional, Any


def detect_rna_edit_variant(variant) -> bool:
    """
    Detect if a variant has RNA editing annotation.
    
    Checks for REDIportal or DARNED annotations in INFO fields.
    
    Args:
        variant: cyvcf2.Variant object
        
    Returns:
        True if variant has RNA editing annotation, False otherwise
    """
    try:
        if hasattr(variant, 'INFO'):
            info = variant.INFO
            # Check for explicit RNA editing field (from rna_editing stage)
            if info.get('RNA_EDIT', False):
                return True
            # Check for REDIportal/DARNED directly
            if info.get('REDIportal') or info.get('DARNED'):
                return True
            # Check for REDI_* fields (explicit REDIportal annotations)
            for key in info.keys():
                if key.startswith('REDI_'):
                    return True
        return False
    except Exception:
        return False


def detect_cosmic_gnomad_annotation(variant) -> bool:
    """
    Detect if a variant has COSMIC or gnomAD annotation.
    
    Checks for COSMIC_CNT, gnomAD_AF, or other database annotations.
    
    Args:
        variant: cyvcf2.Variant object
        
    Returns:
        True if variant has database annotation, False otherwise
    """
    try:
        if hasattr(variant, 'INFO'):
            info = variant.INFO
            # Check for COSMIC annotation
            if info.get('COSMIC_CNT', 0) > 0:
                return True
            if info.get('cosmic_count', 0) > 0:
                return True
            if info.get('COSMIC_ID'):
                return True
            # Check for gnomAD annotation
            gnomad_af = info.get('gnomAD_AF')
            if gnomad_af is not None:
                try:
                    if float(gnomad_af) > 0:
                        return True
                except (ValueError, TypeError):
                    pass
            # Check for other database fields
            if info.get('AF_gnomad'):
                return True
        return False
    except Exception:
        return False


def detect_no_consensus_variant(variant) -> bool:
    """
    Detect if a variant lacks consensus (DNA/RNA disagree or no consensus).
    
    Args:
        variant: cyvcf2.Variant object
        
    Returns:
        True if variant is NoConsensus, False otherwise
    """
    try:
        if hasattr(variant, 'FILTER'):
            filter_val = variant.FILTER
            if filter_val and 'NoConsensus' in str(filter_val):
                return True
        if hasattr(variant, 'INFO'):
            info = variant.INFO
            if info.get('NO_CONSENSUS'):
                return True
            if info.get('NoConsensus'):
                return True
        return False
    except Exception:
        return False


def classify_annotated_variant(variant) -> str:
    """
    Classify a variant based on annotation, returning special categories.
    
    Returns special categories for annotated VCFs:
    - 'RNA_Edit': If has RNA editing annotation
    - 'Cosmic_gnomAD': If has COSMIC/gnomAD annotation
    - 'NoConsensus': If lacks consensus
    - 'Somatic': Default/other
    
    Args:
        variant: cyvcf2.Variant object
        
    Returns:
        Annotation category (RNA_Edit, Cosmic_gnomAD, NoConsensus, Somatic)
    """
    # Check annotations in priority order
    if detect_no_consensus_variant(variant):
        return 'NoConsensus'
    
    if detect_rna_edit_variant(variant):
        return 'RNA_Edit'
    
    if detect_cosmic_gnomad_annotation(variant):
        return 'Cosmic_gnomAD'
    
    # Default to Somatic for anything else
    return 'Somatic'
