#!/usr/bin/env python3
"""
Database Evidence Checker Module

This module implements database evidence detection for the hybrid tiering system.
It checks whether variants are present in external databases:
- gnomAD: Population frequency database (germline variants)
- COSMIC: Cancer somatic mutation database  
- REDIportal: RNA editing database
- DARNED: RNA editing database

Usage:
    from database_checker import check_database_support, compute_database_tier
    
    # Check if variant has database support
    has_support = check_database_support(
        gnomad_af=0.002,
        cosmic_cnt=5,
        rediportal=True
    )
    # Returns: True (gnomAD AF > 0.001 threshold)
    
    # Get database tier
    tier = compute_database_tier(
        gnomad_af=0.002,
        cosmic_cnt=5
    )
    # Returns: "D1" (has database support)
"""

from typing import Dict, Optional, Any, List
import re

# Import centralized annotation configuration (fallback-safe)
try:
    from annotation_config import (
        GNOMAD_THRESHOLDS,
        COSMIC_THRESHOLDS,
        REDIPORTAL_THRESHOLDS,
    )
except Exception:
    # Fallback defaults to preserve runtime if config is unavailable
    GNOMAD_THRESHOLDS = {"germline_frequency": 0.001}
    COSMIC_THRESHOLDS = {"recurrence_minimum": 0}
    REDIPORTAL_THRESHOLDS = {"min_rna_support": 2}


# ============================================================================
# DATABASE FIELD DEFINITIONS
# ============================================================================

# INFO field names for database annotations
DATABASE_INFO_FIELDS = {
    "gnomAD": [
        "gnomAD_AF",           # Allele frequency in gnomAD
        "gnomAD_AC",           # Allele count in gnomAD
        "AF_gnomad",           # Alternative field name
        "gnomad_af",           # Case-insensitive variant
    ],
    "COSMIC": [
        "COSMIC_CNT",          # Count in COSMIC database
        "COSMIC_ID",           # COSMIC identifier
        "cosmic_count",        # Alternative field name
    ],
    "REDIportal": [
        # Explicit REDIportal INFO fields produced by the pipeline
        "REDI_ACCESSION",
        "REDI_DB",
        "REDI_TYPE",
        "REDI_REPEAT",
        "REDI_FUNC",
        "REDI_STRAND",
        # Legacy flags retained for backward compatibility (optional)
        "REDIportal",
        "REDI",
        "RNA_EDIT_REDI",
    ],
}

# Thresholds for database significance (config-driven, with safe fallbacks)
DATABASE_THRESHOLDS = {
    # gnomAD: use configured germline frequency threshold
    "gnomAD_AF": GNOMAD_THRESHOLDS.get("germline_frequency", 0.001),
    # COSMIC: use configured recurrence minimum (any count above this is supportive)
    "COSMIC_CNT": COSMIC_THRESHOLDS.get("recurrence_minimum", 0),
    # REDIportal: presence-driven evidence
    "REDIportal": True,
}


# ============================================================================
# DATABASE SUPPORT DETECTION
# ============================================================================

def check_gnomad_support(info_dict: Dict[str, Any]) -> bool:
    """
    Check if variant has significant gnomAD support.
    
    Args:
        info_dict: Dictionary of INFO fields from VCF record
    
    Returns:
        True if gnomAD AF > threshold (0.001)
    """
    threshold = DATABASE_THRESHOLDS["gnomAD_AF"]
    
    for field in DATABASE_INFO_FIELDS["gnomAD"]:
        value = info_dict.get(field)
        if value is not None:
            try:
                af = float(value)
                if af > threshold:
                    return True
            except (ValueError, TypeError):
                continue
    
    return False


def check_cosmic_support(info_dict: Dict[str, Any]) -> bool:
    """
    Check if variant has COSMIC database support.
    
    Args:
        info_dict: Dictionary of INFO fields from VCF record
    
    Returns:
        True if COSMIC count > 0 or COSMIC ID present
    """
    threshold = DATABASE_THRESHOLDS["COSMIC_CNT"]
    
    for field in DATABASE_INFO_FIELDS["COSMIC"]:
        value = info_dict.get(field)
        if value is not None:
            # Check for COSMIC ID (string)
            if isinstance(value, str) and value not in [".", "", "NA"]:
                return True
            
            # Check for COSMIC count (numeric)
            try:
                count = int(value)
                if count > threshold:
                    return True
            except (ValueError, TypeError):
                continue
    
    return False


def check_rediportal_support(info_dict: Dict[str, Any]) -> bool:
    """
    Check if variant is in REDIportal RNA editing database.
    
    Args:
        info_dict: Dictionary of INFO fields from VCF record
    
    Returns:
        True if present in REDIportal
    """
    for field in DATABASE_INFO_FIELDS["REDIportal"]:
        value = info_dict.get(field)
        if value is not None:
            # Check boolean/flag fields
            if isinstance(value, bool):
                return value
            
            # Check string fields (True, 1, yes, etc.)
            if isinstance(value, str):
                value_lower = value.lower()
                if value_lower in ["true", "1", "yes", "t", "y"]:
                    return True
                if value_lower not in ["false", "0", "no", "f", "n", ".", "na"]:
                    # Non-empty string that's not explicitly false
                    return True
            
            # Check numeric fields
            if isinstance(value, (int, float)) and value > 0:
                return True
    
    return False


def check_database_support(
    gnomad_af: Optional[float] = None,
    cosmic_cnt: Optional[int] = None,
    rediportal: Optional[bool] = None,
    info_dict: Optional[Dict[str, Any]] = None,
) -> bool:
    """
    Check if variant has support from any database.
    
    Can be called with individual parameters or with an info_dict containing
    all INFO fields from VCF record.
    
    Args:
        gnomad_af: gnomAD allele frequency
        cosmic_cnt: COSMIC count
        rediportal: Present in REDIportal
        darned: Present in DARNED
        info_dict: Dictionary of all INFO fields (alternative to individual params)
    
    Returns:
        True if variant has support from at least one database
    
    Example:
        >>> # Using individual parameters
        >>> check_database_support(gnomad_af=0.002)
        True
        
        >>> # Using info_dict
        >>> info = {"gnomAD_AF": 0.002, "COSMIC_CNT": 0}
        >>> check_database_support(info_dict=info)
        True
    """
    # If info_dict provided, use it
    if info_dict is not None:
        return (
            check_gnomad_support(info_dict) or
            check_cosmic_support(info_dict) or
            check_rediportal_support(info_dict)
        )
    
    # Otherwise check individual parameters
    has_gnomad = gnomad_af is not None and gnomad_af > DATABASE_THRESHOLDS["gnomAD_AF"]
    has_cosmic = cosmic_cnt is not None and cosmic_cnt > DATABASE_THRESHOLDS["COSMIC_CNT"]
    has_rediportal = rediportal is True
    
    return has_gnomad or has_cosmic or has_rediportal


def get_supporting_databases(
    info_dict: Dict[str, Any]
) -> List[str]:
    """
    Get list of databases that support this variant.
    
    Args:
        info_dict: Dictionary of INFO fields from VCF record
    
    Returns:
        List of database names with support (e.g., ["gnomAD", "COSMIC"])
    """
    supporting = []
    
    if check_gnomad_support(info_dict):
        supporting.append("gnomAD")
    
    if check_cosmic_support(info_dict):
        supporting.append("COSMIC")
    
    if check_rediportal_support(info_dict):
        supporting.append("REDIportal")
    
    # DARNED not supported in this workflow
    
    return supporting


def get_database_details(
    info_dict: Dict[str, Any]
) -> Dict[str, Any]:
    """
    Get detailed database annotation values.
    
    Args:
        info_dict: Dictionary of INFO fields from VCF record
    
    Returns:
        Dictionary with database values:
        {
            "gnomAD_AF": 0.002,
            "COSMIC_CNT": 5,
            "REDIportal": True,
            "DARNED": False,
            "supporting_databases": ["gnomAD", "COSMIC"]
        }
    """
    details = {
        "gnomAD_AF": None,
        "COSMIC_CNT": None,
        "REDIportal": False,
        # DARNED intentionally unsupported in this workflow
        "DARNED": False,
    }
    
    # Extract gnomAD AF
    for field in DATABASE_INFO_FIELDS["gnomAD"]:
        value = info_dict.get(field)
        if value is not None:
            try:
                details["gnomAD_AF"] = float(value)
                break
            except (ValueError, TypeError):
                continue
    
    # Extract COSMIC count
    for field in DATABASE_INFO_FIELDS["COSMIC"]:
        value = info_dict.get(field)
        if value is not None:
            try:
                details["COSMIC_CNT"] = int(value)
                break
            except (ValueError, TypeError):
                # Could be COSMIC ID string
                if isinstance(value, str) and value not in [".", "", "NA"]:
                    details["COSMIC_CNT"] = 1  # Mark as present
                    break
    
    # Check REDIportal
    details["REDIportal"] = check_rediportal_support(info_dict)
    
    # Add list of supporting databases
    details["supporting_databases"] = get_supporting_databases(info_dict)
    
    return details


# ============================================================================
# DATABASE TIER ASSIGNMENT
# ============================================================================

def compute_database_tier(
    gnomad_af: Optional[float] = None,
    cosmic_cnt: Optional[int] = None,
    rediportal: Optional[bool] = None,
    info_dict: Optional[Dict[str, Any]] = None,
) -> str:
    """
    Compute database evidence tier (D0 or D1).
    
    Args:
        gnomad_af: gnomAD allele frequency
        cosmic_cnt: COSMIC count
        rediportal: Present in REDIportal
        darned: Present in DARNED
        info_dict: Dictionary of all INFO fields (alternative to individual params)
    
    Returns:
        "D1" if has database support, "D0" otherwise
    
    Example:
        >>> compute_database_tier(gnomad_af=0.002)
        "D1"
        
        >>> compute_database_tier(gnomad_af=0.0001)
        "D0"
    """
    has_support = check_database_support(
        gnomad_af=gnomad_af,
        cosmic_cnt=cosmic_cnt,
        rediportal=rediportal,
        info_dict=info_dict
    )
    
    return "D1" if has_support else "D0"


def compute_database_tier_from_vcf_record(vcf_record) -> str:
    """
    Compute database tier directly from a cyvcf2 VCF record.
    
    Args:
        vcf_record: cyvcf2.Variant object with INFO fields
    
    Returns:
        "D1" if has database support, "D0" otherwise
    
    Example:
        >>> import cyvcf2
        >>> vcf = cyvcf2.VCF("rescue.vcf.gz")
        >>> for variant in vcf:
        ...     db_tier = compute_database_tier_from_vcf_record(variant)
        ...     print(f"{variant.CHROM}:{variant.POS} - {db_tier}")
    """
    # Extract INFO dict from VCF record
    if hasattr(vcf_record, 'INFO'):
        info_dict = dict(vcf_record.INFO)
    else:
        info_dict = {}
    
    return compute_database_tier(info_dict=info_dict)


# ============================================================================
# VALIDATION AND DEBUGGING UTILITIES
# ============================================================================

def print_database_summary(info_dict: Dict[str, Any]) -> None:
    """
    Print human-readable summary of database annotations for debugging.
    
    Args:
        info_dict: Dictionary of INFO fields from VCF record
    """
    details = get_database_details(info_dict)
    tier = compute_database_tier(info_dict=info_dict)
    
    print(f"Database Tier: {tier}")
    print(f"\nDatabase Annotations:")
    print(f"  gnomAD AF:     {details['gnomAD_AF']}")
    print(f"  COSMIC Count:  {details['COSMIC_CNT']}")
    print(f"  REDIportal:    {details['REDIportal']}")
    # DARNED removed from workflow
    print(f"\nSupporting Databases: {', '.join(details['supporting_databases']) or 'None'}")
    
    # Show which thresholds are met
    print(f"\nThreshold Check:")
    if details['gnomAD_AF'] is not None:
        status = "✓" if details['gnomAD_AF'] > DATABASE_THRESHOLDS["gnomAD_AF"] else "✗"
        print(f"  {status} gnomAD AF > {DATABASE_THRESHOLDS['gnomAD_AF']}")
    if details['COSMIC_CNT'] is not None:
        status = "✓" if details['COSMIC_CNT'] > DATABASE_THRESHOLDS["COSMIC_CNT"] else "✗"
        print(f"  {status} COSMIC count > {DATABASE_THRESHOLDS['COSMIC_CNT']}")
    if details['REDIportal']:
        print(f"  ✓ Present in REDIportal")
    if details['DARNED']:
        print(f"  ✓ Present in DARNED")


def validate_database_fields(info_dict: Dict[str, Any]) -> Dict[str, Any]:
    """
    Validate and report on database field presence.
    
    Args:
        info_dict: Dictionary of INFO fields from VCF record
    
    Returns:
        Dictionary with validation results:
        {
            "has_gnomad_fields": True/False,
            "has_cosmic_fields": True/False,
            "has_rna_edit_fields": True/False,
            "found_fields": ["gnomAD_AF", "COSMIC_CNT", ...],
        }
    """
    found_fields = []
    
    # Check for any gnomAD field
    has_gnomad = False
    for field in DATABASE_INFO_FIELDS["gnomAD"]:
        if field in info_dict:
            has_gnomad = True
            found_fields.append(field)
    
    # Check for any COSMIC field
    has_cosmic = False
    for field in DATABASE_INFO_FIELDS["COSMIC"]:
        if field in info_dict:
            has_cosmic = True
            found_fields.append(field)
    
    # Check for RNA editing fields (REDIportal only)
    has_rna_edit = False
    for field in DATABASE_INFO_FIELDS["REDIportal"]:
        if field in info_dict:
            has_rna_edit = True
            found_fields.append(field)
    
    return {
        "has_gnomad_fields": has_gnomad,
        "has_cosmic_fields": has_cosmic,
        "has_rna_edit_fields": has_rna_edit,
        "found_fields": found_fields,
    }


if __name__ == "__main__":
    # Example usage and testing
    print("=" * 80)
    print("DATABASE EVIDENCE CHECKER - EXAMPLES")
    print("=" * 80)
    
    # Example 1: Variant with gnomAD support
    print("\nExample 1: Germline variant with gnomAD support")
    info1 = {"gnomAD_AF": 0.002, "COSMIC_CNT": 0}
    print_database_summary(info1)
    
    # Example 2: Somatic variant with COSMIC support
    print("\n" + "=" * 80)
    print("Example 2: Somatic variant with COSMIC support")
    info2 = {"gnomAD_AF": 0.0001, "COSMIC_CNT": 15}
    print_database_summary(info2)
    
    # Example 3: RNA editing site
    print("\n" + "=" * 80)
    print("Example 3: RNA editing site in REDIportal")
    info3 = {"REDIportal": True, "gnomAD_AF": 0.0}
    print_database_summary(info3)
    
    # Example 4: Novel variant (no database support)
    print("\n" + "=" * 80)
    print("Example 4: Novel variant (no database support)")
    info4 = {"gnomAD_AF": 0.0001, "COSMIC_CNT": 0}
    print_database_summary(info4)
    
    print("\n" + "=" * 80)
