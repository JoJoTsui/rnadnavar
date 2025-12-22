#!/usr/bin/env python3
"""
Annotation Parameter Validators

Comprehensive validation framework for annotation parameters with type checking,
range validation, and consistency checks. Prevents invalid parameters from reaching
annotation modules.

This module centralizes parameter validation logic and provides standardized error
messages for debugging.

USAGE:
    from annotation_validators import AnnotationValidator
    
    validator = AnnotationValidator()
    validator.validate_cosmic_parameters(cosmic_vcf, cosmic_threshold)
    validator.validate_gnomad_parameters(gnomad_dir, frequency_threshold)
    validator.validate_rediportal_parameters(rediportal_vcf, min_rna_support)
    validator.validate_classification_parameters(somatic_consensus)

Author: COSMIC/gnomAD/REDIportal Enhancement Pipeline
Date: 2025-12-22
Version: 1.0
"""

import os
import sys
from typing import Dict, Any, List, Tuple, Optional
from pathlib import Path
import gzip

try:
    from annotation_config import (
        COSMIC_THRESHOLDS,
        GNOMAD_THRESHOLDS,
        REDIPORTAL_THRESHOLDS,
        CLASSIFICATION_THRESHOLDS,
    )
except ImportError:
    # Fallback if annotation_config not available
    COSMIC_THRESHOLDS = {"recurrence_minimum": 5, "sample_count_minimum": 2}
    GNOMAD_THRESHOLDS = {"germline_frequency": 0.001, "allele_count_minimum": 2}
    REDIPORTAL_THRESHOLDS = {"min_rna_support": 2, "rna_vaf_minimum": 0.10, "dna_vaf_maximum": 0.05}
    CLASSIFICATION_THRESHOLDS = {"somatic_consensus_threshold": 2}


class ValidationError(Exception):
    """Custom exception for validation failures."""
    pass


class ParameterValidator:
    """Base class for parameter validation with common utility methods."""
    
    @staticmethod
    def is_file_readable(filepath: str) -> bool:
        """Check if file exists and is readable."""
        try:
            path = Path(filepath)
            return path.exists() and os.access(path, os.R_OK)
        except (TypeError, ValueError):
            return False
    
    @staticmethod
    def is_directory_readable(dirpath: str) -> bool:
        """Check if directory exists and is readable."""
        try:
            path = Path(dirpath)
            return path.is_dir() and os.access(path, os.R_OK)
        except (TypeError, ValueError):
            return False
    
    @staticmethod
    def is_vcf_indexed(vcf_path: str) -> bool:
        """Check if VCF file has tabix index."""
        try:
            tbi_path = f"{vcf_path}.tbi"
            csi_path = f"{vcf_path}.csi"
            return os.path.exists(tbi_path) or os.path.exists(csi_path)
        except (TypeError, AttributeError):
            return False
    
    @staticmethod
    def is_vcf_bgzipped(vcf_path: str) -> bool:
        """Check if VCF is bgzip compressed."""
        try:
            with open(vcf_path, "rb") as f:
                return f.read(2) == b"\x1f\x8b"  # gzip magic number
        except (IOError, TypeError):
            return False
    
    @staticmethod
    def in_range(value: float, min_val: float, max_val: float) -> bool:
        """Check if numeric value is in valid range."""
        try:
            return min_val <= float(value) <= max_val
        except (TypeError, ValueError):
            return False
    
    @staticmethod
    def is_positive(value: float) -> bool:
        """Check if value is positive number."""
        try:
            return float(value) > 0
        except (TypeError, ValueError):
            return False


class CosmicValidator(ParameterValidator):
    """Validator for COSMIC annotation parameters."""
    
    def validate_file(self, cosmic_vcf: str) -> Tuple[bool, List[str]]:
        """
        Validate COSMIC VCF file.
        
        Args:
            cosmic_vcf: Path to COSMIC VCF file
        
        Returns:
            (is_valid, list_of_issues)
        """
        issues = []
        
        if not cosmic_vcf:
            issues.append("COSMIC VCF path cannot be empty")
            return False, issues
        
        if not isinstance(cosmic_vcf, str):
            issues.append(f"COSMIC VCF path must be string, got {type(cosmic_vcf)}")
            return False, issues
        
        if not self.is_file_readable(cosmic_vcf):
            issues.append(f"COSMIC VCF file not found or not readable: {cosmic_vcf}")
            return False, issues
        
        if not self.is_vcf_bgzipped(cosmic_vcf):
            issues.append(f"COSMIC VCF must be bgzip compressed (.gz): {cosmic_vcf}")
            return False, issues
        
        if not self.is_vcf_indexed(cosmic_vcf):
            issues.append(f"COSMIC VCF missing tabix index (.tbi or .csi): {cosmic_vcf}")
            return False, issues
        
        return True, issues
    
    def validate_recurrence_threshold(self, threshold: int) -> Tuple[bool, List[str]]:
        """
        Validate recurrence threshold for COSMIC rescue.
        
        Args:
            threshold: Minimum COSMIC recurrence count
        
        Returns:
            (is_valid, list_of_issues)
        """
        issues = []
        
        try:
            threshold = int(threshold)
        except (TypeError, ValueError):
            issues.append(f"Recurrence threshold must be integer, got {type(threshold)}")
            return False, issues
        
        if threshold < 1:
            issues.append(f"Recurrence threshold must be ≥ 1, got {threshold}")
            return False, issues
        
        if threshold > 100:
            issues.append(f"Recurrence threshold {threshold} seems very high (COSMIC_CNT rarely >50)")
        
        default = COSMIC_THRESHOLDS["recurrence_minimum"]
        if threshold != default:
            issues.append(f"Note: Using non-standard threshold {threshold} (default is {default})")
        
        return True, issues
    
    def validate_all_cosmic_parameters(
        self,
        cosmic_vcf: str,
        recurrence_threshold: int = None
    ) -> Dict[str, Any]:
        """
        Validate all COSMIC annotation parameters.
        
        Args:
            cosmic_vcf: Path to COSMIC VCF
            recurrence_threshold: Minimum recurrence count (optional)
        
        Returns:
            Validation results dictionary
        """
        results = {
            "valid": True,
            "issues": [],
            "warnings": [],
        }
        
        # Validate VCF file
        vcf_valid, vcf_issues = self.validate_file(cosmic_vcf)
        if not vcf_valid:
            results["valid"] = False
            results["issues"].extend(vcf_issues)
        else:
            results["warnings"].append(f"✓ COSMIC VCF valid: {cosmic_vcf}")
        
        # Validate threshold if provided
        if recurrence_threshold is not None:
            thresh_valid, thresh_issues = self.validate_recurrence_threshold(recurrence_threshold)
            if not thresh_valid:
                results["valid"] = False
            results["issues"].extend(thresh_issues)
        
        return results


class GnomadValidator(ParameterValidator):
    """Validator for gnomAD annotation parameters."""
    
    def validate_directory(self, gnomad_dir: str) -> Tuple[bool, List[str]]:
        """
        Validate gnomAD database directory.
        
        Args:
            gnomad_dir: Path to gnomAD database directory
        
        Returns:
            (is_valid, list_of_issues)
        """
        issues = []
        
        if not gnomad_dir:
            issues.append("gnomAD directory path cannot be empty")
            return False, issues
        
        if not isinstance(gnomad_dir, str):
            issues.append(f"gnomAD directory path must be string, got {type(gnomad_dir)}")
            return False, issues
        
        if not self.is_directory_readable(gnomad_dir):
            issues.append(f"gnomAD directory not found or not readable: {gnomad_dir}")
            return False, issues
        
        # Check for chromosome VCF files
        has_vcf_files = False
        for i in range(1, 23):
            vcf_file = os.path.join(gnomad_dir, f"chr{i}.vcf.gz")
            if os.path.exists(vcf_file):
                has_vcf_files = True
                break
        
        if not has_vcf_files:
            issues.append(f"No gnomAD chromosome VCF files found in {gnomad_dir}")
            return False, issues
        
        return True, issues
    
    def validate_frequency_threshold(self, threshold: float) -> Tuple[bool, List[str]]:
        """
        Validate germline frequency threshold.
        
        Args:
            threshold: Maximum AF for benign/germline classification
        
        Returns:
            (is_valid, list_of_issues)
        """
        issues = []
        
        try:
            threshold = float(threshold)
        except (TypeError, ValueError):
            issues.append(f"Frequency threshold must be numeric, got {type(threshold)}")
            return False, issues
        
        if not self.in_range(threshold, 0, 1):
            issues.append(f"Frequency threshold must be 0-1, got {threshold}")
            return False, issues
        
        if threshold == 0:
            issues.append("Frequency threshold = 0 will not filter any variants")
        
        if threshold > 0.1:
            issues.append(f"Frequency threshold {threshold} is high (ACMG standard is 0.001 for rare disease)")
        
        default = GNOMAD_THRESHOLDS["germline_frequency"]
        if abs(threshold - default) > 0.0001:
            issues.append(f"Note: Using non-standard threshold {threshold} (default is {default})")
        
        return True, issues
    
    def validate_allele_count_minimum(self, ac_min: int) -> Tuple[bool, List[str]]:
        """
        Validate minimum allele count threshold.
        
        Args:
            ac_min: Minimum allele count
        
        Returns:
            (is_valid, list_of_issues)
        """
        issues = []
        
        try:
            ac_min = int(ac_min)
        except (TypeError, ValueError):
            issues.append(f"Allele count minimum must be integer, got {type(ac_min)}")
            return False, issues
        
        if ac_min < 0:
            issues.append(f"Allele count minimum must be ≥ 0, got {ac_min}")
            return False, issues
        
        if ac_min > 10:
            issues.append(f"Allele count minimum {ac_min} seems high (most variants have AC ≥ 1)")
        
        return True, issues
    
    def validate_all_gnomad_parameters(
        self,
        gnomad_dir: str,
        frequency_threshold: float = None,
        allele_count_minimum: int = None
    ) -> Dict[str, Any]:
        """
        Validate all gnomAD annotation parameters.
        
        Args:
            gnomad_dir: Path to gnomAD database directory
            frequency_threshold: Maximum AF for germline (optional)
            allele_count_minimum: Minimum observation count (optional)
        
        Returns:
            Validation results dictionary
        """
        results = {
            "valid": True,
            "issues": [],
            "warnings": [],
        }
        
        # Validate directory
        dir_valid, dir_issues = self.validate_directory(gnomad_dir)
        if not dir_valid:
            results["valid"] = False
            results["issues"].extend(dir_issues)
        else:
            results["warnings"].append(f"✓ gnomAD directory valid: {gnomad_dir}")
        
        # Validate frequency threshold if provided
        if frequency_threshold is not None:
            freq_valid, freq_issues = self.validate_frequency_threshold(frequency_threshold)
            if not freq_valid:
                results["valid"] = False
            results["issues"].extend(freq_issues)
        
        # Validate AC minimum if provided
        if allele_count_minimum is not None:
            ac_valid, ac_issues = self.validate_allele_count_minimum(allele_count_minimum)
            if not ac_valid:
                results["valid"] = False
            results["issues"].extend(ac_issues)
        
        return results


class RediportalValidator(ParameterValidator):
    """Validator for REDIportal annotation parameters."""
    
    def validate_file(self, rediportal_vcf: str) -> Tuple[bool, List[str]]:
        """
        Validate REDIportal VCF file.
        
        Args:
            rediportal_vcf: Path to REDIportal VCF file
        
        Returns:
            (is_valid, list_of_issues)
        """
        issues = []
        
        if not rediportal_vcf:
            issues.append("REDIportal VCF path cannot be empty")
            return False, issues
        
        if not isinstance(rediportal_vcf, str):
            issues.append(f"REDIportal VCF path must be string, got {type(rediportal_vcf)}")
            return False, issues
        
        if not self.is_file_readable(rediportal_vcf):
            issues.append(f"REDIportal VCF file not found or not readable: {rediportal_vcf}")
            return False, issues
        
        if not self.is_vcf_bgzipped(rediportal_vcf):
            issues.append(f"REDIportal VCF must be bgzip compressed (.gz): {rediportal_vcf}")
            return False, issues
        
        if not self.is_vcf_indexed(rediportal_vcf):
            issues.append(f"REDIportal VCF missing tabix index (.tbi or .csi): {rediportal_vcf}")
            return False, issues
        
        return True, issues
    
    def validate_rna_support_threshold(self, min_rna: int) -> Tuple[bool, List[str]]:
        """
        Validate minimum RNA caller support threshold.
        
        Args:
            min_rna: Minimum number of RNA callers
        
        Returns:
            (is_valid, list_of_issues)
        """
        issues = []
        
        try:
            min_rna = int(min_rna)
        except (TypeError, ValueError):
            issues.append(f"min_rna_support must be integer, got {type(min_rna)}")
            return False, issues
        
        if min_rna < 1:
            issues.append(f"min_rna_support must be ≥ 1, got {min_rna}")
            return False, issues
        
        if min_rna > 3:
            issues.append(f"min_rna_support {min_rna} > 3 (only 3 RNA callers available)")
        
        if min_rna == 1:
            issues.append("min_rna_support = 1 allows single-caller artifacts (consider 2+)")
        
        default = REDIPORTAL_THRESHOLDS["min_rna_support"]
        if min_rna != default:
            issues.append(f"Note: Using non-standard threshold {min_rna} (default is {default})")
        
        return True, issues
    
    def validate_vaf_threshold(self, vaf: float, threshold_type: str = "rna_minimum") -> Tuple[bool, List[str]]:
        """
        Validate VAF threshold (RNA minimum or DNA maximum).
        
        Args:
            vaf: VAF threshold value
            threshold_type: "rna_minimum" or "dna_maximum"
        
        Returns:
            (is_valid, list_of_issues)
        """
        issues = []
        
        try:
            vaf = float(vaf)
        except (TypeError, ValueError):
            issues.append(f"{threshold_type} VAF must be numeric, got {type(vaf)}")
            return False, issues
        
        if not self.in_range(vaf, 0, 1):
            issues.append(f"{threshold_type} VAF must be 0-1, got {vaf}")
            return False, issues
        
        if vaf == 0:
            issues.append(f"{threshold_type} VAF = 0 will include very weak signals")
        
        if vaf > 0.5:
            issues.append(f"{threshold_type} VAF {vaf} seems very high")
        
        return True, issues
    
    def validate_all_rediportal_parameters(
        self,
        rediportal_vcf: str,
        min_rna_support: int = None,
        rna_vaf_minimum: float = None,
        dna_vaf_maximum: float = None
    ) -> Dict[str, Any]:
        """
        Validate all REDIportal annotation parameters.
        
        Args:
            rediportal_vcf: Path to REDIportal VCF
            min_rna_support: Minimum RNA callers (optional)
            rna_vaf_minimum: Minimum RNA VAF (optional)
            dna_vaf_maximum: Maximum DNA VAF (optional)
        
        Returns:
            Validation results dictionary
        """
        results = {
            "valid": True,
            "issues": [],
            "warnings": [],
        }
        
        # Validate VCF file
        vcf_valid, vcf_issues = self.validate_file(rediportal_vcf)
        if not vcf_valid:
            results["valid"] = False
            results["issues"].extend(vcf_issues)
        else:
            results["warnings"].append(f"✓ REDIportal VCF valid: {rediportal_vcf}")
        
        # Validate RNA support threshold if provided
        if min_rna_support is not None:
            rna_valid, rna_issues = self.validate_rna_support_threshold(min_rna_support)
            if not rna_valid:
                results["valid"] = False
            results["issues"].extend(rna_issues)
        
        # Validate RNA VAF minimum if provided
        if rna_vaf_minimum is not None:
            rna_vaf_valid, rna_vaf_issues = self.validate_vaf_threshold(
                rna_vaf_minimum, "rna_minimum"
            )
            if not rna_vaf_valid:
                results["valid"] = False
            results["issues"].extend(rna_vaf_issues)
        
        # Validate DNA VAF maximum if provided
        if dna_vaf_maximum is not None:
            dna_vaf_valid, dna_vaf_issues = self.validate_vaf_threshold(
                dna_vaf_maximum, "dna_maximum"
            )
            if not dna_vaf_valid:
                results["valid"] = False
            results["issues"].extend(dna_vaf_issues)
        
        return results


class AnnotationValidator:
    """Main validator class combining all validators."""
    
    def __init__(self):
        """Initialize all sub-validators."""
        self.cosmic = CosmicValidator()
        self.gnomad = GnomadValidator()
        self.rediportal = RediportalValidator()
    
    def validate_all_parameters(
        self,
        cosmic_vcf: str = None,
        cosmic_threshold: int = None,
        gnomad_dir: str = None,
        gnomad_frequency_threshold: float = None,
        rediportal_vcf: str = None,
        rediportal_min_rna: int = None
    ) -> Dict[str, Any]:
        """
        Validate all annotation parameters at once.
        
        Args:
            cosmic_vcf: COSMIC database VCF
            cosmic_threshold: COSMIC recurrence minimum
            gnomad_dir: gnomAD database directory
            gnomad_frequency_threshold: gnomAD germline AF threshold
            rediportal_vcf: REDIportal database VCF
            rediportal_min_rna: REDIportal minimum RNA caller support
        
        Returns:
            Comprehensive validation results
        """
        all_results = {}
        all_valid = True
        
        if cosmic_vcf or cosmic_threshold is not None:
            cosmic_results = self.cosmic.validate_all_cosmic_parameters(
                cosmic_vcf, cosmic_threshold
            )
            all_results["COSMIC"] = cosmic_results
            if not cosmic_results["valid"]:
                all_valid = False
        
        if gnomad_dir or gnomad_frequency_threshold is not None:
            gnomad_results = self.gnomad.validate_all_gnomad_parameters(
                gnomad_dir, gnomad_frequency_threshold
            )
            all_results["GNOMAD"] = gnomad_results
            if not gnomad_results["valid"]:
                all_valid = False
        
        if rediportal_vcf or rediportal_min_rna is not None:
            rediportal_results = self.rediportal.validate_all_rediportal_parameters(
                rediportal_vcf, rediportal_min_rna
            )
            all_results["REDIPORTAL"] = rediportal_results
            if not rediportal_results["valid"]:
                all_valid = False
        
        return {
            "valid": all_valid,
            "results": all_results,
        }
    
    def print_validation_report(self, validation_results: Dict[str, Any]) -> None:
        """
        Print formatted validation report.
        
        Args:
            validation_results: Results from validation
        """
        print("\n" + "=" * 80)
        print("ANNOTATION PARAMETER VALIDATION REPORT")
        print("=" * 80)
        
        for category, results in validation_results.get("results", {}).items():
            status = "✓ VALID" if results["valid"] else "✗ INVALID"
            print(f"\n{category}: {status}")
            
            if results["warnings"]:
                for warning in results["warnings"]:
                    print(f"  {warning}")
            
            if results["issues"]:
                for issue in results["issues"]:
                    print(f"  ! {issue}")
        
        overall = "✓ ALL PARAMETERS VALID" if validation_results["valid"] else "✗ VALIDATION FAILED"
        print(f"\n{overall}")
        print("=" * 80 + "\n")


if __name__ == "__main__":
    # Example usage
    validator = AnnotationValidator()
    
    results = validator.validate_all_parameters(
        cosmic_threshold=5,
        gnomad_frequency_threshold=0.001,
        rediportal_min_rna=2
    )
    
    validator.print_validation_report(results)
