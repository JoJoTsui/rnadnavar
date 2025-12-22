#!/usr/bin/env python3
"""
Annotation Debugging Utilities

Standardized debugging framework for annotation modules with parameter tracking,
threshold logging, and progress monitoring. Simplifies troubleshooting of annotation
pipeline issues.

This module provides structured logging for all annotation operations, parameter
summaries for validation, and performance tracking.

USAGE:
    from annotation_debug import AnnotationDebugger
    
    debugger = AnnotationDebugger(module_name="cosmic_annotator")
    debugger.log_parameters(cosmic_vcf=vcf_path, threshold=5)
    debugger.log_threshold_application("COSMIC_CNT", value=10, threshold=5, passed=True)
    debugger.print_summary()

Author: COSMIC/gnomAD/REDIportal Enhancement Pipeline
Date: 2025-12-22
Version: 1.0
"""

import sys
import time
import logging
from typing import Dict, Any, List, Tuple, Optional
from datetime import datetime
from pathlib import Path
import json

try:
    from annotation_config import (
        COSMIC_THRESHOLDS,
        GNOMAD_THRESHOLDS,
        REDIPORTAL_THRESHOLDS,
        CLASSIFICATION_THRESHOLDS,
    )
except ImportError:
    # Fallback if annotation_config not available
    COSMIC_THRESHOLDS = {"recurrence_minimum": 5}
    GNOMAD_THRESHOLDS = {"germline_frequency": 0.001}
    REDIPORTAL_THRESHOLDS = {"min_rna_support": 2}
    CLASSIFICATION_THRESHOLDS = {}


# ============================================================================
# LOGGING CONFIGURATION
# ============================================================================

class ColoredFormatter(logging.Formatter):
    """Custom formatter with color support for different log levels."""
    
    # ANSI color codes
    COLORS = {
        "DEBUG": "\033[36m",      # Cyan
        "INFO": "\033[32m",       # Green
        "WARNING": "\033[33m",    # Yellow
        "ERROR": "\033[31m",      # Red
        "CRITICAL": "\033[35m",   # Magenta
        "RESET": "\033[0m",       # Reset
    }
    
    def format(self, record):
        """Format log record with colors."""
        levelname = record.levelname
        if levelname in self.COLORS:
            record.levelname = f"{self.COLORS[levelname]}{levelname}{self.COLORS['RESET']}"
        return super().format(record)


def setup_logger(
    name: str,
    level: int = logging.INFO,
    log_file: Optional[str] = None,
    use_color: bool = True
) -> logging.Logger:
    """
    Set up a logger with file and console handlers.
    
    Args:
        name: Logger name
        level: Logging level (DEBUG, INFO, WARNING, ERROR)
        log_file: Optional file path for logging
        use_color: Whether to use colored output
    
    Returns:
        Configured logger instance
    """
    logger = logging.getLogger(name)
    logger.setLevel(level)
    
    # Console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(level)
    
    if use_color:
        formatter = ColoredFormatter(
            fmt="[%(levelname)s] %(asctime)s - %(name)s - %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S"
        )
    else:
        formatter = logging.Formatter(
            fmt="[%(levelname)s] %(asctime)s - %(name)s - %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S"
        )
    
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)
    
    # File handler (if specified)
    if log_file:
        try:
            file_handler = logging.FileHandler(log_file)
            file_handler.setLevel(level)
            formatter = logging.Formatter(
                fmt="[%(levelname)s] %(asctime)s - %(name)s - %(message)s",
                datefmt="%Y-%m-%d %H:%M:%S"
            )
            file_handler.setFormatter(formatter)
            logger.addHandler(file_handler)
        except IOError as e:
            logger.warning(f"Could not create log file {log_file}: {e}")
    
    return logger


# ============================================================================
# ANNOTATION DEBUGGER CLASS
# ============================================================================

class AnnotationDebugger:
    """Main debugging class for annotation pipeline operations."""
    
    def __init__(
        self,
        module_name: str,
        log_level: int = logging.INFO,
        log_file: Optional[str] = None,
        use_color: bool = True
    ):
        """
        Initialize annotation debugger.
        
        Args:
            module_name: Name of annotation module (cosmic, gnomad, rediportal, etc.)
            log_level: Logging level
            log_file: Optional log file path
            use_color: Whether to use colored output
        """
        self.module_name = module_name
        self.logger = setup_logger(
            f"AnnotationDebug.{module_name}",
            level=log_level,
            log_file=log_file,
            use_color=use_color
        )
        
        # Debug tracking dictionaries
        self.parameters: Dict[str, Any] = {}
        self.thresholds_applied: List[Dict[str, Any]] = []
        self.performance_metrics: Dict[str, Any] = {}
        self.variant_stats: Dict[str, int] = {}
        self.warnings: List[str] = []
        self.errors: List[str] = []
        
        # Timing
        self.start_time = time.time()
        self.logger.debug(f"AnnotationDebugger initialized for module: {module_name}")
    
    # ========================================================================
    # PARAMETER LOGGING
    # ========================================================================
    
    def log_parameters(self, **kwargs) -> None:
        """
        Log annotation parameters for debugging.
        
        Args:
            **kwargs: Parameter name-value pairs
        
        Example:
            debugger.log_parameters(
                cosmic_vcf="/path/to/cosmic.vcf.gz",
                threshold=5,
                sample_id="SAMPLE001"
            )
        """
        for key, value in kwargs.items():
            self.parameters[key] = value
            # Mask file paths if too long
            display_value = value
            if isinstance(value, str) and len(value) > 80:
                display_value = f"...{value[-77:]}"
            self.logger.debug(f"Parameter: {key} = {display_value}")
    
    def get_parameters(self) -> Dict[str, Any]:
        """Get all logged parameters."""
        return self.parameters.copy()
    
    # ========================================================================
    # THRESHOLD LOGGING
    # ========================================================================
    
    def log_threshold_application(
        self,
        threshold_name: str,
        value: float,
        threshold: float,
        passed: bool,
        variant_id: str = None,
        extra_info: Dict[str, Any] = None
    ) -> None:
        """
        Log threshold application for individual variants.
        
        Args:
            threshold_name: Name of threshold (e.g., "COSMIC_CNT")
            value: Actual value from variant
            threshold: Threshold value being compared against
            passed: Whether variant passed threshold
            variant_id: Optional variant identifier
            extra_info: Additional context information
        
        Example:
            debugger.log_threshold_application(
                threshold_name="COSMIC_CNT",
                value=10,
                threshold=5,
                passed=True,
                variant_id="chr1:1000-1001"
            )
        """
        status = "PASS" if passed else "FAIL"
        
        entry = {
            "threshold": threshold_name,
            "value": value,
            "threshold_cutoff": threshold,
            "status": status,
            "variant_id": variant_id,
            "timestamp": datetime.now().isoformat(),
            "extra_info": extra_info or {}
        }
        
        self.thresholds_applied.append(entry)
        
        log_msg = f"Threshold [{threshold_name}]: {value} vs {threshold} → {status}"
        if variant_id:
            log_msg += f" ({variant_id})"
        
        if passed:
            self.logger.debug(log_msg)
        else:
            self.logger.debug(log_msg)
    
    def get_threshold_statistics(self) -> Dict[str, Any]:
        """Get statistics about threshold applications."""
        if not self.thresholds_applied:
            return {}
        
        passed = sum(1 for t in self.thresholds_applied if t["status"] == "PASS")
        failed = sum(1 for t in self.thresholds_applied if t["status"] == "FAIL")
        
        return {
            "total_applications": len(self.thresholds_applied),
            "passed": passed,
            "failed": failed,
            "pass_rate": passed / len(self.thresholds_applied) if self.thresholds_applied else 0,
        }
    
    # ========================================================================
    # VARIANT STATISTICS
    # ========================================================================
    
    def increment_variant_count(self, category: str, count: int = 1) -> None:
        """
        Increment variant count in a category.
        
        Args:
            category: Variant category (e.g., "cosmic_hit", "gnomad_common", "rediportal_match")
            count: Number to increment by
        """
        if category not in self.variant_stats:
            self.variant_stats[category] = 0
        self.variant_stats[category] += count
        self.logger.debug(f"Variant category [{category}]: {self.variant_stats[category]}")
    
    def get_variant_statistics(self) -> Dict[str, int]:
        """Get all variant statistics."""
        return self.variant_stats.copy()
    
    # ========================================================================
    # PERFORMANCE TRACKING
    # ========================================================================
    
    def track_performance(
        self,
        metric_name: str,
        value: float,
        unit: str = "s",
        context: Dict[str, Any] = None
    ) -> None:
        """
        Track performance metrics.
        
        Args:
            metric_name: Name of metric (e.g., "annotation_time")
            value: Metric value
            unit: Unit of measurement
            context: Additional context
        """
        entry = {
            "metric": metric_name,
            "value": value,
            "unit": unit,
            "context": context or {},
            "timestamp": datetime.now().isoformat(),
        }
        
        if metric_name not in self.performance_metrics:
            self.performance_metrics[metric_name] = []
        
        self.performance_metrics[metric_name].append(entry)
        
        self.logger.debug(f"Performance [{metric_name}]: {value} {unit}")
    
    def get_performance_summary(self) -> Dict[str, Dict[str, float]]:
        """Get performance summary statistics."""
        summary = {}
        
        for metric_name, entries in self.performance_metrics.items():
            values = [e["value"] for e in entries]
            summary[metric_name] = {
                "count": len(values),
                "total": sum(values),
                "average": sum(values) / len(values) if values else 0,
                "min": min(values) if values else 0,
                "max": max(values) if values else 0,
            }
        
        return summary
    
    # ========================================================================
    # WARNING AND ERROR TRACKING
    # ========================================================================
    
    def log_warning(self, message: str, variant_id: str = None, context: Dict[str, Any] = None) -> None:
        """
        Log a warning message.
        
        Args:
            message: Warning message
            variant_id: Optional variant identifier
            context: Additional context
        """
        full_msg = message
        if variant_id:
            full_msg = f"{message} ({variant_id})"
        
        self.warnings.append({
            "message": message,
            "variant_id": variant_id,
            "context": context or {},
            "timestamp": datetime.now().isoformat(),
        })
        
        self.logger.warning(full_msg)
    
    def log_error(self, message: str, exception: Exception = None, context: Dict[str, Any] = None) -> None:
        """
        Log an error message.
        
        Args:
            message: Error message
            exception: Optional exception object
            context: Additional context
        """
        self.errors.append({
            "message": message,
            "exception": str(exception) if exception else None,
            "context": context or {},
            "timestamp": datetime.now().isoformat(),
        })
        
        if exception:
            self.logger.error(f"{message}: {exception}")
        else:
            self.logger.error(message)
    
    def get_warnings(self) -> List[Dict[str, Any]]:
        """Get all logged warnings."""
        return self.warnings.copy()
    
    def get_errors(self) -> List[Dict[str, Any]]:
        """Get all logged errors."""
        return self.errors.copy()
    
    # ========================================================================
    # SUMMARY AND REPORTING
    # ========================================================================
    
    def print_parameter_summary(self) -> None:
        """Print all logged parameters in a formatted table."""
        print("\n" + "=" * 80)
        print(f"ANNOTATION PARAMETERS - {self.module_name}")
        print("=" * 80)
        
        if not self.parameters:
            print("No parameters logged")
        else:
            for key, value in self.parameters.items():
                display_value = value
                if isinstance(value, str) and len(str(value)) > 60:
                    display_value = f"...{str(value)[-57:]}"
                print(f"  {key:<30} = {display_value}")
        
        print("=" * 80 + "\n")
    
    def print_threshold_summary(self) -> None:
        """Print threshold application summary."""
        stats = self.get_threshold_statistics()
        
        if not stats:
            return
        
        print("\n" + "=" * 80)
        print(f"THRESHOLD APPLICATION SUMMARY - {self.module_name}")
        print("=" * 80)
        print(f"Total Applications: {stats['total_applications']}")
        print(f"  Passed: {stats['passed']}")
        print(f"  Failed: {stats['failed']}")
        print(f"  Pass Rate: {stats['pass_rate']:.1%}")
        print("=" * 80 + "\n")
    
    def print_variant_summary(self) -> None:
        """Print variant statistics summary."""
        if not self.variant_stats:
            return
        
        print("\n" + "=" * 80)
        print(f"VARIANT STATISTICS - {self.module_name}")
        print("=" * 80)
        
        total_variants = sum(self.variant_stats.values())
        print(f"Total Variants: {total_variants}")
        
        for category, count in sorted(self.variant_stats.items(), key=lambda x: x[1], reverse=True):
            percentage = (count / total_variants * 100) if total_variants > 0 else 0
            print(f"  {category:<30} : {count:>6} ({percentage:>5.1f}%)")
        
        print("=" * 80 + "\n")
    
    def print_performance_summary(self) -> None:
        """Print performance metrics summary."""
        perf_summary = self.get_performance_summary()
        
        if not perf_summary:
            return
        
        print("\n" + "=" * 80)
        print(f"PERFORMANCE METRICS - {self.module_name}")
        print("=" * 80)
        
        for metric_name, stats in perf_summary.items():
            print(f"\n{metric_name}:")
            print(f"  Count:   {stats['count']}")
            print(f"  Total:   {stats['total']:.4f} s")
            print(f"  Average: {stats['average']:.4f} s")
            print(f"  Min:     {stats['min']:.4f} s")
            print(f"  Max:     {stats['max']:.4f} s")
        
        print("\n" + "=" * 80 + "\n")
    
    def print_issues_summary(self) -> None:
        """Print warnings and errors summary."""
        if not self.warnings and not self.errors:
            return
        
        print("\n" + "=" * 80)
        print(f"ISSUES SUMMARY - {self.module_name}")
        print("=" * 80)
        
        if self.warnings:
            print(f"\nWarnings ({len(self.warnings)}):")
            for warning in self.warnings[:10]:  # Show first 10
                print(f"  • {warning['message']}")
            if len(self.warnings) > 10:
                print(f"  ... and {len(self.warnings) - 10} more warnings")
        
        if self.errors:
            print(f"\nErrors ({len(self.errors)}):")
            for error in self.errors[:10]:  # Show first 10
                print(f"  • {error['message']}")
            if len(self.errors) > 10:
                print(f"  ... and {len(self.errors) - 10} more errors")
        
        print("\n" + "=" * 80 + "\n")
    
    def print_summary(self) -> None:
        """Print comprehensive debug summary."""
        elapsed_time = time.time() - self.start_time
        
        print("\n" + "=" * 80)
        print(f"ANNOTATION DEBUG SUMMARY - {self.module_name}")
        print(f"Elapsed Time: {elapsed_time:.2f} seconds")
        print("=" * 80)
        
        self.print_parameter_summary()
        self.print_threshold_summary()
        self.print_variant_summary()
        self.print_performance_summary()
        self.print_issues_summary()
        
        print("\n" + "=" * 80)
        print("END DEBUG SUMMARY")
        print("=" * 80 + "\n")
    
    # ========================================================================
    # JSON EXPORT FOR ANALYSIS
    # ========================================================================
    
    def export_to_json(self, output_file: str) -> None:
        """
        Export all debug information to JSON file.
        
        Args:
            output_file: Path to output JSON file
        """
        export_data = {
            "module_name": self.module_name,
            "elapsed_time_seconds": time.time() - self.start_time,
            "timestamp": datetime.now().isoformat(),
            "parameters": self.parameters,
            "threshold_applications": self.thresholds_applied,
            "variant_statistics": self.variant_stats,
            "performance_metrics": self.performance_metrics,
            "warnings": self.warnings,
            "errors": self.errors,
            "summary": {
                "thresholds": self.get_threshold_statistics(),
                "performance": self.get_performance_summary(),
            }
        }
        
        try:
            with open(output_file, "w") as f:
                json.dump(export_data, f, indent=2, default=str)
            self.logger.info(f"Debug information exported to {output_file}")
        except IOError as e:
            self.logger.error(f"Failed to export debug information: {e}")


# ============================================================================
# CONVENIENCE FUNCTION FOR QUICK SETUP
# ============================================================================

def create_debugger(
    module_name: str,
    log_file: Optional[str] = None,
    log_level: int = logging.INFO
) -> AnnotationDebugger:
    """
    Create and return a configured AnnotationDebugger instance.
    
    Args:
        module_name: Name of annotation module
        log_file: Optional log file path
        log_level: Logging level
    
    Returns:
        Configured AnnotationDebugger instance
    
    Example:
        debugger = create_debugger("cosmic_annotator", log_file="debug.log")
    """
    return AnnotationDebugger(
        module_name=module_name,
        log_level=log_level,
        log_file=log_file,
        use_color=True
    )


if __name__ == "__main__":
    # Example usage
    debugger = AnnotationDebugger("example_module")
    
    # Log parameters
    debugger.log_parameters(
        cosmic_vcf="/path/to/cosmic.vcf.gz",
        threshold=5,
        output_dir="/path/to/output"
    )
    
    # Log some threshold applications
    for i in range(5):
        debugger.log_threshold_application(
            threshold_name="COSMIC_CNT",
            value=10 + i,
            threshold=5,
            passed=True,
            variant_id=f"chr1:{1000+i*100}-{1001+i*100}"
        )
    
    # Log some variant statistics
    debugger.increment_variant_count("cosmic_hit", 5)
    debugger.increment_variant_count("no_match", 3)
    
    # Log performance
    debugger.track_performance("annotation_time", 2.345, unit="s")
    debugger.track_performance("annotation_time", 1.234, unit="s")
    
    # Print summary
    debugger.print_summary()
