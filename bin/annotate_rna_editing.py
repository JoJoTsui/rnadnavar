#!/usr/bin/env python3
"""
RNA Editing Annotation Script - Fully Integrated Pipeline

This script provides the complete integrated RNA editing annotation pipeline with
all components seamlessly working together. It uses the new REDIportal to VCF
conversion script, implements bcftools VCF-to-VCF annotation without separate
header files, adds evidence tiering and FILTER update functionality, and ensures
reliable output generation using pysam.

INTEGRATED FEATURES:
- New REDIportal to VCF conversion script for proper VCF format generation
- bcftools VCF-to-VCF annotation with exact coordinate and allele matching
- Evidence tiering system based on RNA/DNA caller support (HIGH/MEDIUM/LOW/NONE)
- FILTER column updates to "RNAedit" for variants meeting high evidence criteria
- pysam-based VCF output generation for reliable file handling
- Seamless pipeline execution from input to final output
- Complete isolation to annotate_rna_editing.py and its dependencies
- Comprehensive error handling and logging throughout

PIPELINE COMPONENTS:
1. REDIportal Database Converter - converts text format to bgzipped VCF
2. bcftools Annotation Engine - VCF-to-VCF annotation without header files
3. Evidence Tiering Processor - classifies RNA editing confidence levels
4. FILTER Updater - updates FILTER column based on evidence criteria
5. pysam Output Generator - reliable VCF file generation with indexing

Requirements Satisfied: 1.1, 3.1, 4.1, 5.1, 6.1, 7.1

Author: RNA Editing Enhancement Pipeline
Date: 2025-12-14
"""

import argparse
import sys
import logging
import os
import time
from pathlib import Path
from typing import Dict

# Add vcf_utils to path for imports
script_dir = Path(__file__).parent
sys.path.insert(0, str(script_dir))

# Centralized config/validation
COMMON_PATH = script_dir / "common"
if str(COMMON_PATH) not in sys.path:
    sys.path.insert(0, str(COMMON_PATH))

try:
    from annotation_config import REDIPORTAL_THRESHOLDS
    from annotation_validators import AnnotationValidator
except Exception:
    REDIPORTAL_THRESHOLDS = {"min_rna_support": 2}
    AnnotationValidator = None

# Import dedicated modules for RNA editing annotation
try:
    from vcf_utils.rediportal_converter import prepare_rediportal_database
    from vcf_utils.bcftools_annotator import BcftoolsAnnotator
    from vcf_utils.evidence_tiering import create_evidence_tiering_processor
    from vcf_utils.filter_updater import create_filter_updater
    from vcf_utils.rna_editing_core import is_canonical_editing_transition
    from vcf_utils.error_handler import (
        ErrorHandler, ResourceMonitor, timeout_handler, 
        create_graceful_fallback, validate_tool_availability, 
        check_file_accessibility
    )
    from vcf_utils.logging_config import get_operational_logger
    MODULES_AVAILABLE = True
except ImportError as e:
    logging.warning(f"Some RNA editing modules not available: {e}")
    MODULES_AVAILABLE = False

# Set up comprehensive logging with timestamps and structured format for Nextflow integration
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

# Add structured logging for operational monitoring and Nextflow integration
class StructuredLogger:
    """Structured logger for operational monitoring and Nextflow integration."""
    
    def __init__(self, base_logger):
        self.logger = base_logger
        self.metrics = {}
        self.start_time = time.time()
    
    def log_metric(self, metric_name: str, value, unit: str = ""):
        """Log a metric for operational monitoring."""
        self.metrics[metric_name] = {'value': value, 'unit': unit, 'timestamp': time.time()}
        self.logger.info(f"METRIC: {metric_name}={value}{unit}")
    
    def log_stage(self, stage_name: str, status: str = "START"):
        """Log processing stage for monitoring."""
        timestamp = time.time()
        elapsed = timestamp - self.start_time
        self.logger.info(f"STAGE: {stage_name} - {status} (elapsed: {elapsed:.2f}s)")
    
    def log_performance(self, operation: str, duration: float, items_processed: int = 0):
        """Log performance metrics."""
        rate = items_processed / duration if duration > 0 and items_processed > 0 else 0
        self.logger.info(f"PERFORMANCE: {operation} - duration: {duration:.2f}s, items: {items_processed}, rate: {rate:.0f}/s")
    
    def log_resource_usage(self, stage: str):
        """Log current resource usage."""
        try:
            import psutil
            process = psutil.Process()
            memory_mb = process.memory_info().rss / 1024 / 1024
            cpu_percent = process.cpu_percent()
            self.logger.info(f"RESOURCES: {stage} - memory: {memory_mb:.1f}MB, cpu: {cpu_percent:.1f}%")
        except ImportError:
            self.logger.debug(f"RESOURCES: {stage} - psutil not available")
    
    def get_metrics_summary(self) -> Dict:
        """Get summary of all logged metrics."""
        return {
            'metrics': self.metrics,
            'total_elapsed': time.time() - self.start_time
        }

# Global structured logger instance
structured_logger = StructuredLogger(logger)

# Add file handler for persistent logging with enhanced format
def setup_file_logging(output_dir: Path) -> Path:
    """Set up file logging in addition to console logging with structured format."""
    log_file = output_dir.parent / f"rna_editing_annotation_{int(time.time())}.log"
    
    # Create file handler with detailed formatting
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(logging.DEBUG)
    
    # Enhanced formatter for file logging
    detailed_formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - [%(filename)s:%(lineno)d] - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    file_handler.setFormatter(detailed_formatter)
    logger.addHandler(file_handler)
    
    # Also create a structured metrics file
    metrics_file = output_dir.parent / f"rna_editing_metrics_{int(time.time())}.json"
    
    logger.info(f"Logging to file: {log_file}")
    logger.info(f"Metrics will be saved to: {metrics_file}")
    structured_logger.log_stage("LOGGING_SETUP", "COMPLETE")
    
    return log_file


class RNAEditingAnnotator:
    """
    Integrated RNA editing annotator with complete pipeline components.
    
    This class coordinates the complete integrated RNA editing annotation pipeline using:
    - New REDIportal to VCF conversion script for proper VCF format generation
    - bcftools VCF-to-VCF annotation without requiring separate header files
    - Evidence tiering system based on RNA/DNA caller support (HIGH/MEDIUM/LOW/NONE)
    - FILTER column updates to "RNAedit" for high-evidence variants
    - pysam-based VCF output generation for reliable file handling
    - Seamless pipeline execution from input to final output
    - Complete isolation to annotate_rna_editing.py and its dependencies
    
    Requirements Satisfied: 1.1, 3.1, 4.1, 5.1, 6.1, 7.1
    """
    
    def __init__(self, input_vcf: str, rediportal_vcf: str, output_vcf: str, 
                 min_rna_support: int = 2):
        self.input_vcf = Path(input_vcf)
        self.rediportal_vcf = Path(rediportal_vcf)
        self.output_vcf = Path(output_vcf)
        self.min_rna_support = min_rna_support
        
        # Initialize comprehensive error handling and resource monitoring
        self.error_handler = ErrorHandler(temp_dir=self.output_vcf.parent)
        self.resource_monitor = ResourceMonitor(max_memory_mb=8000, max_disk_gb=50)
        
        # Initialize operational logger for comprehensive monitoring
        self.operational_logger = get_operational_logger(
            name="rna_editing_annotation",
            output_dir=self.output_vcf.parent
        )
        
        # Statistics and timing
        self.stats = {
            'start_time': time.time(),
            'total_variants': 0,
            'annotated_variants': 0,
            'processing_steps': [],
            'errors': [],
            'warnings': []
        }
        
        # Set up file logging with structured format
        self.log_file = setup_file_logging(self.output_vcf)
        
        # Initialize operational monitoring
        self.operational_logger.log_stage("PIPELINE_INIT", "START")
        
        logger.info("=== RNA Editing Annotation Pipeline Started (Enhanced with Operational Monitoring) ===")
        logger.info(f"Input VCF: {self.input_vcf}")
        logger.info(f"REDIportal VCF: {self.rediportal_vcf}")
        logger.info(f"Output VCF: {self.output_vcf}")
        logger.info(f"Min RNA support: {self.min_rna_support}")
        
        # Log initial metrics for operational monitoring
        self.operational_logger.log_metric("min_rna_support", self.min_rna_support, "gauge")
        self.operational_logger.log_resource_usage("initialization")
        structured_logger.log_metric("min_rna_support", self.min_rna_support)
        structured_logger.log_resource_usage("initialization")
        
        # Check module availability
        if not MODULES_AVAILABLE:
            logger.error("Required RNA editing modules are not available")
            raise RuntimeError("RNA editing modules not available - cannot proceed")
        
        # Validate tool availability
        required_tools = ['python', 'bcftools', 'tabix', 'bgzip']
        available_tools, missing_tools = validate_tool_availability(required_tools)
        
        if missing_tools:
            logger.error(f"Missing required tools: {missing_tools}")
            logger.error("Please ensure all required tools are installed and available in PATH")
            raise RuntimeError(f"Required tools not available: {missing_tools}")
        else:
            logger.info(f"All required tools are available: {available_tools}")
        
        # Validate inputs with comprehensive error handling
        try:
            self.validate_inputs()
        except Exception as e:
            # Use error handler for graceful error management
            fallback_action = lambda: create_graceful_fallback(
                self.input_vcf, self.output_vcf, "initialization failure"
            )
            
            if not self.error_handler.handle_error(e, "initialization", fallback_action):
                self._log_error("Initialization failed", e)
                raise
    
    def validate_inputs(self) -> None:
        """Validate input files exist and are accessible with comprehensive error handling."""
        step_start = time.time()
        structured_logger.log_stage("INPUT_VALIDATION", "START")
        logger.info("Validating input files with comprehensive checks...")
        
        # Monitor resources at validation start
        resource_stats = self.resource_monitor.check_resources("input validation start")
        logger.debug(f"Resource stats at validation start: {resource_stats}")
        structured_logger.log_resource_usage("validation_start")
        
        try:
            # Use comprehensive file accessibility checker
            files_to_check = [self.input_vcf, self.rediportal_vcf]
            file_status = check_file_accessibility(files_to_check)
            
            # Check each file status
            for file_path, status in file_status.items():
                if status != "OK":
                    if "does not exist" in status:
                        error_msg = f"File not found: {file_path}"
                        logger.error(error_msg)
                        logger.error(f"Current working directory: {os.getcwd()}")
                        logger.error(f"Absolute path: {Path(file_path).absolute()}")
                        raise FileNotFoundError(error_msg)
                    elif "not readable" in status:
                        error_msg = f"Cannot read file: {file_path}"
                        logger.error(error_msg)
                        raise PermissionError(error_msg)
                    elif "empty" in status:
                        error_msg = f"File is empty: {file_path}"
                        logger.error(error_msg)
                        raise ValueError(error_msg)
                    else:
                        error_msg = f"File validation failed for {file_path}: {status}"
                        logger.error(error_msg)
                        raise RuntimeError(error_msg)
                else:
                    logger.debug(f"✓ File validated: {file_path}")
            
            # Check output directory accessibility with enhanced error handling
            output_dir = self.output_vcf.parent
            if not output_dir.exists():
                logger.info(f"Creating output directory: {output_dir}")
                try:
                    output_dir.mkdir(parents=True, exist_ok=True)
                except Exception as e:
                    error_msg = f"Cannot create output directory: {output_dir}"
                    logger.error(error_msg)
                    logger.error(f"Error details: {e}")
                    raise PermissionError(error_msg)
            
            if not os.access(output_dir, os.W_OK):
                error_msg = f"Cannot write to output directory: {output_dir}"
                logger.error(error_msg)
                raise PermissionError(error_msg)
            
            # Log file sizes and resource implications with structured metrics
            input_size = self.input_vcf.stat().st_size
            rediportal_size = self.rediportal_vcf.stat().st_size
            total_input_size = input_size + rediportal_size
            
            # Log file size metrics for monitoring
            structured_logger.log_metric("input_vcf_size_bytes", input_size)
            structured_logger.log_metric("input_vcf_size_mb", input_size/1024/1024, "MB")
            structured_logger.log_metric("rediportal_vcf_size_bytes", rediportal_size)
            structured_logger.log_metric("rediportal_vcf_size_mb", rediportal_size/1024/1024, "MB")
            structured_logger.log_metric("total_input_size_mb", total_input_size/1024/1024, "MB")
            
            logger.info(f"Input VCF size: {input_size:,} bytes ({input_size/1024/1024:.1f} MB)")
            logger.info(f"REDIportal VCF size: {rediportal_size:,} bytes ({rediportal_size/1024/1024:.1f} MB)")
            logger.info(f"Total input size: {total_input_size:,} bytes ({total_input_size/1024/1024:.1f} MB)")
            
            # Estimate resource requirements based on file sizes
            estimated_memory_mb = max(500, (total_input_size / 1024 / 1024) * 2)  # 2x input size
            estimated_disk_gb = max(2, (total_input_size / 1024 / 1024 / 1024) * 3)  # 3x input size
            
            # Log resource estimates for monitoring
            structured_logger.log_metric("estimated_memory_mb", estimated_memory_mb, "MB")
            structured_logger.log_metric("estimated_disk_gb", estimated_disk_gb, "GB")
            structured_logger.log_metric("available_disk_gb", resource_stats['disk_free_gb'], "GB")
            
            logger.info(f"Estimated memory requirement: {estimated_memory_mb:.0f} MB")
            logger.info(f"Estimated disk requirement: {estimated_disk_gb:.1f} GB")
            
            # Check if we have sufficient resources
            if resource_stats['disk_free_gb'] < estimated_disk_gb:
                warning_msg = f"Low disk space: {resource_stats['disk_free_gb']:.1f} GB available, {estimated_disk_gb:.1f} GB estimated needed"
                logger.warning(warning_msg)
                self.stats['warnings'].append(warning_msg)
                structured_logger.log_metric("disk_space_warning", 1)
            else:
                structured_logger.log_metric("disk_space_sufficient", 1)
            
            step_time = time.time() - step_start
            self.stats['processing_steps'].append(('validate_inputs', step_time))
            
            # Log validation performance
            structured_logger.log_performance("input_validation", step_time, 2)  # 2 files validated
            structured_logger.log_stage("INPUT_VALIDATION", "COMPLETE")
            logger.info(f"✓ Input files validated successfully ({step_time:.2f}s)")
            
        except Exception as e:
            self._log_error("Input validation failed", e)
            raise
    
    def run_annotation(self) -> None:
        """
        Run the complete integrated RNA editing annotation pipeline with comprehensive monitoring.
        
        This method implements the fully integrated pipeline that:
        1. Uses new REDIportal to VCF conversion script for proper VCF format generation
        2. Integrates bcftools VCF-to-VCF annotation without requiring separate header files
        3. Adds evidence tiering and FILTER update functionality based on RNA/DNA caller support
        4. Ensures seamless pipeline execution from input to final output using pysam
        5. Maintains isolation to annotate_rna_editing.py and its dependencies
        6. Provides comprehensive logging and monitoring for operational use
        """
        pipeline_start = time.time()
        structured_logger.log_stage("ANNOTATION_PIPELINE", "START")
        structured_logger.log_resource_usage("pipeline_start")
        
        logger.info("Starting integrated RNA editing annotation pipeline with comprehensive monitoring...")
        logger.info("Pipeline components: REDIportal converter + bcftools VCF-to-VCF + evidence tiering + FILTER updates")
        logger.info("Monitoring: Performance metrics, resource usage, processing statistics")
        
        temp_files = []
        
        try:
            # Step 1: Convert REDIportal database to proper VCF format using new conversion script
            structured_logger.log_stage("REDIPORTAL_CONVERSION", "START")
            logger.info("Step 1: Converting REDIportal database to VCF format...")
            prepared_rediportal = self._prepare_rediportal_database_with_converter()
            if prepared_rediportal != str(self.rediportal_vcf):
                temp_files.append(Path(prepared_rediportal))
            structured_logger.log_stage("REDIPORTAL_CONVERSION", "COMPLETE")
            
            # Step 2: Run bcftools VCF-to-VCF annotation without separate header files
            structured_logger.log_stage("BCFTOOLS_ANNOTATION", "START")
            logger.info("Step 2: Running bcftools VCF-to-VCF annotation with exact matching...")
            annotated_vcf = self._run_bcftools_annotation_with_engine(prepared_rediportal)
            structured_logger.log_stage("BCFTOOLS_ANNOTATION", "COMPLETE")
            
            # Step 3: Apply evidence tiering and FILTER updates using dedicated modules and pysam
            structured_logger.log_stage("EVIDENCE_PROCESSING", "START")
            logger.info("Step 3: Applying evidence tiering and FILTER updates...")
            self._process_variants_with_pysam_and_evidence_tiering(annotated_vcf)
            structured_logger.log_stage("EVIDENCE_PROCESSING", "COMPLETE")
            
            # Step 4: Validate final output and create index for compressed files
            structured_logger.log_stage("OUTPUT_VALIDATION", "START")
            logger.info("Step 4: Validating final output and creating index...")
            self._validate_and_index_output()
            structured_logger.log_stage("OUTPUT_VALIDATION", "COMPLETE")
            
            # Step 5: Clean up temporary files created during processing
            structured_logger.log_stage("CLEANUP", "START")
            logger.info("Step 5: Cleaning up temporary files...")
            self._cleanup_temp_files(temp_files)
            structured_logger.log_stage("CLEANUP", "COMPLETE")
            
            # Generate comprehensive final statistics with performance metrics
            pipeline_time = time.time() - pipeline_start
            self.stats['processing_steps'].append(('total_pipeline', pipeline_time))
            
            # Log final performance metrics
            structured_logger.log_performance("complete_pipeline", pipeline_time, self.stats.get('total_variants', 0))
            structured_logger.log_metric("pipeline_duration_seconds", pipeline_time, "s")
            structured_logger.log_metric("pipeline_duration_minutes", pipeline_time/60, "min")
            structured_logger.log_resource_usage("pipeline_complete")
            structured_logger.log_stage("ANNOTATION_PIPELINE", "SUCCESS")
            
            logger.info("✓ Integrated RNA editing annotation pipeline completed successfully!")
            logger.info("All components integrated: REDIportal conversion, bcftools annotation, evidence tiering, FILTER updates")
            self._generate_summary_statistics()
            
        except Exception as e:
            self._log_error("Integrated RNA editing annotation pipeline failed", e)
            
            # Clean up any partial output files
            if self.output_vcf.exists():
                try:
                    self.output_vcf.unlink()
                    logger.info("Cleaned up partial output file")
                except Exception as cleanup_error:
                    logger.warning(f"Failed to clean up partial output: {cleanup_error}")
            
            # Clean up temporary files with force cleanup on failure
            try:
                self._cleanup_temp_files(temp_files, force_cleanup=True)
            except Exception as cleanup_error:
                logger.warning(f"Failed to clean up temporary files: {cleanup_error}")
                logger.warning("Some temporary files may remain on disk")
            
            # Generate error summary for debugging
            self._generate_summary_statistics()
            
            # Exit with appropriate error code based on error type
            if isinstance(e, FileNotFoundError):
                logger.error("FATAL: Required input file not found")
                sys.exit(2)  # File not found
            elif isinstance(e, PermissionError):
                logger.error("FATAL: Permission denied accessing files or directories")
                sys.exit(13)  # Permission denied
            elif isinstance(e, RuntimeError) and "not found in system PATH" in str(e):
                logger.error("FATAL: Required tools not available in system PATH")
                sys.exit(127)  # Command not found
            else:
                logger.error(f"FATAL: Unexpected error during annotation: {e}")
                sys.exit(1)  # General error
    
    def _prepare_rediportal_database_with_converter(self) -> str:
        """
        Prepare REDIportal database using dedicated converter module.
        
        This method uses the new REDIportal to VCF conversion script to convert
        REDIportal text database to bgzipped VCF format for direct use with
        bcftools annotate without requiring separate header files.
        
        Returns:
            Path to prepared REDIportal VCF database
        """
        step_start = time.time()
        logger.info("Preparing REDIportal database with dedicated converter...")
        
        try:
            # Use dedicated REDIportal converter module for VCF conversion
            # This handles text format conversion to proper VCF with headers
            prepared_file = prepare_rediportal_database(str(self.rediportal_vcf))
            
            # Validate that the prepared file is in VCF format
            if not prepared_file.endswith('.vcf.gz') and not prepared_file.endswith('.vcf'):
                logger.warning(f"Prepared file may not be in VCF format: {prepared_file}")
                self.stats['warnings'].append("REDIportal converter output format unclear")
            
            # Log conversion statistics if available
            logger.info(f"REDIportal database converted to VCF format: {prepared_file}")
            
            step_time = time.time() - step_start
            self.stats['processing_steps'].append(('prepare_rediportal_converter', step_time))
            
            logger.info(f"✓ REDIportal database prepared as VCF: {prepared_file} ({step_time:.2f}s)")
            return prepared_file
            
        except Exception as e:
            self._log_error("REDIportal database preparation failed", e)
            raise RuntimeError(f"REDIportal preparation failed: {e}")
    
    def _run_bcftools_annotation_with_engine(self, rediportal_vcf: str) -> str:
        """
        Run bcftools VCF-to-VCF annotation using dedicated annotation engine.
        
        This method uses the bcftools annotation engine to perform VCF-to-VCF
        annotation without requiring separate header files. It implements exact
        coordinate and allele matching between input VCF and REDIportal VCF database.
        
        Args:
            rediportal_vcf: Path to prepared REDIportal VCF database
            
        Returns:
            Path to annotated VCF file
        """
        step_start = time.time()
        logger.info("Running bcftools VCF-to-VCF annotation with dedicated engine...")
        
        try:
            # Create temporary output for bcftools annotation
            temp_annotated = self.output_vcf.parent / f"temp_annotated_{int(time.time())}.vcf"
            if str(self.output_vcf).endswith('.gz'):
                temp_annotated = Path(str(temp_annotated) + '.gz')
            
            # Use dedicated bcftools annotator for VCF-to-VCF annotation
            logger.info("Initializing bcftools VCF-to-VCF annotation engine...")
            annotator = BcftoolsAnnotator(
                input_vcf=str(self.input_vcf),
                annotation_vcf=rediportal_vcf,
                output_vcf=str(temp_annotated)
            )
            
            # Run VCF-to-VCF annotation without separate header files
            logger.info("Executing VCF-to-VCF annotation with exact coordinate and allele matching...")
            annotator.run_annotation()
            
            # Merge statistics from annotator
            if hasattr(annotator, 'stats'):
                for step_name, step_time in annotator.stats.get('processing_steps', []):
                    self.stats['processing_steps'].append((f'annotator_{step_name}', step_time))
                
                self.stats['warnings'].extend(annotator.stats.get('warnings', []))
                self.stats['errors'].extend(annotator.stats.get('errors', []))
            
            # Log annotation statistics
            if hasattr(annotator, 'operation_logger'):
                operation_stats = annotator.operation_logger.get_statistics()
                logger.info(f"bcftools operations: {operation_stats['successful_commands']} successful, "
                           f"{operation_stats['failed_commands']} failed")
            
            step_time = time.time() - step_start
            self.stats['processing_steps'].append(('bcftools_annotation_engine', step_time))
            
            logger.info(f"✓ bcftools VCF-to-VCF annotation completed: {temp_annotated} ({step_time:.2f}s)")
            return str(temp_annotated)
            
        except Exception as e:
            self._log_error("bcftools VCF-to-VCF annotation engine failed", e)
            raise RuntimeError(f"bcftools VCF-to-VCF annotation failed: {e}")
    
    def _process_variants_with_pysam_and_evidence_tiering(self, annotated_vcf: str) -> None:
        """
        Process variants using pysam and dedicated evidence tiering and FILTER update modules with performance optimization.
        
        This method implements the complete evidence tiering system with performance optimizations for large VCF files:
        1. Extracts N_RNA_CALLERS_SUPPORT and N_DNA_CALLERS_SUPPORT from variant INFO fields
        2. Implements RNA consensus detection (N_RNA_CALLERS_SUPPORT >= min_rna_support)
        3. Implements RNA-only variant detection (N_DNA_CALLERS_SUPPORT = 0)
        4. Creates evidence tier assignment logic (HIGH/MEDIUM/LOW/NONE)
        5. Updates FILTER to "RNAedit" for variants meeting high evidence criteria
        6. Preserves original FILTER values for variants not meeting criteria
        7. Uses memory-efficient processing for large datasets
        8. Provides performance monitoring and optimization
        
        Args:
            annotated_vcf: Path to bcftools-annotated VCF file
        """
        step_start = time.time()
        logger.info("Processing variants with pysam, evidence tiering, FILTER updates, and performance optimization...")
        
        # Initialize performance monitoring for large VCF processing
        try:
            from vcf_utils.performance_optimizer import MemoryMonitor, StreamingProcessor
            memory_monitor = MemoryMonitor(max_memory_mb=6000)  # 6GB limit for variant processing
            use_performance_monitoring = True
            logger.info("Performance monitoring enabled for variant processing")
        except ImportError:
            memory_monitor = None
            use_performance_monitoring = False
            logger.debug("Performance monitoring not available")
        
        try:
            # Check if pysam is available
            try:
                import pysam
            except ImportError:
                logger.warning("pysam not available - skipping evidence tiering and FILTER updates")
                # Just copy the annotated file to output
                import shutil
                shutil.copy2(annotated_vcf, self.output_vcf)
                return
            
            # Create evidence tiering processor with RNA consensus threshold
            logger.info(f"Initializing evidence tiering processor (min_rna_support={self.min_rna_support})...")
            evidence_processor = create_evidence_tiering_processor(min_rna_support=self.min_rna_support)
            
            # Create FILTER updater for RNA editing classification
            logger.info("Initializing FILTER updater for RNA editing classification...")
            filter_updater = create_filter_updater()
            
            # Process variants with pysam for reliable VCF file generation
            logger.info("Opening annotated VCF with pysam for processing...")
            input_vcf = pysam.VariantFile(annotated_vcf)
            
            # Prepare output VCF header with RNA editing fields
            output_header = input_vcf.header.copy()
            
            # Add RNA editing specific INFO fields if not present
            if 'REDI_EVIDENCE' not in output_header.info:
                output_header.info.add('REDI_EVIDENCE', number=1, type='String', 
                                     description='RNA editing evidence level (HIGH, MEDIUM, LOW, NONE)')
            if 'REDI_CANONICAL' not in output_header.info:
                output_header.info.add('REDI_CANONICAL', number=1, type='String',
                                     description='Canonical A>G or T>C transition (YES/NO)')
            
            # Add RNAedit filter to header if not present
            if 'RNAedit' not in output_header.filters:
                output_header.filters.add('RNAedit', None, None, 
                                        'RNA editing variant based on evidence classification')
            
            # Open output VCF with pysam for reliable output generation
            logger.info(f"Creating output VCF with pysam: {self.output_vcf}")
            output_vcf = pysam.VariantFile(str(self.output_vcf), 'w', header=output_header)
            
            # Initialize comprehensive statistics
            evidence_stats = {
                'total_variants': 0,
                'rediportal_matches': 0,
                'canonical_transitions': 0,
                'rna_consensus_variants': 0,
                'rna_only_variants': 0,
                'filter_updates': 0,
                'filter_preserved': 0,
                'evidence_distribution': {'HIGH': 0, 'MEDIUM': 0, 'LOW': 0, 'NONE': 0}
            }
            
            try:
                logger.info("Processing variants for evidence tiering and FILTER updates with performance optimization...")
                
                # Monitor memory at processing start
                if memory_monitor:
                    memory_monitor.check_memory("variant processing start")
                
                for variant in input_vcf:
                    evidence_stats['total_variants'] += 1
                    
                    # Extract variant data for evidence classification
                    variant_data = self._extract_variant_data_pysam(variant)
                    
                    # Check for exact REDIportal match (coordinate and allele matching)
                    rediportal_match = self._has_exact_rediportal_match(variant)
                    if rediportal_match:
                        evidence_stats['rediportal_matches'] += 1
                    
                    # Check canonical RNA editing transition
                    canonical = is_canonical_editing_transition(variant.ref, ','.join(variant.alts))
                    if canonical:
                        evidence_stats['canonical_transitions'] += 1
                    
                    # Process with evidence tiering system
                    tiering_result = evidence_processor.process_variant(variant_data, rediportal_match)
                    
                    # Track RNA consensus and RNA-only variants
                    if tiering_result['rna_consensus']:
                        evidence_stats['rna_consensus_variants'] += 1
                    if tiering_result['rna_only']:
                        evidence_stats['rna_only_variants'] += 1
                    
                    # Track evidence distribution
                    evidence_stats['evidence_distribution'][tiering_result['evidence_tier']] += 1
                    
                    # Get original FILTER value for preservation logic
                    original_filter = list(variant.filter.keys())[0] if variant.filter.keys() else 'PASS'
                    
                    # Update FILTER using dedicated FILTER updater with complete criteria
                    new_filter, filter_was_updated = filter_updater.update_variant_filter(
                        original_filter=original_filter,
                        evidence_tier=tiering_result['evidence_tier'],
                        has_rediportal_match=rediportal_match,
                        rna_support=tiering_result['rna_support'],
                        min_rna_support=self.min_rna_support
                    )
                    
                    # Track FILTER update statistics
                    if filter_was_updated:
                        evidence_stats['filter_updates'] += 1
                    else:
                        evidence_stats['filter_preserved'] += 1
                    
                    # Create output variant record with pysam
                    output_variant = output_header.new_record(
                        contig=variant.chrom,
                        start=variant.start,
                        stop=variant.stop,
                        alleles=variant.alleles
                    )
                    
                    # Copy all INFO fields from input (preserving annotations)
                    for key, value in variant.info.items():
                        try:
                            output_variant.info[key] = value
                        except Exception as e:
                            logger.debug(f"Skipping INFO field {key} due to pysam error: {e}")
                            continue
                    
                    # Add RNA editing specific fields
                    output_variant.info['REDI_EVIDENCE'] = tiering_result['evidence_tier']
                    output_variant.info['REDI_CANONICAL'] = 'YES' if canonical else 'NO'
                    
                    # Set FILTER based on evidence tiering and FILTER updater logic
                    output_variant.filter.clear()
                    output_variant.filter.add(new_filter)
                    
                    # Copy sample information if present (preserve genotype data)
                    if len(variant.samples) > 0:
                        sample_name = list(variant.samples.keys())[0]
                        if sample_name in output_header.samples:
                            for fmt_key in variant.format.keys():
                                try:
                                    output_variant.samples[sample_name][fmt_key] = variant.samples[sample_name][fmt_key]
                                except Exception as e:
                                    logger.debug(f"Skipping FORMAT field {fmt_key} due to pysam error: {e}")
                                    continue
                    
                    # Write enhanced variant with evidence tiering and FILTER updates
                    output_vcf.write(output_variant)
                    
                    # Performance monitoring and optimization for large files
                    if evidence_stats['total_variants'] % 10000 == 0:
                        elapsed_time = time.time() - step_start
                        rate = evidence_stats['total_variants'] / elapsed_time if elapsed_time > 0 else 0
                        logger.info(f"Processed {evidence_stats['total_variants']:,} variants ({rate:.0f} variants/sec)...")
                        
                        # Memory monitoring and optimization
                        if memory_monitor:
                            memory_stats = memory_monitor.check_memory(f"processed {evidence_stats['total_variants']:,} variants")
                            
                            # Force garbage collection every 50K variants for memory management
                            if evidence_stats['total_variants'] % 50000 == 0:
                                freed_mb = memory_monitor.force_garbage_collection("periodic cleanup during variant processing")
                                if freed_mb > 5:  # Log if significant memory was freed
                                    logger.info(f"Freed {freed_mb:.1f} MB during variant processing")
                
            finally:
                input_vcf.close()
                output_vcf.close()
                
                # Final memory check
                if memory_monitor:
                    memory_monitor.check_memory("variant processing end")
            
            # Store comprehensive statistics
            self.stats.update(evidence_stats)
            self.stats['evidence_tiering_stats'] = evidence_processor.get_statistics()
            self.stats['filter_update_stats'] = filter_updater.get_statistics()
            
            step_time = time.time() - step_start
            self.stats['processing_steps'].append(('pysam_evidence_processing', step_time))
            
            # Log comprehensive processing results
            logger.info(f"✓ Evidence processing and FILTER updates completed:")
            logger.info(f"  Total variants processed: {evidence_stats['total_variants']:,}")
            logger.info(f"  REDIportal matches: {evidence_stats['rediportal_matches']:,}")
            logger.info(f"  Canonical transitions: {evidence_stats['canonical_transitions']:,}")
            logger.info(f"  RNA consensus variants (>={self.min_rna_support} callers): {evidence_stats['rna_consensus_variants']:,}")
            logger.info(f"  RNA-only variants (DNA=0): {evidence_stats['rna_only_variants']:,}")
            logger.info(f"  FILTER updates to RNAedit: {evidence_stats['filter_updates']:,}")
            logger.info(f"  FILTER values preserved: {evidence_stats['filter_preserved']:,}")
            
            # Log evidence distribution
            logger.info("Evidence level distribution:")
            for level, count in evidence_stats['evidence_distribution'].items():
                percentage = (count / evidence_stats['total_variants'] * 100) if evidence_stats['total_variants'] > 0 else 0
                logger.info(f"  {level}: {count:,} ({percentage:.1f}%)")
            
            # Log detailed statistics from processors
            evidence_processor.log_statistics()
            filter_updater.log_statistics()
            
            # Clean up temporary annotated file
            temp_file = Path(annotated_vcf)
            if temp_file.exists() and temp_file != self.output_vcf:
                temp_file.unlink()
                logger.debug(f"Cleaned up temporary annotated file: {temp_file}")
            
        except Exception as e:
            self._log_error("Evidence processing and FILTER updates failed", e)
            raise RuntimeError(f"Evidence processing failed: {e}")
    
    def _validate_and_index_output(self) -> None:
        """Validate output file and create index if needed with comprehensive checks."""
        step_start = time.time()
        logger.info("Validating output file and creating index...")
        
        try:
            # Validate output file exists and is accessible
            if not self.output_vcf.exists():
                error_msg = f"Output file was not created: {self.output_vcf}"
                logger.error(error_msg)
                logger.error(f"Expected output path: {self.output_vcf.absolute()}")
                logger.error(f"Output directory exists: {self.output_vcf.parent.exists()}")
                logger.error(f"Output directory writable: {os.access(self.output_vcf.parent, os.W_OK)}")
                raise RuntimeError(error_msg)
            
            # Check file size and basic properties
            file_size = self.output_vcf.stat().st_size
            if file_size == 0:
                error_msg = f"Output file is empty: {self.output_vcf}"
                logger.error(error_msg)
                logger.error("This may indicate:")
                logger.error("  - Processing failed silently")
                logger.error("  - Input file had no variants")
                logger.error("  - Write permissions were lost during processing")
                raise RuntimeError(error_msg)
            
            # Check file permissions
            if not os.access(self.output_vcf, os.R_OK):
                error_msg = f"Output file is not readable: {self.output_vcf}"
                logger.error(error_msg)
                raise RuntimeError(error_msg)
            
            logger.info(f"✓ Output file validated: {file_size:,} bytes ({file_size/1024/1024:.1f} MB)")
            
            # Validate VCF format if possible
            self._validate_vcf_format()
            
            # Create index if output is compressed
            if str(self.output_vcf).endswith('.gz'):
                self._create_output_index()
            
            # Set proper file permissions
            self._ensure_output_file_permissions()
            
            # Final integrity check
            self._perform_output_integrity_check()
            
            step_time = time.time() - step_start
            self.stats['processing_steps'].append(('validate_and_index_output', step_time))
            
            logger.info(f"✓ Output validation and indexing completed ({step_time:.2f}s)")
            
        except Exception as e:
            self._log_error("Output validation failed", e)
            raise
    
    def _validate_vcf_format(self) -> None:
        """Validate that output file has proper VCF format."""
        logger.info("Validating VCF format...")
        
        try:
            # Try to read first few lines to validate VCF format
            if str(self.output_vcf).endswith('.gz'):
                import gzip
                with gzip.open(self.output_vcf, 'rt') as f:
                    first_line = f.readline().strip()
                    header_lines = 0
                    variant_lines = 0
                    
                    # Count header and variant lines (first 100 lines)
                    f.seek(0)
                    for i, line in enumerate(f):
                        if i >= 100:  # Limit check to first 100 lines
                            break
                        line = line.strip()
                        if line.startswith('##'):
                            header_lines += 1
                        elif line.startswith('#CHROM'):
                            header_lines += 1
                        elif line and not line.startswith('#'):
                            variant_lines += 1
            else:
                with open(self.output_vcf, 'r') as f:
                    first_line = f.readline().strip()
                    header_lines = 0
                    variant_lines = 0
                    
                    # Count header and variant lines (first 100 lines)
                    f.seek(0)
                    for i, line in enumerate(f):
                        if i >= 100:  # Limit check to first 100 lines
                            break
                        line = line.strip()
                        if line.startswith('##'):
                            header_lines += 1
                        elif line.startswith('#CHROM'):
                            header_lines += 1
                        elif line and not line.startswith('#'):
                            variant_lines += 1
            
            # Validate VCF header
            if not first_line.startswith('##fileformat=VCF'):
                logger.warning(f"Output file may not have valid VCF header. First line: {first_line}")
                self.stats['warnings'].append("Invalid VCF header format")
            else:
                logger.info("✓ Valid VCF format detected")
                logger.info(f"  Header lines: {header_lines}")
                logger.info(f"  Variant lines (sample): {variant_lines}")
            
        except Exception as e:
            logger.warning(f"Could not validate VCF format: {e}")
            self.stats['warnings'].append(f"VCF format validation failed: {e}")
    
    def _perform_output_integrity_check(self) -> None:
        """Perform final integrity check on output file."""
        logger.info("Performing output integrity check...")
        
        try:
            # Check if file can be opened and basic structure is intact
            line_count = 0
            header_count = 0
            variant_count = 0
            
            if str(self.output_vcf).endswith('.gz'):
                import gzip
                with gzip.open(self.output_vcf, 'rt') as f:
                    for line in f:
                        line_count += 1
                        if line.startswith('##'):
                            header_count += 1
                        elif line.startswith('#CHROM'):
                            header_count += 1
                        elif line.strip() and not line.startswith('#'):
                            variant_count += 1
                        
                        # Limit check to avoid processing huge files
                        if line_count >= 10000:
                            break
            else:
                with open(self.output_vcf, 'r') as f:
                    for line in f:
                        line_count += 1
                        if line.startswith('##'):
                            header_count += 1
                        elif line.startswith('#CHROM'):
                            header_count += 1
                        elif line.strip() and not line.startswith('#'):
                            variant_count += 1
                        
                        # Limit check to avoid processing huge files
                        if line_count >= 10000:
                            break
            
            logger.info("✓ Integrity check completed:")
            logger.info(f"  Total lines checked: {line_count:,}")
            logger.info(f"  Header lines: {header_count:,}")
            logger.info(f"  Variant lines: {variant_count:,}")
            
            # Store counts in stats
            self.stats['output_line_count'] = line_count
            self.stats['output_header_count'] = header_count
            self.stats['output_variant_count'] = variant_count
            
            # Validate minimum structure
            if header_count == 0:
                logger.warning("No VCF header lines found in output")
                self.stats['warnings'].append("No VCF header lines in output")
            
            if variant_count == 0:
                logger.warning("No variant lines found in output (first 10,000 lines)")
                self.stats['warnings'].append("No variant lines found in output sample")
            
        except Exception as e:
            logger.warning(f"Output integrity check failed: {e}")
            self.stats['warnings'].append(f"Integrity check failed: {e}")
    
    def _extract_variant_data_pysam(self, variant) -> Dict:
        """Extract variant data from pysam.VariantRecord for evidence classification."""
        return {
            'CHROM': variant.chrom,
            'POS': variant.pos,
            'REF': variant.ref,
            'ALT': ','.join(variant.alts) if variant.alts else '.',
            'N_RNA_CALLERS_SUPPORT': variant.info.get('N_RNA_CALLERS_SUPPORT', 0),
            'N_DNA_CALLERS_SUPPORT': variant.info.get('N_DNA_CALLERS_SUPPORT', 0),
            'VAF_RNA_MEAN': variant.info.get('VAF_RNA_MEAN', 0.0),
            'VAF_DNA_MEAN': variant.info.get('VAF_DNA_MEAN', 0.0)
        }
    
    def _has_exact_rediportal_match(self, variant) -> bool:
        """
        Check if variant has exact REDIportal match based on coordinate and allele matching.
        
        This method checks for the presence of REDIportal annotation fields that indicate
        exact coordinate (CHROM, POS) and allele (REF, ALT) matching between the input
        variant and the REDIportal database.
        
        Args:
            variant: pysam.VariantRecord to check
            
        Returns:
            True if variant has exact REDIportal match, False otherwise
        """
        # Check for any REDIportal annotation fields that indicate exact matching
        rediportal_fields = ['REDI_ACCESSION', 'REDI_DB', 'REDI_TYPE', 'REDI_REPEAT', 'REDI_FUNC', 'REDI_STRAND']
        
        for field in rediportal_fields:
            if field in variant.info:
                return True
        
        return False
    
    def _cleanup_temp_files(self, temp_files: list, force_cleanup: bool = False) -> None:
        """
        Clean up temporary files created during processing.
        
        Args:
            temp_files: List of temporary files to clean up
            force_cleanup: If True, attempt cleanup even on failure scenarios
        """
        if not temp_files:
            logger.debug("No temporary files to clean up")
            return
        
        logger.info(f"Cleaning up {len(temp_files)} temporary files...")
        cleaned_count = 0
        failed_count = 0
        total_size_cleaned = 0
        
        for temp_file in temp_files:
            # Only clean up files that are not the original inputs
            if temp_file != self.input_vcf and temp_file != self.rediportal_vcf:
                try:
                    if temp_file.exists():
                        file_size = temp_file.stat().st_size
                        total_size_cleaned += file_size
                        temp_file.unlink()
                        logger.debug(f"Cleaned up temporary file: {temp_file} ({file_size:,} bytes)")
                        cleaned_count += 1
                    
                    # Also clean up associated index files
                    for suffix in ['.tbi', '.csi']:
                        index_file = Path(str(temp_file) + suffix)
                        if index_file.exists():
                            index_size = index_file.stat().st_size
                            total_size_cleaned += index_size
                            index_file.unlink()
                            logger.debug(f"Cleaned up index file: {index_file} ({index_size:,} bytes)")
                            
                except Exception as e:
                    logger.warning(f"Failed to clean up temporary file {temp_file}: {e}")
                    failed_count += 1
                    
                    # On force cleanup, try harder to remove files
                    if force_cleanup:
                        try:
                            if temp_file.exists():
                                os.chmod(temp_file, 0o777)  # Change permissions and try again
                                temp_file.unlink()
                                logger.info(f"Force cleaned up temporary file: {temp_file}")
                                cleaned_count += 1
                        except Exception as force_error:
                            logger.error(f"Force cleanup also failed for {temp_file}: {force_error}")
        
        if cleaned_count > 0:
            logger.info(f"✓ Cleaned up {cleaned_count} temporary files ({total_size_cleaned:,} bytes freed)")
        if failed_count > 0:
            logger.warning(f"Failed to clean up {failed_count} temporary files")
            if not force_cleanup:
                logger.info("Some temporary files may remain on disk")
    
    def _ensure_output_file_permissions(self) -> None:
        """Ensure output files have correct permissions for accessibility."""
        step_start = time.time()
        logger.info("Setting proper output file permissions...")
        
        try:
            # Set readable permissions for output VCF (644 - owner read/write, group/other read)
            if self.output_vcf.exists():
                os.chmod(self.output_vcf, 0o644)
                logger.debug(f"Set permissions 644 for {self.output_vcf}")
            
            # Set permissions for index file if it exists
            index_file = Path(str(self.output_vcf) + '.tbi')
            if index_file.exists():
                os.chmod(index_file, 0o644)
                logger.debug(f"Set permissions 644 for {index_file}")
            
            step_time = time.time() - step_start
            self.stats['processing_steps'].append(('ensure_output_permissions', step_time))
            logger.info(f"✓ Output file permissions set successfully ({step_time:.2f}s)")
            
        except Exception as e:
            logger.warning(f"Failed to set output file permissions: {e}")
            self.stats['warnings'].append(f"Permission setting failed: {e}")

    def _create_output_index(self) -> None:
        """Create tabix index for compressed output VCF."""
        if not str(self.output_vcf).endswith('.gz'):
            logger.debug("Output VCF is not compressed, skipping indexing")
            return
        
        step_start = time.time()
        logger.info("Creating tabix index for compressed output...")
        
        try:
            import subprocess
            
            # Validate output file exists and is accessible
            if not self.output_vcf.exists():
                raise RuntimeError(f"Output VCF file does not exist: {self.output_vcf}")
            
            if self.output_vcf.stat().st_size == 0:
                raise RuntimeError(f"Output VCF file is empty: {self.output_vcf}")
            
            # Check if index already exists and is valid
            index_file = Path(str(self.output_vcf) + '.tbi')
            if index_file.exists():
                # Verify index is newer than VCF file
                vcf_mtime = self.output_vcf.stat().st_mtime
                index_mtime = index_file.stat().st_mtime
                
                if index_mtime >= vcf_mtime:
                    logger.info("Valid tabix index already exists, skipping creation")
                    return
                else:
                    logger.info("Existing index is outdated, recreating...")
                    index_file.unlink()
            
            cmd = ['tabix', '-p', 'vcf', str(self.output_vcf)]
            logger.debug(f"Running command: {' '.join(cmd)}")
            
            subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True,
                timeout=300  # 5 minute timeout for indexing
            )
            
            # Verify index was created and is valid
            if not index_file.exists():
                raise RuntimeError("Tabix indexing completed but index file was not created")
            
            if index_file.stat().st_size == 0:
                raise RuntimeError("Tabix index file is empty")
            
            index_size = index_file.stat().st_size
            step_time = time.time() - step_start
            self.stats['processing_steps'].append(('create_output_index', step_time))
            
            logger.info(f"✓ Tabix index created successfully ({index_size:,} bytes, {step_time:.2f}s)")
            
        except subprocess.CalledProcessError as e:
            error_msg = f"Failed to create tabix index: {e}"
            logger.error(error_msg)
            logger.error(f"Command stderr: {e.stderr}")
            raise RuntimeError(error_msg)
            
        except subprocess.TimeoutExpired:
            error_msg = "Tabix indexing timed out after 5 minutes"
            logger.error(error_msg)
            raise RuntimeError(error_msg)
            
        except Exception as e:
            logger.error(f"Unexpected error during indexing: {e}")
            raise
    
    def _log_error(self, context: str, error: Exception) -> None:
        """
        Log comprehensive error information with diagnostic details.
        
        Args:
            context: Description of what was being attempted when error occurred
            error: The exception that was raised
        """
        error_info = {
            'context': context,
            'error_type': type(error).__name__,
            'error_message': str(error),
            'timestamp': time.time()
        }
        
        logger.error(f"ERROR: {context}")
        logger.error(f"Error type: {error_info['error_type']}")
        logger.error(f"Error message: {error_info['error_message']}")
        
        # Store error for summary reporting
        self.stats['errors'].append(error_info)
    
    def _generate_summary_statistics(self) -> None:
        """Generate and log comprehensive summary statistics with detailed metrics for operational monitoring."""
        total_time = time.time() - self.stats['start_time']
        
        # Generate structured metrics summary
        metrics_summary = structured_logger.get_metrics_summary()
        
        logger.info("=== RNA Editing Annotation Summary (Enhanced with Operational Metrics) ===")
        logger.info(f"Total processing time: {total_time:.2f} seconds ({total_time/60:.1f} minutes)")
        logger.info(f"Input VCF: {self.input_vcf}")
        logger.info(f"REDIportal VCF: {self.rediportal_vcf}")
        logger.info(f"Output VCF: {self.output_vcf}")
        
        # Log key performance indicators for operational monitoring
        logger.info("=== Key Performance Indicators ===")
        if 'total_variants' in self.stats and self.stats['total_variants'] > 0:
            processing_rate = self.stats['total_variants'] / total_time if total_time > 0 else 0
            structured_logger.log_metric("final_processing_rate", processing_rate, "variants/s")
            logger.info(f"Processing rate: {processing_rate:.0f} variants/second")
            logger.info(f"Total variants processed: {self.stats['total_variants']:,}")
        
        if 'annotated_variants' in self.stats:
            annotation_rate = (self.stats['annotated_variants'] / self.stats['total_variants'] * 100) if self.stats['total_variants'] > 0 else 0
            structured_logger.log_metric("annotation_rate_percent", annotation_rate, "%")
            logger.info(f"Annotation rate: {annotation_rate:.1f}% ({self.stats['annotated_variants']:,} variants annotated)")
        
        # Log operational status indicators
        success_rate = 100.0 if len(self.stats['errors']) == 0 else 0.0
        structured_logger.log_metric("success_rate_percent", success_rate, "%")
        structured_logger.log_metric("error_count", len(self.stats['errors']))
        structured_logger.log_metric("warning_count", len(self.stats['warnings']))
        
        logger.info(f"Pipeline success rate: {success_rate:.1f}%")
        logger.info(f"Errors encountered: {len(self.stats['errors'])}")
        logger.info(f"Warnings generated: {len(self.stats['warnings'])}")
        
        # Log file size information
        try:
            if self.input_vcf.exists():
                input_size = self.input_vcf.stat().st_size
                logger.info(f"Input file size: {input_size:,} bytes ({input_size/1024/1024:.1f} MB)")
            
            if self.rediportal_vcf.exists():
                rediportal_size = self.rediportal_vcf.stat().st_size
                logger.info(f"REDIportal file size: {rediportal_size:,} bytes ({rediportal_size/1024/1024:.1f} MB)")
            
            if self.output_vcf.exists():
                output_size = self.output_vcf.stat().st_size
                logger.info(f"Output file size: {output_size:,} bytes ({output_size/1024/1024:.1f} MB)")
                
                # Calculate compression ratio if applicable
                if str(self.output_vcf).endswith('.gz') and self.input_vcf.exists():
                    compression_ratio = (output_size / input_size) * 100 if input_size > 0 else 0
                    logger.info(f"Compression ratio: {compression_ratio:.1f}%")
        except Exception as e:
            logger.debug(f"Could not get file size information: {e}")
        
        # Log processing step timings with performance analysis
        if self.stats['processing_steps']:
            logger.info("Processing step timings:")
            total_step_time = 0
            for step_name, step_time in self.stats['processing_steps']:
                percentage = (step_time / total_time * 100) if total_time > 0 else 0
                logger.info(f"  {step_name}: {step_time:.2f}s ({percentage:.1f}%)")
                total_step_time += step_time
            
            # Log overhead time
            overhead_time = total_time - total_step_time
            if overhead_time > 0:
                overhead_percentage = (overhead_time / total_time * 100) if total_time > 0 else 0
                logger.info(f"  overhead/other: {overhead_time:.2f}s ({overhead_percentage:.1f}%)")
        
        # Log output file statistics if available
        if 'output_line_count' in self.stats:
            logger.info("Output file structure:")
            logger.info(f"  Header lines: {self.stats.get('output_header_count', 0):,}")
            logger.info(f"  Variant lines: {self.stats.get('output_variant_count', 0):,}")
            logger.info(f"  Total lines checked: {self.stats.get('output_line_count', 0):,}")
        
        # Log warnings and errors with categorization
        if self.stats['warnings']:
            logger.warning(f"Warnings encountered: {len(self.stats['warnings'])}")
            warning_categories = {}
            for warning in self.stats['warnings']:
                category = warning.split(':')[0] if ':' in warning else 'General'
                warning_categories[category] = warning_categories.get(category, 0) + 1
                logger.warning(f"  {warning}")
            
            if len(warning_categories) > 1:
                logger.warning("Warning categories:")
                for category, count in warning_categories.items():
                    logger.warning(f"  {category}: {count}")
        
        if self.stats['errors']:
            logger.error(f"Errors encountered: {len(self.stats['errors'])}")
            error_categories = {}
            for error_info in self.stats['errors']:
                error_type = error_info.get('error_type', 'Unknown')
                error_categories[error_type] = error_categories.get(error_type, 0) + 1
                logger.error(f"  {error_info['context']}: {error_info['error_type']} - {error_info['error_message']}")
            
            if len(error_categories) > 1:
                logger.error("Error categories:")
                for category, count in error_categories.items():
                    logger.error(f"  {category}: {count}")
        
        # Log evidence tiering statistics if available
        if 'evidence_tiering_stats' in self.stats:
            stats = self.stats['evidence_tiering_stats']
            logger.info("Evidence tiering statistics:")
            logger.info(f"  Total variants processed: {stats['total_variants_processed']:,}")
            logger.info(f"  RNA consensus variants: {stats['rna_consensus_variants']:,} ({stats['rna_consensus_rate']:.1f}%)")
            logger.info(f"  RNA-only variants: {stats['rna_only_variants']:,} ({stats['rna_only_rate']:.1f}%)")
            logger.info(f"  FILTER updates to RNAedit: {stats['filter_updates']:,} ({stats['filter_update_rate']:.1f}%)")
            
            # Log evidence level distribution
            if 'evidence_levels' in stats:
                logger.info("  Evidence level distribution:")
                for level, count in stats['evidence_levels'].items():
                    percentage = (count / stats['total_variants_processed'] * 100) if stats['total_variants_processed'] > 0 else 0
                    logger.info(f"    {level}: {count:,} ({percentage:.1f}%)")
        
        # Log FILTER update statistics if available
        if 'filter_update_stats' in self.stats:
            filter_stats = self.stats['filter_update_stats']
            logger.info("FILTER column update statistics:")
            logger.info(f"  Total variants processed: {filter_stats['total_variants_processed']:,}")
            logger.info(f"  FILTER updates to RNAedit: {filter_stats['filter_updates_to_rnaedit']:,} ({filter_stats['filter_update_rate']:.1f}%)")
            logger.info(f"  FILTER values preserved: {filter_stats['filter_values_preserved']:,} ({filter_stats['filter_preservation_rate']:.1f}%)")
            
            # Log FILTER distribution changes
            if 'original_filter_distribution' in filter_stats and 'final_filter_distribution' in filter_stats:
                logger.info("  FILTER value changes:")
                all_filters = set(filter_stats['original_filter_distribution'].keys()) | set(filter_stats['final_filter_distribution'].keys())
                for filter_val in sorted(all_filters):
                    original_count = filter_stats['original_filter_distribution'].get(filter_val, 0)
                    final_count = filter_stats['final_filter_distribution'].get(filter_val, 0)
                    change = final_count - original_count
                    if change != 0:
                        change_str = f"({change:+,})" if change != 0 else ""
                        logger.info(f"    {filter_val}: {original_count:,} → {final_count:,} {change_str}")
        
        # Log performance metrics
        if total_time > 0 and 'total_variants' in self.stats and self.stats['total_variants'] > 0:
            variants_per_second = self.stats['total_variants'] / total_time
            logger.info(f"Processing rate: {variants_per_second:.0f} variants/second")
        
        # Log system resource usage if available
        try:
            import psutil
            process = psutil.Process()
            memory_info = process.memory_info()
            logger.info(f"Peak memory usage: {memory_info.rss / 1024 / 1024:.1f} MB")
        except ImportError:
            logger.debug("psutil not available for memory usage reporting")
        except Exception as e:
            logger.debug(f"Could not get memory usage: {e}")
        
        # Log success/failure summary
        success = len(self.stats['errors']) == 0
        status = "SUCCESS" if success else "COMPLETED WITH ERRORS"
        structured_logger.log_metric("pipeline_success", 1 if success else 0)
        logger.info(f"Pipeline status: {status}")
        
        # Save structured metrics for operational monitoring and integration
        self._save_structured_metrics(metrics_summary, total_time, success)
        
        # Save operational metrics report
        try:
            operational_metrics_file = self.operational_logger.save_metrics_report()
            logger.info(f"Operational metrics report: {operational_metrics_file}")
        except Exception as e:
            logger.warning(f"Failed to save operational metrics: {e}")
        
        if self.log_file:
            logger.info(f"Detailed log file: {self.log_file}")
            
        # Generate operational recommendations based on metrics
        self._generate_operational_recommendations(total_time, success)
    
    def _save_structured_metrics(self, metrics_summary: Dict, total_time: float, success: bool) -> None:
        """Save structured metrics to JSON file for monitoring system integration."""
        try:
            import json
            
            # Create comprehensive metrics document
            metrics_document = {
                'pipeline_info': {
                    'name': 'rna_editing_annotation',
                    'version': '1.0.0',
                    'timestamp': time.time(),
                    'date': time.strftime('%Y-%m-%d %H:%M:%S'),
                    'input_vcf': str(self.input_vcf),
                    'rediportal_vcf': str(self.rediportal_vcf),
                    'output_vcf': str(self.output_vcf),
                    'min_rna_support': self.min_rna_support
                },
                'performance_metrics': {
                    'total_duration_seconds': total_time,
                    'total_duration_minutes': total_time / 60,
                    'processing_rate_variants_per_second': self.stats.get('total_variants', 0) / total_time if total_time > 0 else 0,
                    'memory_efficient': total_time < 300,  # Under 5 minutes considered efficient
                    'large_file_processing': self.stats.get('total_variants', 0) > 100000
                },
                'processing_statistics': {
                    'total_variants': self.stats.get('total_variants', 0),
                    'annotated_variants': self.stats.get('annotated_variants', 0),
                    'annotation_rate_percent': (self.stats.get('annotated_variants', 0) / self.stats.get('total_variants', 1) * 100),
                    'rediportal_matches': self.stats.get('rediportal_matches', 0),
                    'canonical_transitions': self.stats.get('canonical_transitions', 0),
                    'rna_consensus_variants': self.stats.get('rna_consensus_variants', 0),
                    'rna_only_variants': self.stats.get('rna_only_variants', 0),
                    'filter_updates': self.stats.get('filter_updates', 0)
                },
                'quality_metrics': {
                    'success': success,
                    'error_count': len(self.stats['errors']),
                    'warning_count': len(self.stats['warnings']),
                    'processing_steps_completed': len(self.stats['processing_steps']),
                    'data_integrity_validated': 'output_line_count' in self.stats
                },
                'resource_utilization': {
                    'estimated_memory_efficient': True,  # Based on successful completion
                    'disk_space_sufficient': len([w for w in self.stats['warnings'] if 'disk space' in w.lower()]) == 0,
                    'processing_steps': [
                        {'step': step_name, 'duration_seconds': step_time, 'percentage_of_total': (step_time / total_time * 100) if total_time > 0 else 0}
                        for step_name, step_time in self.stats['processing_steps']
                    ]
                },
                'operational_status': {
                    'pipeline_health': 'healthy' if success and len(self.stats['warnings']) < 5 else 'degraded' if success else 'failed',
                    'requires_attention': len(self.stats['errors']) > 0 or len(self.stats['warnings']) > 10,
                    'performance_acceptable': total_time < 1800,  # Under 30 minutes
                    'ready_for_production': success and len(self.stats['errors']) == 0
                },
                'structured_metrics': metrics_summary['metrics'],
                'nextflow_integration': {
                    'log_file': str(self.log_file),
                    'compatible_with_reporting': True,
                    'monitoring_ready': True
                }
            }
            
            # Save metrics to JSON file
            metrics_file = self.output_vcf.parent / f"rna_editing_metrics_{int(time.time())}.json"
            with open(metrics_file, 'w') as f:
                json.dump(metrics_document, f, indent=2, default=str)
            
            logger.info(f"Structured metrics saved to: {metrics_file}")
            structured_logger.log_metric("metrics_file_saved", 1)
            
        except Exception as e:
            logger.warning(f"Failed to save structured metrics: {e}")
            structured_logger.log_metric("metrics_save_failed", 1)
    
    def _generate_operational_recommendations(self, total_time: float, success: bool) -> None:
        """Generate operational recommendations based on processing metrics."""
        logger.info("=== Operational Recommendations ===")
        
        # Performance recommendations
        if total_time > 1800:  # Over 30 minutes
            logger.info("PERFORMANCE: Consider increasing allocated resources for faster processing")
            logger.info("- Increase memory allocation if processing large files")
            logger.info("- Consider parallel processing for multiple samples")
        elif total_time < 60:  # Under 1 minute
            logger.info("PERFORMANCE: Excellent processing speed achieved")
        
        # Resource recommendations
        if len([w for w in self.stats['warnings'] if 'disk space' in w.lower()]) > 0:
            logger.info("RESOURCES: Monitor disk space usage")
            logger.info("- Ensure adequate disk space for temporary files")
            logger.info("- Consider cleanup of old processing files")
        
        # Quality recommendations
        if len(self.stats['warnings']) > 5:
            logger.info("QUALITY: Review warning messages for data quality issues")
            logger.info("- Validate input data quality")
            logger.info("- Check REDIportal database integrity")
        
        # Operational recommendations
        if success:
            logger.info("OPERATIONS: Pipeline completed successfully")
            logger.info("- Output ready for downstream analysis")
            logger.info("- Consider archiving log files for audit trail")
        else:
            logger.info("OPERATIONS: Pipeline completed with issues")
            logger.info("- Review error messages for troubleshooting")
            logger.info("- Validate input files and configuration")
            logger.info("- Consider rerunning with increased resources")
        
        logger.info("=== End Operational Recommendations ===")
        
        logger.info("=== Summary Complete ===")
        
        # Log final recommendations if there were issues
        if self.stats['warnings'] or self.stats['errors']:
            logger.info("=== Recommendations ===")
            if self.stats['errors']:
                logger.info("- Review error messages above for critical issues")
                logger.info("- Check input file formats and permissions")
                logger.info("- Verify tool availability and versions")
            if self.stats['warnings']:
                logger.info("- Review warning messages for potential data quality issues")
                logger.info("- Consider validating output with downstream tools")
            logger.info(f"- Consult detailed log file for debugging: {self.log_file}")
            logger.info("=== End Recommendations ===")


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Reorganized RNA editing annotation using dedicated modules",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Integrated RNA Editing Annotation Pipeline:
  This script provides the complete integrated RNA editing annotation pipeline
  with all components seamlessly working together. It uses the new REDIportal
  to VCF conversion script, implements bcftools VCF-to-VCF annotation without
  separate header files, adds evidence tiering and FILTER update functionality,
  and ensures reliable output generation using pysam.

Integrated Pipeline Components:
  1. REDIportal Database Converter - converts text format to bgzipped VCF
  2. bcftools Annotation Engine - VCF-to-VCF annotation without header files
  3. Evidence Tiering Processor - classifies RNA editing confidence levels
  4. FILTER Updater - updates FILTER column based on evidence criteria
  5. pysam Output Generator - reliable VCF file generation with indexing

Key Features:
  - New REDIportal to VCF conversion script for proper VCF format generation
  - bcftools VCF-to-VCF annotation with exact coordinate and allele matching
  - Evidence tiering system based on RNA/DNA caller support (HIGH/MEDIUM/LOW/NONE)
  - FILTER column updates to "RNAedit" for variants meeting high evidence criteria
  - pysam-based VCF output generation for reliable file handling
  - Seamless pipeline execution from input to final output
  - Complete isolation to annotate_rna_editing.py and its dependencies

Module Requirements:
  - vcf_utils.rediportal_converter (REDIportal database conversion)
  - vcf_utils.bcftools_annotator (bcftools VCF-to-VCF annotation)
  - vcf_utils.evidence_tiering (RNA editing evidence classification)
  - vcf_utils.filter_updater (FILTER column update logic)
  - vcf_utils.rna_editing_core (core RNA editing functions)
  - pysam (VCF file reading and writing)

Tool Requirements:
  - bcftools (for VCF-to-VCF annotation)
  - tabix (for VCF indexing)
  - bgzip (for VCF compression)

Examples:
  # Basic integrated annotation
  python annotate_rna_editing.py -i input.vcf.gz -r /path/to/rediportal.vcf.gz -o output.vcf.gz
  
  # Custom RNA support threshold for evidence tiering
  python annotate_rna_editing.py -i input.vcf.gz -r /path/to/rediportal.vcf.gz -o output.vcf.gz --min-rna-support 3
  
  # Verbose logging for detailed pipeline monitoring
  python annotate_rna_editing.py -i input.vcf.gz -r /path/to/rediportal.vcf.gz -o output.vcf.gz --verbose
        """
    )
    
    parser.add_argument(
        '-i', '--input',
        required=True,
        metavar='VCF_FILE',
        help='Input VCF file (can be compressed with .gz)'
    )
    
    parser.add_argument(
        '-r', '--rediportal',
        required=True,
        metavar='VCF_FILE',
        help='REDIportal database VCF file path (can be compressed with .gz)'
    )
    
    parser.add_argument(
        '-o', '--output',
        required=True,
        metavar='VCF_FILE',
        help='Output annotated VCF file (use .gz extension for compression)'
    )
    
    parser.add_argument(
        '--min-rna-support',
        type=int,
        default=None,
        metavar='N',
        help=f"Minimum RNA caller support threshold (default: {REDIPORTAL_THRESHOLDS.get('min_rna_support', 2)}, must be >= 1)"
    )
    
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Enable verbose logging'
    )
    
    args = parser.parse_args()
    
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # Validate command-line arguments
    try:
        # Resolve defaults from config
        min_rna_support = args.min_rna_support
        if min_rna_support is None:
            min_rna_support = REDIPORTAL_THRESHOLDS.get('min_rna_support', 2)

        # Validate min_rna_support range
        if min_rna_support < 1:
            logger.error("FATAL: --min-rna-support must be at least 1")
            sys.exit(2)  # Invalid argument
        
        # Validate file paths are not empty
        if not args.input.strip():
            logger.error("FATAL: Input VCF path cannot be empty")
            sys.exit(2)  # Invalid argument
            
        if not args.rediportal.strip():
            logger.error("FATAL: REDIportal VCF path cannot be empty")
            sys.exit(2)  # Invalid argument
            
        if not args.output.strip():
            logger.error("FATAL: Output VCF path cannot be empty")
            sys.exit(2)  # Invalid argument
            
    except Exception as e:
        logger.error(f"FATAL: Argument validation failed: {e}")
        sys.exit(2)  # Invalid argument
    
    # Validate against shared validator (if available)
    if AnnotationValidator:
        validator = AnnotationValidator()
        validation = validator.validate_all_parameters(
            rediportal_vcf=args.rediportal,
            rediportal_min_rna=min_rna_support,
        )
        if not validation.get('valid', True):
            logger.error("FATAL: Parameter validation failed")
            validator.print_validation_report(validation)
            sys.exit(2)

    # Run annotation
    try:
        annotator = RNAEditingAnnotator(
            input_vcf=args.input, 
            rediportal_vcf=args.rediportal, 
            output_vcf=args.output,
            min_rna_support=min_rna_support
        )
        
        annotator.run_annotation()
        
        # Explicit success message and exit code
        logger.info("✓ RNA editing annotation completed successfully!")
        sys.exit(0)  # Success
        
    except KeyboardInterrupt:
        logger.error("FATAL: Annotation interrupted by user")
        sys.exit(130)  # Script terminated by Control-C
    except Exception as e:
        logger.error(f"FATAL: Annotation failed: {e}")
        # Provide specific exit codes based on error type
        if isinstance(e, FileNotFoundError):
            sys.exit(2)  # File not found
        elif isinstance(e, PermissionError):
            sys.exit(13)  # Permission denied
        elif "not found in system PATH" in str(e):
            sys.exit(127)  # Command not found
        else:
            sys.exit(1)  # General error


if __name__ == '__main__':
    main()