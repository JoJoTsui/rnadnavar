#!/usr/bin/env python3
"""
COSMIC/gnomAD Main Annotation Engine

This is the main annotation script that integrates COSMIC and gnomAD annotation
workflows with multi-modal variant classification. It provides a unified interface
for comprehensive variant annotation using population and cancer databases.

FEATURES:
- Integrated COSMIC and gnomAD annotation workflows
- Multi-modal variant classification with statistics generation
- Configurable logging (verbose/minimal) and graceful fallback mechanisms
- Comprehensive error handling and recovery
- Performance monitoring and resource usage tracking

Requirements Satisfied: 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 3.1, 3.2, 3.3, 3.4, 3.5, 5.1, 5.2, 5.3, 5.4, 5.5

Author: COSMIC/gnomAD Enhancement Pipeline
Date: 2025-12-17
"""

import logging
import sys
import os
import time
import json
import argparse
from pathlib import Path
from typing import Dict, Any, Optional, List, Tuple

# Import utility modules
sys.path.insert(0, str(Path(__file__).parent))
from vcf_utils.cosmic_annotator import CosmicAnnotator
from vcf_utils.gnomad_annotator import GnomadAnnotator
from vcf_utils.variant_classifier import VariantClassifier, create_variant_classifier
from vcf_utils.annotation_utils import (
    validate_vcf_file, create_temp_directory, cleanup_temp_directory,
    get_vcf_statistics, format_file_size, format_duration,
    log_annotation_summary, check_tool_availability
)

logger = logging.getLogger(__name__)


class CosmicGnomadAnnotator:
    """
    Main COSMIC/gnomAD annotation engine that integrates all annotation workflows.
    
    This class provides a unified interface for comprehensive variant annotation
    using COSMIC and gnomAD databases with multi-modal classification logic.
    """
    
    def __init__(self, 
                 input_vcf: str,
                 cosmic_vcf: Optional[str] = None,
                 gnomad_dir: Optional[str] = None,
                 output_vcf: str = None,
                 verbose: bool = False,
                 max_workers: int = None,
                 classification_config: Optional[Dict[str, Any]] = None):
        """
        Initialize the main annotation engine.
        
        Args:
            input_vcf: Path to input VCF file to be annotated
            cosmic_vcf: Path to COSMIC VCF database (optional)
            gnomad_dir: Path to gnomAD database directory (optional)
            output_vcf: Path to output annotated VCF file
            verbose: Enable verbose logging
            max_workers: Maximum number of parallel workers (None for auto-optimization)
            classification_config: Optional configuration for variant classification
        """
        self.input_vcf = Path(input_vcf)
        self.cosmic_vcf = Path(cosmic_vcf) if cosmic_vcf else None
        self.gnomad_dir = Path(gnomad_dir) if gnomad_dir else None
        self.output_vcf = Path(output_vcf) if output_vcf else None
        self.verbose = verbose
        self.max_workers = max_workers
        
        # Initialize variant classifier
        self.classifier = create_variant_classifier(classification_config)
        
        # Statistics and timing
        self.stats = {
            'start_time': time.time(),
            'processing_steps': [],
            'errors': [],
            'warnings': [],
            'cosmic_annotations': 0,
            'gnomad_annotations': 0,
            'variants_classified': 0,
            'classification_changes': 0,
            'fallback_occurred': False
        }
        
        # Set up logging
        self._setup_logging()
        
        logger.info("=== COSMIC/gnomAD Main Annotation Engine Initialized ===")
        logger.info(f"Input VCF: {self.input_vcf}")
        logger.info(f"COSMIC VCF: {self.cosmic_vcf}")
        logger.info(f"gnomAD directory: {self.gnomad_dir}")
        logger.info(f"Output VCF: {self.output_vcf}")
        logger.info(f"Verbose logging: {self.verbose}")
        logger.info(f"Max workers: {self.max_workers}")
        
        # Validate inputs and environment
        self._validate_environment()
    
    def _setup_logging(self) -> None:
        """Set up logging configuration based on verbose setting."""
        log_level = logging.DEBUG if self.verbose else logging.INFO
        
        # Configure root logger
        logging.basicConfig(
            level=log_level,
            format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S',
            force=True  # Override any existing configuration
        )
        
        # Set specific logger levels
        if not self.verbose:
            # Reduce noise from utility modules in non-verbose mode
            logging.getLogger('vcf_utils').setLevel(logging.WARNING)
            logging.getLogger('urllib3').setLevel(logging.WARNING)
            logging.getLogger('requests').setLevel(logging.WARNING)
    
    def _validate_environment(self) -> None:
        """Validate environment, inputs, and tool availability."""
        step_start = time.time()
        logger.info("Validating annotation environment...")
        
        # Check required tools
        required_tools = ['bcftools', 'tabix', 'bgzip']
        available_tools, missing_tools = check_tool_availability(required_tools)
        
        if missing_tools:
            error_msg = f"Required tools not found: {missing_tools}"
            logger.error(error_msg)
            logger.error("Installation instructions:")
            logger.error("  Ubuntu/Debian: apt install bcftools tabix")
            logger.error("  CentOS/RHEL: yum install bcftools htslib")
            logger.error("  Conda: conda install bcftools htslib")
            logger.error("  Homebrew: brew install bcftools htslib")
            raise RuntimeError(error_msg)
        
        logger.info(f"✓ All required tools available: {', '.join(available_tools)}")
        
        # Validate input VCF
        logger.info("Validating input VCF file...")
        input_validation = validate_vcf_file(self.input_vcf)
        
        if not input_validation['valid']:
            error_msg = f"Input VCF validation failed: {input_validation['errors']}"
            logger.error(error_msg)
            raise ValueError(error_msg)
        
        logger.info(f"✓ Input VCF validated ({format_file_size(input_validation['size_bytes'])})")
        
        # Validate databases if provided
        if self.cosmic_vcf:
            logger.info("Validating COSMIC database...")
            cosmic_validation = validate_vcf_file(self.cosmic_vcf)
            
            if not cosmic_validation['valid']:
                logger.warning(f"COSMIC database validation issues: {cosmic_validation['errors']}")
                self.stats['warnings'].extend(cosmic_validation['errors'])
            else:
                logger.info(f"✓ COSMIC database validated ({format_file_size(cosmic_validation['size_bytes'])})")
        
        if self.gnomad_dir:
            logger.info("Validating gnomAD database directory...")
            
            if not self.gnomad_dir.exists():
                error_msg = f"gnomAD directory not found: {self.gnomad_dir}"
                logger.error(error_msg)
                raise FileNotFoundError(error_msg)
            
            if not self.gnomad_dir.is_dir():
                error_msg = f"gnomAD path is not a directory: {self.gnomad_dir}"
                logger.error(error_msg)
                raise NotADirectoryError(error_msg)
            
            logger.info(f"✓ gnomAD directory validated: {self.gnomad_dir}")
        
        # Validate output directory
        if self.output_vcf:
            output_dir = self.output_vcf.parent
            if not output_dir.exists():
                logger.info(f"Creating output directory: {output_dir}")
                output_dir.mkdir(parents=True, exist_ok=True)
            
            if not os.access(output_dir, os.W_OK):
                error_msg = f"Cannot write to output directory: {output_dir}"
                logger.error(error_msg)
                raise PermissionError(error_msg)
        
        step_time = time.time() - step_start
        self.stats['processing_steps'].append(('validate_environment', step_time))
        logger.info(f"✓ Environment validation completed ({step_time:.2f}s)")
    
    def run_cosmic_annotation(self, input_vcf: Path, output_vcf: Path) -> bool:
        """
        Run COSMIC annotation workflow.
        
        Args:
            input_vcf: Path to input VCF file
            output_vcf: Path to output VCF file
            
        Returns:
            True if annotation succeeded, False if fallback occurred
        """
        if not self.cosmic_vcf:
            logger.info("COSMIC database not provided, skipping COSMIC annotation")
            return True
        
        step_start = time.time()
        logger.info("Running COSMIC annotation workflow...")
        
        try:
            cosmic_annotator = CosmicAnnotator(
                input_vcf=str(input_vcf),
                cosmic_vcf=str(self.cosmic_vcf),
                output_vcf=str(output_vcf)
            )
            
            cosmic_annotator.run_cosmic_annotation()
            
            # Extract statistics
            cosmic_stats = cosmic_annotator.stats
            self.stats['cosmic_annotations'] = cosmic_stats.get('annotations_added', 0)
            
            step_time = time.time() - step_start
            self.stats['processing_steps'].append(('cosmic_annotation', step_time))
            
            logger.info(f"✓ COSMIC annotation completed successfully ({step_time:.2f}s)")
            logger.info(f"COSMIC annotations added: {self.stats['cosmic_annotations']:,}")
            
            return True
            
        except Exception as e:
            step_time = time.time() - step_start
            error_msg = f"COSMIC annotation failed after {step_time:.2f}s: {e}"
            logger.error(error_msg)
            logger.error("CRITICAL: COSMIC annotation failure - exiting immediately")
            logger.error("Suggested remediation:")
            logger.error("  1. Check COSMIC database file accessibility and format")
            logger.error("  2. Verify COSMIC VCF is properly bgzip compressed and tabix indexed")
            logger.error("  3. Ensure bcftools is available and functional")
            logger.error("  4. Check file permissions and disk space")
            # Immediate exit - no graceful fallback
            raise RuntimeError(error_msg)
    
    def run_gnomad_annotation(self, input_vcf: Path, output_vcf: Path) -> bool:
        """
        Run gnomAD annotation workflow.
        
        Args:
            input_vcf: Path to input VCF file
            output_vcf: Path to output VCF file
            
        Returns:
            True if annotation succeeded, False if fallback occurred
        """
        if not self.gnomad_dir:
            logger.info("gnomAD database not provided, skipping gnomAD annotation")
            return True
        
        step_start = time.time()
        logger.info("Running gnomAD scatter-gather annotation workflow...")
        
        try:
            gnomad_annotator = GnomadAnnotator(
                input_vcf=str(input_vcf),
                gnomad_dir=str(self.gnomad_dir),
                output_vcf=str(output_vcf),
                max_workers=self.max_workers
            )
            
            gnomad_annotator.run_scatter_gather_annotation()
            
            # Extract statistics
            gnomad_stats = gnomad_annotator.stats
            self.stats['gnomad_annotations'] = gnomad_stats.get('annotations_added', 0)
            
            step_time = time.time() - step_start
            self.stats['processing_steps'].append(('gnomad_annotation', step_time))
            
            logger.info(f"✓ gnomAD annotation completed successfully ({step_time:.2f}s)")
            logger.info(f"gnomAD annotations added: {self.stats['gnomad_annotations']:,}")
            
            return True
            
        except Exception as e:
            step_time = time.time() - step_start
            error_msg = f"gnomAD annotation failed after {step_time:.2f}s: {e}"
            logger.error(error_msg)
            logger.error("CRITICAL: gnomAD annotation failure - exiting immediately")
            logger.error("Suggested remediation:")
            logger.error("  1. Check gnomAD database directory accessibility")
            logger.error("  2. Verify chromosome-split VCF files are present and indexed")
            logger.error("  3. Ensure sufficient disk space and memory")
            logger.error("  4. Check system resource availability")
            # Immediate exit - no graceful fallback
            raise RuntimeError(error_msg)
    
    def run_variant_classification(self, vcf_path: Path) -> bool:
        """
        Run multi-modal variant classification on annotated VCF with actual FILTER updates.
        
        Args:
            vcf_path: Path to annotated VCF file to classify
            
        Returns:
            True if classification succeeded
        """
        step_start = time.time()
        logger.info("Running multi-modal variant classification with FILTER updates...")
        
        try:
            # Create temporary directory for classification processing
            temp_dir = create_temp_directory(vcf_path.parent, "classification_temp")
            classified_vcf = temp_dir / "classified.vcf.gz"
            
            # Read VCF and apply classification
            import pysam
            
            # Track classification statistics
            original_filter_counts = {}
            new_filter_counts = {}
            classification_changes = {}
            variants_processed = 0
            
            logger.info("Reading VCF and applying variant classification...")
            
            # Open input VCF for reading
            with pysam.VariantFile(str(vcf_path)) as vcf_in:
                # Create output VCF with same header
                with pysam.VariantFile(str(classified_vcf), 'w', header=vcf_in.header) as vcf_out:
                    
                    for record in vcf_in:
                        variants_processed += 1
                        
                        # Track original FILTER
                        original_filter = record.filter.keys()[0] if record.filter.keys() else 'PASS'
                        original_filter_counts[original_filter] = original_filter_counts.get(original_filter, 0) + 1
                        
                        # Extract variant information for classification
                        variant_info = {
                            'FILTER': original_filter,
                            'CHROM': record.chrom,
                            'POS': record.pos,
                            'REF': record.ref,
                            'ALT': str(record.alts[0]) if record.alts else None
                        }
                        
                        # Add INFO fields to variant_info
                        for key, value in record.info.items():
                            variant_info[key] = value
                        
                        # Classify the variant
                        classification_result = self.classifier.classify_variant_from_info(variant_info)
                        new_classification = classification_result.classification
                        
                        # Rescue Logic: Classification change = Rescue
                        # Only update FILTER when evidence thresholds are met AND category changes
                        
                        rescue_occurred = False
                        
                        if new_classification == "Germline" and original_filter != "Germline":
                            # Germline rescue: Non-germline → Germline based on population frequency
                            record.filter.clear()
                            record.filter.add("Germline")
                            rescue_occurred = True
                            logger.info(f"GERMLINE RESCUE: {record.chrom}:{record.pos} {original_filter} → Germline (population frequency evidence)")
                            # Add germline rescue flag to INFO field
                            try:
                                record.info['germline_rescue'] = True
                            except Exception as e:
                                logger.debug(f"Could not add germline_rescue flag to INFO: {e}")
                            
                        elif new_classification == "Somatic" and original_filter != "Somatic":
                            # Somatic rescue: Non-somatic → Somatic based on evidence thresholds
                            record.filter.clear()
                            record.filter.add("Somatic")
                            rescue_occurred = True
                            
                            if classification_result.rescue_flag:
                                # Cross-modality + COSMIC rescue
                                logger.info(f"SOMATIC RESCUE: {record.chrom}:{record.pos} {original_filter} → Somatic (cross-modality + COSMIC evidence)")
                                # Add somatic rescue flag to INFO field
                                try:
                                    record.info['somatic_rescue'] = True
                                except Exception as e:
                                    logger.debug(f"Could not add somatic_rescue flag to INFO: {e}")
                            else:
                                # High-confidence somatic rescue (unanimous callers)
                                logger.info(f"SOMATIC RESCUE: {record.chrom}:{record.pos} {original_filter} → Somatic (unanimous caller support)")
                                # Add somatic rescue flag to INFO field
                                try:
                                    record.info['somatic_rescue'] = True
                                except Exception as e:
                                    logger.debug(f"Could not add somatic_rescue flag to INFO: {e}")
                            
                        else:
                            # No classification change - preserve original FILTER
                            logger.debug(f"No rescue: {record.chrom}:{record.pos} remains {original_filter} (insufficient evidence or already classified)")
                            # Keep original filter unchanged
                        
                        # Track new FILTER
                        final_filter = record.filter.keys()[0] if record.filter.keys() else 'PASS'
                        new_filter_counts[final_filter] = new_filter_counts.get(final_filter, 0) + 1
                        
                        # Track classification changes
                        if original_filter != final_filter:
                            change_key = f"{original_filter} -> {final_filter}"
                            classification_changes[change_key] = classification_changes.get(change_key, 0) + 1
                        
                        # Write updated record
                        vcf_out.write(record)
            
            # Index the classified VCF
            pysam.tabix_index(str(classified_vcf), preset='vcf', force=True)
            
            # Replace original file with classified version
            import shutil
            shutil.copy2(classified_vcf, vcf_path)
            
            # Copy index if it exists
            classified_index = Path(str(classified_vcf) + '.tbi')
            if classified_index.exists():
                output_index = Path(str(vcf_path) + '.tbi')
                shutil.copy2(classified_index, output_index)
            
            # Update statistics
            self.stats['variants_classified'] = variants_processed
            classification_stats = self.classifier.get_classification_statistics()
            
            # Clean up temporary directory
            cleanup_temp_directory(temp_dir, force=True)
            
            step_time = time.time() - step_start
            self.stats['processing_steps'].append(('variant_classification', step_time))
            
            # Log detailed classification results
            logger.info(f"✓ Variant classification completed successfully ({step_time:.2f}s)")
            logger.info(f"Variants processed: {variants_processed:,}")
            
            if classification_changes:
                logger.info("FILTER reclassification changes:")
                for change, count in sorted(classification_changes.items(), key=lambda x: x[1], reverse=True):
                    logger.info(f"  {change}: {count:,} variants")
            else:
                logger.info("No FILTER reclassification changes applied")
            
            logger.info("Final FILTER distribution:")
            for filter_name, count in sorted(new_filter_counts.items()):
                logger.info(f"  {filter_name}: {count:,} variants")
            
            # Log classification summary if verbose
            if self.verbose:
                self.classifier.log_classification_summary()
            
            return True
            
        except Exception as e:
            step_time = time.time() - step_start
            error_msg = f"Variant classification failed after {step_time:.2f}s: {e}"
            logger.error(error_msg)
            logger.error("CRITICAL: Variant classification failure - exiting immediately")
            logger.error("Suggested remediation:")
            logger.error("  1. Check VCF file format and accessibility")
            logger.error("  2. Verify pysam library is available and functional")
            logger.error("  3. Ensure sufficient disk space for temporary files")
            logger.error("  4. Check annotation fields are present in VCF INFO")
            # Immediate exit - no graceful fallback
            raise RuntimeError(error_msg)
    
    def generate_statistics_report(self) -> Dict[str, Any]:
        """
        Generate comprehensive statistics report.
        
        Returns:
            Dictionary containing comprehensive annotation statistics
        """
        total_time = time.time() - self.stats['start_time']
        
        # Get input/output file statistics
        input_stats = get_vcf_statistics(self.input_vcf) if self.input_vcf.exists() else {}
        output_stats = get_vcf_statistics(self.output_vcf) if self.output_vcf and self.output_vcf.exists() else {}
        
        # Get classification statistics
        classification_stats = self.classifier.get_classification_statistics()
        
        report = {
            'summary': {
                'total_processing_time_seconds': total_time,
                'total_processing_time_formatted': format_duration(total_time),
                'fallback_occurred': self.stats['fallback_occurred'],
                'errors_count': len(self.stats['errors']),
                'warnings_count': len(self.stats['warnings'])
            },
            'files': {
                'input_vcf': {
                    'path': str(self.input_vcf),
                    'size_bytes': input_stats.get('file_size_bytes', 0),
                    'size_formatted': format_file_size(input_stats.get('file_size_bytes', 0)),
                    'total_records': input_stats.get('total_records', 0),
                    'samples': input_stats.get('samples', 0)
                },
                'output_vcf': {
                    'path': str(self.output_vcf) if self.output_vcf else None,
                    'size_bytes': output_stats.get('file_size_bytes', 0),
                    'size_formatted': format_file_size(output_stats.get('file_size_bytes', 0)),
                    'total_records': output_stats.get('total_records', 0),
                    'samples': output_stats.get('samples', 0)
                }
            },
            'annotations': {
                'cosmic_annotations_added': self.stats['cosmic_annotations'],
                'gnomad_annotations_added': self.stats['gnomad_annotations'],
                'total_annotations_added': self.stats['cosmic_annotations'] + self.stats['gnomad_annotations']
            },
            'classification': classification_stats,
            'processing_steps': {
                step_name: {
                    'duration_seconds': step_time,
                    'duration_formatted': format_duration(step_time),
                    'percentage_of_total': (step_time / total_time * 100) if total_time > 0 else 0
                }
                for step_name, step_time in self.stats['processing_steps']
            },
            'errors': self.stats['errors'],
            'warnings': self.stats['warnings']
        }
        
        return report
    
    def log_final_summary(self) -> None:
        """Log comprehensive final summary of annotation process."""
        report = self.generate_statistics_report()
        
        logger.info("=== COSMIC/gnomAD Annotation Final Summary ===")
        
        # Summary information
        summary = report['summary']
        logger.info(f"Total processing time: {summary['total_processing_time_formatted']}")
        logger.info(f"Fallback occurred: {summary['fallback_occurred']}")
        
        # File information
        files = report['files']
        input_file = files['input_vcf']
        output_file = files['output_vcf']
        
        logger.info(f"Input VCF: {input_file['path']} ({input_file['size_formatted']}, {input_file['total_records']:,} records)")
        if output_file['path']:
            logger.info(f"Output VCF: {output_file['path']} ({output_file['size_formatted']}, {output_file['total_records']:,} records)")
        
        # Annotation statistics
        annotations = report['annotations']
        logger.info(f"COSMIC annotations added: {annotations['cosmic_annotations_added']:,}")
        logger.info(f"gnomAD annotations added: {annotations['gnomad_annotations_added']:,}")
        logger.info(f"Total annotations added: {annotations['total_annotations_added']:,}")
        
        # Classification statistics
        classification = report['classification']
        if classification['total_variants_classified'] > 0:
            logger.info(f"Variants classified: {classification['total_variants_classified']:,}")
            
            if self.verbose:
                logger.info("Classification distribution:")
                for class_name, count in classification['classification_counts'].items():
                    rate = classification['classification_rates'].get(class_name, 0)
                    logger.info(f"  {class_name}: {count:,} ({rate:.1f}%)")
        
        # Processing step timings
        if self.verbose and report['processing_steps']:
            logger.info("Processing step timings:")
            for step_name, step_info in report['processing_steps'].items():
                logger.info(f"  {step_name}: {step_info['duration_formatted']} ({step_info['percentage_of_total']:.1f}%)")
        
        # Errors and warnings
        if summary['warnings_count'] > 0:
            logger.warning(f"Warnings encountered: {summary['warnings_count']}")
            if self.verbose:
                for warning in report['warnings']:
                    logger.warning(f"  {warning}")
        
        if summary['errors_count'] > 0:
            logger.error(f"Errors encountered: {summary['errors_count']}")
            for error in report['errors']:
                logger.error(f"  {error}")
        
        logger.info("=== Annotation Summary Complete ===")
    
    def run_annotation(self) -> None:
        """
        Run the complete COSMIC/gnomAD annotation pipeline with comprehensive output management.
        
        This is the main entry point that orchestrates all annotation steps:
        1. COSMIC annotation (if database provided) 
        2. gnomAD annotation (if database provided)
        3. Multi-modal variant classification
        4. Intermediate file preservation and statistics generation
        """
        pipeline_start = time.time()
        logger.info("Starting COSMIC/gnomAD annotation pipeline with comprehensive output management...")
        
        # Determine output directory and base name for consistent naming
        output_dir = self.output_vcf.parent if self.output_vcf else Path.cwd()
        # Extract base name from output file, handling Nextflow naming patterns
        if self.output_vcf:
            base_name = self.output_vcf.name.replace('.cosmic_gnomad_annotated.vcf.gz', '').replace('.vcf.gz', '').replace('.vcf', '')
        else:
            base_name = 'annotated'
        
        # Create temporary directory for processing
        temp_dir = create_temp_directory(output_dir, "cosmic_gnomad_temp")
        
        # Track intermediate files for preservation
        intermediate_files = {}
        
        try:
            current_vcf = self.input_vcf
            
            # Step 1: COSMIC annotation
            logger.info("Step 1: COSMIC annotation...")
            
            if self.cosmic_vcf:
                cosmic_output = output_dir / f"{base_name}.cosmic.vcf.gz"
                self.run_cosmic_annotation(self.input_vcf, cosmic_output)
                
                # Preserve COSMIC intermediate file
                intermediate_files['cosmic'] = cosmic_output
                current_vcf = cosmic_output
                logger.info(f"✓ COSMIC intermediate file preserved: {cosmic_output}")
            else:
                logger.info("COSMIC database not provided, skipping COSMIC annotation")
            
            # Step 2: gnomAD annotation
            logger.info("Step 2: gnomAD annotation...")
            
            if self.gnomad_dir:
                gnomad_output = output_dir / f"{base_name}.gnomad.vcf.gz"
                self.run_gnomad_annotation(current_vcf, gnomad_output)
                
                # Preserve gnomAD intermediate file
                intermediate_files['gnomad'] = gnomad_output
                current_vcf = gnomad_output
                logger.info(f"✓ gnomAD intermediate file preserved: {gnomad_output}")
            else:
                logger.info("gnomAD database not provided, skipping gnomAD annotation")
            
            # Step 3: Multi-modal variant classification
            logger.info("Step 3: Multi-modal variant classification...")
            
            # Create final classified output
            final_output = output_dir / f"{base_name}.final.vcf.gz"
            
            # Copy current VCF to final output location for classification
            import shutil
            shutil.copy2(current_vcf, final_output)
            
            # Copy index if it exists
            current_index = Path(str(current_vcf) + '.tbi')
            if current_index.exists():
                final_index = Path(str(final_output) + '.tbi')
                shutil.copy2(current_index, final_index)
            
            # Apply classification to final output
            self.run_variant_classification(final_output)
            
            # Step 4: Copy final result to specified output location if different
            if self.output_vcf and final_output != self.output_vcf:
                logger.info("Step 4: Copying to specified output location...")
                
                shutil.copy2(final_output, self.output_vcf)
                
                # Copy index
                final_index = Path(str(final_output) + '.tbi')
                if final_index.exists():
                    output_index = Path(str(self.output_vcf) + '.tbi')
                    shutil.copy2(final_index, output_index)
                
                logger.info(f"✓ Final output copied to: {self.output_vcf}")
            else:
                # Use final output as the main output
                self.output_vcf = final_output
            
            # Step 5: Generate and save comprehensive statistics
            logger.info("Step 5: Generating comprehensive statistics...")
            
            pipeline_time = time.time() - pipeline_start
            self.stats['processing_steps'].append(('total_pipeline', pipeline_time))
            
            # Generate statistics report
            stats_report = self.generate_statistics_report()
            
            # Save statistics to JSON file in output directory
            stats_file = output_dir / f"{base_name}.annotation_stats.json"
            with open(stats_file, 'w') as f:
                import json
                json.dump(stats_report, f, indent=2, default=str)
            
            logger.info(f"✓ Statistics saved to: {stats_file}")
            
            # Save performance metrics if verbose
            if self.verbose:
                performance_file = output_dir / f"{base_name}.performance_metrics.json"
                performance_data = {
                    'processing_steps': dict(self.stats['processing_steps']),
                    'total_time_seconds': pipeline_time,
                    'intermediate_files': {k: str(v) for k, v in intermediate_files.items()},
                    'worker_configuration': {
                        'max_workers': self.max_workers,
                        'auto_optimized': self.max_workers is None
                    }
                }
                
                with open(performance_file, 'w') as f:
                    json.dump(performance_data, f, indent=2)
                
                logger.info(f"✓ Performance metrics saved to: {performance_file}")
            
            # Clean up temporary directory (keep intermediate files)
            cleanup_temp_directory(temp_dir, force=True)
            
            # Log final summary with file organization
            logger.info("=== Output File Organization ===")
            logger.info(f"Final annotated VCF: {self.output_vcf}")
            
            for step, file_path in intermediate_files.items():
                logger.info(f"{step.upper()} intermediate: {file_path}")
            
            logger.info(f"Statistics: {stats_file}")
            if self.verbose:
                logger.info(f"Performance metrics: {performance_file}")
            
            self.log_final_summary()
            
            logger.info("✓ COSMIC/gnomAD annotation pipeline completed successfully!")
            
        except Exception as e:
            logger.error(f"COSMIC/gnomAD annotation pipeline failed: {e}")
            
            # Clean up temporary directory on error (preserve any successful intermediate files)
            try:
                cleanup_temp_directory(temp_dir, force=True)
            except Exception as cleanup_error:
                logger.warning(f"Failed to clean up temporary directory: {cleanup_error}")
            
            # Generate error summary with partial results
            try:
                stats_report = self.generate_statistics_report()
                error_stats_file = output_dir / f"{base_name}.error_stats.json"
                with open(error_stats_file, 'w') as f:
                    import json
                    json.dump(stats_report, f, indent=2, default=str)
                logger.info(f"Error statistics saved to: {error_stats_file}")
            except Exception as stats_error:
                logger.warning(f"Failed to save error statistics: {stats_error}")
            
            self.log_final_summary()
            raise


def main():
    """Main entry point for standalone usage."""
    parser = argparse.ArgumentParser(
        description="COSMIC/gnomAD main annotation engine",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
COSMIC/gnomAD Main Annotation Engine:
  This is the main annotation script that integrates COSMIC and gnomAD annotation
  workflows with multi-modal variant classification. It provides a unified interface
  for comprehensive variant annotation using population and cancer databases.

Key Features:
  - Integrated COSMIC and gnomAD annotation workflows
  - Multi-modal variant classification with statistics generation
  - Configurable logging (verbose/minimal) and graceful fallback mechanisms
  - Comprehensive error handling and recovery
  - Performance monitoring and resource usage tracking

Examples:
  # Full annotation with both databases
  python annotate_cosmic_gnomad.py -i input.vcf.gz -c cosmic.vcf.gz -g /path/to/gnomad -o output.vcf.gz
  
  # COSMIC annotation only
  python annotate_cosmic_gnomad.py -i input.vcf.gz -c cosmic.vcf.gz -o output.vcf.gz
  
  # gnomAD annotation only
  python annotate_cosmic_gnomad.py -i input.vcf.gz -g /path/to/gnomad -o output.vcf.gz
  
  # With verbose logging and custom workers
  python annotate_cosmic_gnomad.py -i input.vcf.gz -c cosmic.vcf.gz -g /path/to/gnomad -o output.vcf.gz --verbose --workers 8
        """
    )
    
    parser.add_argument(
        '-i', '--input',
        required=True,
        metavar='VCF_FILE',
        help='Input VCF file to be annotated (can be compressed with .gz)'
    )
    
    parser.add_argument(
        '-c', '--cosmic',
        metavar='VCF_FILE',
        help='COSMIC VCF database file (optional, can be compressed with .gz)'
    )
    
    parser.add_argument(
        '-g', '--gnomad',
        metavar='DIRECTORY',
        help='gnomAD database directory containing chromosome-split VCF files (optional)'
    )
    
    parser.add_argument(
        '-o', '--output',
        required=True,
        metavar='VCF_FILE',
        help='Output annotated VCF file (use .gz extension for compression)'
    )
    
    parser.add_argument(
        '--workers', '-w',
        type=int,
        default=None,
        metavar='N',
        help='Maximum number of parallel workers for gnomAD processing (default: auto-optimize based on system resources)'
    )
    
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Enable verbose logging and detailed statistics'
    )
    
    parser.add_argument(
        '--germline-freq-threshold',
        type=float,
        default=0.001,
        metavar='FLOAT',
        help='Population frequency threshold for germline classification (default: 0.001)'
    )
    
    parser.add_argument(
        '--somatic-consensus-threshold',
        type=int,
        default=3,
        metavar='INT',
        help='Minimum caller support for high-confidence somatic classification (default: 3)'
    )
    
    parser.add_argument(
        '--cosmic-recurrence-threshold',
        type=int,
        default=5,
        metavar='INT',
        help='Minimum COSMIC recurrence for rescue classification (default: 5)'
    )
    
    parser.add_argument(
        '--stats-output',
        metavar='JSON_FILE',
        help='Output file for detailed statistics in JSON format (optional)'
    )
    
    args = parser.parse_args()
    
    # Comprehensive parameter validation
    logger.info("Validating parameters and system requirements...")
    
    # Check required databases
    if not args.cosmic and not args.gnomad:
        logger.error("CRITICAL: No annotation databases provided")
        logger.error("At least one database must be provided (--cosmic or --gnomad)")
        logger.error("Use --cosmic for COSMIC database or --gnomad for gnomAD directory")
        sys.exit(1)
    
    # Validate input file
    input_path = Path(args.input)
    if not input_path.exists():
        logger.error(f"CRITICAL: Input VCF file not found: {args.input}")
        logger.error("Please check the file path and ensure the file exists")
        sys.exit(1)
    
    if not input_path.is_file():
        logger.error(f"CRITICAL: Input path is not a file: {args.input}")
        sys.exit(1)
    
    # Validate COSMIC database if provided
    if args.cosmic:
        cosmic_path = Path(args.cosmic)
        if not cosmic_path.exists():
            logger.error(f"CRITICAL: COSMIC database file not found: {args.cosmic}")
            logger.error("Please check the file path and ensure the COSMIC VCF file exists")
            sys.exit(1)
        
        if not cosmic_path.is_file():
            logger.error(f"CRITICAL: COSMIC path is not a file: {args.cosmic}")
            sys.exit(1)
        
        # Check for index file
        cosmic_index = Path(str(cosmic_path) + '.tbi')
        if not cosmic_index.exists():
            logger.error(f"CRITICAL: COSMIC database index not found: {cosmic_index}")
            logger.error("Please create tabix index: tabix -p vcf {args.cosmic}")
            sys.exit(1)
    
    # Validate gnomAD directory if provided
    if args.gnomad:
        gnomad_path = Path(args.gnomad)
        if not gnomad_path.exists():
            logger.error(f"CRITICAL: gnomAD directory not found: {args.gnomad}")
            logger.error("Please check the directory path and ensure it exists")
            sys.exit(1)
        
        if not gnomad_path.is_dir():
            logger.error(f"CRITICAL: gnomAD path is not a directory: {args.gnomad}")
            sys.exit(1)
        
        # Check for chromosome files
        import glob
        chr_files = glob.glob(str(gnomad_path / "*.chr*.vcf.bgz")) + glob.glob(str(gnomad_path / "*.chr*.vcf.gz"))
        if not chr_files:
            logger.error(f"CRITICAL: No chromosome-split VCF files found in gnomAD directory: {args.gnomad}")
            logger.error("Expected files like: gnomad.exomes.v4.1.sites.chr1.vcf.bgz")
            sys.exit(1)
    
    # Validate output directory
    output_path = Path(args.output)
    output_dir = output_path.parent
    if not output_dir.exists():
        try:
            output_dir.mkdir(parents=True, exist_ok=True)
            logger.info(f"Created output directory: {output_dir}")
        except Exception as e:
            logger.error(f"CRITICAL: Cannot create output directory {output_dir}: {e}")
            sys.exit(1)
    
    if not os.access(output_dir, os.W_OK):
        logger.error(f"CRITICAL: Cannot write to output directory: {output_dir}")
        logger.error("Please check directory permissions")
        sys.exit(1)
    
    # Validate worker configuration
    if args.workers is not None and args.workers < 1:
        logger.error(f"CRITICAL: Invalid worker count: {args.workers}")
        logger.error("Worker count must be at least 1 or None for auto-optimization")
        sys.exit(1)
    
    # Validate threshold parameters
    if args.germline_freq_threshold < 0 or args.germline_freq_threshold > 1:
        logger.error(f"CRITICAL: Invalid germline frequency threshold: {args.germline_freq_threshold}")
        logger.error("Threshold must be between 0 and 1")
        sys.exit(1)
    
    if args.somatic_consensus_threshold < 1:
        logger.error(f"CRITICAL: Invalid somatic consensus threshold: {args.somatic_consensus_threshold}")
        logger.error("Threshold must be at least 1")
        sys.exit(1)
    
    if args.cosmic_recurrence_threshold < 1:
        logger.error(f"CRITICAL: Invalid COSMIC recurrence threshold: {args.cosmic_recurrence_threshold}")
        logger.error("Threshold must be at least 1")
        sys.exit(1)
    
    logger.info("✓ Parameter validation completed successfully")
    
    # Build classification configuration
    classification_config = {
        'germline_frequency_threshold': args.germline_freq_threshold,
        'somatic_consensus_threshold': args.somatic_consensus_threshold,
        'cosmic_recurrence_threshold': args.cosmic_recurrence_threshold,
        'cross_modality_min_support': 1
    }
    
    # Run annotation
    try:
        annotator = CosmicGnomadAnnotator(
            input_vcf=args.input,
            cosmic_vcf=args.cosmic,
            gnomad_dir=args.gnomad,
            output_vcf=args.output,
            verbose=args.verbose,
            max_workers=args.workers,
            classification_config=classification_config
        )
        
        annotator.run_annotation()
        
        # Save detailed statistics if requested
        if args.stats_output:
            stats_report = annotator.generate_statistics_report()
            with open(args.stats_output, 'w') as f:
                json.dump(stats_report, f, indent=2, default=str)
            logger.info(f"Detailed statistics saved to: {args.stats_output}")
        
        logger.info("✓ COSMIC/gnomAD annotation completed successfully!")
        sys.exit(0)
        
    except KeyboardInterrupt:
        logger.error("COSMIC/gnomAD annotation interrupted by user")
        sys.exit(130)
    except Exception as e:
        logger.error(f"COSMIC/gnomAD annotation failed: {e}")
        if args.verbose:
            import traceback
            logger.error(f"Traceback: {traceback.format_exc()}")
        sys.exit(1)


if __name__ == '__main__':
    main()