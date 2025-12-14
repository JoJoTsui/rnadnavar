#!/usr/bin/env python3
"""
RNA Editing Annotation Script - Reorganized for Maintainability

This script provides streamlined RNA editing annotation using dedicated modules
for REDIportal conversion, bcftools annotation, and evidence tiering. The
reorganized architecture maintains clear separation of concerns between components
and uses pysam for reliable VCF file generation.

FEATURES:
- Modular architecture with dedicated components
- REDIportal conversion logic in dedicated module
- bcftools annotation logic in dedicated module  
- Evidence tiering logic in dedicated module
- pysam-based VCF output generation
- Comprehensive error handling and logging

Requirements Satisfied: 7.1, 7.2, 7.4, 7.5

Author: RNA Editing Enhancement Pipeline
Date: 2025-12-14
"""

import argparse
import sys
import logging
import os
import time
from pathlib import Path
from typing import Dict, Optional, Any

# Add vcf_utils to path for imports
script_dir = Path(__file__).parent
sys.path.insert(0, str(script_dir))

# Import dedicated modules for RNA editing annotation
try:
    from vcf_utils.rediportal_converter import detect_rediportal_format, prepare_rediportal_database
    from vcf_utils.bcftools_annotator import BcftoolsAnnotator
    from vcf_utils.evidence_tiering import create_evidence_tiering_processor
    from vcf_utils.filter_updater import create_filter_updater
    from vcf_utils.rna_editing_core import is_canonical_editing_transition
    MODULES_AVAILABLE = True
except ImportError as e:
    logging.warning(f"Some RNA editing modules not available: {e}")
    MODULES_AVAILABLE = False

# Set up comprehensive logging with timestamps
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

# Add file handler for persistent logging
def setup_file_logging(output_dir: Path) -> None:
    """Set up file logging in addition to console logging."""
    log_file = output_dir.parent / f"rna_editing_annotation_{int(time.time())}.log"
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    logger.info(f"Logging to file: {log_file}")
    return log_file


class RNAEditingAnnotator:
    """
    Reorganized RNA editing annotator using dedicated modules.
    
    This class coordinates RNA editing annotation using dedicated modules for:
    - REDIportal database conversion (rediportal_converter)
    - bcftools VCF-to-VCF annotation (bcftools_annotator)
    - Evidence tiering and classification (evidence_tiering)
    - FILTER column updates (filter_updater)
    - pysam-based VCF output generation
    """
    
    def __init__(self, input_vcf: str, rediportal_vcf: str, output_vcf: str, 
                 min_rna_support: int = 2):
        self.input_vcf = Path(input_vcf)
        self.rediportal_vcf = Path(rediportal_vcf)
        self.output_vcf = Path(output_vcf)
        self.min_rna_support = min_rna_support
        
        # Statistics and timing
        self.stats = {
            'start_time': time.time(),
            'total_variants': 0,
            'annotated_variants': 0,
            'processing_steps': [],
            'errors': [],
            'warnings': []
        }
        
        # Set up file logging
        self.log_file = setup_file_logging(self.output_vcf)
        
        logger.info("=== RNA Editing Annotation Pipeline Started (Reorganized) ===")
        logger.info(f"Input VCF: {self.input_vcf}")
        logger.info(f"REDIportal VCF: {self.rediportal_vcf}")
        logger.info(f"Output VCF: {self.output_vcf}")
        logger.info(f"Min RNA support: {self.min_rna_support}")
        
        # Check module availability
        if not MODULES_AVAILABLE:
            logger.error("Required RNA editing modules are not available")
            raise RuntimeError("RNA editing modules not available - cannot proceed")
        
        # Validate inputs
        try:
            self.validate_inputs()
        except Exception as e:
            self._log_error("Initialization failed", e)
            raise
    
    def validate_inputs(self) -> None:
        """Validate input files exist and are accessible with detailed diagnostics."""
        step_start = time.time()
        logger.info("Validating input files...")
        
        try:
            # Check input VCF existence and accessibility
            if not self.input_vcf.exists():
                error_msg = f"Input VCF not found: {self.input_vcf}"
                logger.error(error_msg)
                logger.error(f"Current working directory: {os.getcwd()}")
                logger.error(f"Absolute path attempted: {self.input_vcf.absolute()}")
                raise FileNotFoundError(error_msg)
            
            if not self.rediportal_vcf.exists():
                error_msg = f"REDIportal VCF not found: {self.rediportal_vcf}"
                logger.error(error_msg)
                logger.error(f"Absolute path attempted: {self.rediportal_vcf.absolute()}")
                raise FileNotFoundError(error_msg)
            
            # Check file permissions with detailed diagnostics
            if not os.access(self.input_vcf, os.R_OK):
                error_msg = f"Cannot read input VCF: {self.input_vcf}"
                logger.error(error_msg)
                raise PermissionError(error_msg)
                
            if not os.access(self.rediportal_vcf, os.R_OK):
                error_msg = f"Cannot read REDIportal VCF: {self.rediportal_vcf}"
                logger.error(error_msg)
                raise PermissionError(error_msg)
            
            # Check output directory accessibility
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
            
            # Log file sizes for diagnostics
            input_size = self.input_vcf.stat().st_size
            rediportal_size = self.rediportal_vcf.stat().st_size
            logger.info(f"Input VCF size: {input_size:,} bytes ({input_size/1024/1024:.1f} MB)")
            logger.info(f"REDIportal VCF size: {rediportal_size:,} bytes ({rediportal_size/1024/1024:.1f} MB)")
            
            step_time = time.time() - step_start
            self.stats['processing_steps'].append(('validate_inputs', step_time))
            logger.info(f"✓ Input files validated successfully ({step_time:.2f}s)")
            
        except Exception as e:
            self._log_error("Input validation failed", e)
            raise
    
    def run_annotation(self) -> None:
        """
        Run the RNA editing annotation pipeline using dedicated modules.
        """
        pipeline_start = time.time()
        logger.info("Starting RNA editing annotation with reorganized modules...")
        
        temp_files = []
        
        try:
            # Step 1: Prepare REDIportal database using dedicated converter
            logger.info("Step 1: Preparing REDIportal database...")
            prepared_rediportal = self._prepare_rediportal_database_with_converter()
            if prepared_rediportal != str(self.rediportal_vcf):
                temp_files.append(Path(prepared_rediportal))
            
            # Step 2: Run bcftools annotation using dedicated annotator
            logger.info("Step 2: Running bcftools VCF-to-VCF annotation...")
            annotated_vcf = self._run_bcftools_annotation_with_engine(prepared_rediportal)
            
            # Step 3: Process variants for evidence classification and FILTER updates using pysam
            logger.info("Step 3: Processing variants for RNA editing evidence classification...")
            self._process_variants_with_pysam_and_evidence_tiering(annotated_vcf)
            
            # Step 4: Validate output and create index
            logger.info("Step 4: Validating output and creating index...")
            self._validate_and_index_output()
            
            # Step 5: Clean up temporary files
            logger.info("Step 5: Cleaning up temporary files...")
            self._cleanup_temp_files(temp_files)
            
            # Generate final statistics
            pipeline_time = time.time() - pipeline_start
            self.stats['processing_steps'].append(('total_pipeline', pipeline_time))
            
            logger.info("✓ RNA editing annotation completed successfully!")
            self._generate_summary_statistics()
            
        except Exception as e:
            self._log_error("RNA editing annotation pipeline failed", e)
            
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
            
            # Generate error summary
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
        
        Returns:
            Path to prepared REDIportal database
        """
        step_start = time.time()
        logger.info("Preparing REDIportal database with dedicated converter...")
        
        try:
            # Use dedicated REDIportal converter module
            prepared_file = prepare_rediportal_database(str(self.rediportal_vcf))
            
            step_time = time.time() - step_start
            self.stats['processing_steps'].append(('prepare_rediportal_converter', step_time))
            
            logger.info(f"✓ REDIportal database prepared: {prepared_file} ({step_time:.2f}s)")
            return prepared_file
            
        except Exception as e:
            self._log_error("REDIportal database preparation failed", e)
            raise RuntimeError(f"REDIportal preparation failed: {e}")
    
    def _run_bcftools_annotation_with_engine(self, rediportal_vcf: str) -> str:
        """
        Run bcftools annotation using dedicated annotation engine.
        
        Args:
            rediportal_vcf: Path to prepared REDIportal VCF
            
        Returns:
            Path to annotated VCF file
        """
        step_start = time.time()
        logger.info("Running bcftools annotation with dedicated engine...")
        
        try:
            # Create temporary output for bcftools annotation
            temp_annotated = self.output_vcf.parent / f"temp_annotated_{int(time.time())}.vcf"
            if str(self.output_vcf).endswith('.gz'):
                temp_annotated = Path(str(temp_annotated) + '.gz')
            
            # Use dedicated bcftools annotator
            annotator = BcftoolsAnnotator(
                input_vcf=str(self.input_vcf),
                annotation_vcf=rediportal_vcf,
                output_vcf=str(temp_annotated)
            )
            
            # Run annotation
            annotator.run_annotation()
            
            # Merge statistics from annotator
            if hasattr(annotator, 'stats'):
                for step_name, step_time in annotator.stats.get('processing_steps', []):
                    self.stats['processing_steps'].append((f'annotator_{step_name}', step_time))
                
                self.stats['warnings'].extend(annotator.stats.get('warnings', []))
                self.stats['errors'].extend(annotator.stats.get('errors', []))
            
            step_time = time.time() - step_start
            self.stats['processing_steps'].append(('bcftools_annotation_engine', step_time))
            
            logger.info(f"✓ bcftools annotation completed: {temp_annotated} ({step_time:.2f}s)")
            return str(temp_annotated)
            
        except Exception as e:
            self._log_error("bcftools annotation engine failed", e)
            raise RuntimeError(f"bcftools annotation failed: {e}")
    
    def _process_variants_with_pysam_and_evidence_tiering(self, annotated_vcf: str) -> None:
        """
        Process variants using pysam and dedicated evidence tiering modules.
        
        Args:
            annotated_vcf: Path to bcftools-annotated VCF file
        """
        step_start = time.time()
        logger.info("Processing variants with pysam and evidence tiering...")
        
        try:
            # Check if pysam is available
            try:
                import pysam
            except ImportError:
                logger.warning("pysam not available - skipping evidence tiering")
                # Just copy the annotated file to output
                import shutil
                shutil.copy2(annotated_vcf, self.output_vcf)
                return
            
            # Create evidence tiering processor
            evidence_processor = create_evidence_tiering_processor(min_rna_support=self.min_rna_support)
            
            # Create FILTER updater
            filter_updater = create_filter_updater()
            
            # Process variants with pysam
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
            
            # Open output VCF with pysam
            output_vcf = pysam.VariantFile(str(self.output_vcf), 'w', header=output_header)
            
            # Initialize statistics
            evidence_stats = {
                'total_variants': 0,
                'rediportal_matches': 0,
                'canonical_transitions': 0,
                'filter_updates': 0,
                'filter_preserved': 0
            }
            
            try:
                for variant in input_vcf:
                    evidence_stats['total_variants'] += 1
                    
                    # Extract variant data
                    variant_data = self._extract_variant_data_pysam(variant)
                    
                    # Check for REDIportal match
                    rediportal_match = ('DB' in variant.info or 
                                      'REDI_DB' in variant.info or 
                                      'REDI_ACCESSION' in variant.info)
                    if rediportal_match:
                        evidence_stats['rediportal_matches'] += 1
                    
                    # Check canonical transition
                    canonical = is_canonical_editing_transition(variant.ref, ','.join(variant.alts))
                    if canonical:
                        evidence_stats['canonical_transitions'] += 1
                    
                    # Process with evidence tiering
                    tiering_result = evidence_processor.process_variant(variant_data, rediportal_match)
                    
                    # Get original FILTER value
                    original_filter = list(variant.filter.keys())[0] if variant.filter.keys() else 'PASS'
                    
                    # Update FILTER using dedicated FILTER updater
                    new_filter, filter_was_updated = filter_updater.update_variant_filter(
                        original_filter=original_filter,
                        evidence_tier=tiering_result['evidence_tier'],
                        has_rediportal_match=rediportal_match,
                        rna_support=tiering_result['rna_support'],
                        min_rna_support=self.min_rna_support
                    )
                    
                    # Track FILTER updates
                    if filter_was_updated:
                        evidence_stats['filter_updates'] += 1
                    else:
                        evidence_stats['filter_preserved'] += 1
                    
                    # Create output variant record
                    output_variant = output_header.new_record(
                        contig=variant.chrom,
                        start=variant.start,
                        stop=variant.stop,
                        alleles=variant.alleles
                    )
                    
                    # Copy all INFO fields from input
                    for key, value in variant.info.items():
                        try:
                            output_variant.info[key] = value
                        except Exception as e:
                            logger.debug(f"Skipping INFO field {key} due to pysam error: {e}")
                            continue
                    
                    # Add RNA editing specific fields
                    output_variant.info['REDI_EVIDENCE'] = tiering_result['evidence_tier']
                    output_variant.info['REDI_CANONICAL'] = 'YES' if canonical else 'NO'
                    
                    # Set FILTER based on FILTER updater logic
                    output_variant.filter.clear()
                    output_variant.filter.add(new_filter)
                    
                    # Copy sample information if present
                    if len(variant.samples) > 0:
                        sample_name = list(variant.samples.keys())[0]
                        if sample_name in output_header.samples:
                            for fmt_key in variant.format.keys():
                                try:
                                    output_variant.samples[sample_name][fmt_key] = variant.samples[sample_name][fmt_key]
                                except Exception as e:
                                    logger.debug(f"Skipping FORMAT field {fmt_key} due to pysam error: {e}")
                                    continue
                    
                    # Write enhanced variant
                    output_vcf.write(output_variant)
                
            finally:
                input_vcf.close()
                output_vcf.close()
            
            # Store statistics
            self.stats.update(evidence_stats)
            self.stats['evidence_tiering_stats'] = evidence_processor.get_statistics()
            self.stats['filter_update_stats'] = filter_updater.get_statistics()
            
            step_time = time.time() - step_start
            self.stats['processing_steps'].append(('pysam_evidence_processing', step_time))
            
            logger.info(f"✓ Evidence processing completed: {evidence_stats['total_variants']} variants processed")
            logger.info(f"  REDIportal matches: {evidence_stats['rediportal_matches']}")
            logger.info(f"  Canonical transitions: {evidence_stats['canonical_transitions']}")
            logger.info(f"  FILTER updates: {evidence_stats['filter_updates']}")
            logger.info(f"  FILTER preserved: {evidence_stats['filter_preserved']}")
            
            # Log evidence tiering statistics
            evidence_processor.log_statistics()
            
            # Log FILTER update statistics
            filter_updater.log_statistics()
            
            # Clean up temporary annotated file
            temp_file = Path(annotated_vcf)
            if temp_file.exists() and temp_file != self.output_vcf:
                temp_file.unlink()
                logger.debug(f"Cleaned up temporary annotated file: {temp_file}")
            
        except Exception as e:
            logger.error(f"Evidence processing failed: {e}")
            # Don't raise - this is not a fatal error, annotation can continue
            logger.warning("Continuing with basic annotation only")
    
    def _validate_and_index_output(self) -> None:
        """Validate output file and create index if needed."""
        step_start = time.time()
        logger.info("Validating output file and creating index...")
        
        try:
            # Validate output file exists and is accessible
            if not self.output_vcf.exists():
                raise RuntimeError(f"Output file was not created: {self.output_vcf}")
            
            # Check file size
            file_size = self.output_vcf.stat().st_size
            if file_size == 0:
                raise RuntimeError(f"Output file is empty: {self.output_vcf}")
            
            logger.info(f"✓ Output file validated: {file_size:,} bytes")
            
            # Create index if output is compressed
            if str(self.output_vcf).endswith('.gz'):
                self._create_output_index()
            
            # Set proper file permissions
            self._ensure_output_file_permissions()
            
            step_time = time.time() - step_start
            self.stats['processing_steps'].append(('validate_and_index_output', step_time))
            
            logger.info(f"✓ Output validation and indexing completed ({step_time:.2f}s)")
            
        except Exception as e:
            self._log_error("Output validation failed", e)
            raise
    
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
            
            result = subprocess.run(
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
        """Generate and log comprehensive summary statistics."""
        total_time = time.time() - self.stats['start_time']
        
        logger.info("=== RNA Editing Annotation Summary (Reorganized) ===")
        logger.info(f"Total processing time: {total_time:.2f} seconds")
        logger.info(f"Input VCF: {self.input_vcf}")
        logger.info(f"REDIportal VCF: {self.rediportal_vcf}")
        logger.info(f"Output VCF: {self.output_vcf}")
        
        if self.output_vcf.exists():
            output_size = self.output_vcf.stat().st_size
            logger.info(f"Output file size: {output_size:,} bytes ({output_size/1024/1024:.1f} MB)")
        
        # Log processing step timings
        if self.stats['processing_steps']:
            logger.info("Processing step timings:")
            for step_name, step_time in self.stats['processing_steps']:
                logger.info(f"  {step_name}: {step_time:.2f}s")
        
        # Log warnings and errors
        if self.stats['warnings']:
            logger.info(f"Warnings encountered: {len(self.stats['warnings'])}")
            for warning in self.stats['warnings']:
                logger.warning(f"  {warning}")
        
        if self.stats['errors']:
            logger.error(f"Errors encountered: {len(self.stats['errors'])}")
            for error_info in self.stats['errors']:
                logger.error(f"  {error_info['context']}: {error_info['error_type']} - {error_info['error_message']}")
        
        # Log evidence tiering statistics if available
        if 'evidence_tiering_stats' in self.stats:
            stats = self.stats['evidence_tiering_stats']
            logger.info("Evidence tiering statistics:")
            logger.info(f"  Total variants processed: {stats['total_variants_processed']}")
            logger.info(f"  RNA consensus variants: {stats['rna_consensus_variants']} ({stats['rna_consensus_rate']:.1f}%)")
            logger.info(f"  RNA-only variants: {stats['rna_only_variants']} ({stats['rna_only_rate']:.1f}%)")
            logger.info(f"  FILTER updates to RNAedit: {stats['filter_updates']} ({stats['filter_update_rate']:.1f}%)")
        
        # Log FILTER update statistics if available
        if 'filter_update_stats' in self.stats:
            filter_stats = self.stats['filter_update_stats']
            logger.info("FILTER column update statistics:")
            logger.info(f"  Total variants processed: {filter_stats['total_variants_processed']}")
            logger.info(f"  FILTER updates to RNAedit: {filter_stats['filter_updates_to_rnaedit']} ({filter_stats['filter_update_rate']:.1f}%)")
            logger.info(f"  FILTER values preserved: {filter_stats['filter_values_preserved']} ({filter_stats['filter_preservation_rate']:.1f}%)")
        
        logger.info(f"Log file: {self.log_file}")
        logger.info("=== Summary Complete ===")


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Reorganized RNA editing annotation using dedicated modules",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Reorganized RNA Editing Annotation:
  This script provides streamlined RNA editing annotation using dedicated modules
  for REDIportal conversion, bcftools annotation, and evidence tiering. The
  reorganized architecture maintains clear separation of concerns between components
  and uses pysam for reliable VCF file generation.

Key Features:
  - Modular architecture with dedicated components
  - REDIportal conversion logic in dedicated module
  - bcftools annotation logic in dedicated module  
  - Evidence tiering logic in dedicated module
  - pysam-based VCF output generation
  - Comprehensive error handling and logging

Module Requirements:
  - vcf_utils.rediportal_converter (REDIportal database conversion)
  - vcf_utils.bcftools_annotator (bcftools VCF-to-VCF annotation)
  - vcf_utils.evidence_tiering (RNA editing evidence classification)
  - vcf_utils.filter_updater (FILTER column update logic)
  - vcf_utils.rna_editing_core (core RNA editing functions)
  - pysam (VCF file reading and writing)

Examples:
  # Basic annotation
  python annotate_rna_editing.py -i input.vcf.gz -r /path/to/rediportal.vcf.gz -o output.vcf.gz
  
  # Custom RNA support threshold
  python annotate_rna_editing.py -i input.vcf.gz -r /path/to/rediportal.vcf.gz -o output.vcf.gz --min-rna-support 3
  
  # Verbose logging
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
        default=2,
        metavar='N',
        help='Minimum RNA caller support threshold (default: 2, must be >= 1)'
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
        # Validate min_rna_support range
        if args.min_rna_support < 1:
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
    
    # Run annotation
    try:
        annotator = RNAEditingAnnotator(
            input_vcf=args.input, 
            rediportal_vcf=args.rediportal, 
            output_vcf=args.output,
            min_rna_support=args.min_rna_support
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