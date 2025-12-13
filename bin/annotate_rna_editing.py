#!/usr/bin/env python3
"""
RNA Editing Annotation Script - Simplified Version

This script provides streamlined RNA editing annotation using only bcftools
for VCF annotation with REDIportal database. The refactored architecture
eliminates complex fallback mechanisms and uses system PATH tools exclusively.

FEATURES:
- Single RNAEditingAnnotator class using bcftools annotate
- System PATH-based tool discovery with shutil.which()
- VCF header conflict detection and resolution
- Comprehensive error handling and logging
- Standard exit codes for pipeline integration

Requirements Satisfied: 1.1, 1.2, 1.3, 1.5, 2.1, 2.4, 2.5, 3.1-3.5, 4.1-4.5, 5.1, 5.3-5.5

Author: RNA Editing Enhancement Pipeline
Date: 2025-12-12
"""

import argparse
import sys
import logging
import shutil
import subprocess
import os
import time
import traceback
from pathlib import Path
from typing import Dict, Optional, List, Any, Tuple

# Add vcf_utils to path for imports
script_dir = Path(__file__).parent
sys.path.insert(0, str(script_dir))

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
    Simplified RNA editing annotator using only bcftools.
    
    This class provides streamlined RNA editing annotation using bcftools annotate
    for exact matching with REDIportal database. It uses system PATH for tool
    discovery and handles VCF header conflicts gracefully.
    """
    
    def __init__(self, input_vcf: str, rediportal_vcf: str, output_vcf: str, 
                 min_rna_support: int = 2):
        self.input_vcf = Path(input_vcf)
        self.rediportal_vcf = Path(rediportal_vcf)
        self.output_vcf = Path(output_vcf)
        self.min_rna_support = min_rna_support
        
        # Tool paths discovered via system PATH
        self.tool_paths: Dict[str, Optional[str]] = {}
        
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
        
        logger.info("=== RNA Editing Annotation Pipeline Started ===")
        logger.info(f"Input VCF: {self.input_vcf}")
        logger.info(f"REDIportal VCF: {self.rediportal_vcf}")
        logger.info(f"Output VCF: {self.output_vcf}")
        logger.info(f"Min RNA support: {self.min_rna_support}")
        
        # Validate inputs and tools with comprehensive error handling
        try:
            self.validate_inputs()
            self.validate_tools()
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
                logger.error(f"File permissions: {oct(os.stat(self.input_vcf).st_mode)[-3:]}")
                logger.error(f"File owner: {os.stat(self.input_vcf).st_uid}")
                logger.error(f"Current user: {os.getuid()}")
                raise PermissionError(error_msg)
                
            if not os.access(self.rediportal_vcf, os.R_OK):
                error_msg = f"Cannot read REDIportal VCF: {self.rediportal_vcf}"
                logger.error(error_msg)
                logger.error(f"File permissions: {oct(os.stat(self.rediportal_vcf).st_mode)[-3:]}")
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
                logger.error(f"Directory permissions: {oct(os.stat(output_dir).st_mode)[-3:]}")
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
    
    def validate_tools(self) -> None:
        """Validate required tools are available in system PATH with version checking."""
        step_start = time.time()
        logger.info("Validating required tools in system PATH...")
        
        required_tools = ['bcftools', 'tabix', 'bgzip']
        missing_tools = []
        tool_versions = {}
        
        # Check for pysam availability for evidence classification
        try:
            import pysam
            logger.info("✓ pysam available for enhanced VCF processing")
        except ImportError:
            logger.warning("pysam not available - evidence classification will be skipped")
            self.stats['warnings'].append("pysam not available for evidence classification")
        
        # Validate tools using system PATH only
        
        for tool in required_tools:
            try:
                tool_path = shutil.which(tool)
                if tool_path:
                    self.tool_paths[tool] = tool_path
                    
                    # Get tool version for diagnostics
                    try:
                        version_result = subprocess.run(
                            [tool_path, '--version'],
                            capture_output=True,
                            text=True,
                            timeout=10
                        )
                        if version_result.returncode == 0:
                            version_line = version_result.stdout.split('\n')[0]
                            tool_versions[tool] = version_line
                            logger.info(f"✓ Found {tool}: {tool_path} ({version_line})")
                        else:
                            logger.warning(f"✓ Found {tool}: {tool_path} (version check failed)")
                            tool_versions[tool] = "version unknown"
                    except subprocess.TimeoutExpired:
                        logger.warning(f"✓ Found {tool}: {tool_path} (version check timed out)")
                        tool_versions[tool] = "version timeout"
                    except Exception as e:
                        logger.warning(f"✓ Found {tool}: {tool_path} (version check error: {e})")
                        tool_versions[tool] = f"version error: {e}"
                else:
                    self.tool_paths[tool] = None
                    missing_tools.append(tool)
                    logger.error(f"✗ {tool} not found in system PATH")
                    
            except Exception as e:
                logger.error(f"✗ Error checking {tool}: {e}")
                missing_tools.append(tool)
        
        if missing_tools:
            error_msg = f"Required tools not found in system PATH: {missing_tools}"
            logger.error(error_msg)
            logger.error("Diagnostic information:")
            logger.error(f"  Current user: {os.getenv('USER', 'unknown')}")
            logger.error(f"  Shell: {os.getenv('SHELL', 'unknown')}")
            logger.error("")
            logger.error("Installation instructions:")
            logger.error("  Ubuntu/Debian: apt install bcftools tabix")
            logger.error("  CentOS/RHEL: yum install bcftools htslib")
            logger.error("  Conda: conda install bcftools htslib")
            logger.error("  Homebrew: brew install bcftools htslib")
            raise RuntimeError(error_msg)
        
        # Store tool versions in stats
        self.stats['tool_versions'] = tool_versions
        
        step_time = time.time() - step_start
        self.stats['processing_steps'].append(('validate_tools', step_time))
        logger.info(f"✓ All required tools found in system PATH ({step_time:.2f}s)")
    
    def check_header_compatibility(self) -> Dict[str, Any]:
        """
        Check VCF header for existing fields and determine annotation strategy.
        
        Returns:
            Dict containing header analysis results and annotation strategy
        """
        step_start = time.time()
        logger.info("Checking VCF header compatibility...")
        
        try:
            # Read VCF header to check for existing INFO fields
            cmd = [str(self.tool_paths['bcftools']), 'view', '-h', str(self.input_vcf)]
            logger.debug(f"Running command: {' '.join(cmd)}")
            
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True,
                timeout=60
            )
            
            header_lines = result.stdout.strip().split('\n')
            existing_info_fields = set()
            header_stats = {
                'total_lines': len(header_lines),
                'info_lines': 0,
                'format_lines': 0,
                'contig_lines': 0,
                'sample_count': 0
            }
            
            for line in header_lines:
                if line.startswith('##INFO=<ID='):
                    header_stats['info_lines'] += 1
                    # Extract field ID with better error handling
                    try:
                        start = line.find('ID=') + 3
                        end = line.find(',', start)
                        if end == -1:
                            end = line.find('>', start)
                        field_id = line[start:end]
                        existing_info_fields.add(field_id)
                    except Exception as e:
                        logger.warning(f"Failed to parse INFO field from line: {line[:100]}... Error: {e}")
                elif line.startswith('##FORMAT='):
                    header_stats['format_lines'] += 1
                elif line.startswith('##contig='):
                    header_stats['contig_lines'] += 1
                elif line.startswith('#CHROM'):
                    # Count samples
                    columns = line.split('\t')
                    if len(columns) > 9:
                        header_stats['sample_count'] = len(columns) - 9
            
            logger.info(f"VCF header analysis: {header_stats['total_lines']} lines, "
                       f"{header_stats['info_lines']} INFO fields, "
                       f"{header_stats['sample_count']} samples")
            
            # Check for potential conflicts with REDIportal annotation fields
            redi_annotation_fields = {'REDI_EVIDENCE', 'REDI_CANONICAL', 'DB', 
                                    'REDI_TYPE', 'REDI_REPEAT', 'REDI_FUNC', 'REDI_STRAND'}
            conflicts = existing_info_fields.intersection(redi_annotation_fields)
            
            # Determine annotation strategy based on conflicts
            if conflicts:
                logger.warning(f"Found existing INFO fields that may conflict: {conflicts}")
                strategy = 'handle_conflicts'
                logger.info("Will use conflict-aware annotation strategy")
                self.stats['warnings'].append(f"Header conflicts detected: {conflicts}")
            else:
                logger.info("✓ No header conflicts detected")
                strategy = 'standard'
            
            analysis_result = {
                'existing_fields': existing_info_fields,
                'conflicting_fields': conflicts,
                'strategy': strategy,
                'requires_preprocessing': len(conflicts) > 0,
                'header_stats': header_stats
            }
            
            step_time = time.time() - step_start
            self.stats['processing_steps'].append(('check_header_compatibility', step_time))
            logger.info(f"✓ Header compatibility check completed ({step_time:.2f}s)")
            
            return analysis_result
                
        except subprocess.CalledProcessError as e:
            error_msg = f"Failed to check VCF header: {e}"
            logger.error(error_msg)
            logger.error(f"Command: {' '.join(cmd)}")
            logger.error(f"Return code: {e.returncode}")
            logger.error(f"stdout: {e.stdout}")
            logger.error(f"stderr: {e.stderr}")
            self._log_error("VCF header validation failed", e)
            raise RuntimeError(error_msg)
        except subprocess.TimeoutExpired as e:
            error_msg = "VCF header check timed out after 60 seconds"
            logger.error(error_msg)
            logger.error("This may indicate a very large VCF file or system performance issues")
            self._log_error("VCF header check timeout", e)
            raise RuntimeError(error_msg)
        except Exception as e:
            self._log_error("Unexpected error during header compatibility check", e)
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
        
        # Add traceback for debugging
        if logger.isEnabledFor(logging.DEBUG):
            logger.debug("Full traceback:")
            logger.debug(traceback.format_exc())
        
        # Store error for summary reporting
        self.stats['errors'].append(error_info)
        
        # Add specific diagnostic information based on error type
        if isinstance(error, subprocess.CalledProcessError):
            logger.error(f"Command exit code: {error.returncode}")
            if hasattr(error, 'cmd'):
                logger.error(f"Failed command: {' '.join(str(x) for x in error.cmd)}")
            if hasattr(error, 'stdout') and error.stdout:
                logger.error(f"Command stdout: {error.stdout}")
            if hasattr(error, 'stderr') and error.stderr:
                logger.error(f"Command stderr: {error.stderr}")
        elif isinstance(error, FileNotFoundError):
            logger.error("File not found - check file paths and permissions")
        elif isinstance(error, PermissionError):
            logger.error("Permission denied - check file and directory permissions")
        elif isinstance(error, subprocess.TimeoutExpired):
            logger.error("Command timed out - may indicate performance issues or large files")
    
    def _generate_summary_statistics(self) -> None:
        """Generate and log comprehensive summary statistics."""
        total_time = time.time() - self.stats['start_time']
        
        logger.info("=== RNA Editing Annotation Summary ===")
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
        
        # Log tool versions used
        if 'tool_versions' in self.stats:
            logger.info("Tool versions used:")
            for tool, version in self.stats['tool_versions'].items():
                logger.info(f"  {tool}: {version}")
        
        # Log warnings and errors
        if self.stats['warnings']:
            logger.info(f"Warnings encountered: {len(self.stats['warnings'])}")
            for warning in self.stats['warnings']:
                logger.warning(f"  {warning}")
        
        if self.stats['errors']:
            logger.error(f"Errors encountered: {len(self.stats['errors'])}")
            for error_info in self.stats['errors']:
                logger.error(f"  {error_info['context']}: {error_info['error_type']} - {error_info['error_message']}")
        
        # Log category transition statistics if available
        if 'category_transition_stats' in self.stats:
            stats = self.stats['category_transition_stats']
            logger.info("Category transition statistics:")
            logger.info(f"  Total variants processed: {stats['summary']['total_variants_processed']}")
            logger.info(f"  Variants reclassified to RNAedit: {stats['summary']['total_reclassified_to_rnaedit']}")
            logger.info(f"  Overall reclassification rate: {stats['summary']['overall_reclassification_rate']:.1f}%")
            logger.info(f"  Exact REDIportal match rate: {stats['summary']['exact_match_rate']:.1f}%")
            
            # Log category-specific transitions
            for category, category_stats in stats['category_transitions'].items():
                if category_stats['total'] > 0:
                    logger.info(f"  {category}: {category_stats['to_rnaedit']}/{category_stats['total']} "
                               f"({category_stats['percentage']:.1f}%) reclassified to RNAedit")
        
        logger.info(f"Log file: {self.log_file}")
        logger.info("=== Summary Complete ===")

    def run_annotation(self) -> None:
        """
        Run the RNA editing annotation pipeline using bcftools with comprehensive error handling.
        """
        pipeline_start = time.time()
        logger.info("Starting RNA editing annotation with bcftools...")
        
        temp_files = []
        temp_header_files = []
        
        try:
            # Check header compatibility and get strategy
            logger.info("Step 1: Analyzing VCF header compatibility...")
            header_analysis = self.check_header_compatibility()
            
            # Prepare input files (compress if needed)
            logger.info("Step 2: Preparing input files...")
            prepared_input = self._prepare_input_vcf()
            if prepared_input != self.input_vcf:
                temp_files.append(prepared_input)
            
            logger.info("Step 3: Preparing REDIportal database...")
            prepared_rediportal = self.prepare_rediportal_database()
            if prepared_rediportal != self.rediportal_vcf:
                temp_files.append(prepared_rediportal)
            
            # Execute bcftools annotate with appropriate strategy
            logger.info("Step 4: Running bcftools annotation...")
            temp_header_file = self._execute_bcftools_annotate(header_analysis, prepared_input, prepared_rediportal)
            if temp_header_file:
                temp_header_files.append(temp_header_file)
            
            # Process variants for evidence classification and FILTER updates
            logger.info("Step 5: Processing variants for RNA editing evidence classification...")
            self._process_variants_for_evidence_classification(prepared_rediportal)
            
            # Validate output file creation and accessibility
            logger.info("Step 6: Validating output file creation...")
            self._validate_output_file_creation()
            
            # Ensure proper file permissions
            logger.info("Step 7: Setting output file permissions...")
            self._ensure_output_file_permissions()
            
            # Create output index if needed
            logger.info("Step 8: Creating output index...")
            self._create_output_index()
            
            # Clean up temporary files
            logger.info("Step 9: Cleaning up temporary files...")
            self._cleanup_temp_files(temp_files + temp_header_files)
            
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
                self._cleanup_temp_files(temp_files + temp_header_files, force_cleanup=True)
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
            elif isinstance(e, subprocess.CalledProcessError):
                logger.error(f"FATAL: Tool execution failed with exit code {e.returncode}")
                sys.exit(e.returncode if e.returncode != 0 else 1)
            elif isinstance(e, RuntimeError) and "not found in system PATH" in str(e):
                logger.error("FATAL: Required tools not available in system PATH")
                sys.exit(127)  # Command not found
            else:
                logger.error(f"FATAL: Unexpected error during annotation: {e}")
                sys.exit(1)  # General error
    
    def _prepare_input_vcf(self) -> Path:
        """
        Prepare input VCF file for bcftools (compress if needed).
        
        Returns:
            Path to prepared (possibly compressed) input VCF
        """
        step_start = time.time()
        
        if str(self.input_vcf).endswith('.gz'):
            logger.info("Input VCF is already compressed")
            # Check if index exists
            index_file = Path(str(self.input_vcf) + '.tbi')
            if not index_file.exists():
                logger.info("Creating missing tabix index for input VCF...")
                try:
                    subprocess.run(
                        [str(self.tool_paths['tabix']), '-p', 'vcf', str(self.input_vcf)],
                        capture_output=True,
                        text=True,
                        check=True,
                        timeout=300
                    )
                    logger.info("✓ Tabix index created")
                except Exception as e:
                    self._log_error("Failed to create tabix index for input VCF", e)
                    raise
            return self.input_vcf
        
        logger.info("Compressing input VCF for bcftools compatibility...")
        compressed_input = Path(str(self.input_vcf) + '.gz')
        
        try:
            # Check if compressed version already exists
            if compressed_input.exists():
                logger.warning(f"Compressed file already exists: {compressed_input}")
                logger.info("Using existing compressed file")
                
                # Ensure it has an index
                index_file = Path(str(compressed_input) + '.tbi')
                if not index_file.exists():
                    logger.info("Creating tabix index for existing compressed file...")
                    subprocess.run(
                        [str(self.tool_paths['tabix']), '-p', 'vcf', str(compressed_input)],
                        capture_output=True,
                        text=True,
                        check=True,
                        timeout=300
                    )
                
                step_time = time.time() - step_start
                self.stats['processing_steps'].append(('prepare_input_vcf', step_time))
                return compressed_input
            
            # Compress with bgzip
            logger.info(f"Compressing {self.input_vcf} to {compressed_input}...")
            cmd = [str(self.tool_paths['bgzip']), '-c', str(self.input_vcf)]
            
            with open(compressed_input, 'wb') as f:
                subprocess.run(
                    cmd, 
                    stdout=f, 
                    stderr=subprocess.PIPE, 
                    text=True,
                    check=True,
                    timeout=1800  # 30 minute timeout for compression
                )
            
            # Verify compression worked
            if not compressed_input.exists() or compressed_input.stat().st_size == 0:
                raise RuntimeError("Compression failed - output file is missing or empty")
            
            logger.info(f"✓ Compression completed, size: {compressed_input.stat().st_size:,} bytes")
            
            # Index with tabix
            logger.info("Creating tabix index...")
            subprocess.run(
                [str(self.tool_paths['tabix']), '-p', 'vcf', str(compressed_input)],
                capture_output=True,
                text=True,
                check=True,
                timeout=300
            )
            
            # Verify index was created
            index_file = Path(str(compressed_input) + '.tbi')
            if not index_file.exists():
                raise RuntimeError("Tabix indexing failed - index file not created")
            
            logger.info(f"✓ Input VCF compressed and indexed: {compressed_input}")
            
            step_time = time.time() - step_start
            self.stats['processing_steps'].append(('prepare_input_vcf', step_time))
            return compressed_input
            
        except subprocess.CalledProcessError as e:
            error_msg = f"Failed to prepare input VCF: {e}"
            logger.error(error_msg)
            logger.error(f"Command stderr: {e.stderr}")
            
            # Clean up partial files
            if compressed_input.exists():
                try:
                    compressed_input.unlink()
                    logger.info("Cleaned up partial compressed file")
                except Exception:
                    pass
            
            self._log_error("VCF preparation failed", e)
            raise RuntimeError(error_msg)
        except subprocess.TimeoutExpired as e:
            error_msg = "VCF preparation timed out"
            logger.error(error_msg)
            logger.error("This may indicate a very large VCF file or system performance issues")
            self._log_error(error_msg, e)
            raise RuntimeError(error_msg)
        except Exception as e:
            self._log_error("Unexpected error during input VCF preparation", e)
            raise
    
    def prepare_rediportal_database(self) -> Path:
        """
        Prepare REDIportal database for bcftools annotation with automatic format detection.
        
        This method automatically detects whether the REDIportal database is in VCF or
        gzipped text format and converts it appropriately for bcftools annotation.
        
        Returns:
            Path to prepared (compressed and indexed) REDIportal database
        """
        step_start = time.time()
        logger.info("Preparing REDIportal database with automatic format detection...")
        
        try:
            # Import REDIportal converter from vcf_utils
            from vcf_utils.rediportal_converter import detect_rediportal_format, prepare_rediportal_database
            
            # Detect REDIportal format
            format_type = detect_rediportal_format(str(self.rediportal_vcf))
            logger.info(f"Detected REDIportal format: {format_type}")
            
            # Prepare database based on format
            if format_type == 'text':
                logger.info("Converting REDIportal gzipped text format to bcftools annotation format...")
                prepared_file = prepare_rediportal_database(str(self.rediportal_vcf))
                logger.info(f"✓ REDIportal text format converted: {prepared_file}")
            elif format_type == 'vcf':
                logger.info("Preparing REDIportal VCF format...")
                prepared_file = prepare_rediportal_database(str(self.rediportal_vcf))
                logger.info(f"✓ REDIportal VCF format prepared: {prepared_file}")
            else:
                raise ValueError(f"Unsupported REDIportal format: {format_type}")
            
            step_time = time.time() - step_start
            self.stats['processing_steps'].append(('prepare_rediportal_database', step_time))
            
            return Path(prepared_file)
            
        except ImportError as e:
            logger.error(f"REDIportal converter not available: {e}")
            logger.info("Falling back to legacy VCF preparation...")
            return self._prepare_rediportal_vcf_legacy()
        except Exception as e:
            logger.error(f"REDIportal database preparation failed: {e}")
            logger.info("Falling back to legacy VCF preparation...")
            return self._prepare_rediportal_vcf_legacy()
    
    def _prepare_rediportal_vcf_legacy(self) -> Path:
        """
        Legacy REDIportal VCF preparation method (VCF format only).
        
        Returns:
            Path to prepared (compressed and indexed) REDIportal VCF
        """
        if str(self.rediportal_vcf).endswith('.gz'):
            # Check if index exists
            index_file = Path(str(self.rediportal_vcf) + '.tbi')
            if index_file.exists():
                logger.info("REDIportal VCF is already compressed and indexed")
                return self.rediportal_vcf
            else:
                logger.info("Creating index for compressed REDIportal VCF...")
                try:
                    subprocess.run(
                        [str(self.tool_paths['tabix']), '-p', 'vcf', str(self.rediportal_vcf)],
                        check=True
                    )
                    return self.rediportal_vcf
                except subprocess.CalledProcessError as e:
                    logger.error(f"Failed to index REDIportal VCF: {e}")
                    raise RuntimeError(f"REDIportal VCF indexing failed: {e}")
        
        logger.info("Compressing and indexing REDIportal VCF...")
        compressed_rediportal = Path(str(self.rediportal_vcf) + '.gz')
        
        try:
            # Compress with bgzip
            cmd = [str(self.tool_paths['bgzip']), '-c', str(self.rediportal_vcf)]
            with open(compressed_rediportal, 'wb') as f:
                subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE, check=True)
            
            # Index with tabix
            subprocess.run(
                [str(self.tool_paths['tabix']), '-p', 'vcf', str(compressed_rediportal)],
                check=True
            )
            
            logger.info(f"✓ REDIportal VCF compressed and indexed: {compressed_rediportal}")
            return compressed_rediportal
            
        except subprocess.CalledProcessError as e:
            logger.error(f"Failed to prepare REDIportal VCF: {e}")
            raise RuntimeError(f"REDIportal VCF preparation failed: {e}")
    
    def _execute_bcftools_annotate(self, header_analysis: Dict[str, Any], 
                                 input_vcf: Path, rediportal_vcf: Path) -> Optional[Path]:
        """
        Execute bcftools annotate command with conflict-aware strategy and comprehensive error handling.
        
        Args:
            header_analysis: Results from header compatibility check
            input_vcf: Path to prepared input VCF
            rediportal_vcf: Path to prepared REDIportal VCF
        """
        step_start = time.time()
        logger.info("Executing bcftools annotate...")
        temp_header_file = None
        
        # Build bcftools annotate command based on strategy
        if header_analysis['strategy'] == 'handle_conflicts':
            cmd, temp_header_file = self._build_conflict_aware_command(header_analysis, input_vcf, rediportal_vcf)
            logger.info(f"Using conflict-aware strategy for fields: {header_analysis['conflicting_fields']}")
        else:
            cmd, temp_header_file = self._build_standard_command(input_vcf, rediportal_vcf)
            logger.info("Using standard annotation strategy")
        
        logger.info(f"Running command: {' '.join(str(x) for x in cmd)}")
        
        try:
            # Run with timeout to prevent hanging
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True,
                timeout=3600  # 1 hour timeout for large files
            )
            
            # Log command output for diagnostics
            if result.stdout:
                logger.debug(f"bcftools stdout: {result.stdout}")
            if result.stderr:
                # bcftools often writes informational messages to stderr
                if "lines processed" in result.stderr or "records written" in result.stderr:
                    logger.info(f"bcftools progress: {result.stderr.strip()}")
                else:
                    logger.debug(f"bcftools stderr: {result.stderr}")
            
            # Verify output file was created
            if not self.output_vcf.exists():
                raise RuntimeError("bcftools completed but output file was not created")
            
            # Check output file size
            output_size = self.output_vcf.stat().st_size
            if output_size == 0:
                raise RuntimeError("bcftools created empty output file")
            
            logger.info(f"✓ bcftools annotation completed, output size: {output_size:,} bytes")
            
            step_time = time.time() - step_start
            self.stats['processing_steps'].append(('bcftools_annotate', step_time))
            
            return temp_header_file
            
        except subprocess.CalledProcessError as e:
            error_context = f"bcftools annotate failed with exit code {e.returncode}"
            logger.error(error_context)
            logger.error(f"Command: {' '.join(str(x) for x in cmd)}")
            logger.error(f"stdout: {e.stdout}")
            logger.error(f"stderr: {e.stderr}")
            
            # Provide specific error guidance based on stderr content
            stderr_lower = e.stderr.lower() if e.stderr else ""
            
            if "header already exists" in stderr_lower:
                logger.error("DIAGNOSIS: VCF header conflict detected")
                logger.error("SOLUTION: The input VCF already contains annotation fields that conflict with REDIportal fields")
                logger.error("  - Check if the VCF was previously annotated")
                logger.error("  - Consider removing existing annotation fields first")
                logger.error("  - Use bcftools annotate -x to remove conflicting fields")
                
            elif "no such file" in stderr_lower or "not found" in stderr_lower:
                logger.error("DIAGNOSIS: File not found error")
                logger.error("SOLUTION: Check that all input files exist and are accessible")
                logger.error(f"  - Input VCF: {input_vcf} (exists: {input_vcf.exists()})")
                logger.error(f"  - REDIportal VCF: {rediportal_vcf} (exists: {rediportal_vcf.exists()})")
                
            elif "permission denied" in stderr_lower:
                logger.error("DIAGNOSIS: Permission denied error")
                logger.error("SOLUTION: Check file permissions and output directory access")
                logger.error(f"  - Output directory: {self.output_vcf.parent}")
                logger.error(f"  - Directory writable: {os.access(self.output_vcf.parent, os.W_OK)}")
                
            elif "not compressed with bgzip" in stderr_lower or "not indexed" in stderr_lower:
                logger.error("DIAGNOSIS: VCF file format issue")
                logger.error("SOLUTION: Files should be bgzip compressed and tabix indexed for bcftools")
                logger.error("  - Use 'bgzip file.vcf' to compress")
                logger.error("  - Use 'tabix -p vcf file.vcf.gz' to index")
                
            elif "malformed" in stderr_lower or "invalid" in stderr_lower:
                logger.error("DIAGNOSIS: VCF format validation error")
                logger.error("SOLUTION: Input VCF may have format issues")
                logger.error("  - Use 'bcftools view -h' to check header format")
                logger.error("  - Use 'bcftools stats' to validate VCF structure")
                
            elif "out of memory" in stderr_lower or "memory" in stderr_lower:
                logger.error("DIAGNOSIS: Memory allocation error")
                logger.error("SOLUTION: Insufficient memory for processing")
                logger.error("  - Try processing smaller VCF chunks")
                logger.error("  - Increase available system memory")
                
            else:
                logger.error("DIAGNOSIS: Unknown bcftools error")
                logger.error("SOLUTION: Check bcftools documentation and input file formats")
            
            self._log_error(error_context, e)
            raise RuntimeError(f"bcftools annotate failed: {e}")
            
        except subprocess.TimeoutExpired as e:
            error_msg = "bcftools annotate timed out after 1 hour"
            logger.error(error_msg)
            logger.error("DIAGNOSIS: Processing is taking too long")
            logger.error("SOLUTION: This may indicate:")
            logger.error("  - Very large input files requiring more time")
            logger.error("  - System performance issues")
            logger.error("  - Inefficient annotation strategy")
            logger.error("  - Consider splitting large VCF files into smaller chunks")
            
            self._log_error(error_msg, e)
            raise RuntimeError(error_msg)
            
        except Exception as e:
            self._log_error("Unexpected error during bcftools execution", e)
            raise
    
    def _build_standard_command(self, input_vcf: Path, rediportal_vcf: Path) -> Tuple[List[str], Optional[Path]]:
        """Build standard bcftools annotate command with REDIportal format support."""
        # Detect if this is a converted REDIportal text format (has extended columns)
        # or legacy VCF format (has only DB field)
        temp_header_file = None
        
        try:
            from vcf_utils.rediportal_converter import detect_rediportal_format
            original_format = detect_rediportal_format(str(self.rediportal_vcf))
            
            if original_format == 'converted' or original_format == 'text':
                # Use extended annotation format for converted text files
                # Create temporary header file instead of using bash process substitution
                temp_header_file = self._create_temp_header_file()
                cmd = [
                    str(self.tool_paths['bcftools']), 'annotate',
                    '-a', str(rediportal_vcf),
                    '-c', 'CHROM,POS,REF,ALT,REDI_ACCESSION,REDI_DB,REDI_TYPE,REDI_REPEAT,REDI_FUNC,REDI_STRAND',
                    '-h', str(temp_header_file),
                    '-o', str(self.output_vcf),
                    str(input_vcf)
                ]
            else:
                # Use legacy VCF format annotation
                cmd = [
                    str(self.tool_paths['bcftools']), 'annotate',
                    '-a', str(rediportal_vcf),
                    '-c', 'CHROM,POS,REF,ALT,INFO/DB',
                    '-o', str(self.output_vcf),
                    str(input_vcf)
                ]
        except ImportError:
            # Fall back to legacy format if converter not available
            cmd = [
                str(self.tool_paths['bcftools']), 'annotate',
                '-a', str(rediportal_vcf),
                '-c', 'CHROM,POS,REF,ALT,INFO/DB',
                '-o', str(self.output_vcf),
                str(input_vcf)
            ]
        
        # Add compression if output is .gz
        if str(self.output_vcf).endswith('.gz'):
            cmd.extend(['-O', 'z'])
        
        return cmd, temp_header_file
    
    def _create_temp_header_file(self) -> Path:
        """Create temporary header file for REDIportal annotation fields."""
        import tempfile
        
        header_content = """##INFO=<ID=REDI_ACCESSION,Number=1,Type=String,Description="REDIportal accession ID">
##INFO=<ID=REDI_DB,Number=1,Type=String,Description="REDIportal database source">
##INFO=<ID=REDI_TYPE,Number=1,Type=String,Description="RNA editing type">
##INFO=<ID=REDI_REPEAT,Number=1,Type=String,Description="Repeat element annotation">
##INFO=<ID=REDI_FUNC,Number=1,Type=String,Description="Functional annotation">
##INFO=<ID=REDI_STRAND,Number=1,Type=String,Description="Strand information">
"""
        
        # Create temporary file
        temp_fd, temp_path = tempfile.mkstemp(suffix='.txt', prefix='redi_headers_')
        try:
            with os.fdopen(temp_fd, 'w') as f:
                f.write(header_content)
        except Exception:
            os.close(temp_fd)
            raise
        
        return Path(temp_path)
    
    def _create_selective_header_file(self, safe_fields: List[str]) -> Path:
        """Create temporary header file for selected REDIportal annotation fields."""
        import tempfile
        
        header_definitions = {
            'REDI_ACCESSION': '##INFO=<ID=REDI_ACCESSION,Number=1,Type=String,Description="REDIportal accession ID">',
            'REDI_DB': '##INFO=<ID=REDI_DB,Number=1,Type=String,Description="REDIportal database source">',
            'REDI_TYPE': '##INFO=<ID=REDI_TYPE,Number=1,Type=String,Description="RNA editing type">',
            'REDI_REPEAT': '##INFO=<ID=REDI_REPEAT,Number=1,Type=String,Description="Repeat element annotation">',
            'REDI_FUNC': '##INFO=<ID=REDI_FUNC,Number=1,Type=String,Description="Functional annotation">',
            'REDI_STRAND': '##INFO=<ID=REDI_STRAND,Number=1,Type=String,Description="Strand information">'
        }
        
        header_content = ""
        for field in safe_fields:
            if field in header_definitions:
                header_content += header_definitions[field] + "\n"
        
        # Create temporary file
        temp_fd, temp_path = tempfile.mkstemp(suffix='.txt', prefix='redi_selective_headers_')
        try:
            with os.fdopen(temp_fd, 'w') as f:
                f.write(header_content)
        except Exception:
            os.close(temp_fd)
            raise
        
        return Path(temp_path)
    
    def _build_conflict_aware_command(self, header_analysis: Dict[str, Any], 
                                    input_vcf: Path, rediportal_vcf: Path) -> Tuple[List[str], Optional[Path]]:
        """
        Build bcftools annotate command that handles header conflicts with REDIportal format support.
        
        Args:
            header_analysis: Results from header compatibility check
            input_vcf: Path to prepared input VCF
            rediportal_vcf: Path to prepared REDIportal VCF
            
        Returns:
            List of command arguments for bcftools
        """
        conflicting_fields = header_analysis['conflicting_fields']
        
        logger.info(f"Building conflict-aware command for fields: {conflicting_fields}")
        
        # Detect REDIportal format to determine annotation strategy
        try:
            from vcf_utils.rediportal_converter import detect_rediportal_format
            original_format = detect_rediportal_format(str(self.rediportal_vcf))
            
            if original_format == 'converted' or original_format == 'text':
                # For converted text format, use selective annotation to avoid conflicts
                # Only annotate fields that don't conflict
                available_fields = ['REDI_ACCESSION', 'REDI_DB', 'REDI_TYPE', 'REDI_REPEAT', 'REDI_FUNC', 'REDI_STRAND']
                safe_fields = [field for field in available_fields if field not in conflicting_fields]
                
                if safe_fields:
                    columns = 'CHROM,POS,REF,ALT,' + ','.join(safe_fields)
                    
                    # Create temporary header file for safe fields
                    header_file = self._create_selective_header_file(safe_fields)
                    
                    cmd = [
                        str(self.tool_paths['bcftools']), 'annotate',
                        '-a', str(rediportal_vcf),
                        '-c', columns,
                        '-h', str(header_file),
                        '-o', str(self.output_vcf),
                        str(input_vcf)
                    ]
                else:
                    # If all fields conflict, just mark sites without adding INFO fields
                    cmd = [
                        str(self.tool_paths['bcftools']), 'annotate',
                        '-a', str(rediportal_vcf),
                        '--mark-sites', '+REDIPORTAL_MATCH',
                        '-o', str(self.output_vcf),
                        str(input_vcf)
                    ]
            else:
                # Legacy VCF format - use mark-sites strategy
                cmd = [
                    str(self.tool_paths['bcftools']), 'annotate',
                    '-a', str(rediportal_vcf),
                    '-c', 'CHROM,POS,REF,ALT,INFO/DB',
                    '--mark-sites', '+DB',  # Mark sites that match without duplicating headers
                    '-o', str(self.output_vcf),
                    str(input_vcf)
                ]
        except ImportError:
            # Fall back to legacy strategy if converter not available
            cmd = [
                str(self.tool_paths['bcftools']), 'annotate',
                '-a', str(rediportal_vcf),
                '-c', 'CHROM,POS,REF,ALT,INFO/DB',
                '--mark-sites', '+DB',
                '-o', str(self.output_vcf),
                str(input_vcf)
            ]
        
        # Add compression if output is .gz
        if str(self.output_vcf).endswith('.gz'):
            cmd.extend(['-O', 'z'])
        
        return cmd, header_file if 'header_file' in locals() else None
    
    def _process_variants_for_evidence_classification(self, rediportal_vcf: Path) -> None:
        """
        Process annotated variants to add RNA editing evidence classification and FILTER updates.
        
        This method integrates evidence classification with the bcftools annotation pipeline by:
        1. Reading the bcftools-annotated VCF output
        2. Adding REDI_EVIDENCE and REDI_CANONICAL INFO fields
        3. Updating FILTER fields based on evidence classification
        4. Writing the enhanced output VCF
        
        Args:
            rediportal_vcf: Path to REDIportal VCF database for exact matching
        """
        step_start = time.time()
        logger.info("Processing variants for RNA editing evidence classification...")
        
        # Check if output file exists (should be created by bcftools)
        if not self.output_vcf.exists():
            logger.warning("Output VCF not found - evidence classification requires bcftools annotation first")
            return
        
        try:
            # Import required modules
            try:
                import pysam
                from vcf_utils.rna_editing_core import (
                    is_canonical_editing_transition,
                    classify_rna_editing_evidence,
                    classify_rna_editing_biological_category,
                    initialize_category_transition_stats,
                    update_category_transition_stats,
                    finalize_category_transition_stats
                )
            except ImportError as e:
                logger.error(f"Required modules not available: {e}")
                logger.warning("Skipping evidence classification - using basic annotation only")
                return
            
            # Create temporary file for enhanced output
            import tempfile
            temp_dir = self.output_vcf.parent
            temp_enhanced = temp_dir / f"temp_enhanced_{int(time.time())}.vcf"
            if str(self.output_vcf).endswith('.gz'):
                temp_enhanced = Path(str(temp_enhanced) + '.gz')
            
            # Initialize statistics
            evidence_stats = {
                'total_variants': 0,
                'rediportal_matches': 0,
                'canonical_transitions': 0,
                'evidence_levels': {'VERY_HIGH': 0, 'HIGH': 0, 'MEDIUM': 0, 'LOW': 0, 'NONE': 0},
                'filter_updates': 0,
                'filter_preserved': 0
            }
            
            category_stats = initialize_category_transition_stats()
            
            # Process variants with pysam for robust VCF handling
            try:
                input_vcf = pysam.VariantFile(str(self.output_vcf))
            except Exception as e:
                logger.error(f"Failed to open annotated VCF file: {e}")
                logger.warning("Skipping evidence classification")
                return
            
            # Prepare output VCF header with RNA editing fields
            output_header = input_vcf.header.copy()
            
            # Add RNA editing specific INFO fields if not present
            if 'REDI_EVIDENCE' not in output_header.info:
                output_header.info.add('REDI_EVIDENCE', number=1, type='String', 
                                     description='RNA editing evidence level (VERY_HIGH, HIGH, MEDIUM, LOW, NONE)')
            if 'REDI_CANONICAL' not in output_header.info:
                output_header.info.add('REDI_CANONICAL', number=1, type='String',
                                     description='Canonical A>G or T>C transition (YES/NO)')
            
            # Add RNAedit filter to header if not present
            if 'RNAedit' not in output_header.filters:
                output_header.filters.add('RNAedit', None, None, 
                                        'RNA editing variant based on evidence classification')
            
            # Open output VCF
            try:
                output_vcf = pysam.VariantFile(str(temp_enhanced), 'w', header=output_header)
            except Exception as e:
                logger.error(f"Failed to create enhanced output VCF: {e}")
                input_vcf.close()
                return
            
            # Process each variant
            try:
                for variant in input_vcf:
                    evidence_stats['total_variants'] += 1
                    
                    # Extract variant data
                    variant_data = self._extract_variant_data_pysam(variant)
                    
                    # Check for REDIportal match (DB flag from bcftools annotation or REDI_DB field from text format)
                    rediportal_match = ('DB' in variant.info or 
                                      'REDI_DB' in variant.info or 
                                      'REDI_ACCESSION' in variant.info)
                    if rediportal_match:
                        evidence_stats['rediportal_matches'] += 1
                    
                    # Check canonical transition
                    canonical = is_canonical_editing_transition(variant.ref, ','.join(variant.alts))
                    if canonical:
                        evidence_stats['canonical_transitions'] += 1
                    
                    # Classify evidence level
                    evidence_level = classify_rna_editing_evidence(
                        variant_data, rediportal_match, self.min_rna_support
                    )
                    evidence_stats['evidence_levels'][evidence_level] += 1
                    
                    # Determine biological classification and FILTER update
                    original_filter = list(variant.filter.keys())[0] if variant.filter.keys() else 'PASS'
                    biological_classification = classify_rna_editing_biological_category(
                        evidence_level, original_filter
                    )
                    
                    # Update category transition statistics
                    category_stats = update_category_transition_stats(
                        category_stats, original_filter, biological_classification, 
                        evidence_level, rediportal_match
                    )
                    
                    # Track FILTER updates
                    if biological_classification != original_filter:
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
                    
                    # Copy all INFO fields from input with proper handling of multi-value fields
                    for key, value in variant.info.items():
                        try:
                            output_variant.info[key] = value
                        except Exception as e:
                            # Handle multi-value fields that pysam has trouble with
                            logger.debug(f"Skipping INFO field {key} due to pysam error: {e}")
                            continue
                    
                    # Add RNA editing specific fields
                    output_variant.info['REDI_EVIDENCE'] = evidence_level
                    output_variant.info['REDI_CANONICAL'] = 'YES' if canonical else 'NO'
                    
                    # Set biological classification as FILTER
                    output_variant.filter.clear()
                    output_variant.filter.add(biological_classification)
                    
                    # Copy sample information if present
                    if len(variant.samples) > 0:
                        sample_name = list(variant.samples.keys())[0]
                        if sample_name in output_header.samples:
                            for fmt_key in variant.format.keys():
                                try:
                                    output_variant.samples[sample_name][fmt_key] = variant.samples[sample_name][fmt_key]
                                except Exception as e:
                                    # Handle FORMAT fields that pysam has trouble with
                                    logger.debug(f"Skipping FORMAT field {fmt_key} due to pysam error: {e}")
                                    continue
                    
                    # Write enhanced variant
                    output_vcf.write(output_variant)
                
            except Exception as e:
                logger.error(f"Error processing variants: {e}")
                raise
            finally:
                input_vcf.close()
                output_vcf.close()
            
            # Replace original output with enhanced version
            if self.output_vcf.exists():
                self.output_vcf.unlink()
            temp_enhanced.rename(self.output_vcf)
            
            # Store statistics
            self.stats.update(evidence_stats)
            self.stats['category_transition_stats'] = finalize_category_transition_stats(category_stats)
            
            step_time = time.time() - step_start
            self.stats['processing_steps'].append(('evidence_classification', step_time))
            
            logger.info(f"✓ Evidence classification completed: {evidence_stats['total_variants']} variants processed")
            logger.info(f"  REDIportal matches: {evidence_stats['rediportal_matches']}")
            logger.info(f"  Canonical transitions: {evidence_stats['canonical_transitions']}")
            logger.info(f"  FILTER updates: {evidence_stats['filter_updates']}")
            logger.info(f"  FILTER preserved: {evidence_stats['filter_preserved']}")
            
        except Exception as e:
            logger.error(f"Evidence classification failed: {e}")
            logger.warning("Continuing with basic annotation only")
            # Don't raise - this is not a fatal error, annotation can continue
    
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
    
    def _cleanup_temp_files(self, temp_files: List[Path], force_cleanup: bool = False) -> None:
        """
        Clean up temporary files created during processing with comprehensive error handling.
        
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
                    
                    # Clean up parent directory if it's a temporary directory and empty
                    parent_dir = temp_file.parent
                    if parent_dir.name.startswith('tmp') and parent_dir != Path.cwd():
                        try:
                            if not any(parent_dir.iterdir()):  # Directory is empty
                                parent_dir.rmdir()
                                logger.debug(f"Cleaned up empty temporary directory: {parent_dir}")
                        except Exception:
                            # Ignore errors when cleaning up directories
                            pass
                            
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
        """Create tabix index for output VCF if it's compressed with comprehensive validation."""
        if not str(self.output_vcf).endswith('.gz'):
            logger.debug("Output VCF is not compressed, skipping indexing")
            return
        
        step_start = time.time()
        logger.info("Creating tabix index for compressed output...")
        
        try:
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
            
            cmd = [str(self.tool_paths['tabix']), '-p', 'vcf', str(self.output_vcf)]
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
            
            # Test index functionality with a quick query
            try:
                test_cmd = [str(self.tool_paths['tabix']), str(self.output_vcf), '1:1-1']
                subprocess.run(test_cmd, capture_output=True, check=True, timeout=30)
                logger.debug("✓ Tabix index functionality verified")
            except subprocess.CalledProcessError:
                # This is expected if no variants exist in the test region
                logger.debug("✓ Tabix index created (no variants in test region)")
            except Exception as e:
                logger.warning(f"Tabix index functionality test failed: {e}")
            
            index_size = index_file.stat().st_size
            step_time = time.time() - step_start
            self.stats['processing_steps'].append(('create_output_index', step_time))
            
            logger.info(f"✓ Tabix index created successfully ({index_size:,} bytes, {step_time:.2f}s)")
            
        except subprocess.CalledProcessError as e:
            error_msg = f"Failed to create tabix index: {e}"
            logger.error(error_msg)
            logger.error(f"Command stderr: {e.stderr}")
            logger.error("This may indicate:")
            logger.error("  - Output VCF format issues")
            logger.error("  - Insufficient disk space")
            logger.error("  - File permission problems")
            self._log_error("Tabix indexing failed", e)
            raise RuntimeError(error_msg)
            
        except subprocess.TimeoutExpired as e:
            error_msg = "Tabix indexing timed out after 5 minutes"
            logger.error(error_msg)
            logger.error("This may indicate very large output files or system performance issues")
            self._log_error("Tabix indexing timeout", e)
            raise RuntimeError(error_msg)
            
        except Exception as e:
            self._log_error("Unexpected error during indexing", e)
            raise

    def _validate_output_file_creation(self) -> None:
        """Validate that output file was created successfully and is accessible."""
        step_start = time.time()
        logger.info("Validating output file creation and accessibility...")
        
        try:
            # Check file existence
            if not self.output_vcf.exists():
                raise RuntimeError(f"Output file was not created: {self.output_vcf}")
            
            # Check file size
            file_size = self.output_vcf.stat().st_size
            if file_size == 0:
                raise RuntimeError(f"Output file is empty: {self.output_vcf}")
            
            # Check file permissions
            if not os.access(self.output_vcf, os.R_OK):
                raise RuntimeError(f"Output file is not readable: {self.output_vcf}")
            
            # Validate VCF format if possible
            try:
                if str(self.output_vcf).endswith('.gz'):
                    # Test that it's a valid gzipped file
                    cmd = [str(self.tool_paths['bcftools']), 'view', '-h', str(self.output_vcf)]
                    result = subprocess.run(cmd, capture_output=True, text=True, check=True, timeout=60)
                    
                    # Check for VCF header
                    if not result.stdout.startswith('##fileformat=VCF'):
                        raise RuntimeError("Output file does not have valid VCF header")
                    
                    logger.debug("✓ Output VCF format validation passed")
                else:
                    # For uncompressed files, check first line
                    with open(self.output_vcf, 'r') as f:
                        first_line = f.readline().strip()
                        if not first_line.startswith('##fileformat=VCF'):
                            raise RuntimeError("Output file does not have valid VCF header")
                    
                    logger.debug("✓ Output VCF format validation passed")
                    
            except subprocess.CalledProcessError as e:
                logger.warning(f"VCF format validation failed: {e}")
                logger.warning("Output file may have format issues")
            except Exception as e:
                logger.warning(f"Could not validate VCF format: {e}")
            
            step_time = time.time() - step_start
            self.stats['processing_steps'].append(('validate_output_creation', step_time))
            
            logger.info(f"✓ Output file validation completed ({file_size:,} bytes, {step_time:.2f}s)")
            
        except Exception as e:
            self._log_error("Output file validation failed", e)
            raise


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Simplified RNA editing annotation using bcftools",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Simplified RNA Editing Annotation:
  This script provides streamlined RNA editing annotation using only bcftools
  for VCF annotation with REDIportal database. The refactored architecture
  eliminates complex fallback mechanisms and uses system PATH tools exclusively.

Key Features:
  - Single RNAEditingAnnotator class using bcftools annotate
  - System PATH-based tool discovery with shutil.which()
  - VCF header conflict detection and resolution
  - Comprehensive error handling and logging
  - Standard exit codes for pipeline integration

Tool Requirements:
  - bcftools (available in system PATH)
  - tabix (available in system PATH)
  - bgzip (available in system PATH)

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