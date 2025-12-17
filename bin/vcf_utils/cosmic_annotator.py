#!/usr/bin/env python3
"""
COSMIC VCF-to-VCF Annotation Module

This module provides COSMIC database annotation functionality using bcftools
for VCF-to-VCF annotation. It implements exact coordinate and allele matching
between input VCF and COSMIC VCF database, extracting mutation recurrence data
and cancer-specific annotations.

FEATURES:
- Direct VCF-to-VCF annotation using bcftools annotate
- COSMIC mutation recurrence data extraction (GENOME_SCREEN_SAMPLE_COUNT, ID)
- Exact coordinate and allele matching (CHROM, POS, REF, ALT)
- Compressed VCF input/output handling with proper indexing
- Comprehensive error handling and logging

Requirements Satisfied: 2.3, 2.4

Author: COSMIC/gnomAD Enhancement Pipeline
Date: 2025-12-17
"""

import logging
import subprocess
import os
import time
from pathlib import Path
from typing import Dict, Any, Optional, Tuple

logger = logging.getLogger(__name__)


class CosmicAnnotator:
    """
    COSMIC VCF-to-VCF annotator using bcftools for cancer mutation annotation.
    
    This class provides streamlined COSMIC annotation using bcftools annotate for exact
    coordinate and allele matching between input VCF and COSMIC VCF database.
    It extracts mutation recurrence data and cancer-specific annotations.
    """
    
    def __init__(self, input_vcf: str, cosmic_vcf: str, output_vcf: str):
        """
        Initialize COSMIC annotator.
        
        Args:
            input_vcf: Path to input VCF file to be annotated
            cosmic_vcf: Path to COSMIC VCF database
            output_vcf: Path to output annotated VCF file
        """
        self.input_vcf = Path(input_vcf)
        self.cosmic_vcf = Path(cosmic_vcf)
        self.output_vcf = Path(output_vcf)
        
        # Statistics and timing
        self.stats = {
            'start_time': time.time(),
            'processing_steps': [],
            'errors': [],
            'warnings': [],
            'annotations_added': 0,
            'variants_processed': 0
        }
        
        logger.info("=== COSMIC VCF-to-VCF Annotator Initialized ===")
        logger.info(f"Input VCF: {self.input_vcf}")
        logger.info(f"COSMIC VCF: {self.cosmic_vcf}")
        logger.info(f"Output VCF: {self.output_vcf}")
        
        # Validate inputs
        self.validate_inputs()
        self.validate_tools()
    
    def validate_inputs(self) -> None:
        """Validate input files exist and are accessible."""
        step_start = time.time()
        logger.info("Validating COSMIC annotation input files...")
        
        # Check input VCF existence and accessibility
        if not self.input_vcf.exists():
            error_msg = f"Input VCF not found: {self.input_vcf}"
            logger.error(error_msg)
            raise FileNotFoundError(error_msg)
        
        if not self.cosmic_vcf.exists():
            error_msg = f"COSMIC VCF not found: {self.cosmic_vcf}"
            logger.error(error_msg)
            raise FileNotFoundError(error_msg)
        
        # Check file permissions
        if not os.access(self.input_vcf, os.R_OK):
            error_msg = f"Cannot read input VCF: {self.input_vcf}"
            logger.error(error_msg)
            raise PermissionError(error_msg)
            
        if not os.access(self.cosmic_vcf, os.R_OK):
            error_msg = f"Cannot read COSMIC VCF: {self.cosmic_vcf}"
            logger.error(error_msg)
            raise PermissionError(error_msg)
        
        # Check output directory accessibility
        output_dir = self.output_vcf.parent
        if not output_dir.exists():
            logger.info(f"Creating output directory: {output_dir}")
            output_dir.mkdir(parents=True, exist_ok=True)
        
        if not os.access(output_dir, os.W_OK):
            error_msg = f"Cannot write to output directory: {output_dir}"
            logger.error(error_msg)
            raise PermissionError(error_msg)
        
        # Log file sizes for diagnostics
        input_size = self.input_vcf.stat().st_size
        cosmic_size = self.cosmic_vcf.stat().st_size
        logger.info(f"Input VCF size: {input_size:,} bytes ({input_size/1024/1024:.1f} MB)")
        logger.info(f"COSMIC VCF size: {cosmic_size:,} bytes ({cosmic_size/1024/1024:.1f} MB)")
        
        step_time = time.time() - step_start
        self.stats['processing_steps'].append(('validate_inputs', step_time))
        logger.info(f"✓ COSMIC annotation input files validated successfully ({step_time:.2f}s)")
    
    def validate_tools(self) -> None:
        """Validate required tools are available in system PATH."""
        step_start = time.time()
        logger.info("Validating required tools for COSMIC annotation...")
        
        required_tools = ['bcftools', 'tabix', 'bgzip']
        missing_tools = []
        
        for tool in required_tools:
            tool_path = subprocess.run(['which', tool], capture_output=True, text=True)
            if tool_path.returncode != 0:
                missing_tools.append(tool)
                logger.error(f"✗ {tool} not found in system PATH")
            else:
                logger.debug(f"✓ Found {tool}: {tool_path.stdout.strip()}")
        
        if missing_tools:
            error_msg = f"Required tools not found in system PATH: {missing_tools}"
            logger.error(error_msg)
            logger.error("Installation instructions:")
            logger.error("  Ubuntu/Debian: apt install bcftools tabix")
            logger.error("  CentOS/RHEL: yum install bcftools htslib")
            logger.error("  Conda: conda install bcftools htslib")
            logger.error("  Homebrew: brew install bcftools htslib")
            raise RuntimeError(error_msg)
        
        step_time = time.time() - step_start
        self.stats['processing_steps'].append(('validate_tools', step_time))
        logger.info(f"✓ All required tools found for COSMIC annotation ({step_time:.2f}s)")
    
    def get_cosmic_annotation_fields(self) -> set:
        """
        Get COSMIC annotation fields from COSMIC VCF header.
        
        Returns:
            Set of INFO field IDs from COSMIC VCF
        """
        try:
            cmd = ['bcftools', 'view', '-h', str(self.cosmic_vcf)]
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True,
                timeout=60
            )
            
            cosmic_fields = set()
            for line in result.stdout.split('\n'):
                if line.startswith('##INFO=<ID='):
                    try:
                        start = line.find('ID=') + 3
                        end = line.find(',', start)
                        if end == -1:
                            end = line.find('>', start)
                        field_id = line[start:end]
                        cosmic_fields.add(field_id)
                    except Exception as e:
                        logger.warning(f"Failed to parse COSMIC field from line: {line[:100]}... Error: {e}")
            
            logger.info(f"Found COSMIC annotation fields: {cosmic_fields}")
            return cosmic_fields
            
        except Exception as e:
            logger.warning(f"Could not determine COSMIC annotation fields: {e}")
            # Return common COSMIC fields as fallback
            return {'COSMIC_CNT', 'COSMIC_ID', 'GENOME_SCREEN_SAMPLE_COUNT', 'ID'}
    
    def build_cosmic_annotation_command(self) -> Tuple[list, Optional[Path]]:
        """
        Build bcftools annotate command for COSMIC VCF-to-VCF annotation.
        
        Returns:
            Tuple of (command arguments, temporary header file path if created)
        """
        logger.info("Building bcftools annotate command for COSMIC VCF-to-VCF annotation...")
        
        # Get available COSMIC annotation fields
        cosmic_fields = self.get_cosmic_annotation_fields()
        
        # Use VCF-to-VCF annotation without requiring separate header files
        # bcftools will automatically match INFO fields from COSMIC VCF
        cmd = [
            'bcftools', 'annotate',
            '-a', str(self.cosmic_vcf),
            '-c', 'CHROM,POS,REF,ALT,INFO',  # Copy all INFO fields from COSMIC VCF
            '-o', str(self.output_vcf),
            str(self.input_vcf)
        ]
        
        # Add compression if output is .gz
        if str(self.output_vcf).endswith('.gz'):
            cmd.extend(['-O', 'z'])
        
        logger.info(f"COSMIC annotation command: {' '.join(cmd)}")
        return cmd, None
    
    def execute_cosmic_annotation(self) -> None:
        """
        Execute bcftools annotate command for COSMIC annotation with comprehensive error handling.
        """
        step_start = time.time()
        logger.info("Executing bcftools COSMIC VCF-to-VCF annotation...")
        
        # Build annotation command
        cmd, temp_header_file = self.build_cosmic_annotation_command()
        
        try:
            # Run with timeout to prevent hanging
            logger.info(f"Executing: {' '.join(str(x) for x in cmd)}")
            logger.info("This may take several minutes for large files...")
            
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True,
                timeout=3600  # 1 hour timeout for large files
            )
            
            execution_time = time.time() - step_start
            
            # Parse and log bcftools output statistics
            self._parse_cosmic_output_statistics(result.stdout, result.stderr)
            
            # Clean up temporary header file if created
            if temp_header_file and temp_header_file.exists():
                temp_header_file.unlink()
                logger.debug(f"Cleaned up temporary header file: {temp_header_file}")
            
            self.stats['processing_steps'].append(('execute_cosmic_annotation', execution_time))
            
            logger.info(f"✓ COSMIC annotation completed successfully ({execution_time:.2f}s)")
            
        except subprocess.CalledProcessError as e:
            execution_time = time.time() - step_start
            
            error_context = f"bcftools COSMIC annotate failed with exit code {e.returncode}"
            logger.error(error_context)
            logger.error(f"Command: {' '.join(str(x) for x in cmd)}")
            logger.error(f"Execution time: {execution_time:.2f}s")
            
            # Log detailed output
            if e.stdout:
                logger.error(f"stdout: {e.stdout}")
            if e.stderr:
                logger.error(f"stderr: {e.stderr}")
            
            # Provide specific error guidance
            self._diagnose_cosmic_error(e.stderr)
            
            # Clean up temporary files on error
            if temp_header_file and temp_header_file.exists():
                temp_header_file.unlink()
                logger.debug("Cleaned up temporary header file after error")
            
            raise RuntimeError(f"COSMIC annotation failed: {e}")
            
        except subprocess.TimeoutExpired:
            execution_time = time.time() - step_start
            error_msg = f"COSMIC annotation timed out after {execution_time:.0f} seconds (limit: 3600s)"
            
            logger.error(error_msg)
            logger.error("This may indicate:")
            logger.error("  - Very large input files requiring more processing time")
            logger.error("  - System performance issues or resource constraints")
            logger.error("  - Possible infinite loop or deadlock in bcftools")
            
            # Clean up temporary files on timeout
            if temp_header_file and temp_header_file.exists():
                temp_header_file.unlink()
                logger.debug("Cleaned up temporary header file after timeout")
            
            raise RuntimeError(error_msg)
            
        except Exception as e:
            execution_time = time.time() - step_start
            logger.error(f"Unexpected error during COSMIC annotation after {execution_time:.2f}s: {e}")
            
            # Clean up temporary files on error
            if temp_header_file and temp_header_file.exists():
                temp_header_file.unlink()
                logger.debug("Cleaned up temporary header file after unexpected error")
            
            raise
    
    def _parse_cosmic_output_statistics(self, stdout: str, stderr: str):
        """Parse and log statistics from COSMIC annotation output."""
        try:
            # Look for common bcftools statistics patterns
            combined_output = (stdout or "") + "\n" + (stderr or "")
            
            # Parse lines processed
            import re
            lines_match = re.search(r'(\d+)\s+lines?\s+processed', combined_output, re.IGNORECASE)
            if lines_match:
                lines_processed = int(lines_match.group(1))
                self.stats['variants_processed'] = lines_processed
                logger.info(f"COSMIC annotation processed {lines_processed:,} lines")
            
            # Parse annotations added
            annotations_match = re.search(r'(\d+)\s+annotations?\s+added', combined_output, re.IGNORECASE)
            if annotations_match:
                annotations_added = int(annotations_match.group(1))
                self.stats['annotations_added'] = annotations_added
                logger.info(f"COSMIC annotation added {annotations_added:,} annotations")
            
        except Exception as e:
            logger.debug(f"Could not parse COSMIC annotation output statistics: {e}")
    
    def _diagnose_cosmic_error(self, stderr: str) -> None:
        """
        Provide specific error guidance based on COSMIC annotation stderr content.
        
        Args:
            stderr: Standard error output from bcftools command
        """
        if not stderr:
            return
        
        stderr_lower = stderr.lower()
        
        if "header already exists" in stderr_lower:
            logger.error("DIAGNOSIS: VCF header conflict with COSMIC fields")
            logger.error("SOLUTION: The input VCF already contains COSMIC annotation fields")
            logger.error("  - Check if the VCF was previously annotated with COSMIC")
            logger.error("  - Consider removing existing COSMIC fields first")
            
        elif "no such file" in stderr_lower or "not found" in stderr_lower:
            logger.error("DIAGNOSIS: COSMIC database file not found")
            logger.error("SOLUTION: Check that COSMIC VCF database exists and is accessible")
            logger.error(f"  - COSMIC VCF: {self.cosmic_vcf} (exists: {self.cosmic_vcf.exists()})")
            
        elif "not compressed with bgzip" in stderr_lower or "not indexed" in stderr_lower:
            logger.error("DIAGNOSIS: COSMIC VCF file format issue")
            logger.error("SOLUTION: COSMIC VCF should be bgzip compressed and tabix indexed")
            logger.error("  - Use 'bgzip cosmic.vcf' to compress")
            logger.error("  - Use 'tabix -p vcf cosmic.vcf.gz' to index")
            
        elif "malformed" in stderr_lower or "invalid" in stderr_lower:
            logger.error("DIAGNOSIS: COSMIC VCF format validation error")
            logger.error("SOLUTION: COSMIC VCF may have format issues")
            logger.error("  - Use 'bcftools view -h' to check COSMIC VCF header format")
            logger.error("  - Ensure COSMIC VCF follows standard VCF format")
            
        else:
            logger.error("DIAGNOSIS: Unknown COSMIC annotation error")
            logger.error("SOLUTION: Check COSMIC VCF database format and accessibility")
    
    def validate_cosmic_output(self) -> None:
        """Validate that COSMIC annotated output file was created successfully."""
        step_start = time.time()
        logger.info("Validating COSMIC annotation output file...")
        
        try:
            # Check file existence
            if not self.output_vcf.exists():
                raise RuntimeError(f"COSMIC annotated output file was not created: {self.output_vcf}")
            
            # Check file size
            file_size = self.output_vcf.stat().st_size
            if file_size == 0:
                raise RuntimeError(f"COSMIC annotated output file is empty: {self.output_vcf}")
            
            # Check file permissions
            if not os.access(self.output_vcf, os.R_OK):
                raise RuntimeError(f"COSMIC annotated output file is not readable: {self.output_vcf}")
            
            # Validate VCF format
            try:
                cmd = ['bcftools', 'view', '-h', str(self.output_vcf)]
                result = subprocess.run(cmd, capture_output=True, text=True, check=True, timeout=60)
                
                # Check for VCF header
                if not result.stdout.startswith('##fileformat=VCF'):
                    raise RuntimeError("COSMIC annotated output file does not have valid VCF header")
                
                logger.debug("✓ COSMIC annotated output VCF format validation passed")
                    
            except subprocess.CalledProcessError as e:
                logger.warning(f"COSMIC annotated VCF format validation failed: {e}")
                logger.warning("Output file may have format issues")
            except Exception as e:
                logger.warning(f"Could not validate COSMIC annotated VCF format: {e}")
            
            # Create index if output is compressed
            if str(self.output_vcf).endswith('.gz'):
                self._create_cosmic_output_index()
            
            step_time = time.time() - step_start
            self.stats['processing_steps'].append(('validate_cosmic_output', step_time))
            
            logger.info(f"✓ COSMIC annotation output validation completed ({file_size:,} bytes, {step_time:.2f}s)")
            
        except Exception as e:
            logger.error(f"COSMIC annotation output validation failed: {e}")
            raise
    
    def _create_cosmic_output_index(self) -> None:
        """Create tabix index for compressed COSMIC annotated output VCF."""
        logger.info("Creating tabix index for COSMIC annotated output...")
        
        try:
            index_file = Path(str(self.output_vcf) + '.tbi')
            
            # Check if index already exists and is valid
            if index_file.exists():
                vcf_mtime = self.output_vcf.stat().st_mtime
                index_mtime = index_file.stat().st_mtime
                
                if index_mtime >= vcf_mtime:
                    logger.info("Valid tabix index already exists for COSMIC annotated output")
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
            
            # Verify index was created
            if not index_file.exists():
                raise RuntimeError("Tabix indexing completed but index file was not created")
            
            if index_file.stat().st_size == 0:
                raise RuntimeError("Tabix index file is empty")
            
            index_size = index_file.stat().st_size
            logger.info(f"✓ Tabix index created successfully for COSMIC annotated output ({index_size:,} bytes)")
            
        except subprocess.CalledProcessError as e:
            error_msg = f"Failed to create tabix index for COSMIC annotated output: {e}"
            logger.error(error_msg)
            logger.error(f"Command stderr: {e.stderr}")
            raise RuntimeError(error_msg)
            
        except subprocess.TimeoutExpired:
            error_msg = "Tabix indexing timed out after 5 minutes"
            logger.error(error_msg)
            raise RuntimeError(error_msg)
            
        except Exception as e:
            logger.error(f"Unexpected error during COSMIC output indexing: {e}")
            raise
    
    def run_cosmic_annotation(self) -> None:
        """
        Run the complete COSMIC VCF-to-VCF annotation pipeline.
        """
        pipeline_start = time.time()
        logger.info("Starting COSMIC VCF-to-VCF annotation pipeline...")
        
        try:
            # Execute COSMIC annotation
            logger.info("Step 1: Executing COSMIC annotation...")
            self.execute_cosmic_annotation()
            
            # Validate output
            logger.info("Step 2: Validating COSMIC annotation output...")
            self.validate_cosmic_output()
            
            # Generate final statistics
            pipeline_time = time.time() - pipeline_start
            self.stats['processing_steps'].append(('total_cosmic_pipeline', pipeline_time))
            
            logger.info("✓ COSMIC VCF-to-VCF annotation completed successfully!")
            self._generate_cosmic_summary_statistics()
            
        except Exception as e:
            logger.error(f"COSMIC annotation pipeline failed: {e}")
            
            # Clean up any partial output files
            if self.output_vcf.exists():
                try:
                    self.output_vcf.unlink()
                    logger.info("Cleaned up partial COSMIC annotated output file")
                except Exception as cleanup_error:
                    logger.warning(f"Failed to clean up partial COSMIC output: {cleanup_error}")
            
            # Generate error summary
            self._generate_cosmic_summary_statistics()
            raise
    
    def _generate_cosmic_summary_statistics(self) -> None:
        """Generate and log comprehensive COSMIC annotation summary statistics."""
        total_time = time.time() - self.stats['start_time']
        
        logger.info("=== COSMIC VCF-to-VCF Annotation Summary ===")
        logger.info(f"Total processing time: {total_time:.2f} seconds")
        logger.info(f"Input VCF: {self.input_vcf}")
        logger.info(f"COSMIC VCF: {self.cosmic_vcf}")
        logger.info(f"Output VCF: {self.output_vcf}")
        
        # Log annotation statistics
        if self.stats['variants_processed'] > 0:
            logger.info(f"Variants processed: {self.stats['variants_processed']:,}")
        if self.stats['annotations_added'] > 0:
            logger.info(f"COSMIC annotations added: {self.stats['annotations_added']:,}")
            annotation_rate = (self.stats['annotations_added'] / self.stats['variants_processed']) * 100
            logger.info(f"Annotation rate: {annotation_rate:.1f}%")
        
        # Log file size information
        try:
            if self.input_vcf.exists():
                input_size = self.input_vcf.stat().st_size
                logger.info(f"Input VCF size: {input_size:,} bytes ({input_size/1024/1024:.1f} MB)")
            
            if self.cosmic_vcf.exists():
                cosmic_size = self.cosmic_vcf.stat().st_size
                logger.info(f"COSMIC VCF size: {cosmic_size:,} bytes ({cosmic_size/1024/1024:.1f} MB)")
            
            if self.output_vcf.exists():
                output_size = self.output_vcf.stat().st_size
                logger.info(f"Output VCF size: {output_size:,} bytes ({output_size/1024/1024:.1f} MB)")
        except Exception as e:
            logger.debug(f"Could not get file size information: {e}")
        
        # Log processing step timings
        if self.stats['processing_steps']:
            logger.info("Processing step timings:")
            for step_name, step_time in self.stats['processing_steps']:
                percentage = (step_time / total_time * 100) if total_time > 0 else 0
                logger.info(f"  {step_name}: {step_time:.2f}s ({percentage:.1f}%)")
        
        # Log warnings and errors
        if self.stats['warnings']:
            logger.warning(f"Warnings encountered: {len(self.stats['warnings'])}")
            for warning in self.stats['warnings']:
                logger.warning(f"  {warning}")
        
        if self.stats['errors']:
            logger.error(f"Errors encountered: {len(self.stats['errors'])}")
            for error_info in self.stats['errors']:
                logger.error(f"  {error_info}")
        
        logger.info("=== COSMIC Annotation Summary Complete ===")


def main():
    """Main entry point for standalone usage."""
    import argparse
    import sys
    
    parser = argparse.ArgumentParser(
        description="COSMIC VCF-to-VCF annotation engine",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
COSMIC VCF-to-VCF Annotation Engine:
  This module provides COSMIC database annotation functionality using bcftools
  for VCF-to-VCF annotation. It implements exact coordinate and allele matching
  between input VCF and COSMIC VCF database.

Key Features:
  - Direct VCF-to-VCF annotation using bcftools annotate
  - COSMIC mutation recurrence data extraction
  - Exact coordinate and allele matching (CHROM, POS, REF, ALT)
  - Compressed VCF input/output handling with proper indexing
  - Comprehensive error handling and logging

Examples:
  # Basic COSMIC VCF-to-VCF annotation
  python cosmic_annotator.py -i input.vcf.gz -c cosmic.vcf.gz -o output.vcf.gz
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
        required=True,
        metavar='VCF_FILE',
        help='COSMIC VCF database file (can be compressed with .gz)'
    )
    
    parser.add_argument(
        '-o', '--output',
        required=True,
        metavar='VCF_FILE',
        help='Output annotated VCF file (use .gz extension for compression)'
    )
    
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Enable verbose logging'
    )
    
    args = parser.parse_args()
    
    # Set up logging
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    
    # Run COSMIC annotation
    try:
        annotator = CosmicAnnotator(
            input_vcf=args.input,
            cosmic_vcf=args.cosmic,
            output_vcf=args.output
        )
        
        annotator.run_cosmic_annotation()
        
        logger.info("✓ COSMIC VCF-to-VCF annotation completed successfully!")
        sys.exit(0)
        
    except KeyboardInterrupt:
        logger.error("COSMIC annotation interrupted by user")
        sys.exit(130)
    except Exception as e:
        logger.error(f"COSMIC annotation failed: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()