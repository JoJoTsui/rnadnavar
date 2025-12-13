#!/usr/bin/env python3
"""
Standalone REDIportal Text Format Converter

This script converts REDIportal TABLE1_hg38_v3.txt.gz format to bcftools-compatible
annotation format. It can be used independently or as part of the annotation pipeline.

SAFETY FEATURES:
- Comprehensive input validation
- No destructive operations on original files
- Extensive error handling and logging
- Atomic file operations with temporary files
- Detailed progress reporting

Usage:
    python scripts/convert_rediportal_text.py -i TABLE1_hg38_v3.txt.gz -o rediportal_annotations

Author: RNA Editing Enhancement Pipeline
Date: 2025-12-13
"""

import argparse
import gzip
import logging
import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Dict, Optional

# Add vcf_utils to path for imports
script_dir = Path(__file__).parent.parent / 'bin'
sys.path.insert(0, str(script_dir))

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

# REDIportal TABLE1 format columns (0-based indexing)
REDIPORTAL_COLUMNS = {
    'Accession': 0,
    'Region': 1, 
    'Position': 2,
    'Ref': 3,
    'Ed': 4,
    'Strand': 5,
    'db': 6,
    'type': 7,
    'repeat': 8,
    'Func': 9
}

def validate_input_file(file_path: Path) -> None:
    """
    Validate input REDIportal file.
    
    Args:
        file_path: Path to input file
        
    Raises:
        FileNotFoundError: If file doesn't exist
        PermissionError: If file is not readable
        ValueError: If file format is invalid
    """
    if not file_path.exists():
        raise FileNotFoundError(f"Input file not found: {file_path}")
    
    if not os.access(file_path, os.R_OK):
        raise PermissionError(f"Cannot read input file: {file_path}")
    
    # Check file size
    file_size = file_path.stat().st_size
    if file_size == 0:
        raise ValueError(f"Input file is empty: {file_path}")
    
    logger.info(f"Input file validation passed: {file_path}")
    logger.info(f"File size: {file_size:,} bytes ({file_size/1024/1024:.1f} MB)")
    
    # Try to read first few lines to validate format
    try:
        if str(file_path).endswith('.gz'):
            with gzip.open(file_path, 'rt') as f:
                first_line = f.readline().strip()
                second_line = f.readline().strip()
        else:
            with open(file_path, 'r') as f:
                first_line = f.readline().strip()
                second_line = f.readline().strip()
        
        # Check if it looks like REDIportal format
        if first_line.startswith('Accession') or 'Region' in first_line:
            logger.info("Detected REDIportal text format with header")
        else:
            # Check if first line has expected number of fields
            fields = first_line.split('\t')
            if len(fields) >= 10:
                logger.info("Detected REDIportal text format without header")
            else:
                raise ValueError(f"Invalid REDIportal format. First line has {len(fields)} fields, expected ≥10")
        
        # Validate second line has proper format
        if second_line:
            fields = second_line.split('\t')
            if len(fields) < 10:
                logger.warning(f"Second line has only {len(fields)} fields, may indicate format issues")
        
    except Exception as e:
        raise ValueError(f"Cannot validate REDIportal format: {e}")

# Tool validation and conversion functions are now imported from vcf_utils

# This function has been moved to vcf_utils.rediportal_converter
def convert_rediportal_text_legacy(input_file: Path, output_prefix: str, tool_paths: Dict[str, str]) -> Path:
    """
    Convert REDIportal text format to bcftools annotation format.
    
    Args:
        input_file: Path to input REDIportal file
        output_prefix: Output file prefix
        tool_paths: Dictionary of tool paths
        
    Returns:
        Path to compressed annotation file
        
    Raises:
        RuntimeError: If conversion fails
    """
    logger.info("Starting REDIportal text format conversion...")
    start_time = time.time()
    
    # Determine output paths
    output_file = Path(f"{output_prefix}_annotations.txt")
    compressed_output = Path(str(output_file) + '.gz')
    
    # Safety check - don't overwrite existing files
    if compressed_output.exists():
        logger.warning(f"Output file already exists: {compressed_output}")
        response = input("Overwrite existing file? (y/N): ").strip().lower()
        if response != 'y':
            logger.info("Conversion cancelled by user")
            return compressed_output
        else:
            logger.info("Proceeding with overwrite...")
    
    logger.info(f"Converting to: {compressed_output}")
    
    try:
        # Create temporary file for atomic operation
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt', prefix='rediportal_') as temp_file:
            temp_path = Path(temp_file.name)
            logger.info(f"Using temporary file: {temp_path}")
            
            # Process REDIportal text file
            processed_count = 0
            skipped_count = 0
            error_count = 0
            
            # Open input file
            if str(input_file).endswith('.gz'):
                input_handle = gzip.open(input_file, 'rt')
            else:
                input_handle = open(input_file, 'r')
            
            try:
                logger.info("Processing REDIportal entries...")
                
                for line_num, line in enumerate(input_handle, 1):
                    line = line.strip()
                    
                    # Skip header line
                    if line_num == 1 and (line.startswith('Accession') or 'Region' in line):
                        logger.info("Skipping header line")
                        continue
                    
                    if not line:
                        continue
                    
                    try:
                        fields = line.split('\t')
                        if len(fields) < 10:
                            logger.debug(f"Line {line_num}: insufficient fields ({len(fields)}), skipping")
                            skipped_count += 1
                            continue
                        
                        # Extract required fields with validation
                        accession = fields[REDIPORTAL_COLUMNS['Accession']].strip()
                        region = fields[REDIPORTAL_COLUMNS['Region']].strip()
                        position = fields[REDIPORTAL_COLUMNS['Position']].strip()
                        ref = fields[REDIPORTAL_COLUMNS['Ref']].strip()
                        ed = fields[REDIPORTAL_COLUMNS['Ed']].strip()
                        strand = fields[REDIPORTAL_COLUMNS['Strand']].strip()
                        db = fields[REDIPORTAL_COLUMNS['db']].strip()
                        rna_type = fields[REDIPORTAL_COLUMNS['type']].strip()
                        repeat = fields[REDIPORTAL_COLUMNS['repeat']].strip()
                        func = fields[REDIPORTAL_COLUMNS['Func']].strip()
                        
                        # Validate essential fields
                        if not all([region, position, ref, ed]):
                            logger.debug(f"Line {line_num}: missing essential fields, skipping")
                            skipped_count += 1
                            continue
                        
                        # Parse chromosome and position
                        # In REDIportal format, Region contains chromosome (e.g., "chr1") 
                        # and Position contains the position separately
                        chrom = region.strip()
                        pos = position.strip()
                        
                        # Validate chromosome format
                        if not chrom.startswith('chr'):
                            chrom = f"chr{chrom}"
                        
                        # Validate position is numeric
                        try:
                            int(pos)
                        except ValueError:
                            logger.debug(f"Line {line_num}: invalid position format: {pos}")
                            skipped_count += 1
                            continue
                        
                        # Validate nucleotides
                        if ref not in ['A', 'T', 'G', 'C'] or ed not in ['A', 'T', 'G', 'C']:
                            logger.debug(f"Line {line_num}: invalid nucleotides (REF={ref}, ALT={ed})")
                            skipped_count += 1
                            continue
                        
                        # Clean up fields (replace empty with '.')
                        accession = accession if accession else '.'
                        db = db if db else '.'
                        rna_type = rna_type if rna_type else '.'
                        repeat = repeat if repeat else '.'
                        func = func if func else '.'
                        strand = strand if strand else '.'
                        
                        # Create bcftools annotation format line
                        # Format: CHROM POS REF ALT REDI_ACCESSION REDI_DB REDI_TYPE REDI_REPEAT REDI_FUNC REDI_STRAND
                        annotation_line = f"{chrom}\t{pos}\t{ref}\t{ed}\t{accession}\t{db}\t{rna_type}\t{repeat}\t{func}\t{strand}\n"
                        temp_file.write(annotation_line)
                        processed_count += 1
                        
                        # Progress reporting
                        if processed_count % 100000 == 0:
                            elapsed = time.time() - start_time
                            rate = processed_count / elapsed
                            logger.info(f"Processed {processed_count:,} entries ({rate:.0f} entries/sec)...")
                            
                    except Exception as e:
                        logger.debug(f"Line {line_num}: parsing error: {e}")
                        error_count += 1
                        if error_count > 1000:  # Stop if too many errors
                            logger.error("Too many parsing errors, stopping conversion")
                            raise RuntimeError("Excessive parsing errors in input file")
                        continue
                        
            finally:
                input_handle.close()
        
        conversion_time = time.time() - start_time
        logger.info(f"✓ Conversion completed in {conversion_time:.1f} seconds")
        logger.info(f"  Processed: {processed_count:,} entries")
        logger.info(f"  Skipped: {skipped_count:,} entries")
        logger.info(f"  Errors: {error_count:,} entries")
        
        if processed_count == 0:
            raise RuntimeError("No valid entries were processed")
        
        # Compress with bgzip
        logger.info("Compressing annotation file with bgzip...")
        compress_start = time.time()
        
        with open(compressed_output, 'wb') as f:
            result = subprocess.run(
                [tool_paths['bgzip'], '-c', str(temp_path)], 
                stdout=f, 
                stderr=subprocess.PIPE,
                text=True,
                check=True
            )
        
        compress_time = time.time() - compress_start
        compressed_size = compressed_output.stat().st_size
        logger.info(f"✓ Compression completed in {compress_time:.1f} seconds")
        logger.info(f"  Compressed size: {compressed_size:,} bytes ({compressed_size/1024/1024:.1f} MB)")
        
        # Index with tabix
        logger.info("Creating tabix index...")
        index_start = time.time()
        
        result = subprocess.run(
            [tool_paths['tabix'], '-s1', '-b2', '-e2', str(compressed_output)], 
            capture_output=True,
            text=True,
            check=True
        )
        
        index_time = time.time() - index_start
        index_file = Path(str(compressed_output) + '.tbi')
        index_size = index_file.stat().st_size
        logger.info(f"✓ Indexing completed in {index_time:.1f} seconds")
        logger.info(f"  Index size: {index_size:,} bytes")
        
        # Clean up temporary file
        temp_path.unlink()
        logger.info("✓ Temporary files cleaned up")
        
        # Final validation
        if not compressed_output.exists():
            raise RuntimeError("Compressed output file was not created")
        
        if not index_file.exists():
            raise RuntimeError("Index file was not created")
        
        total_time = time.time() - start_time
        logger.info(f"✓ REDIportal conversion completed successfully in {total_time:.1f} seconds")
        logger.info(f"  Output file: {compressed_output}")
        logger.info(f"  Index file: {index_file}")
        logger.info(f"  Final entries: {processed_count:,}")
        
        return compressed_output
        
    except subprocess.CalledProcessError as e:
        logger.error(f"Command failed during conversion: {e}")
        logger.error(f"Command stderr: {e.stderr}")
        raise RuntimeError(f"REDIportal conversion failed: {e}")
    except Exception as e:
        logger.error(f"Unexpected error during conversion: {e}")
        # Clean up partial files
        if compressed_output.exists():
            try:
                compressed_output.unlink()
                logger.info("Cleaned up partial output file")
            except Exception:
                pass
        if temp_path.exists():
            try:
                temp_path.unlink()
                logger.info("Cleaned up temporary file")
            except Exception:
                pass
        raise RuntimeError(f"REDIportal conversion failed: {e}")

# This function has been moved to vcf_utils.rediportal_converter
def validate_output_legacy(output_file: Path) -> None:
    """
    Validate converted output file.
    
    Args:
        output_file: Path to output file
        
    Raises:
        RuntimeError: If validation fails
    """
    logger.info("Validating converted output...")
    
    # Check files exist
    if not output_file.exists():
        raise RuntimeError(f"Output file does not exist: {output_file}")
    
    index_file = Path(str(output_file) + '.tbi')
    if not index_file.exists():
        raise RuntimeError(f"Index file does not exist: {index_file}")
    
    # Test file can be read
    try:
        result = subprocess.run(['zcat', str(output_file)], capture_output=True, text=True, timeout=30)
        if result.returncode != 0:
            raise RuntimeError("Cannot read compressed file with zcat")
        
        lines = result.stdout.strip().split('\n')
        valid_lines = [line for line in lines if line.strip()]
        
        if len(valid_lines) == 0:
            raise RuntimeError("Output file contains no valid lines")
        
        # Validate format of first few lines
        for i, line in enumerate(valid_lines[:5]):
            fields = line.split('\t')
            if len(fields) < 6:
                raise RuntimeError(f"Line {i+1} has insufficient fields: {len(fields)}")
        
        logger.info(f"✓ Output validation passed: {len(valid_lines):,} lines")
        
    except subprocess.TimeoutExpired:
        logger.warning("Output validation timed out (file may be very large)")
    except Exception as e:
        raise RuntimeError(f"Output validation failed: {e}")

def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Convert REDIportal text format to bcftools annotation format",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic conversion
  python scripts/convert_rediportal_text.py -i TABLE1_hg38_v3.txt.gz -o rediportal_annotations
  
  # With verbose output
  python scripts/convert_rediportal_text.py -i TABLE1_hg38_v3.txt.gz -o rediportal_annotations -v
  
Output files:
  - {output_prefix}_annotations.txt.gz (compressed annotation file)
  - {output_prefix}_annotations.txt.gz.tbi (tabix index)
        """
    )
    
    parser.add_argument(
        '-i', '--input',
        required=True,
        metavar='FILE',
        help='Input REDIportal text file (can be gzipped)'
    )
    
    parser.add_argument(
        '-o', '--output',
        required=True,
        metavar='PREFIX',
        help='Output file prefix'
    )
    
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Enable verbose logging'
    )
    
    args = parser.parse_args()
    
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    try:
        # Validate inputs
        input_file = Path(args.input)
        validate_input_file(input_file)
        
        # Import converter functions from vcf_utils
        from vcf_utils.rediportal_converter import prepare_rediportal_database, validate_converted_file
        
        # Convert file using vcf_utils
        output_file = prepare_rediportal_database(str(input_file), args.output)
        
        # Validate output
        validation = validate_converted_file(output_file)
        if not validation['valid']:
            raise RuntimeError(f"Validation failed: {validation['error']}")
        
        print(f"✓ Conversion completed successfully!")
        print(f"  Output: {output_file}")
        print(f"  Index: {output_file}.tbi")
        
        sys.exit(0)
        
    except KeyboardInterrupt:
        logger.error("Conversion interrupted by user")
        sys.exit(130)
    except Exception as e:
        logger.error(f"Conversion failed: {e}")
        sys.exit(1)

if __name__ == '__main__':
    main()