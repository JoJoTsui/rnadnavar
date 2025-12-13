#!/usr/bin/env python3
"""
REDIportal Database to VCF Converter

This script converts REDIportal TABLE1_hg38_v3.txt.gz format to complete VCF format
with proper headers and variant records. The output can be used directly with
bcftools annotate as a VCF database without additional header files.

Key Features:
- Converts REDIportal text format to complete VCF format
- Generates proper VCF headers with INFO field definitions
- Creates bgzip-compressed VCF with tabix indexing
- Handles large files efficiently with streaming processing
- Comprehensive error handling and validation

Usage:
    python scripts/convert_rediportal_to_vcf.py -i TABLE1_hg38_v3.txt.gz -o rediportal.vcf.gz

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
from typing import Dict, Optional, TextIO

# Add vcf_utils to path for imports
script_dir = Path(__file__).parent.parent / 'bin'
sys.path.insert(0, str(script_dir))

from vcf_utils.field_mapping import REDIportalFieldMapper

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

def validate_tools() -> Dict[str, str]:
    """
    Validate required tools are available.
    
    Returns:
        Dictionary mapping tool names to their paths
        
    Raises:
        RuntimeError: If required tools are missing
    """
    required_tools = ['bgzip', 'tabix']
    tool_paths = {}
    
    for tool in required_tools:
        try:
            result = subprocess.run(['which', tool], capture_output=True, text=True)
            if result.returncode == 0:
                tool_paths[tool] = result.stdout.strip()
                logger.debug(f"Found {tool}: {tool_paths[tool]}")
            else:
                raise RuntimeError(f"Required tool not found: {tool}")
        except Exception as e:
            raise RuntimeError(f"Error checking for {tool}: {e}")
    
    logger.info("✓ All required tools are available")
    return tool_paths

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

def write_vcf_header(output_file: TextIO) -> None:
    """
    Write complete VCF header with INFO field definitions.
    
    Args:
        output_file: Output file handle
    """
    logger.info("Writing VCF header with INFO field definitions")
    
    # VCF format version
    output_file.write("##fileformat=VCFv4.2\n")
    
    # Source information
    output_file.write("##source=REDIportal_Database_Converter\n")
    output_file.write(f"##fileDate={time.strftime('%Y%m%d')}\n")
    
    # Reference genome
    output_file.write("##reference=hg38\n")
    
    # INFO field definitions for REDIportal data
    info_definitions = REDIportalFieldMapper.get_info_field_definitions()
    for field_name in ['REDI_ACCESSION', 'REDI_DB', 'REDI_TYPE', 'REDI_REPEAT', 'REDI_STRAND']:
        output_file.write(f"{info_definitions[field_name]}\n")
    
    # Add REDI_FUNC which is not part of the main mapping but still used
    output_file.write('##INFO=<ID=REDI_FUNC,Number=1,Type=String,Description="REDIportal functional annotation">\n')
    
    # Column header line
    output_file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

# Field mapping functions are now provided by REDIportalFieldMapper

def convert_rediportal_to_vcf(input_file: Path, output_file: Path, tool_paths: Dict[str, str]) -> None:
    """
    Convert REDIportal text format to complete VCF format.
    
    Args:
        input_file: Path to input REDIportal file
        output_file: Path to output VCF file
        tool_paths: Dictionary of tool paths
        
    Raises:
        RuntimeError: If conversion fails
    """
    logger.info("Starting REDIportal to VCF conversion...")
    start_time = time.time()
    
    # Safety check - don't overwrite existing files
    if output_file.exists():
        logger.warning(f"Output file already exists: {output_file}")
        response = input("Overwrite existing file? (y/N): ").strip().lower()
        if response != 'y':
            logger.info("Conversion cancelled by user")
            return
        else:
            logger.info("Proceeding with overwrite...")
    
    logger.info(f"Converting to: {output_file}")
    
    try:
        # Create temporary file for atomic operation
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.vcf', prefix='rediportal_') as temp_file:
            temp_path = Path(temp_file.name)
            logger.info(f"Using temporary file: {temp_path}")
            
            # Write VCF header
            write_vcf_header(temp_file)
            
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
                        chrom = region.strip()
                        pos = position.strip()
                        
                        # Validate chromosome format
                        if not chrom.startswith('chr'):
                            chrom = f"chr{chrom}"
                        
                        # Validate position is numeric
                        try:
                            pos_int = int(pos)
                            if pos_int <= 0:
                                raise ValueError("Position must be positive")
                        except ValueError:
                            logger.debug(f"Line {line_num}: invalid position format: {pos}")
                            skipped_count += 1
                            continue
                        
                        # Validate nucleotides
                        if ref not in ['A', 'T', 'G', 'C'] or ed not in ['A', 'T', 'G', 'C']:
                            logger.debug(f"Line {line_num}: invalid nucleotides (REF={ref}, ALT={ed})")
                            skipped_count += 1
                            continue
                        
                        # Map fields using the field mapping system
                        field_mappings = REDIportalFieldMapper.map_all_fields(
                            accession, db, rna_type, repeat, strand
                        )
                        
                        # Add REDI_FUNC which is not part of the main mapping
                        func_escaped = REDIportalFieldMapper.escape_info_value(func) if func else '.'
                        if func_escaped != '.':
                            field_mappings['REDI_FUNC'] = func_escaped
                        
                        # Build INFO field string
                        info_field = REDIportalFieldMapper.build_info_string(field_mappings)
                        
                        # Create VCF record
                        # Format: CHROM POS ID REF ALT QUAL FILTER INFO
                        vcf_record = f"{chrom}\t{pos}\t.\t{ref}\t{ed}\t.\tPASS\t{info_field}\n"
                        temp_file.write(vcf_record)
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
        logger.info(f"✓ VCF conversion completed in {conversion_time:.1f} seconds")
        logger.info(f"  Processed: {processed_count:,} entries")
        logger.info(f"  Skipped: {skipped_count:,} entries")
        logger.info(f"  Errors: {error_count:,} entries")
        
        if processed_count == 0:
            raise RuntimeError("No valid entries were processed")
        
        # Compress with bgzip
        logger.info("Compressing VCF file with bgzip...")
        compress_start = time.time()
        
        with open(output_file, 'wb') as f:
            result = subprocess.run(
                [tool_paths['bgzip'], '-c', str(temp_path)], 
                stdout=f, 
                stderr=subprocess.PIPE,
                check=True
            )
        
        compress_time = time.time() - compress_start
        compressed_size = output_file.stat().st_size
        logger.info(f"✓ Compression completed in {compress_time:.1f} seconds")
        logger.info(f"  Compressed size: {compressed_size:,} bytes ({compressed_size/1024/1024:.1f} MB)")
        
        # Index with tabix
        logger.info("Creating tabix index...")
        index_start = time.time()
        
        result = subprocess.run(
            [tool_paths['tabix'], '-p', 'vcf', str(output_file)], 
            capture_output=True,
            text=True,
            check=True
        )
        
        index_time = time.time() - index_start
        index_file = Path(str(output_file) + '.tbi')
        index_size = index_file.stat().st_size
        logger.info(f"✓ Indexing completed in {index_time:.1f} seconds")
        logger.info(f"  Index size: {index_size:,} bytes")
        
        # Clean up temporary file
        temp_path.unlink()
        logger.info("✓ Temporary files cleaned up")
        
        # Final validation
        if not output_file.exists():
            raise RuntimeError("Compressed VCF output file was not created")
        
        if not index_file.exists():
            raise RuntimeError("Index file was not created")
        
        total_time = time.time() - start_time
        logger.info(f"✓ REDIportal to VCF conversion completed successfully in {total_time:.1f} seconds")
        logger.info(f"  Output VCF: {output_file}")
        logger.info(f"  Index file: {index_file}")
        logger.info(f"  Final entries: {processed_count:,}")
        
    except subprocess.CalledProcessError as e:
        logger.error(f"Command failed during conversion: {e}")
        logger.error(f"Command stderr: {e.stderr}")
        raise RuntimeError(f"REDIportal to VCF conversion failed: {e}")
    except Exception as e:
        logger.error(f"Unexpected error during conversion: {e}")
        # Clean up partial files
        if output_file.exists():
            try:
                output_file.unlink()
                logger.info("Cleaned up partial output file")
            except Exception:
                pass
        if 'temp_path' in locals() and temp_path.exists():
            try:
                temp_path.unlink()
                logger.info("Cleaned up temporary file")
            except Exception:
                pass
        raise RuntimeError(f"REDIportal to VCF conversion failed: {e}")

def validate_vcf_output(output_file: Path) -> None:
    """
    Validate converted VCF output file.
    
    Args:
        output_file: Path to output VCF file
        
    Raises:
        RuntimeError: If validation fails
    """
    logger.info("Validating converted VCF output...")
    
    # Check files exist
    if not output_file.exists():
        raise RuntimeError(f"Output VCF file does not exist: {output_file}")
    
    index_file = Path(str(output_file) + '.tbi')
    if not index_file.exists():
        raise RuntimeError(f"Index file does not exist: {index_file}")
    
    # Test file can be read and has proper VCF format
    try:
        result = subprocess.run(['zcat', str(output_file)], capture_output=True, text=True, timeout=30)
        if result.returncode != 0:
            raise RuntimeError("Cannot read compressed VCF file with zcat")
        
        lines = result.stdout.strip().split('\n')
        
        # Check VCF header
        if not lines[0].startswith('##fileformat=VCF'):
            raise RuntimeError("Missing VCF format header")
        
        # Find column header line
        header_line = None
        variant_lines = []
        for line in lines:
            if line.startswith('#CHROM'):
                header_line = line
            elif not line.startswith('#') and line.strip():
                variant_lines.append(line)
        
        if not header_line:
            raise RuntimeError("Missing VCF column header line")
        
        if len(variant_lines) == 0:
            raise RuntimeError("VCF file contains no variant records")
        
        # Validate format of first few variant lines
        for i, line in enumerate(variant_lines[:5]):
            fields = line.split('\t')
            if len(fields) < 8:
                raise RuntimeError(f"Variant line {i+1} has insufficient fields: {len(fields)}")
            
            # Check chromosome format
            if not fields[0].startswith('chr'):
                raise RuntimeError(f"Invalid chromosome format in line {i+1}: {fields[0]}")
            
            # Check position is numeric
            try:
                int(fields[1])
            except ValueError:
                raise RuntimeError(f"Invalid position in line {i+1}: {fields[1]}")
            
            # Check REF and ALT are valid nucleotides
            if fields[3] not in ['A', 'T', 'G', 'C']:
                raise RuntimeError(f"Invalid REF nucleotide in line {i+1}: {fields[3]}")
            if fields[4] not in ['A', 'T', 'G', 'C']:
                raise RuntimeError(f"Invalid ALT nucleotide in line {i+1}: {fields[4]}")
        
        logger.info(f"✓ VCF validation passed: {len(variant_lines):,} variant records")
        
    except subprocess.TimeoutExpired:
        logger.warning("VCF validation timed out (file may be very large)")
    except Exception as e:
        raise RuntimeError(f"VCF validation failed: {e}")

def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Convert REDIportal text format to complete VCF format",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic conversion
  python scripts/convert_rediportal_to_vcf.py -i TABLE1_hg38_v3.txt.gz -o rediportal.vcf.gz
  
  # With verbose output
  python scripts/convert_rediportal_to_vcf.py -i TABLE1_hg38_v3.txt.gz -o rediportal.vcf.gz -v
  
Output files:
  - {output}.vcf.gz (compressed VCF file)
  - {output}.vcf.gz.tbi (tabix index)
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
        metavar='FILE',
        help='Output VCF file (will be bgzip compressed)'
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
        # Validate tools
        tool_paths = validate_tools()
        
        # Validate inputs
        input_file = Path(args.input)
        validate_input_file(input_file)
        
        # Ensure output has .gz extension
        output_file = Path(args.output)
        if not str(output_file).endswith('.gz'):
            output_file = Path(str(output_file) + '.gz')
        
        # Convert file
        convert_rediportal_to_vcf(input_file, output_file, tool_paths)
        
        # Validate output
        validate_vcf_output(output_file)
        
        print("✓ Conversion completed successfully!")
        print(f"  Output VCF: {output_file}")
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