#!/usr/bin/env python3
"""
COSMIC Database Preparation Script

This script prepares COSMIC VCF database files for annotation by normalizing
chromosome naming conventions and ensuring proper VCF format, bgzip compression,
and tabix indexing for efficient coordinate-based queries with bcftools.

Key Features:
- Chromosome name normalization (1 → chr1, X → chrX, MT → chrM)
- VCF format validation and error reporting
- Proper bgzip compression and tabix indexing
- Field preservation during chromosome mapping
- Comprehensive error handling and logging

Usage:
    python scripts/prepare_cosmic_database.py -i cosmic_input.vcf -o cosmic_prepared.vcf.gz

Requirements Satisfied: 1.1, 1.2, 1.3, 1.4, 1.5

Author: COSMIC/gnomAD Annotation Enhancement
Date: 2025-12-16
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
from typing import Dict, Optional, TextIO, Set

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

# Chromosome mapping for normalization (1 → chr1, X → chrX, MT → chrM)
CHROMOSOME_MAPPING = {
    '1': 'chr1', '2': 'chr2', '3': 'chr3', '4': 'chr4', '5': 'chr5',
    '6': 'chr6', '7': 'chr7', '8': 'chr8', '9': 'chr9', '10': 'chr10',
    '11': 'chr11', '12': 'chr12', '13': 'chr13', '14': 'chr14', '15': 'chr15',
    '16': 'chr16', '17': 'chr17', '18': 'chr18', '19': 'chr19', '20': 'chr20',
    '21': 'chr21', '22': 'chr22', 'X': 'chrX', 'Y': 'chrY', 
    'MT': 'chrM', 'M': 'chrM'
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
    Validate input COSMIC VCF file.
    
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
    
    # Try to read first few lines to validate VCF format
    try:
        if str(file_path).endswith('.gz'):
            with gzip.open(file_path, 'rt') as f:
                first_line = f.readline().strip()
                second_line = f.readline().strip()
        else:
            with open(file_path, 'r') as f:
                first_line = f.readline().strip()
                second_line = f.readline().strip()
        
        # Check if it looks like VCF format
        if not first_line.startswith('##fileformat=VCF'):
            raise ValueError(f"Invalid VCF format. First line should start with '##fileformat=VCF', got: {first_line[:50]}")
        
        logger.info("Detected valid VCF format")
        
        # Check if we can find the header line
        found_header = False
        line_count = 0
        
        if str(file_path).endswith('.gz'):
            with gzip.open(file_path, 'rt') as f:
                for line in f:
                    line_count += 1
                    if line.startswith('#CHROM'):
                        found_header = True
                        break
                    if line_count > 1000:  # Don't read too many lines
                        break
        else:
            with open(file_path, 'r') as f:
                for line in f:
                    line_count += 1
                    if line.startswith('#CHROM'):
                        found_header = True
                        break
                    if line_count > 1000:  # Don't read too many lines
                        break
        
        if not found_header:
            raise ValueError("Could not find VCF column header line (#CHROM)")
        
        logger.info("✓ VCF header structure validated")
        
    except Exception as e:
        raise ValueError(f"Cannot validate VCF format: {e}")

def normalize_chromosome(chrom: str) -> str:
    """
    Normalize chromosome name to standard VCF convention.
    
    Args:
        chrom: Original chromosome name
        
    Returns:
        Normalized chromosome name with 'chr' prefix
    """
    # Remove any existing 'chr' prefix for consistent mapping
    clean_chrom = chrom.replace('chr', '').replace('Chr', '')
    
    # Apply mapping
    if clean_chrom in CHROMOSOME_MAPPING:
        return CHROMOSOME_MAPPING[clean_chrom]
    else:
        # For unmapped chromosomes, add chr prefix if not already present
        if not chrom.startswith('chr'):
            return f"chr{chrom}"
        return chrom

def write_vcf_header(input_file: Path, output_file: TextIO) -> None:
    """
    Copy VCF header from input to output, updating chromosome contigs.
    
    Args:
        input_file: Path to input VCF file
        output_file: Output file handle
    """
    logger.info("Processing VCF header with chromosome normalization")
    
    # Open input file
    if str(input_file).endswith('.gz'):
        input_handle = gzip.open(input_file, 'rt')
    else:
        input_handle = open(input_file, 'r')
    
    try:
        header_lines = 0
        contig_updates = 0
        
        for line in input_handle:
            line = line.strip()
            
            if line.startswith('#CHROM'):
                # This is the column header line - write it and stop
                output_file.write(line + '\n')
                header_lines += 1
                break
            elif line.startswith('##contig=<ID='):
                # Update contig chromosome names
                try:
                    start = line.find('ID=') + 3
                    end = line.find(',', start)
                    if end == -1:
                        end = line.find('>', start)
                    
                    original_chrom = line[start:end]
                    normalized_chrom = normalize_chromosome(original_chrom)
                    
                    if original_chrom != normalized_chrom:
                        updated_line = line.replace(f'ID={original_chrom}', f'ID={normalized_chrom}')
                        output_file.write(updated_line + '\n')
                        contig_updates += 1
                        logger.debug(f"Updated contig: {original_chrom} → {normalized_chrom}")
                    else:
                        output_file.write(line + '\n')
                    
                    header_lines += 1
                except Exception as e:
                    logger.warning(f"Failed to parse contig line: {line[:100]}... Error: {e}")
                    # Write original line if parsing fails
                    output_file.write(line + '\n')
                    header_lines += 1
            else:
                # Regular header line - copy as-is
                output_file.write(line + '\n')
                header_lines += 1
        
        logger.info(f"✓ Processed {header_lines} header lines with {contig_updates} contig updates")
        
    finally:
        input_handle.close()

def convert_cosmic_vcf(input_file: Path, output_file: Path, tool_paths: Dict[str, str]) -> None:
    """
    Convert COSMIC VCF with chromosome normalization and field preservation.
    
    Args:
        input_file: Path to input COSMIC VCF file
        output_file: Path to output VCF file
        tool_paths: Dictionary of tool paths
        
    Raises:
        RuntimeError: If conversion fails
    """
    logger.info("Starting COSMIC VCF conversion with chromosome normalization...")
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
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.vcf', prefix='cosmic_') as temp_file:
            temp_path = Path(temp_file.name)
            logger.info(f"Using temporary file: {temp_path}")
            
            # Write VCF header with chromosome normalization
            write_vcf_header(input_file, temp_file)
            
            # Process VCF records
            processed_count = 0
            skipped_count = 0
            error_count = 0
            chromosome_changes = {}
            
            # Open input file
            if str(input_file).endswith('.gz'):
                input_handle = gzip.open(input_file, 'rt')
            else:
                input_handle = open(input_file, 'r')
            
            try:
                logger.info("Processing COSMIC variant records...")
                
                # Skip header lines
                for line in input_handle:
                    if line.startswith('#CHROM'):
                        break
                
                for line_num, line in enumerate(input_handle, 1):
                    line = line.strip()
                    
                    if not line or line.startswith('#'):
                        continue
                    
                    try:
                        fields = line.split('\t')
                        if len(fields) < 8:
                            logger.debug(f"Line {line_num}: insufficient fields ({len(fields)}), skipping")
                            skipped_count += 1
                            continue
                        
                        # Extract and normalize chromosome
                        original_chrom = fields[0]
                        normalized_chrom = normalize_chromosome(original_chrom)
                        
                        # Track chromosome changes for statistics
                        if original_chrom != normalized_chrom:
                            if original_chrom not in chromosome_changes:
                                chromosome_changes[original_chrom] = normalized_chrom
                        
                        # Update chromosome field
                        fields[0] = normalized_chrom
                        
                        # Validate essential fields
                        pos = fields[1]
                        ref = fields[3]
                        alt = fields[4]
                        
                        # Validate position is numeric
                        try:
                            pos_int = int(pos)
                            if pos_int <= 0:
                                raise ValueError("Position must be positive")
                        except ValueError:
                            logger.debug(f"Line {line_num}: invalid position format: {pos}")
                            skipped_count += 1
                            continue
                        
                        # Validate nucleotides (allow complex variants)
                        if not ref or not alt or ref == '.' or alt == '.':
                            logger.debug(f"Line {line_num}: invalid REF/ALT fields (REF={ref}, ALT={alt})")
                            skipped_count += 1
                            continue
                        
                        # Write normalized VCF record
                        normalized_record = '\t'.join(fields) + '\n'
                        temp_file.write(normalized_record)
                        processed_count += 1
                        
                        # Progress reporting
                        if processed_count % 100000 == 0:
                            elapsed = time.time() - start_time
                            rate = processed_count / elapsed
                            logger.info(f"Processed {processed_count:,} records ({rate:.0f} records/sec)...")
                            
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
        logger.info(f"  Processed: {processed_count:,} records")
        logger.info(f"  Skipped: {skipped_count:,} records")
        logger.info(f"  Errors: {error_count:,} records")
        
        if chromosome_changes:
            logger.info(f"  Chromosome mappings applied: {len(chromosome_changes)}")
            for orig, norm in sorted(chromosome_changes.items()):
                logger.info(f"    {orig} → {norm}")
        
        if processed_count == 0:
            raise RuntimeError("No valid records were processed")
        
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
        logger.info(f"✓ COSMIC database preparation completed successfully in {total_time:.1f} seconds")
        logger.info(f"  Output VCF: {output_file}")
        logger.info(f"  Index file: {index_file}")
        logger.info(f"  Final records: {processed_count:,}")
        
    except subprocess.CalledProcessError as e:
        logger.error(f"Command failed during conversion: {e}")
        logger.error(f"Command stderr: {e.stderr}")
        raise RuntimeError(f"COSMIC database preparation failed: {e}")
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
        raise RuntimeError(f"COSMIC database preparation failed: {e}")

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
        normalized_chromosomes = set()
        for i, line in enumerate(variant_lines[:5]):
            fields = line.split('\t')
            if len(fields) < 8:
                raise RuntimeError(f"Variant line {i+1} has insufficient fields: {len(fields)}")
            
            # Check chromosome format (should be normalized)
            chrom = fields[0]
            normalized_chromosomes.add(chrom)
            if not chrom.startswith('chr'):
                raise RuntimeError(f"Invalid chromosome format in line {i+1}: {chrom} (should start with 'chr')")
            
            # Check position is numeric
            try:
                int(fields[1])
            except ValueError:
                raise RuntimeError(f"Invalid position in line {i+1}: {fields[1]}")
            
            # Check REF and ALT are not empty
            if not fields[3] or fields[3] == '.':
                raise RuntimeError(f"Invalid REF field in line {i+1}: {fields[3]}")
            if not fields[4] or fields[4] == '.':
                raise RuntimeError(f"Invalid ALT field in line {i+1}: {fields[4]}")
        
        logger.info(f"✓ VCF validation passed: {len(variant_lines):,} variant records")
        logger.info(f"  Normalized chromosomes found: {sorted(normalized_chromosomes)}")
        
    except subprocess.TimeoutExpired:
        logger.warning("VCF validation timed out (file may be very large)")
    except Exception as e:
        raise RuntimeError(f"VCF validation failed: {e}")

def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Prepare COSMIC VCF database for annotation with chromosome normalization",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic conversion
  python scripts/prepare_cosmic_database.py -i cosmic_input.vcf -o cosmic_prepared.vcf.gz
  
  # With verbose output
  python scripts/prepare_cosmic_database.py -i cosmic_input.vcf.gz -o cosmic_prepared.vcf.gz -v
  
Chromosome Normalization:
  - 1 → chr1, 2 → chr2, ..., 22 → chr22
  - X → chrX, Y → chrY
  - MT → chrM, M → chrM
  
Output files:
  - {output}.vcf.gz (compressed VCF file)
  - {output}.vcf.gz.tbi (tabix index)
        """
    )
    
    parser.add_argument(
        '-i', '--input',
        required=True,
        metavar='FILE',
        help='Input COSMIC VCF file (can be gzipped)'
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
        convert_cosmic_vcf(input_file, output_file, tool_paths)
        
        # Validate output
        validate_vcf_output(output_file)
        
        print("✓ COSMIC database preparation completed successfully!")
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