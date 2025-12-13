#!/usr/bin/env python3
"""
REDIportal Database Converter

This module provides functionality to detect and convert REDIportal database formats
for use with bcftools annotation. It supports both VCF and gzipped text formats
from the REDIportal database.

SAFETY FEATURES:
- Comprehensive input validation
- No destructive operations on original files
- Extensive error handling and logging
- Atomic file operations with temporary files

Author: RNA Editing Enhancement Pipeline
Date: 2025-12-13
"""

import gzip
import logging
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, Optional, Union

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

def detect_rediportal_format(file_path: str) -> str:
    """
    Detect REDIportal database format (VCF or gzipped text).
    
    Args:
        file_path: Path to REDIportal database file
        
    Returns:
        'vcf' for VCF format, 'text' for gzipped text format
        
    Raises:
        FileNotFoundError: If file doesn't exist
        ValueError: If format cannot be determined
    """
    file_path = Path(file_path)
    
    if not file_path.exists():
        raise FileNotFoundError(f"REDIportal file not found: {file_path}")
    
    logger.info(f"Detecting format of REDIportal file: {file_path}")
    
    try:
        # Try to read first few lines to determine format
        if str(file_path).endswith('.gz'):
            with gzip.open(file_path, 'rt') as f:
                first_line = f.readline().strip()
        else:
            with open(file_path, 'r') as f:
                first_line = f.readline().strip()
        
        if first_line.startswith('##fileformat=VCF'):
            logger.info("Detected VCF format")
            return 'vcf'
        elif first_line.startswith('Accession') or 'Region' in first_line:
            logger.info("Detected REDIportal text format")
            return 'text'
        else:
            # Try to parse as tab-delimited to see if it matches expected columns
            fields = first_line.split('\t')
            if len(fields) >= 10:
                logger.info("Detected REDIportal text format (no header)")
                return 'text'
            else:
                raise ValueError(f"Unknown REDIportal format. First line: {first_line[:100]}...")
                
    except Exception as e:
        logger.error(f"Failed to detect REDIportal format: {e}")
        raise ValueError(f"Cannot determine REDIportal format for {file_path}: {e}")

def prepare_rediportal_database(file_path: str, output_prefix: Optional[str] = None) -> str:
    """
    Prepare REDIportal database for bcftools annotation.
    
    Args:
        file_path: Path to REDIportal database file
        output_prefix: Optional output prefix (default: auto-generated)
        
    Returns:
        Path to prepared annotation file
        
    Raises:
        FileNotFoundError: If input file doesn't exist
        RuntimeError: If conversion fails
    """
    file_path = Path(file_path)
    
    if not file_path.exists():
        raise FileNotFoundError(f"REDIportal file not found: {file_path}")
    
    format_type = detect_rediportal_format(str(file_path))
    
    if format_type == 'vcf':
        return _prepare_vcf_format(file_path, output_prefix)
    elif format_type == 'text':
        return _prepare_text_format(file_path, output_prefix)
    else:
        raise ValueError(f"Unsupported format: {format_type}")

def _prepare_vcf_format(file_path: Path, output_prefix: Optional[str] = None) -> str:
    """
    Prepare VCF format REDIportal database.
    
    Args:
        file_path: Path to VCF REDIportal file
        output_prefix: Optional output prefix
        
    Returns:
        Path to prepared VCF file
    """
    logger.info(f"Preparing VCF format REDIportal database: {file_path}")
    
    # For VCF format, we just need to ensure it's compressed and indexed
    if str(file_path).endswith('.gz'):
        # Check if index exists
        index_file = Path(str(file_path) + '.tbi')
        if index_file.exists():
            logger.info("VCF is already compressed and indexed")
            return str(file_path)
        else:
            logger.info("Creating index for compressed VCF...")
            try:
                subprocess.run(['tabix', '-p', 'vcf', str(file_path)], check=True)
                return str(file_path)
            except subprocess.CalledProcessError as e:
                raise RuntimeError(f"Failed to index VCF: {e}")
    else:
        # Need to compress and index
        logger.info("Compressing and indexing VCF...")
        if output_prefix:
            compressed_file = Path(f"{output_prefix}.vcf.gz")
        else:
            compressed_file = Path(str(file_path) + '.gz')
        
        try:
            # Compress with bgzip
            with open(compressed_file, 'wb') as f:
                subprocess.run(['bgzip', '-c', str(file_path)], stdout=f, check=True)
            
            # Index with tabix
            subprocess.run(['tabix', '-p', 'vcf', str(compressed_file)], check=True)
            
            logger.info(f"✓ VCF prepared: {compressed_file}")
            return str(compressed_file)
            
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"Failed to prepare VCF: {e}")

def _prepare_text_format(file_path: Path, output_prefix: Optional[str] = None) -> str:
    """
    Convert REDIportal text format to bcftools annotation format.
    
    Args:
        file_path: Path to gzipped text REDIportal file
        output_prefix: Optional output prefix
        
    Returns:
        Path to prepared annotation file
    """
    logger.info(f"Converting REDIportal text format: {file_path}")
    
    if output_prefix:
        output_file = Path(f"{output_prefix}_annotations.txt")
    else:
        output_file = file_path.parent / f"{file_path.stem}_annotations.txt"
    
    compressed_output = Path(str(output_file) + '.gz')
    
    # Safety check - don't overwrite existing files without explicit confirmation
    if compressed_output.exists():
        logger.warning(f"Output file already exists: {compressed_output}")
        logger.info("Using existing converted file")
        return str(compressed_output)
    
    logger.info(f"Converting to: {compressed_output}")
    
    try:
        # Create temporary file for atomic operation
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as temp_file:
            temp_path = Path(temp_file.name)
            
            # Process REDIportal text file
            processed_count = 0
            skipped_count = 0
            
            if str(file_path).endswith('.gz'):
                input_file = gzip.open(file_path, 'rt')
            else:
                input_file = open(file_path, 'r')
            
            try:
                for line_num, line in enumerate(input_file, 1):
                    line = line.strip()
                    
                    # Skip header line
                    if line_num == 1 and (line.startswith('Accession') or 'Region' in line):
                        logger.debug("Skipping header line")
                        continue
                    
                    if not line:
                        continue
                    
                    try:
                        fields = line.split('\t')
                        if len(fields) < 10:
                            logger.warning(f"Line {line_num}: insufficient fields ({len(fields)}), skipping")
                            skipped_count += 1
                            continue
                        
                        # Extract required fields
                        accession = fields[REDIPORTAL_COLUMNS['Accession']]
                        region = fields[REDIPORTAL_COLUMNS['Region']]
                        position = fields[REDIPORTAL_COLUMNS['Position']]
                        ref = fields[REDIPORTAL_COLUMNS['Ref']]
                        ed = fields[REDIPORTAL_COLUMNS['Ed']]
                        strand = fields[REDIPORTAL_COLUMNS['Strand']]
                        db = fields[REDIPORTAL_COLUMNS['db']]
                        rna_type = fields[REDIPORTAL_COLUMNS['type']]
                        repeat = fields[REDIPORTAL_COLUMNS['repeat']]
                        func = fields[REDIPORTAL_COLUMNS['Func']]
                        
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
                            logger.warning(f"Line {line_num}: invalid position format: {pos}")
                            skipped_count += 1
                            continue
                        
                        # Validate nucleotides
                        if ref not in ['A', 'T', 'G', 'C'] or ed not in ['A', 'T', 'G', 'C']:
                            logger.debug(f"Line {line_num}: invalid nucleotides (REF={ref}, ALT={ed})")
                            skipped_count += 1
                            continue
                        
                        # Clean up fields (replace empty with '.')
                        accession = accession.strip() if accession.strip() else '.'
                        db = db.strip() if db.strip() else '.'
                        rna_type = rna_type.strip() if rna_type.strip() else '.'
                        repeat = repeat.strip() if repeat.strip() else '.'
                        func = func.strip() if func.strip() else '.'
                        strand = strand.strip() if strand.strip() else '.'
                        
                        # Create bcftools annotation format line
                        # Format: CHROM POS REF ALT REDI_ACCESSION REDI_DB REDI_TYPE REDI_REPEAT REDI_FUNC REDI_STRAND
                        annotation_line = f"{chrom}\t{pos}\t{ref}\t{ed}\t{accession}\t{db}\t{rna_type}\t{repeat}\t{func}\t{strand}\n"
                        temp_file.write(annotation_line)
                        processed_count += 1
                        
                        if processed_count % 100000 == 0:
                            logger.info(f"Processed {processed_count:,} entries...")
                            
                    except Exception as e:
                        logger.warning(f"Line {line_num}: parsing error: {e}")
                        skipped_count += 1
                        continue
                        
            finally:
                input_file.close()
        
        logger.info(f"✓ Conversion completed: {processed_count:,} entries processed, {skipped_count:,} skipped")
        
        # Compress with bgzip
        logger.info("Compressing annotation file...")
        with open(compressed_output, 'wb') as f:
            subprocess.run(['bgzip', '-c', str(temp_path)], stdout=f, check=True)
        
        # Index with tabix
        logger.info("Creating tabix index...")
        subprocess.run(['tabix', '-s1', '-b2', '-e2', str(compressed_output)], check=True)
        
        # Clean up temporary file
        temp_path.unlink()
        
        # Verify output
        if not compressed_output.exists():
            raise RuntimeError("Compressed output file was not created")
        
        index_file = Path(str(compressed_output) + '.tbi')
        if not index_file.exists():
            raise RuntimeError("Index file was not created")
        
        output_size = compressed_output.stat().st_size
        logger.info("✓ REDIportal text format converted successfully")
        logger.info(f"  Output: {compressed_output} ({output_size:,} bytes)")
        logger.info(f"  Index: {index_file}")
        logger.info(f"  Entries: {processed_count:,} processed, {skipped_count:,} skipped")
        
        return str(compressed_output)
        
    except subprocess.CalledProcessError as e:
        logger.error(f"Command failed during conversion: {e}")
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
        raise RuntimeError(f"REDIportal conversion failed: {e}")

def validate_converted_file(file_path: str) -> Dict[str, Union[int, bool]]:
    """
    Validate converted REDIportal annotation file.
    
    Args:
        file_path: Path to converted annotation file
        
    Returns:
        Dictionary with validation results
    """
    file_path = Path(file_path)
    
    if not file_path.exists():
        return {'valid': False, 'error': 'File does not exist'}
    
    # Check index exists
    index_file = Path(str(file_path) + '.tbi')
    if not index_file.exists():
        return {'valid': False, 'error': 'Index file missing'}
    
    try:
        # Count lines in compressed file
        result = subprocess.run(['zcat', str(file_path)], capture_output=True, text=True)
        if result.returncode != 0:
            return {'valid': False, 'error': 'Cannot read compressed file'}
        
        lines = result.stdout.strip().split('\n')
        line_count = len([line for line in lines if line.strip()])
        
        # Validate format of first few lines
        valid_lines = 0
        for i, line in enumerate(lines[:10]):
            if not line.strip():
                continue
            fields = line.split('\t')
            if len(fields) >= 6:  # Minimum required fields
                valid_lines += 1
        
        return {
            'valid': True,
            'line_count': line_count,
            'valid_sample_lines': valid_lines,
            'file_size': file_path.stat().st_size,
            'index_size': index_file.stat().st_size
        }
        
    except Exception as e:
        return {'valid': False, 'error': str(e)}

if __name__ == '__main__':
    # Simple CLI for testing
    import argparse
    
    parser = argparse.ArgumentParser(description='Convert REDIportal database for bcftools')
    parser.add_argument('-i', '--input', required=True, help='Input REDIportal file')
    parser.add_argument('-o', '--output', help='Output prefix')
    parser.add_argument('-v', '--verbose', action='store_true', help='Verbose logging')
    
    args = parser.parse_args()
    
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    
    try:
        output_file = prepare_rediportal_database(args.input, args.output)
        print(f"✓ Conversion completed: {output_file}")
        
        # Validate output
        validation = validate_converted_file(output_file)
        if validation['valid']:
            print(f"✓ Validation passed: {validation['line_count']:,} entries")
        else:
            print(f"✗ Validation failed: {validation['error']}")
            
    except Exception as e:
        print(f"✗ Conversion failed: {e}")
        exit(1)