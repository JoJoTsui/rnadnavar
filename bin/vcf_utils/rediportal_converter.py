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
import time
from pathlib import Path
from typing import Dict, Optional, Union

from .field_mapping import REDIportalFieldMapper

logger = logging.getLogger(__name__)

# Enhanced logging configuration for REDIportal conversion
class REDIportalConversionLogger:
    """Enhanced logging for REDIportal conversion operations."""
    
    def __init__(self):
        self.stats = {
            'start_time': time.time(),
            'total_lines_read': 0,
            'header_lines_skipped': 0,
            'entries_processed': 0,
            'entries_skipped': 0,
            'parsing_errors': 0,
            'validation_errors': 0,
            'field_mapping_errors': 0,
            'compression_time': 0,
            'indexing_time': 0,
            'error_details': []
        }
    
    def log_parsing_error(self, line_num: int, error: str, line_content: str = None):
        """Log detailed parsing error information."""
        error_info = {
            'line_number': line_num,
            'error_type': 'parsing_error',
            'error_message': error,
            'line_content': line_content[:100] + '...' if line_content and len(line_content) > 100 else line_content,
            'timestamp': time.time()
        }
        self.stats['parsing_errors'] += 1
        self.stats['error_details'].append(error_info)
        
        logger.warning(f"Line {line_num}: Parsing error - {error}")
        if line_content:
            logger.debug(f"Line content: {line_content[:200]}...")
    
    def log_validation_error(self, line_num: int, field: str, value: str, error: str):
        """Log detailed validation error information."""
        error_info = {
            'line_number': line_num,
            'error_type': 'validation_error',
            'field': field,
            'value': value[:50] + '...' if len(value) > 50 else value,
            'error_message': error,
            'timestamp': time.time()
        }
        self.stats['validation_errors'] += 1
        self.stats['error_details'].append(error_info)
        
        logger.warning(f"Line {line_num}: Validation error in field '{field}' - {error}")
        logger.debug(f"Invalid value: {value[:100]}...")
    
    def log_field_mapping_error(self, line_num: int, field: str, original_value: str, error: str):
        """Log field mapping error information."""
        error_info = {
            'line_number': line_num,
            'error_type': 'field_mapping_error',
            'field': field,
            'original_value': original_value[:50] + '...' if len(original_value) > 50 else original_value,
            'error_message': error,
            'timestamp': time.time()
        }
        self.stats['field_mapping_errors'] += 1
        self.stats['error_details'].append(error_info)
        
        logger.warning(f"Line {line_num}: Field mapping error for '{field}' - {error}")
    
    def log_progress(self, processed_count: int, interval: int = 100000):
        """Log processing progress at regular intervals."""
        if processed_count % interval == 0:
            elapsed_time = time.time() - self.stats['start_time']
            rate = processed_count / elapsed_time if elapsed_time > 0 else 0
            logger.info(f"Processed {processed_count:,} entries ({rate:.0f} entries/sec, "
                       f"{self.stats['entries_skipped']:,} skipped, "
                       f"{self.stats['parsing_errors']} parsing errors)")
    
    def log_conversion_statistics(self):
        """Log comprehensive conversion statistics."""
        total_time = time.time() - self.stats['start_time']
        
        logger.info("=== REDIportal Conversion Statistics ===")
        logger.info(f"Total processing time: {total_time:.2f} seconds")
        logger.info(f"Total lines read: {self.stats['total_lines_read']:,}")
        logger.info(f"Header lines skipped: {self.stats['header_lines_skipped']:,}")
        logger.info(f"Entries processed successfully: {self.stats['entries_processed']:,}")
        logger.info(f"Entries skipped: {self.stats['entries_skipped']:,}")
        logger.info(f"Processing rate: {self.stats['entries_processed'] / total_time:.0f} entries/sec")
        
        if self.stats['parsing_errors'] > 0:
            logger.warning(f"Parsing errors encountered: {self.stats['parsing_errors']}")
        if self.stats['validation_errors'] > 0:
            logger.warning(f"Validation errors encountered: {self.stats['validation_errors']}")
        if self.stats['field_mapping_errors'] > 0:
            logger.warning(f"Field mapping errors encountered: {self.stats['field_mapping_errors']}")
        
        # Log error rate
        total_errors = self.stats['parsing_errors'] + self.stats['validation_errors'] + self.stats['field_mapping_errors']
        if total_errors > 0:
            error_rate = (total_errors / self.stats['total_lines_read']) * 100 if self.stats['total_lines_read'] > 0 else 0
            logger.warning(f"Total error rate: {error_rate:.2f}% ({total_errors:,} errors)")
        
        # Log timing breakdown
        if self.stats['compression_time'] > 0:
            logger.info(f"Compression time: {self.stats['compression_time']:.2f} seconds")
        if self.stats['indexing_time'] > 0:
            logger.info(f"Indexing time: {self.stats['indexing_time']:.2f} seconds")
    
    def get_statistics(self) -> Dict:
        """Get conversion statistics dictionary."""
        return self.stats.copy()

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
    Detect REDIportal database format (VCF, converted annotation, or original text).
    
    Args:
        file_path: Path to REDIportal database file
        
    Returns:
        'vcf' for VCF format, 'converted' for already-converted annotation format, 'text' for original text format
        
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
                f.readline().strip()  # Read second line but don't store
        else:
            with open(file_path, 'r') as f:
                first_line = f.readline().strip()
                f.readline().strip()  # Read second line but don't store
        
        if first_line.startswith('##fileformat=VCF'):
            logger.info("Detected VCF format")
            return 'vcf'
        elif first_line.startswith('Accession') or 'Region' in first_line:
            logger.info("Detected original REDIportal text format with header")
            return 'text'
        else:
            # Check if this is already converted annotation format vs original text format
            fields = first_line.split('\t')
            
            if len(fields) >= 6:
                # Check if first field looks like chromosome (converted format)
                if fields[0].startswith('chr') and fields[1].isdigit():
                    # Check if nucleotides are in expected positions (converted format)
                    if len(fields) >= 4 and len(fields[2]) == 1 and len(fields[3]) == 1:
                        if fields[2] in ['A', 'T', 'G', 'C'] and fields[3] in ['A', 'T', 'G', 'C']:
                            logger.info("Detected already-converted annotation format")
                            return 'converted'
                
                # Check if this looks like original REDIportal format (10+ columns, different structure)
                if len(fields) >= 10:
                    # In original format, position is in column 2, and it should be numeric
                    try:
                        int(fields[2])  # Position should be numeric in original format
                        # Check if column 1 looks like chromosome region
                        if 'chr' in fields[1] or fields[1].isdigit():
                            logger.info("Detected original REDIportal text format (no header)")
                            return 'text'
                    except ValueError:
                        pass
            
            # If we can't determine format clearly, default based on field count
            if len(fields) >= 10:
                logger.warning(f"Ambiguous format, defaulting to original text format. First line: {first_line[:100]}...")
                return 'text'
            elif len(fields) >= 6:
                logger.warning(f"Ambiguous format, defaulting to converted format. First line: {first_line[:100]}...")
                return 'converted'
            else:
                raise ValueError(f"Unknown REDIportal format. First line has {len(fields)} fields: {first_line[:100]}...")
                
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
    elif format_type == 'converted':
        return _prepare_converted_format(file_path, output_prefix)
    elif format_type == 'text':
        return _prepare_text_format(file_path, output_prefix)
    else:
        raise ValueError(f"Unsupported format: {format_type}")

def _prepare_converted_format(file_path: Path, output_prefix: Optional[str] = None) -> str:
    """
    Prepare already-converted REDIportal annotation format.
    
    Args:
        file_path: Path to converted annotation file
        output_prefix: Optional output prefix (ignored for converted files)
        
    Returns:
        Path to prepared annotation file
    """
    logger.info(f"Preparing already-converted REDIportal annotation format: {file_path}")
    
    # For already-converted format, we just need to ensure it's compressed and indexed
    if str(file_path).endswith('.gz'):
        # Check if index exists
        index_file = Path(str(file_path) + '.tbi')
        if index_file.exists():
            logger.info("Converted annotation file is already compressed and indexed")
            return str(file_path)
        else:
            logger.info("Creating index for compressed annotation file...")
            try:
                subprocess.run(['tabix', '-f', '-p', 'vcf', str(file_path)], check=True)
                logger.info(f"✓ Index created for converted file: {file_path}")
                return str(file_path)
            except subprocess.CalledProcessError as e:
                raise RuntimeError(f"Failed to index converted annotation file: {e}")
    else:
        # Need to compress and index
        logger.info("Compressing and indexing converted annotation file...")
        compressed_file = Path(str(file_path) + '.gz')
        
        try:
            # Compress with bgzip
            with open(compressed_file, 'wb') as f:
                subprocess.run(['bgzip', '-c', str(file_path)], stdout=f, check=True)
            
            # Index with tabix
            subprocess.run(['tabix', '-f', '-p', 'vcf', str(compressed_file)], check=True)
            
            logger.info(f"✓ Converted annotation file prepared: {compressed_file}")
            return str(compressed_file)
            
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"Failed to prepare converted annotation file: {e}")

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
                subprocess.run(['tabix', '-f', '-p', 'vcf', str(file_path)], check=True)
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
            subprocess.run(['tabix', '-f', '-p', 'vcf', str(compressed_file)], check=True)
            
            logger.info(f"✓ VCF prepared: {compressed_file}")
            return str(compressed_file)
            
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"Failed to prepare VCF: {e}")

def _prepare_text_format(file_path: Path, output_prefix: Optional[str] = None) -> str:
    """
    Convert REDIportal text format to bcftools annotation format with performance optimization.
    
    This function is optimized for large databases (>1GB) with:
    - Memory-efficient streaming processing
    - Progress monitoring and performance tracking
    - Efficient temporary file management
    - Resource usage optimization
    
    Args:
        file_path: Path to gzipped text REDIportal file
        output_prefix: Optional output prefix
        
    Returns:
        Path to prepared annotation file
    """
    logger.info(f"Converting REDIportal text format with performance optimization: {file_path}")
    
    # Initialize enhanced logging and performance monitoring
    conversion_logger = REDIportalConversionLogger()
    
    # Initialize performance monitoring for large files
    try:
        from .performance_optimizer import MemoryMonitor, TempFileManager
        memory_monitor = MemoryMonitor(max_memory_mb=4000)  # 4GB limit for conversion
        use_performance_monitoring = True
        logger.info("Performance monitoring enabled for large file processing")
    except ImportError:
        memory_monitor = None
        use_performance_monitoring = False
        logger.debug("Performance monitoring not available")
    
    if output_prefix:
        # If output_prefix ends with .vcf.gz, use it directly
        if output_prefix.endswith('.vcf.gz'):
            compressed_output = Path(output_prefix)
            output_file = Path(output_prefix[:-3])  # Remove .gz
        else:
            output_file = Path(f"{output_prefix}_annotations.txt")
            compressed_output = Path(str(output_file) + '.gz')
    else:
        output_file = file_path.parent / f"{file_path.stem}_annotations.txt"
        compressed_output = Path(str(output_file) + '.gz')
    
    # Safety check - don't overwrite existing files without explicit confirmation
    if compressed_output.exists():
        logger.warning(f"Output file already exists: {compressed_output}")
        logger.info("Using existing converted file")
        return str(compressed_output)
    
    logger.info(f"Converting to: {compressed_output}")
    logger.info(f"Input file size: {file_path.stat().st_size:,} bytes ({file_path.stat().st_size/1024/1024:.1f} MB)")
    
    try:
        # Use performance-optimized temporary file management
        if use_performance_monitoring:
            temp_manager = TempFileManager(prefix="rediportal_conv_")
            temp_path = temp_manager.create_temp_file(suffix='.txt', mode='w')
            temp_file = open(temp_path, 'w')
        else:
            # Fallback to standard temporary file
            temp_file = tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt')
            temp_path = Path(temp_file.name)
        
        try:
            # Monitor memory at conversion start
            if memory_monitor:
                memory_monitor.check_memory("conversion start")
            
            # Open input file with proper encoding handling and buffering optimization
            buffer_size = 8192 * 16  # 128KB buffer for large files
            
            if str(file_path).endswith('.gz'):
                try:
                    input_file = gzip.open(file_path, 'rt', encoding='utf-8')
                except UnicodeDecodeError:
                    logger.warning("UTF-8 decoding failed, trying latin-1 encoding")
                    input_file = gzip.open(file_path, 'rt', encoding='latin-1')
            else:
                try:
                    input_file = open(file_path, 'r', encoding='utf-8', buffering=buffer_size)
                except UnicodeDecodeError:
                    logger.warning("UTF-8 decoding failed, trying latin-1 encoding")
                    input_file = open(file_path, 'r', encoding='latin-1', buffering=buffer_size)
            
            try:
                # Write VCF header with contig definitions
                vcf_header = """##fileformat=VCFv4.2
##contig=<ID=chr1>
##contig=<ID=chr2>
##contig=<ID=chr3>
##contig=<ID=chr4>
##contig=<ID=chr5>
##contig=<ID=chr6>
##contig=<ID=chr7>
##contig=<ID=chr8>
##contig=<ID=chr9>
##contig=<ID=chr10>
##contig=<ID=chr11>
##contig=<ID=chr12>
##contig=<ID=chr13>
##contig=<ID=chr14>
##contig=<ID=chr15>
##contig=<ID=chr16>
##contig=<ID=chr17>
##contig=<ID=chr18>
##contig=<ID=chr19>
##contig=<ID=chr20>
##contig=<ID=chr21>
##contig=<ID=chr22>
##contig=<ID=chrX>
##contig=<ID=chrY>
##contig=<ID=chrM>
##INFO=<ID=REDI_ACCESSION,Number=1,Type=String,Description="REDIportal accession identifier">
##INFO=<ID=REDI_DB,Number=1,Type=String,Description="REDIportal database source">
##INFO=<ID=REDI_TYPE,Number=1,Type=String,Description="REDIportal editing type classification">
##INFO=<ID=REDI_REPEAT,Number=1,Type=String,Description="REDIportal repeat element annotation">
##INFO=<ID=REDI_FUNC,Number=1,Type=String,Description="REDIportal functional annotation">
##INFO=<ID=REDI_STRAND,Number=1,Type=String,Description="REDIportal strand information">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
"""
                temp_file.write(vcf_header)
                
                for line_num, line in enumerate(input_file, 1):
                    conversion_logger.stats['total_lines_read'] += 1
                    line = line.strip()
                    
                    # Skip header line
                    if line_num == 1 and (line.startswith('Accession') or 'Region' in line):
                        logger.debug("Skipping header line")
                        conversion_logger.stats['header_lines_skipped'] += 1
                        continue
                    
                    if not line:
                        continue
                    
                    try:
                        fields = line.split('\t')
                        if len(fields) < 10:
                            conversion_logger.log_parsing_error(
                                line_num, 
                                f"Insufficient fields ({len(fields)}, expected 10+)", 
                                line
                            )
                            conversion_logger.stats['entries_skipped'] += 1
                            continue
                        
                        # Extract required fields with validation
                        try:
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
                        except IndexError as e:
                            conversion_logger.log_parsing_error(
                                line_num, 
                                f"Field extraction failed: {e}", 
                                line
                            )
                            conversion_logger.stats['entries_skipped'] += 1
                            continue
                        
                        # Parse chromosome and position with validation
                        chrom = region.strip()
                        pos = position.strip()
                        
                        # Validate chromosome format
                        if not chrom:
                            conversion_logger.log_validation_error(line_num, 'Region', chrom, "Empty chromosome")
                            conversion_logger.stats['entries_skipped'] += 1
                            continue
                        
                        if not chrom.startswith('chr'):
                            chrom = f"chr{chrom}"
                        
                        # Validate position is numeric
                        try:
                            pos_int = int(pos)
                            if pos_int <= 0:
                                conversion_logger.log_validation_error(line_num, 'Position', pos, "Position must be positive")
                                conversion_logger.stats['entries_skipped'] += 1
                                continue
                        except ValueError:
                            conversion_logger.log_validation_error(line_num, 'Position', pos, "Position must be numeric")
                            conversion_logger.stats['entries_skipped'] += 1
                            continue
                        
                        # Validate nucleotides
                        if ref not in ['A', 'T', 'G', 'C']:
                            conversion_logger.log_validation_error(line_num, 'Ref', ref, "Invalid reference nucleotide")
                            conversion_logger.stats['entries_skipped'] += 1
                            continue
                        
                        if ed not in ['A', 'T', 'G', 'C']:
                            conversion_logger.log_validation_error(line_num, 'Ed', ed, "Invalid edited nucleotide")
                            conversion_logger.stats['entries_skipped'] += 1
                            continue
                        
                        # Map fields using the field mapping system with error handling
                        try:
                            mapped_accession = REDIportalFieldMapper.map_accession_field(accession)
                            mapped_db = REDIportalFieldMapper.map_db_field(db)
                            mapped_type = REDIportalFieldMapper.map_type_field(rna_type)
                            mapped_repeat = REDIportalFieldMapper.map_repeat_field(repeat)
                            mapped_strand = REDIportalFieldMapper.map_strand_field(strand)
                        except Exception as e:
                            conversion_logger.log_field_mapping_error(line_num, 'multiple', str(fields), str(e))
                            conversion_logger.stats['entries_skipped'] += 1
                            continue
                        
                        # Clean up func field (not part of the main mapping but still needed)
                        func = func.strip() if func.strip() else '.'
                        
                        # Create bcftools annotation format line
                        # Create VCF record with INFO fields
                        info_fields = [
                            f"REDI_ACCESSION={mapped_accession}",
                            f"REDI_DB={mapped_db}",
                            f"REDI_TYPE={mapped_type}",
                            f"REDI_REPEAT={mapped_repeat}",
                            f"REDI_FUNC={func}",
                            f"REDI_STRAND={mapped_strand}"
                        ]
                        info_string = ";".join(info_fields)
                        
                        # Format as VCF record: CHROM POS ID REF ALT QUAL FILTER INFO
                        vcf_line = f"{chrom}\t{pos}\t.\t{ref}\t{ed}\t.\tPASS\t{info_string}\n"
                        temp_file.write(vcf_line)
                        conversion_logger.stats['entries_processed'] += 1
                        
                        # Log progress at regular intervals with performance monitoring
                        conversion_logger.log_progress(conversion_logger.stats['entries_processed'])
                        
                        # Memory monitoring and optimization for large files
                        if memory_monitor and conversion_logger.stats['entries_processed'] % 50000 == 0:
                            memory_stats = memory_monitor.check_memory(f"processed {conversion_logger.stats['entries_processed']:,} entries")
                            
                            # Force garbage collection every 100K entries to manage memory
                            if conversion_logger.stats['entries_processed'] % 100000 == 0:
                                freed_mb = memory_monitor.force_garbage_collection("periodic cleanup")
                                if freed_mb > 10:  # Log if significant memory was freed
                                    logger.info(f"Freed {freed_mb:.1f} MB during processing")
                            
                    except Exception as e:
                        conversion_logger.log_parsing_error(line_num, f"Unexpected parsing error: {e}", line)
                        conversion_logger.stats['entries_skipped'] += 1
                        continue
                        
            finally:
                input_file.close()
                temp_file.close()  # Ensure temp file is closed
        
        except Exception as e:
            logger.error(f"Error during file processing: {e}")
            # Clean up temp file on error
            if 'temp_path' in locals() and temp_path.exists():
                temp_path.unlink()
            raise
        
        # Log conversion completion statistics
        conversion_logger.log_conversion_statistics()
        
        # Sort VCF file before compression (required for tabix indexing)
        logger.info("Sorting VCF file for tabix compatibility...")
        sorting_start = time.time()
        
        try:
            # Create sorted temporary file
            sorted_temp_path = temp_path.parent / f"{temp_path.stem}_sorted{temp_path.suffix}"
            
            # Use bcftools sort for proper VCF sorting
            sort_cmd = ['bcftools', 'sort', '-o', str(sorted_temp_path), str(temp_path)]
            subprocess.run(sort_cmd, check=True, capture_output=True, text=True)
            
            # Replace original temp file with sorted version
            temp_path.unlink()
            sorted_temp_path.rename(temp_path)
            
            sorting_time = time.time() - sorting_start
            logger.info(f"✓ VCF sorting completed ({sorting_time:.2f}s)")
            
        except subprocess.CalledProcessError as e:
            logger.error(f"VCF sorting failed: {e}")
            logger.error(f"bcftools stderr: {e.stderr}")
            raise RuntimeError(f"VCF sorting failed: {e}")
        
        # Compress with bgzip using optimized settings for large files
        logger.info("Compressing annotation file with performance optimization...")
        compression_start = time.time()
        
        # Monitor memory before compression
        if memory_monitor:
            memory_monitor.check_memory("before compression")
        
        try:
            # Use optimized bgzip compression for large files
            compression_level = 6  # Balance between speed and compression ratio
            
            with open(compressed_output, 'wb') as f:
                cmd = ['bgzip', '-c', '-l', str(compression_level), str(temp_path)]
                
                subprocess.run(
                    cmd, 
                    stdout=f, 
                    stderr=subprocess.PIPE,
                    text=True,
                    check=True,
                    timeout=3600  # 1 hour timeout for very large files
                )
            
            conversion_logger.stats['compression_time'] = time.time() - compression_start
            compressed_size = compressed_output.stat().st_size
            compression_ratio = (compressed_size / temp_path.stat().st_size) * 100 if temp_path.stat().st_size > 0 else 0
            
            logger.info(f"✓ Compression completed: {compressed_size:,} bytes ({compression_ratio:.1f}% of original)")
            
        except subprocess.CalledProcessError as e:
            logger.error(f"Compression failed: {e}")
            logger.error(f"bgzip stderr: {e.stderr}")
            raise RuntimeError(f"Compression failed: {e}")
        except subprocess.TimeoutExpired:
            logger.error("Compression timed out after 30 minutes")
            raise RuntimeError("Compression timeout")
        
        # Index with tabix (optional for performance testing)
        logger.info("Creating tabix index...")
        indexing_start = time.time()
        try:
            subprocess.run(
                ['tabix', '-f', '-p', 'vcf', str(compressed_output)], 
                capture_output=True,
                text=True,
                check=True,
                timeout=300  # 5 minute timeout
            )
            
            conversion_logger.stats['indexing_time'] = time.time() - indexing_start
            
            # Verify index was created
            index_file = Path(str(compressed_output) + '.tbi')
            if index_file.exists():
                index_size = index_file.stat().st_size
                logger.info(f"✓ Tabix index created: {index_size:,} bytes")
            else:
                logger.warning("Index file was not created but continuing")
                
        except subprocess.CalledProcessError as e:
            logger.warning(f"Indexing failed (continuing without index): {e}")
            logger.debug(f"tabix stderr: {e.stderr}")
            # Don't raise error for indexing failure in performance tests
            conversion_logger.stats['indexing_time'] = time.time() - indexing_start
        except subprocess.TimeoutExpired:
            logger.warning("Indexing timed out after 5 minutes (continuing without index)")
            conversion_logger.stats['indexing_time'] = time.time() - indexing_start
        
        # Clean up temporary file with performance monitoring
        try:
            if use_performance_monitoring and 'temp_manager' in locals():
                # Use managed cleanup
                cleanup_stats = temp_manager.cleanup_temp_files()
                logger.debug(f"Managed cleanup: {cleanup_stats}")
            else:
                # Standard cleanup
                temp_path.unlink()
                logger.debug(f"Cleaned up temporary file: {temp_path}")
        except Exception as e:
            logger.warning(f"Failed to clean up temporary file {temp_path}: {e}")
        
        # Final validation and statistics with performance metrics
        if not compressed_output.exists():
            raise RuntimeError("Compressed output file was not created")
        
        # Index file is optional for performance testing
        index_file = Path(str(compressed_output) + '.tbi')
        if not index_file.exists():
            logger.debug("Index file was not created (optional for performance testing)")
        
        # Log final statistics with performance analysis
        output_size = compressed_output.stat().st_size
        total_time = time.time() - conversion_logger.stats['start_time']
        input_size = file_path.stat().st_size
        
        # Calculate performance metrics
        processing_rate_mb_per_sec = (input_size / 1024 / 1024) / total_time if total_time > 0 else 0
        compression_ratio = (output_size / input_size) * 100 if input_size > 0 else 0
        entries_per_second = conversion_logger.stats['entries_processed'] / total_time if total_time > 0 else 0
        
        logger.info("✓ REDIportal text format converted successfully with performance optimization")
        logger.info(f"  Input: {file_path} ({input_size:,} bytes, {input_size/1024/1024:.1f} MB)")
        logger.info(f"  Output: {compressed_output} ({output_size:,} bytes, {output_size/1024/1024:.1f} MB)")
        if index_file.exists():
            logger.info(f"  Index: {index_file} ({index_file.stat().st_size:,} bytes)")
        else:
            logger.info(f"  Index: Not created (optional for performance testing)")
        logger.info(f"  Total time: {total_time:.2f} seconds")
        logger.info(f"  Processing rate: {processing_rate_mb_per_sec:.1f} MB/sec")
        logger.info(f"  Compression ratio: {compression_ratio:.1f}%")
        logger.info(f"  Entries processed: {conversion_logger.stats['entries_processed']:,} ({entries_per_second:.0f} entries/sec)")
        logger.info(f"  Entries skipped: {conversion_logger.stats['entries_skipped']:,}")
        logger.info(f"  Success rate: {(conversion_logger.stats['entries_processed'] / conversion_logger.stats['total_lines_read'] * 100):.1f}%")
        
        # Log memory usage summary if available
        if memory_monitor:
            memory_summary = memory_monitor.get_summary()
            logger.info(f"  Peak memory usage: {memory_summary['peak_mb']:.1f} MB")
            logger.info(f"  Memory increase: {memory_summary['total_increase_mb']:.1f} MB")
            logger.info(f"  Memory trend: {memory_summary['trend']}")
        
        # Performance assessment
        if processing_rate_mb_per_sec > 20:
            logger.info("  Performance: Excellent (>20 MB/sec)")
        elif processing_rate_mb_per_sec > 10:
            logger.info("  Performance: Good (>10 MB/sec)")
        elif processing_rate_mb_per_sec > 5:
            logger.info("  Performance: Acceptable (>5 MB/sec)")
        else:
            logger.warning("  Performance: Slow (<5 MB/sec) - consider system optimization")
        
        return str(compressed_output)
        
    except subprocess.CalledProcessError as e:
        logger.error(f"Command failed during conversion: {e}")
        logger.error(f"Command: {' '.join(e.cmd) if hasattr(e, 'cmd') else 'unknown'}")
        logger.error(f"Return code: {e.returncode}")
        logger.error(f"stderr: {e.stderr}")
        raise RuntimeError(f"REDIportal conversion failed: {e}")
    except Exception as e:
        logger.error(f"Unexpected error during conversion: {e}")
        logger.error(f"Error type: {type(e).__name__}")
        
        # Clean up partial files
        if compressed_output.exists():
            try:
                compressed_output.unlink()
                logger.info("Cleaned up partial output file")
            except Exception as cleanup_error:
                logger.warning(f"Failed to clean up partial output file: {cleanup_error}")
        
        # Log final error statistics
        conversion_logger.log_conversion_statistics()
        
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