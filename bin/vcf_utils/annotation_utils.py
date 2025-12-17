#!/usr/bin/env python3
"""
Annotation Utilities Module

This module provides shared utilities and helper functions for COSMIC and gnomAD
annotation workflows. It includes common functions for file validation, VCF
processing, error handling, and statistics generation.

FEATURES:
- VCF file validation and format checking
- Chromosome name normalization and mapping
- Annotation field extraction and parsing
- Error handling and logging utilities
- Statistics aggregation and reporting
- File compression and indexing helpers

Requirements Satisfied: 2.6

Author: COSMIC/gnomAD Enhancement Pipeline
Date: 2025-12-17
"""

import logging
import subprocess
import os
import time
from pathlib import Path
from typing import Dict, Any, Optional, List, Tuple, Union

logger = logging.getLogger(__name__)


def validate_vcf_file(vcf_path: Union[str, Path], check_index: bool = True) -> Dict[str, Any]:
    """
    Validate VCF file format, accessibility, and indexing.
    
    Args:
        vcf_path: Path to VCF file to validate
        check_index: Whether to check for tabix index (for compressed files)
        
    Returns:
        Dictionary with validation results and file information
    """
    vcf_path = Path(vcf_path)
    
    validation_result = {
        'valid': False,
        'exists': False,
        'readable': False,
        'size_bytes': 0,
        'compressed': False,
        'indexed': False,
        'vcf_format': False,
        'header_lines': 0,
        'sample_count': 0,
        'errors': [],
        'warnings': []
    }
    
    errors = validation_result['errors']
    warnings = validation_result['warnings']
    
    try:
        # Check file existence
        if not vcf_path.exists():
            errors.append(f"File does not exist: {vcf_path}")
            return validation_result
        
        validation_result['exists'] = True
        
        # Check file accessibility
        if not os.access(vcf_path, os.R_OK):
            errors.append(f"File is not readable: {vcf_path}")
            return validation_result
        
        validation_result['readable'] = True
        
        # Get file size
        validation_result['size_bytes'] = vcf_path.stat().st_size
        
        if validation_result['size_bytes'] == 0:
            errors.append("File is empty")
            return validation_result
        
        # Check if compressed
        validation_result['compressed'] = str(vcf_path).endswith(('.gz', '.bgz'))
        
        # Check for index if compressed
        if validation_result['compressed'] and check_index:
            index_extensions = ['.tbi', '.csi']
            for ext in index_extensions:
                index_path = Path(str(vcf_path) + ext)
                if index_path.exists():
                    validation_result['indexed'] = True
                    break
            
            if not validation_result['indexed']:
                warnings.append("Compressed VCF file is not indexed")
        
        # Validate VCF format using bcftools
        try:
            cmd = ['bcftools', 'view', '-h', str(vcf_path)]
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True,
                timeout=60
            )
            
            header_lines = result.stdout.strip().split('\n')
            validation_result['header_lines'] = len(header_lines)
            
            # Check VCF format header
            if header_lines and header_lines[0].startswith('##fileformat=VCF'):
                validation_result['vcf_format'] = True
            else:
                errors.append("Invalid VCF format header")
            
            # Count samples
            for line in header_lines:
                if line.startswith('#CHROM'):
                    columns = line.split('\t')
                    if len(columns) > 9:
                        validation_result['sample_count'] = len(columns) - 9
                    break
            
        except subprocess.CalledProcessError as e:
            errors.append(f"bcftools validation failed: {e.stderr}")
        except subprocess.TimeoutExpired:
            errors.append("bcftools validation timed out")
        except Exception as e:
            errors.append(f"Unexpected error during validation: {e}")
        
        # Overall validation status
        validation_result['valid'] = (
            validation_result['exists'] and
            validation_result['readable'] and
            validation_result['vcf_format'] and
            validation_result['size_bytes'] > 0 and
            len(errors) == 0
        )
        
    except Exception as e:
        errors.append(f"Validation error: {e}")
    
    return validation_result


def normalize_chromosome_name(chrom: str) -> str:
    """
    Normalize chromosome name for consistent processing.
    
    Args:
        chrom: Chromosome name (e.g., 'chr1', 'Chr1', '1', 'chrX')
        
    Returns:
        Normalized chromosome name (e.g., '1', 'X', 'M')
    """
    if not chrom:
        return chrom
    
    # Remove common prefixes
    normalized = str(chrom).replace('chr', '').replace('Chr', '').replace('CHR', '')
    
    # Handle mitochondrial chromosome variants
    if normalized.upper() in ['MT', 'CHRM', 'M']:
        return 'M'
    
    # Handle sex chromosomes
    if normalized.upper() == 'X':
        return 'X'
    if normalized.upper() == 'Y':
        return 'Y'
    
    # Return as-is for numeric and other chromosomes
    return normalized


def add_chromosome_prefix(chrom: str, prefix: str = 'chr') -> str:
    """
    Add chromosome prefix if not already present.
    
    Args:
        chrom: Chromosome name
        prefix: Prefix to add (default: 'chr')
        
    Returns:
        Chromosome name with prefix
    """
    if not chrom:
        return chrom
    
    chrom_str = str(chrom)
    if chrom_str.startswith(prefix):
        return chrom_str
    
    return f"{prefix}{chrom_str}"


def extract_info_fields(vcf_header: str) -> Dict[str, Dict[str, str]]:
    """
    Extract INFO field definitions from VCF header.
    
    Args:
        vcf_header: VCF header string
        
    Returns:
        Dictionary mapping INFO field IDs to their definitions
    """
    info_fields = {}
    
    for line in vcf_header.split('\n'):
        if line.startswith('##INFO=<'):
            try:
                # Parse INFO field definition
                # Example: ##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
                
                # Extract content between < and >
                content = line[8:-1]  # Remove ##INFO=< and >
                
                # Parse key=value pairs
                field_def = {}
                current_key = None
                current_value = ""
                in_quotes = False
                
                i = 0
                while i < len(content):
                    char = content[i]
                    
                    if char == '"':
                        in_quotes = not in_quotes
                        if not in_quotes and current_key:
                            field_def[current_key] = current_value
                            current_key = None
                            current_value = ""
                    elif char == '=' and not in_quotes:
                        if current_value:
                            # This is a new key=value pair
                            if current_key:
                                field_def[current_key] = current_value
                            current_key = current_value
                            current_value = ""
                        else:
                            current_key = current_value
                            current_value = ""
                    elif char == ',' and not in_quotes:
                        if current_key and current_value:
                            field_def[current_key] = current_value
                            current_key = None
                            current_value = ""
                    else:
                        current_value += char
                    
                    i += 1
                
                # Handle last key=value pair
                if current_key and current_value:
                    field_def[current_key] = current_value
                
                # Store field definition
                if 'ID' in field_def:
                    info_fields[field_def['ID']] = field_def
                    
            except Exception as e:
                logger.debug(f"Failed to parse INFO field line: {line[:100]}... Error: {e}")
    
    return info_fields


def parse_info_value(value: str, field_type: str = 'String') -> Any:
    """
    Parse INFO field value according to its type.
    
    Args:
        value: Raw INFO field value string
        field_type: VCF field type (Integer, Float, String, Flag)
        
    Returns:
        Parsed value in appropriate Python type
    """
    if not value or value == '.':
        return None
    
    try:
        if field_type == 'Integer':
            if ',' in value:
                return [int(x) for x in value.split(',') if x.strip()]
            else:
                return int(value)
        elif field_type == 'Float':
            if ',' in value:
                return [float(x) for x in value.split(',') if x.strip()]
            else:
                return float(value)
        elif field_type == 'Flag':
            return True  # Presence of flag indicates True
        else:  # String or unknown type
            return value
    except (ValueError, TypeError):
        # Return as string if parsing fails
        return value


def create_temp_directory(base_dir: Optional[Path] = None, prefix: str = "annotation_temp") -> Path:
    """
    Create a temporary directory for annotation processing.
    
    Args:
        base_dir: Base directory for temporary files (default: current working directory)
        prefix: Prefix for temporary directory name
        
    Returns:
        Path to created temporary directory
    """
    if base_dir is None:
        base_dir = Path.cwd()
    
    timestamp = int(time.time())
    temp_dir = base_dir / f"{prefix}_{timestamp}"
    
    # Ensure unique directory name
    counter = 1
    while temp_dir.exists():
        temp_dir = base_dir / f"{prefix}_{timestamp}_{counter}"
        counter += 1
    
    temp_dir.mkdir(parents=True, exist_ok=True)
    logger.debug(f"Created temporary directory: {temp_dir}")
    
    return temp_dir


def cleanup_temp_directory(temp_dir: Path, force: bool = False) -> bool:
    """
    Clean up temporary directory and all its contents.
    
    Args:
        temp_dir: Path to temporary directory to clean up
        force: If True, ignore errors and force cleanup
        
    Returns:
        True if cleanup was successful, False otherwise
    """
    if not temp_dir.exists():
        logger.debug(f"Temporary directory does not exist: {temp_dir}")
        return True
    
    try:
        import shutil
        shutil.rmtree(temp_dir)
        logger.debug(f"Successfully cleaned up temporary directory: {temp_dir}")
        return True
    except Exception as e:
        if force:
            logger.warning(f"Failed to clean up temporary directory (forced): {temp_dir} - {e}")
            return False
        else:
            logger.error(f"Failed to clean up temporary directory: {temp_dir} - {e}")
            raise


def compress_and_index_vcf(vcf_path: Path, output_path: Optional[Path] = None, 
                          remove_original: bool = False) -> Path:
    """
    Compress VCF file with bgzip and create tabix index.
    
    Args:
        vcf_path: Path to uncompressed VCF file
        output_path: Path for compressed output (default: vcf_path + .gz)
        remove_original: Whether to remove original uncompressed file
        
    Returns:
        Path to compressed and indexed VCF file
    """
    if output_path is None:
        output_path = Path(str(vcf_path) + '.gz')
    
    try:
        # Compress with bgzip
        logger.debug(f"Compressing VCF file: {vcf_path} -> {output_path}")
        
        cmd = ['bgzip', '-c', str(vcf_path)]
        with open(output_path, 'wb') as f:
            subprocess.run(
                cmd,
                stdout=f,
                stderr=subprocess.PIPE,
                text=True,
                check=True,
                timeout=1800  # 30 minute timeout
            )
        
        # Verify compression
        if not output_path.exists() or output_path.stat().st_size == 0:
            raise RuntimeError("Compression failed - output file is missing or empty")
        
        # Create tabix index
        logger.debug(f"Creating tabix index for: {output_path}")
        
        subprocess.run(
            ['tabix', '-p', 'vcf', str(output_path)],
            capture_output=True,
            text=True,
            check=True,
            timeout=300  # 5 minute timeout
        )
        
        # Verify index creation
        index_path = Path(str(output_path) + '.tbi')
        if not index_path.exists():
            raise RuntimeError("Tabix indexing failed - index file not created")
        
        # Remove original file if requested
        if remove_original and vcf_path.exists():
            vcf_path.unlink()
            logger.debug(f"Removed original uncompressed file: {vcf_path}")
        
        logger.debug(f"Successfully compressed and indexed VCF: {output_path}")
        return output_path
        
    except subprocess.CalledProcessError as e:
        error_msg = f"Failed to compress and index VCF: {e}"
        logger.error(error_msg)
        logger.error(f"Command stderr: {e.stderr}")
        raise RuntimeError(error_msg)
    except subprocess.TimeoutExpired:
        error_msg = "VCF compression/indexing timed out"
        logger.error(error_msg)
        raise RuntimeError(error_msg)
    except Exception as e:
        logger.error(f"Unexpected error during VCF compression/indexing: {e}")
        raise


def get_vcf_statistics(vcf_path: Path) -> Dict[str, Any]:
    """
    Get basic statistics about a VCF file using bcftools stats.
    
    Args:
        vcf_path: Path to VCF file
        
    Returns:
        Dictionary with VCF statistics
    """
    stats = {
        'total_records': 0,
        'snps': 0,
        'indels': 0,
        'multiallelic': 0,
        'samples': 0,
        'chromosomes': set(),
        'file_size_bytes': 0
    }
    
    try:
        # Get file size
        stats['file_size_bytes'] = vcf_path.stat().st_size
        
        # Run bcftools stats
        cmd = ['bcftools', 'stats', str(vcf_path)]
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True,
            timeout=300  # 5 minute timeout
        )
        
        # Parse bcftools stats output
        for line in result.stdout.split('\n'):
            if line.startswith('SN\t'):
                # Summary numbers
                parts = line.split('\t')
                if len(parts) >= 4:
                    key = parts[2].strip()
                    value_str = parts[3].strip()
                    
                    try:
                        value = int(value_str)
                        if 'number of records' in key.lower():
                            stats['total_records'] = value
                        elif 'number of snps' in key.lower():
                            stats['snps'] = value
                        elif 'number of indels' in key.lower():
                            stats['indels'] = value
                        elif 'number of multiallelic sites' in key.lower():
                            stats['multiallelic'] = value
                        elif 'number of samples' in key.lower():
                            stats['samples'] = value
                    except ValueError:
                        continue
        
        # Get chromosome information
        cmd = ['bcftools', 'view', '-H', str(vcf_path)]
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True,
            timeout=300
        )
        
        chromosomes = set()
        for line in result.stdout.split('\n'):
            if line.strip():
                chrom = line.split('\t')[0]
                chromosomes.add(normalize_chromosome_name(chrom))
        
        stats['chromosomes'] = chromosomes
        
    except subprocess.CalledProcessError as e:
        logger.warning(f"Failed to get VCF statistics: {e}")
    except subprocess.TimeoutExpired:
        logger.warning("VCF statistics collection timed out")
    except Exception as e:
        logger.warning(f"Unexpected error getting VCF statistics: {e}")
    
    return stats


def format_file_size(size_bytes: int) -> str:
    """
    Format file size in human-readable format.
    
    Args:
        size_bytes: File size in bytes
        
    Returns:
        Formatted file size string
    """
    if size_bytes < 1024:
        return f"{size_bytes} B"
    elif size_bytes < 1024 ** 2:
        return f"{size_bytes / 1024:.1f} KB"
    elif size_bytes < 1024 ** 3:
        return f"{size_bytes / (1024 ** 2):.1f} MB"
    else:
        return f"{size_bytes / (1024 ** 3):.1f} GB"


def format_duration(seconds: float) -> str:
    """
    Format duration in human-readable format.
    
    Args:
        seconds: Duration in seconds
        
    Returns:
        Formatted duration string
    """
    if seconds < 60:
        return f"{seconds:.1f}s"
    elif seconds < 3600:
        minutes = seconds / 60
        return f"{minutes:.1f}m"
    else:
        hours = seconds / 3600
        return f"{hours:.1f}h"


def log_annotation_summary(input_vcf: Path, output_vcf: Path, 
                          processing_time: float, annotations_added: int = 0,
                          annotation_type: str = "annotation") -> None:
    """
    Log comprehensive annotation summary.
    
    Args:
        input_vcf: Path to input VCF file
        output_vcf: Path to output VCF file
        processing_time: Total processing time in seconds
        annotations_added: Number of annotations added
        annotation_type: Type of annotation performed
    """
    logger.info(f"=== {annotation_type.title()} Summary ===")
    
    # File information
    try:
        input_size = input_vcf.stat().st_size if input_vcf.exists() else 0
        output_size = output_vcf.stat().st_size if output_vcf.exists() else 0
        
        logger.info(f"Input file: {input_vcf} ({format_file_size(input_size)})")
        logger.info(f"Output file: {output_vcf} ({format_file_size(output_size)})")
        
        if output_size > input_size:
            size_increase = output_size - input_size
            logger.info(f"Size increase: {format_file_size(size_increase)}")
        
    except Exception as e:
        logger.debug(f"Could not get file size information: {e}")
    
    # Processing information
    logger.info(f"Processing time: {format_duration(processing_time)}")
    
    if annotations_added > 0:
        logger.info(f"Annotations added: {annotations_added:,}")
        
        if processing_time > 0:
            rate = annotations_added / processing_time
            logger.info(f"Annotation rate: {rate:.1f} annotations/second")
    
    logger.info(f"=== {annotation_type.title()} Complete ===")


def check_tool_availability(tools: List[str]) -> Tuple[List[str], List[str]]:
    """
    Check availability of required tools in system PATH.
    
    Args:
        tools: List of tool names to check
        
    Returns:
        Tuple of (available_tools, missing_tools)
    """
    available_tools = []
    missing_tools = []
    
    for tool in tools:
        try:
            result = subprocess.run(
                ['which', tool],
                capture_output=True,
                text=True,
                timeout=10
            )
            
            if result.returncode == 0:
                available_tools.append(tool)
                logger.debug(f"Found tool: {tool} at {result.stdout.strip()}")
            else:
                missing_tools.append(tool)
                logger.debug(f"Tool not found: {tool}")
                
        except Exception as e:
            missing_tools.append(tool)
            logger.debug(f"Error checking tool {tool}: {e}")
    
    return available_tools, missing_tools


def main():
    """Main entry point for standalone usage and testing."""
    import argparse
    import sys
    
    parser = argparse.ArgumentParser(
        description="Annotation utilities for VCF processing",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Annotation Utilities:
  This module provides shared utilities and helper functions for COSMIC and gnomAD
  annotation workflows. It includes functions for file validation, VCF processing,
  error handling, and statistics generation.

Examples:
  # Validate VCF file
  python annotation_utils.py --validate input.vcf.gz
  
  # Get VCF statistics
  python annotation_utils.py --stats input.vcf.gz
  
  # Compress and index VCF
  python annotation_utils.py --compress input.vcf
        """
    )
    
    parser.add_argument(
        '--validate',
        metavar='VCF_FILE',
        help='Validate VCF file format and accessibility'
    )
    
    parser.add_argument(
        '--stats',
        metavar='VCF_FILE',
        help='Get VCF file statistics'
    )
    
    parser.add_argument(
        '--compress',
        metavar='VCF_FILE',
        help='Compress and index VCF file'
    )
    
    parser.add_argument(
        '--check-tools',
        action='store_true',
        help='Check availability of required tools'
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
    
    try:
        if args.validate:
            logger.info(f"Validating VCF file: {args.validate}")
            result = validate_vcf_file(args.validate)
            
            logger.info(f"Validation result: {'PASS' if result['valid'] else 'FAIL'}")
            logger.info(f"File size: {format_file_size(result['size_bytes'])}")
            logger.info(f"Compressed: {result['compressed']}")
            logger.info(f"Indexed: {result['indexed']}")
            logger.info(f"Header lines: {result['header_lines']}")
            logger.info(f"Samples: {result['sample_count']}")
            
            if result['warnings']:
                logger.warning("Warnings:")
                for warning in result['warnings']:
                    logger.warning(f"  {warning}")
            
            if result['errors']:
                logger.error("Errors:")
                for error in result['errors']:
                    logger.error(f"  {error}")
        
        elif args.stats:
            logger.info(f"Getting statistics for VCF file: {args.stats}")
            stats = get_vcf_statistics(Path(args.stats))
            
            logger.info(f"Total records: {stats['total_records']:,}")
            logger.info(f"SNPs: {stats['snps']:,}")
            logger.info(f"Indels: {stats['indels']:,}")
            logger.info(f"Multiallelic sites: {stats['multiallelic']:,}")
            logger.info(f"Samples: {stats['samples']}")
            logger.info(f"Chromosomes: {', '.join(sorted(stats['chromosomes']))}")
            logger.info(f"File size: {format_file_size(stats['file_size_bytes'])}")
        
        elif args.compress:
            logger.info(f"Compressing and indexing VCF file: {args.compress}")
            input_path = Path(args.compress)
            output_path = compress_and_index_vcf(input_path)
            logger.info(f"Compressed file created: {output_path}")
        
        elif args.check_tools:
            logger.info("Checking tool availability...")
            required_tools = ['bcftools', 'tabix', 'bgzip']
            available, missing = check_tool_availability(required_tools)
            
            logger.info(f"Available tools: {', '.join(available)}")
            if missing:
                logger.error(f"Missing tools: {', '.join(missing)}")
                sys.exit(1)
            else:
                logger.info("All required tools are available")
        
        else:
            parser.print_help()
            sys.exit(1)
        
        logger.info("✓ Operation completed successfully!")
        sys.exit(0)
        
    except KeyboardInterrupt:
        logger.error("Operation interrupted by user")
        sys.exit(130)
    except Exception as e:
        logger.error(f"Operation failed: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()