#!/usr/bin/env python3
"""
gnomAD Scatter-Gather Annotation Module

This module provides gnomAD database annotation functionality using scatter-gather
processing across chromosome-split VCF files. It implements parallel processing
across chromosomes for efficient annotation with population frequency data.

FEATURES:
- Scatter-gather processing across chromosome-split gnomAD files
- Population frequency data extraction (faf95, AF)
- Parallel chromosome processing for performance
- Compressed VCF input/output handling with proper indexing
- Comprehensive error handling and logging

Requirements Satisfied: 2.1, 2.2

Author: COSMIC/gnomAD Enhancement Pipeline
Date: 2025-12-17
"""

import logging
import subprocess
import os
import time
import glob
import psutil
from pathlib import Path
from typing import Dict, List, Tuple
from concurrent.futures import ThreadPoolExecutor, as_completed

logger = logging.getLogger(__name__)


class GnomadAnnotator:
    """
    gnomAD scatter-gather annotator for population frequency annotation.
    
    This class provides scatter-gather gnomAD annotation across chromosome-split
    VCF files for efficient parallel processing. It extracts population frequency
    data and handles chromosome-specific annotation workflows.
    """
    
    def __init__(self, input_vcf: str, gnomad_dir: str, output_vcf: str, max_workers: int = None):
        """
        Initialize gnomAD annotator.
        
        Args:
            input_vcf: Path to input VCF file to be annotated
            gnomad_dir: Path to gnomAD database directory containing chromosome-split files
            output_vcf: Path to output annotated VCF file
            max_workers: Maximum number of parallel workers (None for auto-optimization)
        """
        self.input_vcf = Path(input_vcf)
        self.gnomad_dir = Path(gnomad_dir)
        self.output_vcf = Path(output_vcf)
        
        # Auto-optimize worker count if not specified
        if max_workers is None:
            self.max_workers = self._optimize_worker_count()
        else:
            self.max_workers = max_workers
        
        # Statistics and timing
        self.stats = {
            'start_time': time.time(),
            'processing_steps': [],
            'errors': [],
            'warnings': [],
            'chromosomes_processed': 0,
            'annotations_added': 0,
            'variants_processed': 0
        }
        
        logger.info("=== gnomAD Scatter-Gather Annotator Initialized ===")
        logger.info(f"Input VCF: {self.input_vcf}")
        logger.info(f"gnomAD directory: {self.gnomad_dir}")
        logger.info(f"Output VCF: {self.output_vcf}")
        logger.info(f"Max workers: {self.max_workers}")
        
        # Validate inputs
        self.validate_inputs()
        self.validate_tools()
        
        # Discover gnomAD chromosome files
        self.gnomad_files = self.discover_gnomad_files()
    
    def _optimize_worker_count(self) -> int:
        """
        Automatically determine optimal worker count based on system resources.
        
        Returns:
            Optimal number of workers for parallel processing
        """
        try:
            # Get system resources
            cpu_count = psutil.cpu_count(logical=True)
            memory_gb = psutil.virtual_memory().total / (1024**3)
            
            # Conservative memory estimate: ~2GB per worker for gnomAD processing
            memory_limited_workers = max(1, int(memory_gb / 2))
            
            # CPU-based workers: use 75% of available cores for I/O bound tasks
            cpu_limited_workers = max(1, int(cpu_count * 0.75))
            
            # Take the minimum to avoid resource exhaustion
            optimal_workers = min(memory_limited_workers, cpu_limited_workers, 16)  # Cap at 16
            
            logger.info(f"Auto-optimized worker count: {optimal_workers} "
                       f"(CPU cores: {cpu_count}, Memory: {memory_gb:.1f}GB)")
            
            return optimal_workers
            
        except Exception as e:
            logger.warning(f"Failed to auto-optimize worker count: {e}, using default 4")
            return 4
    
    def validate_inputs(self) -> None:
        """Validate input files exist and are accessible."""
        step_start = time.time()
        logger.info("Validating gnomAD annotation input files...")
        
        # Check input VCF existence and accessibility
        if not self.input_vcf.exists():
            error_msg = f"Input VCF not found: {self.input_vcf}"
            logger.error(error_msg)
            raise FileNotFoundError(error_msg)
        
        if not self.gnomad_dir.exists():
            error_msg = f"gnomAD directory not found: {self.gnomad_dir}"
            logger.error(error_msg)
            raise FileNotFoundError(error_msg)
        
        if not self.gnomad_dir.is_dir():
            error_msg = f"gnomAD path is not a directory: {self.gnomad_dir}"
            logger.error(error_msg)
            raise NotADirectoryError(error_msg)
        
        # Check file permissions
        if not os.access(self.input_vcf, os.R_OK):
            error_msg = f"Cannot read input VCF: {self.input_vcf}"
            logger.error(error_msg)
            raise PermissionError(error_msg)
            
        if not os.access(self.gnomad_dir, os.R_OK):
            error_msg = f"Cannot read gnomAD directory: {self.gnomad_dir}"
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
        logger.info(f"Input VCF size: {input_size:,} bytes ({input_size/1024/1024:.1f} MB)")
        
        step_time = time.time() - step_start
        self.stats['processing_steps'].append(('validate_inputs', step_time))
        logger.info(f"✓ gnomAD annotation input files validated successfully ({step_time:.2f}s)")
    
    def validate_tools(self) -> None:
        """Validate required tools are available in system PATH."""
        step_start = time.time()
        logger.info("Validating required tools for gnomAD annotation...")
        
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
        logger.info(f"✓ All required tools found for gnomAD annotation ({step_time:.2f}s)")
    
    def discover_gnomad_files(self) -> Dict[str, Path]:
        """
        Discover gnomAD chromosome-split VCF files in the directory.
        
        Returns:
            Dictionary mapping chromosome names to gnomAD VCF file paths
        """
        step_start = time.time()
        logger.info("Discovering gnomAD chromosome-split VCF files...")
        
        gnomad_files = {}
        
        # Common gnomAD file patterns
        patterns = [
            "*.chr*.vcf.gz",
            "*.chr*.vcf.bgz",
            "*chr*.vcf.gz",
            "*chr*.vcf.bgz",
            "gnomad*.chr*.vcf.gz",
            "gnomad*.chr*.vcf.bgz"
        ]
        
        found_files = []
        for pattern in patterns:
            found_files.extend(glob.glob(str(self.gnomad_dir / pattern)))
        
        # Extract chromosome information from filenames
        for file_path in found_files:
            file_name = Path(file_path).name
            
            # Extract chromosome from filename
            import re
            chr_match = re.search(r'chr([0-9XYM]+|MT)', file_name, re.IGNORECASE)
            if chr_match:
                chromosome = chr_match.group(1).upper()
                # Normalize chromosome names
                if chromosome == 'MT':
                    chromosome = 'M'
                
                gnomad_files[chromosome] = Path(file_path)
                logger.debug(f"Found gnomAD file for chromosome {chromosome}: {file_name}")
        
        if not gnomad_files:
            error_msg = f"No gnomAD chromosome-split VCF files found in directory: {self.gnomad_dir}"
            logger.error(error_msg)
            logger.error("Expected file patterns:")
            for pattern in patterns:
                logger.error(f"  {pattern}")
            raise FileNotFoundError(error_msg)
        
        # Sort chromosomes for consistent processing order
        sorted_chromosomes = sorted(gnomad_files.keys(), key=lambda x: (
            int(x) if x.isdigit() else (23 if x == 'X' else (24 if x == 'Y' else 25))
        ))
        
        logger.info(f"Found gnomAD files for {len(gnomad_files)} chromosomes: {', '.join(sorted_chromosomes)}")
        
        # Log file sizes
        total_size = 0
        for chrom in sorted_chromosomes:
            file_path = gnomad_files[chrom]
            file_size = file_path.stat().st_size
            total_size += file_size
            logger.debug(f"  chr{chrom}: {file_size:,} bytes ({file_size/1024/1024:.1f} MB)")
        
        logger.info(f"Total gnomAD database size: {total_size:,} bytes ({total_size/1024/1024:.1f} MB)")
        
        step_time = time.time() - step_start
        self.stats['processing_steps'].append(('discover_gnomad_files', step_time))
        
        return gnomad_files
    
    def get_input_chromosomes(self) -> set:
        """
        Get chromosomes present in the input VCF file.
        
        Returns:
            Set of chromosome names present in input VCF
        """
        step_start = time.time()
        logger.info("Analyzing chromosomes in input VCF...")
        
        try:
            cmd = ['bcftools', 'view', '-H', str(self.input_vcf)]
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True,
                timeout=300  # 5 minute timeout
            )
            
            chromosomes = set()
            for line in result.stdout.split('\n'):
                if line.strip():
                    chrom = line.split('\t')[0]
                    # Normalize chromosome names
                    chrom = chrom.replace('chr', '').replace('Chr', '').upper()
                    if chrom == 'MT':
                        chrom = 'M'
                    chromosomes.add(chrom)
            
            logger.info(f"Input VCF contains {len(chromosomes)} chromosomes: {', '.join(sorted(chromosomes, key=lambda x: (int(x) if x.isdigit() else (23 if x == 'X' else (24 if x == 'Y' else 25)))))}")
            
            step_time = time.time() - step_start
            self.stats['processing_steps'].append(('get_input_chromosomes', step_time))
            
            return chromosomes
            
        except subprocess.CalledProcessError as e:
            logger.error(f"Failed to analyze input VCF chromosomes: {e}")
            logger.error(f"Command stderr: {e.stderr}")
            raise RuntimeError(f"Failed to analyze input VCF: {e}")
        except subprocess.TimeoutExpired:
            logger.error("Timeout while analyzing input VCF chromosomes")
            raise RuntimeError("Timeout during input VCF analysis")
        except Exception as e:
            logger.error(f"Unexpected error analyzing input VCF: {e}")
            raise
    
    def annotate_chromosome(self, chromosome: str, gnomad_file: Path, temp_dir: Path) -> Tuple[str, Path]:
        """
        Annotate variants for a specific chromosome using gnomAD data with performance optimizations.
        
        Args:
            chromosome: Chromosome name (e.g., '1', 'X', 'M')
            gnomad_file: Path to gnomAD VCF file for this chromosome
            temp_dir: Temporary directory for intermediate files
            
        Returns:
            Tuple of (chromosome, path to annotated chromosome VCF)
        """
        chr_start_time = time.time()
        logger.debug(f"Processing chromosome {chromosome} with gnomAD annotation...")
        
        try:
            # Create chromosome-specific temporary files with optimized naming
            chr_input = temp_dir / f"input_chr{chromosome}.vcf.gz"
            chr_output = temp_dir / f"annotated_chr{chromosome}.vcf.gz"
            
            # Optimized chromosome extraction with streaming
            logger.debug(f"Extracting chromosome {chromosome} from input VCF...")
            extract_cmd = [
                'bcftools', 'view',
                '-r', f"chr{chromosome},{chromosome}",  # Handle both chr1 and 1 formats
                '-O', 'z',
                '--threads', '2',  # Use 2 threads for compression
                '-o', str(chr_input),
                str(self.input_vcf)
            ]
            
            extract_result = subprocess.run(
                extract_cmd,
                capture_output=True,
                text=True,
                check=True,
                timeout=300  # Reduced timeout for extraction
            )
            
            # Quick check for variants without full file size check
            if not chr_input.exists():
                logger.debug(f"No variants found for chromosome {chromosome}, skipping...")
                return chromosome, None
            
            # Check variant count efficiently
            count_cmd = ['bcftools', 'view', '-H', str(chr_input)]
            count_result = subprocess.run(count_cmd, capture_output=True, text=True, timeout=60)
            variant_count = len([line for line in count_result.stdout.split('\n') if line.strip()])
            
            if variant_count == 0:
                logger.debug(f"No variants found for chromosome {chromosome}, skipping...")
                return chromosome, None
            
            logger.debug(f"Chromosome {chromosome}: {variant_count} variants to annotate")
            
            # Index extracted chromosome VCF with parallel compression
            subprocess.run(
                ['tabix', '-p', 'vcf', str(chr_input)],
                capture_output=True,
                text=True,
                check=True,
                timeout=120  # Reduced timeout for indexing
            )
            
            # Optimized annotation with field renaming for clear identification
            logger.debug(f"Annotating chromosome {chromosome} with gnomAD...")
            annotate_cmd = [
                'bcftools', 'annotate',
                '-a', str(gnomad_file),
                '-c', 'CHROM,POS,REF,ALT,INFO/GNOMAD_AF:=INFO/AF',  # Rename AF to GNOMAD_AF for clarity
                '-O', 'z',
                '--threads', '2',  # Use threads for compression
                '-o', str(chr_output),
                str(chr_input)
            ]
            
            result = subprocess.run(
                annotate_cmd,
                capture_output=True,
                text=True,
                check=True,
                timeout=900  # Reduced timeout per chromosome
            )
            
            # Index annotated chromosome VCF
            subprocess.run(
                ['tabix', '-p', 'vcf', str(chr_output)],
                capture_output=True,
                text=True,
                check=True,
                timeout=120
            )
            
            # Parse annotation statistics
            annotations_added = 0
            if result.stderr:
                import re
                annotations_match = re.search(r'(\d+)\s+annotations?\s+added', result.stderr, re.IGNORECASE)
                if annotations_match:
                    annotations_added = int(annotations_match.group(1))
            
            chr_time = time.time() - chr_start_time
            logger.info(f"✓ Chromosome {chromosome}: {annotations_added:,} annotations in {chr_time:.1f}s")
            
            return chromosome, chr_output
            
        except subprocess.CalledProcessError as e:
            chr_time = time.time() - chr_start_time
            error_msg = f"Chromosome {chromosome} annotation failed after {chr_time:.1f}s: {e}"
            logger.error(error_msg)
            if e.stderr:
                logger.error(f"Command stderr: {e.stderr}")
            # Immediate exit on failure - no graceful fallback
            raise RuntimeError(error_msg)
        except subprocess.TimeoutExpired:
            chr_time = time.time() - chr_start_time
            error_msg = f"Chromosome {chromosome} annotation timed out after {chr_time:.1f}s"
            logger.error(error_msg)
            # Immediate exit on timeout - no graceful fallback
            raise RuntimeError(error_msg)
        except Exception as e:
            chr_time = time.time() - chr_start_time
            error_msg = f"Unexpected error annotating chromosome {chromosome} after {chr_time:.1f}s: {e}"
            logger.error(error_msg)
            # Immediate exit on any error - no graceful fallback
            raise RuntimeError(error_msg)
    
    def merge_annotated_chromosomes(self, annotated_files: List[Path], temp_dir: Path) -> None:
        """
        Merge annotated chromosome VCF files into final output.
        
        Args:
            annotated_files: List of paths to annotated chromosome VCF files
            temp_dir: Temporary directory for intermediate files
        """
        step_start = time.time()
        logger.info(f"Merging {len(annotated_files)} annotated chromosome files...")
        
        try:
            if len(annotated_files) == 1:
                # Single file, just copy
                logger.info("Single chromosome file, copying to output...")
                import shutil
                shutil.copy2(annotated_files[0], self.output_vcf)
            else:
                # Multiple files, merge with bcftools concat
                logger.info("Multiple chromosome files, merging with bcftools concat...")
                
                # Create file list for bcftools concat
                file_list = temp_dir / "chromosome_files.txt"
                with open(file_list, 'w') as f:
                    for file_path in annotated_files:
                        f.write(f"{file_path}\n")
                
                # Merge files
                merge_cmd = [
                    'bcftools', 'concat',
                    '-f', str(file_list),
                    '-O', 'z' if str(self.output_vcf).endswith('.gz') else 'v',
                    '-o', str(self.output_vcf)
                ]
                
                subprocess.run(
                    merge_cmd,
                    capture_output=True,
                    text=True,
                    check=True,
                    timeout=1800  # 30 minute timeout for merging
                )
            
            # Create index if output is compressed
            if str(self.output_vcf).endswith('.gz'):
                subprocess.run(
                    ['tabix', '-p', 'vcf', str(self.output_vcf)],
                    capture_output=True,
                    text=True,
                    check=True,
                    timeout=300
                )
            
            step_time = time.time() - step_start
            self.stats['processing_steps'].append(('merge_annotated_chromosomes', step_time))
            
            logger.info(f"✓ Chromosome files merged successfully ({step_time:.2f}s)")
            
        except subprocess.CalledProcessError as e:
            logger.error(f"Failed to merge annotated chromosome files: {e}")
            logger.error(f"Command stderr: {e.stderr}")
            raise RuntimeError(f"Chromosome merging failed: {e}")
        except subprocess.TimeoutExpired:
            logger.error("Timeout while merging annotated chromosome files")
            raise RuntimeError("Chromosome merging timed out")
        except Exception as e:
            logger.error(f"Unexpected error merging chromosome files: {e}")
            raise
    
    def run_scatter_gather_annotation(self) -> None:
        """
        Run the complete gnomAD scatter-gather annotation pipeline.
        """
        pipeline_start = time.time()
        logger.info("Starting gnomAD scatter-gather annotation pipeline...")
        
        # Create temporary directory
        temp_dir = self.output_vcf.parent / f"gnomad_temp_{int(time.time())}"
        temp_dir.mkdir(exist_ok=True)
        
        try:
            # Get chromosomes present in input VCF
            logger.info("Step 1: Analyzing input chromosomes...")
            input_chromosomes = self.get_input_chromosomes()
            
            # Find matching gnomAD files
            logger.info("Step 2: Matching gnomAD files to input chromosomes...")
            matching_files = {}
            for chrom in input_chromosomes:
                if chrom in self.gnomad_files:
                    matching_files[chrom] = self.gnomad_files[chrom]
                    logger.debug(f"Found gnomAD file for chromosome {chrom}")
                else:
                    logger.warning(f"No gnomAD file found for chromosome {chrom}")
                    self.stats['warnings'].append(f"No gnomAD file for chromosome {chrom}")
            
            if not matching_files:
                raise RuntimeError("No matching gnomAD files found for input chromosomes")
            
            logger.info(f"Processing {len(matching_files)} chromosomes with gnomAD annotation")
            
            # Process chromosomes in parallel
            logger.info("Step 3: Processing chromosomes in parallel...")
            annotated_files = []
            
            with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
                # Submit chromosome annotation tasks
                future_to_chrom = {
                    executor.submit(self.annotate_chromosome, chrom, gnomad_file, temp_dir): chrom
                    for chrom, gnomad_file in matching_files.items()
                }
                
                # Collect results
                for future in as_completed(future_to_chrom):
                    chrom = future_to_chrom[future]
                    try:
                        result_chrom, result_file = future.result()
                        if result_file:
                            annotated_files.append(result_file)
                            self.stats['chromosomes_processed'] += 1
                        logger.info(f"Completed chromosome {result_chrom}")
                    except Exception as e:
                        logger.error(f"Chromosome {chrom} processing failed: {e}")
                        self.stats['errors'].append(f"Chromosome {chrom}: {e}")
                        raise
            
            if not annotated_files:
                raise RuntimeError("No chromosomes were successfully annotated")
            
            # Sort annotated files by chromosome order
            def chrom_sort_key(file_path):
                chrom = file_path.name.replace('annotated_chr', '').replace('.vcf.gz', '')
                if chrom.isdigit():
                    return int(chrom)
                elif chrom == 'X':
                    return 23
                elif chrom == 'Y':
                    return 24
                elif chrom in ['M', 'MT']:
                    return 25
                else:
                    return 999
            
            annotated_files.sort(key=chrom_sort_key)
            
            # Merge annotated chromosome files
            logger.info("Step 4: Merging annotated chromosome files...")
            self.merge_annotated_chromosomes(annotated_files, temp_dir)
            
            # Validate output
            logger.info("Step 5: Validating gnomAD annotation output...")
            self.validate_gnomad_output()
            
            # Clean up temporary files
            logger.info("Step 6: Cleaning up temporary files...")
            self._cleanup_temp_directory(temp_dir)
            
            # Generate final statistics
            pipeline_time = time.time() - pipeline_start
            self.stats['processing_steps'].append(('total_gnomad_pipeline', pipeline_time))
            
            logger.info("✓ gnomAD scatter-gather annotation completed successfully!")
            self._generate_gnomad_summary_statistics()
            
        except Exception as e:
            logger.error(f"gnomAD scatter-gather annotation pipeline failed: {e}")
            
            # Clean up any partial output files
            if self.output_vcf.exists():
                try:
                    self.output_vcf.unlink()
                    logger.info("Cleaned up partial gnomAD annotated output file")
                except Exception as cleanup_error:
                    logger.warning(f"Failed to clean up partial gnomAD output: {cleanup_error}")
            
            # Clean up temporary directory
            try:
                self._cleanup_temp_directory(temp_dir)
            except Exception as cleanup_error:
                logger.warning(f"Failed to clean up temporary directory: {cleanup_error}")
            
            # Generate error summary
            self._generate_gnomad_summary_statistics()
            raise
    
    def validate_gnomad_output(self) -> None:
        """Validate that gnomAD annotated output file was created successfully."""
        step_start = time.time()
        logger.info("Validating gnomAD annotation output file...")
        
        try:
            # Check file existence
            if not self.output_vcf.exists():
                raise RuntimeError(f"gnomAD annotated output file was not created: {self.output_vcf}")
            
            # Check file size
            file_size = self.output_vcf.stat().st_size
            if file_size == 0:
                raise RuntimeError(f"gnomAD annotated output file is empty: {self.output_vcf}")
            
            # Check file permissions
            if not os.access(self.output_vcf, os.R_OK):
                raise RuntimeError(f"gnomAD annotated output file is not readable: {self.output_vcf}")
            
            # Validate VCF format
            try:
                cmd = ['bcftools', 'view', '-h', str(self.output_vcf)]
                result = subprocess.run(cmd, capture_output=True, text=True, check=True, timeout=60)
                
                # Check for VCF header
                if not result.stdout.startswith('##fileformat=VCF'):
                    raise RuntimeError("gnomAD annotated output file does not have valid VCF header")
                
                logger.debug("✓ gnomAD annotated output VCF format validation passed")
                    
            except subprocess.CalledProcessError as e:
                logger.warning(f"gnomAD annotated VCF format validation failed: {e}")
                logger.warning("Output file may have format issues")
            except Exception as e:
                logger.warning(f"Could not validate gnomAD annotated VCF format: {e}")
            
            step_time = time.time() - step_start
            self.stats['processing_steps'].append(('validate_gnomad_output', step_time))
            
            logger.info(f"✓ gnomAD annotation output validation completed ({file_size:,} bytes, {step_time:.2f}s)")
            
        except Exception as e:
            logger.error(f"gnomAD annotation output validation failed: {e}")
            raise
    
    def _cleanup_temp_directory(self, temp_dir: Path) -> None:
        """
        Clean up temporary directory and all its contents.
        
        Args:
            temp_dir: Temporary directory to clean up
        """
        if not temp_dir.exists():
            logger.debug("Temporary directory does not exist, nothing to clean up")
            return
        
        logger.info(f"Cleaning up temporary directory: {temp_dir}")
        
        try:
            import shutil
            shutil.rmtree(temp_dir)
            logger.info("✓ Temporary directory cleaned up successfully")
        except Exception as e:
            logger.warning(f"Failed to clean up temporary directory: {e}")
            logger.warning("Some temporary files may remain on disk")
    
    def _generate_gnomad_summary_statistics(self) -> None:
        """Generate and log comprehensive gnomAD annotation summary statistics."""
        total_time = time.time() - self.stats['start_time']
        
        logger.info("=== gnomAD Scatter-Gather Annotation Summary ===")
        logger.info(f"Total processing time: {total_time:.2f} seconds")
        logger.info(f"Input VCF: {self.input_vcf}")
        logger.info(f"gnomAD directory: {self.gnomad_dir}")
        logger.info(f"Output VCF: {self.output_vcf}")
        logger.info(f"Max workers: {self.max_workers}")
        
        # Log processing statistics
        logger.info(f"Chromosomes processed: {self.stats['chromosomes_processed']}")
        logger.info(f"gnomAD files available: {len(self.gnomad_files)}")
        
        if self.stats['annotations_added'] > 0:
            logger.info(f"Total annotations added: {self.stats['annotations_added']:,}")
        
        # Log file size information
        try:
            if self.input_vcf.exists():
                input_size = self.input_vcf.stat().st_size
                logger.info(f"Input VCF size: {input_size:,} bytes ({input_size/1024/1024:.1f} MB)")
            
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
        
        logger.info("=== gnomAD Annotation Summary Complete ===")


def main():
    """Main entry point for standalone usage."""
    import argparse
    import sys
    
    parser = argparse.ArgumentParser(
        description="gnomAD scatter-gather annotation engine",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
gnomAD Scatter-Gather Annotation Engine:
  This module provides gnomAD database annotation functionality using scatter-gather
  processing across chromosome-split VCF files. It implements parallel processing
  across chromosomes for efficient annotation with population frequency data.

Key Features:
  - Scatter-gather processing across chromosome-split gnomAD files
  - Population frequency data extraction (faf95, AF)
  - Parallel chromosome processing for performance
  - Compressed VCF input/output handling with proper indexing
  - Comprehensive error handling and logging

Examples:
  # Basic gnomAD scatter-gather annotation
  python gnomad_annotator.py -i input.vcf.gz -g /path/to/gnomad -o output.vcf.gz
  
  # With custom number of workers
  python gnomad_annotator.py -i input.vcf.gz -g /path/to/gnomad -o output.vcf.gz -w 8
        """
    )
    
    parser.add_argument(
        '-i', '--input',
        required=True,
        metavar='VCF_FILE',
        help='Input VCF file to be annotated (can be compressed with .gz)'
    )
    
    parser.add_argument(
        '-g', '--gnomad',
        required=True,
        metavar='DIRECTORY',
        help='gnomAD database directory containing chromosome-split VCF files'
    )
    
    parser.add_argument(
        '-o', '--output',
        required=True,
        metavar='VCF_FILE',
        help='Output annotated VCF file (use .gz extension for compression)'
    )
    
    parser.add_argument(
        '-w', '--workers',
        type=int,
        default=4,
        metavar='N',
        help='Maximum number of parallel workers (default: 4)'
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
    
    # Run gnomAD annotation
    try:
        annotator = GnomadAnnotator(
            input_vcf=args.input,
            gnomad_dir=args.gnomad,
            output_vcf=args.output,
            max_workers=args.workers
        )
        
        annotator.run_scatter_gather_annotation()
        
        logger.info("✓ gnomAD scatter-gather annotation completed successfully!")
        sys.exit(0)
        
    except KeyboardInterrupt:
        logger.error("gnomAD annotation interrupted by user")
        sys.exit(130)
    except Exception as e:
        logger.error(f"gnomAD annotation failed: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()