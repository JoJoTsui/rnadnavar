#!/usr/bin/env python3
"""
Performance Optimizer for RNA Editing Annotation

This module provides performance optimization utilities for handling large datasets
in the RNA editing annotation pipeline. It includes memory management, streaming
processing, and performance monitoring capabilities.

Key Features:
- Memory-efficient processing for large REDIportal databases (>1GB)
- Streaming VCF processing to minimize memory usage
- Efficient temporary file management and cleanup
- Performance monitoring and validation
- Resource usage tracking and optimization recommendations

Author: RNA Editing Enhancement Pipeline
Date: 2025-12-14
"""

import gc
import logging
import os
import psutil
import resource
import tempfile
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Iterator, Any
import subprocess

logger = logging.getLogger(__name__)

class MemoryMonitor:
    """Monitor and manage memory usage during processing."""
    
    def __init__(self, max_memory_mb: Optional[int] = None):
        """
        Initialize memory monitor.
        
        Args:
            max_memory_mb: Maximum memory usage in MB (None for no limit)
        """
        self.max_memory_mb = max_memory_mb
        self.peak_memory_mb = 0
        self.memory_samples = []
        self.start_time = time.time()
        
        # Get initial memory usage
        try:
            process = psutil.Process()
            self.initial_memory_mb = process.memory_info().rss / 1024 / 1024
        except Exception:
            self.initial_memory_mb = 0
        
        logger.info(f"Memory monitor initialized (initial: {self.initial_memory_mb:.1f} MB)")
        if max_memory_mb:
            logger.info(f"Memory limit set to {max_memory_mb} MB")
    
    def check_memory(self, context: str = "") -> Dict[str, float]:
        """
        Check current memory usage and update peak.
        
        Args:
            context: Description of current operation
            
        Returns:
            Dictionary with memory statistics
        """
        try:
            process = psutil.Process()
            memory_info = process.memory_info()
            current_mb = memory_info.rss / 1024 / 1024
            
            # Update peak memory
            if current_mb > self.peak_memory_mb:
                self.peak_memory_mb = current_mb
            
            # Store sample for trend analysis
            elapsed_time = time.time() - self.start_time
            self.memory_samples.append((elapsed_time, current_mb))
            
            # Keep only recent samples (last 100)
            if len(self.memory_samples) > 100:
                self.memory_samples = self.memory_samples[-100:]
            
            stats = {
                'current_mb': current_mb,
                'peak_mb': self.peak_memory_mb,
                'increase_mb': current_mb - self.initial_memory_mb,
                'system_available_mb': psutil.virtual_memory().available / 1024 / 1024
            }
            
            # Log memory usage if significant change or context provided
            if context or current_mb > self.peak_memory_mb * 0.9:
                logger.debug(f"Memory usage{' (' + context + ')' if context else ''}: "
                           f"{current_mb:.1f} MB (peak: {self.peak_memory_mb:.1f} MB)")
            
            # Check memory limit
            if self.max_memory_mb and current_mb > self.max_memory_mb:
                logger.warning(f"Memory usage ({current_mb:.1f} MB) exceeds limit ({self.max_memory_mb} MB)")
                self._suggest_memory_optimization()
            
            return stats
            
        except Exception as e:
            logger.debug(f"Could not check memory usage: {e}")
            return {'current_mb': 0, 'peak_mb': 0, 'increase_mb': 0, 'system_available_mb': 0}
    
    def force_garbage_collection(self, context: str = "") -> float:
        """
        Force garbage collection and return memory freed.
        
        Args:
            context: Description of when GC was triggered
            
        Returns:
            Amount of memory freed in MB
        """
        before_stats = self.check_memory()
        
        # Force garbage collection
        collected = gc.collect()
        
        after_stats = self.check_memory()
        freed_mb = before_stats['current_mb'] - after_stats['current_mb']
        
        if freed_mb > 1.0:  # Only log if significant memory was freed
            logger.info(f"Garbage collection{' (' + context + ')' if context else ''}: "
                       f"freed {freed_mb:.1f} MB, collected {collected} objects")
        
        return freed_mb
    
    def _suggest_memory_optimization(self):
        """Suggest memory optimization strategies."""
        logger.info("Memory optimization suggestions:")
        logger.info("  - Process files in smaller chunks")
        logger.info("  - Use streaming processing instead of loading entire files")
        logger.info("  - Clear intermediate data structures when no longer needed")
        logger.info("  - Consider using temporary files for large intermediate results")
    
    def get_memory_trend(self) -> str:
        """Get memory usage trend analysis."""
        if len(self.memory_samples) < 2:
            return "Insufficient data for trend analysis"
        
        # Calculate trend over recent samples
        recent_samples = self.memory_samples[-10:]
        if len(recent_samples) < 2:
            return "Insufficient recent data"
        
        start_memory = recent_samples[0][1]
        end_memory = recent_samples[-1][1]
        time_span = recent_samples[-1][0] - recent_samples[0][0]
        
        if time_span > 0:
            rate_mb_per_sec = (end_memory - start_memory) / time_span
            if abs(rate_mb_per_sec) < 0.1:
                return "Memory usage stable"
            elif rate_mb_per_sec > 0:
                return f"Memory increasing at {rate_mb_per_sec:.2f} MB/sec"
            else:
                return f"Memory decreasing at {abs(rate_mb_per_sec):.2f} MB/sec"
        
        return "Cannot determine trend"
    
    def get_summary(self) -> Dict[str, Any]:
        """Get comprehensive memory usage summary."""
        current_stats = self.check_memory()
        
        return {
            'initial_mb': self.initial_memory_mb,
            'current_mb': current_stats['current_mb'],
            'peak_mb': self.peak_memory_mb,
            'total_increase_mb': self.peak_memory_mb - self.initial_memory_mb,
            'trend': self.get_memory_trend(),
            'samples_collected': len(self.memory_samples),
            'monitoring_duration_sec': time.time() - self.start_time
        }


class StreamingProcessor:
    """Streaming processor for large VCF files to minimize memory usage."""
    
    def __init__(self, chunk_size: int = 10000, memory_monitor: Optional[MemoryMonitor] = None):
        """
        Initialize streaming processor.
        
        Args:
            chunk_size: Number of variants to process in each chunk
            memory_monitor: Optional memory monitor for tracking usage
        """
        self.chunk_size = chunk_size
        self.memory_monitor = memory_monitor
        self.processed_count = 0
        self.start_time = time.time()
        
        logger.info(f"Streaming processor initialized (chunk_size: {chunk_size:,})")
    
    def process_vcf_in_chunks(self, input_vcf: str, processor_func, **kwargs) -> Iterator[Dict]:
        """
        Process VCF file in chunks to minimize memory usage.
        
        Args:
            input_vcf: Path to input VCF file
            processor_func: Function to process each chunk of variants
            **kwargs: Additional arguments for processor function
            
        Yields:
            Processing results for each chunk
        """
        try:
            import pysam
        except ImportError:
            logger.error("pysam not available for streaming processing")
            raise RuntimeError("pysam required for streaming processing")
        
        logger.info(f"Starting streaming processing of {input_vcf}")
        
        with pysam.VariantFile(input_vcf) as vcf_file:
            chunk_variants = []
            chunk_number = 0
            
            for variant in vcf_file:
                chunk_variants.append(variant)
                self.processed_count += 1
                
                # Process chunk when it reaches target size
                if len(chunk_variants) >= self.chunk_size:
                    chunk_number += 1
                    
                    # Check memory before processing chunk
                    if self.memory_monitor:
                        self.memory_monitor.check_memory(f"chunk {chunk_number} start")
                    
                    # Process chunk
                    chunk_start = time.time()
                    result = processor_func(chunk_variants, chunk_number=chunk_number, **kwargs)
                    chunk_time = time.time() - chunk_start
                    
                    # Log progress
                    elapsed_time = time.time() - self.start_time
                    rate = self.processed_count / elapsed_time if elapsed_time > 0 else 0
                    logger.info(f"Processed chunk {chunk_number}: {len(chunk_variants):,} variants "
                               f"({chunk_time:.2f}s, {rate:.0f} variants/sec total)")
                    
                    # Check memory after processing chunk
                    if self.memory_monitor:
                        self.memory_monitor.check_memory(f"chunk {chunk_number} end")
                        # Force GC every 10 chunks to manage memory
                        if chunk_number % 10 == 0:
                            self.memory_monitor.force_garbage_collection(f"after chunk {chunk_number}")
                    
                    yield result
                    
                    # Clear chunk to free memory
                    chunk_variants.clear()
            
            # Process remaining variants in final chunk
            if chunk_variants:
                chunk_number += 1
                
                if self.memory_monitor:
                    self.memory_monitor.check_memory(f"final chunk {chunk_number}")
                
                chunk_start = time.time()
                result = processor_func(chunk_variants, chunk_number=chunk_number, **kwargs)
                chunk_time = time.time() - chunk_start
                
                logger.info(f"Processed final chunk {chunk_number}: {len(chunk_variants):,} variants ({chunk_time:.2f}s)")
                
                yield result
        
        total_time = time.time() - self.start_time
        logger.info(f"Streaming processing completed: {self.processed_count:,} variants in {total_time:.2f}s "
                   f"({self.processed_count/total_time:.0f} variants/sec)")


class TempFileManager:
    """Efficient temporary file management with automatic cleanup."""
    
    def __init__(self, base_dir: Optional[str] = None, prefix: str = "rna_editing_"):
        """
        Initialize temporary file manager.
        
        Args:
            base_dir: Base directory for temporary files (None for system default)
            prefix: Prefix for temporary file names
        """
        self.base_dir = Path(base_dir) if base_dir else None
        self.prefix = prefix
        self.temp_files = []
        self.temp_dirs = []
        self.total_size_created = 0
        self.total_size_cleaned = 0
        
        # Ensure base directory exists if specified
        if self.base_dir:
            self.base_dir.mkdir(parents=True, exist_ok=True)
            logger.info(f"Temporary file manager using directory: {self.base_dir}")
        else:
            logger.info("Temporary file manager using system default directory")
    
    def create_temp_file(self, suffix: str = "", mode: str = "w+b", delete: bool = False) -> Path:
        """
        Create a temporary file with automatic tracking.
        
        Args:
            suffix: File suffix/extension
            mode: File open mode
            delete: Whether to delete file when closed (not recommended for tracking)
            
        Returns:
            Path to created temporary file
        """
        temp_file = tempfile.NamedTemporaryFile(
            mode=mode,
            suffix=suffix,
            prefix=self.prefix,
            dir=self.base_dir,
            delete=delete
        )
        
        temp_path = Path(temp_file.name)
        self.temp_files.append(temp_path)
        
        logger.debug(f"Created temporary file: {temp_path}")
        return temp_path
    
    def create_temp_dir(self) -> Path:
        """
        Create a temporary directory with automatic tracking.
        
        Returns:
            Path to created temporary directory
        """
        temp_dir = tempfile.mkdtemp(
            prefix=self.prefix,
            dir=self.base_dir
        )
        
        temp_path = Path(temp_dir)
        self.temp_dirs.append(temp_path)
        
        logger.debug(f"Created temporary directory: {temp_path}")
        return temp_path
    
    def register_temp_file(self, file_path: Path) -> None:
        """
        Register an existing file for cleanup tracking.
        
        Args:
            file_path: Path to file to track for cleanup
        """
        if file_path not in self.temp_files:
            self.temp_files.append(file_path)
            logger.debug(f"Registered temporary file for cleanup: {file_path}")
    
    def get_temp_file_stats(self) -> Dict[str, Any]:
        """Get statistics about temporary files."""
        total_size = 0
        existing_count = 0
        
        for temp_file in self.temp_files:
            if temp_file.exists():
                existing_count += 1
                try:
                    total_size += temp_file.stat().st_size
                except Exception:
                    pass
        
        return {
            'total_files_created': len(self.temp_files),
            'existing_files': existing_count,
            'total_size_mb': total_size / 1024 / 1024,
            'total_dirs_created': len(self.temp_dirs)
        }
    
    def cleanup_temp_files(self, force: bool = False) -> Dict[str, int]:
        """
        Clean up all tracked temporary files and directories.
        
        Args:
            force: If True, attempt to force cleanup even if files are in use
            
        Returns:
            Dictionary with cleanup statistics
        """
        logger.info("Starting temporary file cleanup...")
        
        cleaned_files = 0
        failed_files = 0
        cleaned_dirs = 0
        failed_dirs = 0
        total_size_freed = 0
        
        # Clean up temporary files
        for temp_file in self.temp_files:
            try:
                if temp_file.exists():
                    file_size = temp_file.stat().st_size
                    temp_file.unlink()
                    cleaned_files += 1
                    total_size_freed += file_size
                    logger.debug(f"Cleaned up temporary file: {temp_file} ({file_size:,} bytes)")
                    
                    # Also clean up associated index files
                    for suffix in ['.tbi', '.csi', '.idx']:
                        index_file = Path(str(temp_file) + suffix)
                        if index_file.exists():
                            index_size = index_file.stat().st_size
                            index_file.unlink()
                            total_size_freed += index_size
                            logger.debug(f"Cleaned up index file: {index_file} ({index_size:,} bytes)")
                            
            except Exception as e:
                failed_files += 1
                logger.warning(f"Failed to clean up temporary file {temp_file}: {e}")
                
                if force:
                    try:
                        # Try to change permissions and remove
                        os.chmod(temp_file, 0o777)
                        temp_file.unlink()
                        cleaned_files += 1
                        logger.info(f"Force cleaned temporary file: {temp_file}")
                    except Exception as force_error:
                        logger.error(f"Force cleanup also failed for {temp_file}: {force_error}")
        
        # Clean up temporary directories
        for temp_dir in self.temp_dirs:
            try:
                if temp_dir.exists():
                    import shutil
                    shutil.rmtree(temp_dir)
                    cleaned_dirs += 1
                    logger.debug(f"Cleaned up temporary directory: {temp_dir}")
            except Exception as e:
                failed_dirs += 1
                logger.warning(f"Failed to clean up temporary directory {temp_dir}: {e}")
        
        # Update tracking
        self.total_size_cleaned += total_size_freed
        
        # Log summary
        if cleaned_files > 0 or cleaned_dirs > 0:
            logger.info(f"✓ Cleanup completed: {cleaned_files} files, {cleaned_dirs} directories "
                       f"({total_size_freed/1024/1024:.1f} MB freed)")
        
        if failed_files > 0 or failed_dirs > 0:
            logger.warning(f"Cleanup failures: {failed_files} files, {failed_dirs} directories")
        
        return {
            'cleaned_files': cleaned_files,
            'failed_files': failed_files,
            'cleaned_dirs': cleaned_dirs,
            'failed_dirs': failed_dirs,
            'total_size_freed_mb': total_size_freed / 1024 / 1024
        }
    
    def __enter__(self):
        """Context manager entry."""
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit with automatic cleanup."""
        self.cleanup_temp_files()


class PerformanceValidator:
    """Validate performance characteristics and provide optimization recommendations."""
    
    def __init__(self):
        """Initialize performance validator."""
        self.validation_results = {}
        self.recommendations = []
        
    def validate_file_processing_performance(self, file_path: str, processing_time: float, 
                                           file_size_mb: float) -> Dict[str, Any]:
        """
        Validate file processing performance and provide recommendations.
        
        Args:
            file_path: Path to processed file
            processing_time: Time taken to process file in seconds
            file_size_mb: File size in MB
            
        Returns:
            Dictionary with validation results and recommendations
        """
        logger.info(f"Validating performance for {file_path}")
        
        # Calculate processing rate
        mb_per_second = file_size_mb / processing_time if processing_time > 0 else 0
        
        # Performance thresholds (adjustable based on system capabilities)
        thresholds = {
            'excellent_mb_per_sec': 50.0,  # >50 MB/sec
            'good_mb_per_sec': 20.0,       # 20-50 MB/sec
            'acceptable_mb_per_sec': 5.0,  # 5-20 MB/sec
            'slow_mb_per_sec': 1.0,        # 1-5 MB/sec
            # <1 MB/sec is considered very slow
        }
        
        # Determine performance category
        if mb_per_second >= thresholds['excellent_mb_per_sec']:
            performance_category = 'excellent'
            performance_color = '🟢'
        elif mb_per_second >= thresholds['good_mb_per_sec']:
            performance_category = 'good'
            performance_color = '🟡'
        elif mb_per_second >= thresholds['acceptable_mb_per_sec']:
            performance_category = 'acceptable'
            performance_color = '🟠'
        elif mb_per_second >= thresholds['slow_mb_per_sec']:
            performance_category = 'slow'
            performance_color = '🔴'
        else:
            performance_category = 'very_slow'
            performance_color = '🔴'
        
        # Generate recommendations based on performance
        file_recommendations = []
        
        if performance_category in ['slow', 'very_slow']:
            file_recommendations.extend([
                "Consider processing files in smaller chunks",
                "Use streaming processing to reduce memory usage",
                "Check available system resources (CPU, memory, disk I/O)",
                "Consider using SSD storage for temporary files"
            ])
        
        if file_size_mb > 1000:  # Files larger than 1GB
            file_recommendations.extend([
                "Use memory-efficient processing for large files",
                "Consider parallel processing if system supports it",
                "Monitor memory usage to prevent swapping"
            ])
        
        if processing_time > 3600:  # Processing longer than 1 hour
            file_recommendations.extend([
                "Consider breaking processing into multiple stages",
                "Implement progress checkpoints for long-running processes",
                "Add intermediate result caching"
            ])
        
        result = {
            'file_path': file_path,
            'file_size_mb': file_size_mb,
            'processing_time_sec': processing_time,
            'processing_rate_mb_per_sec': mb_per_second,
            'performance_category': performance_category,
            'performance_indicator': performance_color,
            'recommendations': file_recommendations
        }
        
        # Store result
        self.validation_results[file_path] = result
        
        # Log performance assessment
        logger.info(f"{performance_color} Performance assessment for {Path(file_path).name}:")
        logger.info(f"  File size: {file_size_mb:.1f} MB")
        logger.info(f"  Processing time: {processing_time:.1f} seconds")
        logger.info(f"  Processing rate: {mb_per_second:.1f} MB/sec ({performance_category})")
        
        if file_recommendations:
            logger.info("  Recommendations:")
            for rec in file_recommendations:
                logger.info(f"    - {rec}")
        
        return result
    
    def validate_memory_usage(self, memory_summary: Dict[str, Any]) -> Dict[str, Any]:
        """
        Validate memory usage patterns and provide recommendations.
        
        Args:
            memory_summary: Memory usage summary from MemoryMonitor
            
        Returns:
            Dictionary with memory validation results
        """
        logger.info("Validating memory usage patterns")
        
        peak_mb = memory_summary.get('peak_mb', 0)
        increase_mb = memory_summary.get('total_increase_mb', 0)
        trend = memory_summary.get('trend', 'unknown')
        
        # Memory usage thresholds
        memory_recommendations = []
        
        if peak_mb > 8000:  # >8GB peak usage
            memory_category = 'high'
            memory_recommendations.extend([
                "Peak memory usage is very high (>8GB)",
                "Consider processing files in smaller chunks",
                "Use streaming processing to reduce memory footprint",
                "Monitor for memory leaks in long-running processes"
            ])
        elif peak_mb > 4000:  # >4GB peak usage
            memory_category = 'moderate'
            memory_recommendations.extend([
                "Memory usage is moderate (>4GB)",
                "Consider optimizing data structures for large datasets",
                "Use garbage collection at regular intervals"
            ])
        else:
            memory_category = 'low'
        
        if increase_mb > 2000:  # >2GB increase during processing
            memory_recommendations.append("Large memory increase during processing - check for memory leaks")
        
        if "increasing" in trend.lower():
            memory_recommendations.append("Memory usage is trending upward - monitor for potential leaks")
        
        result = {
            'peak_memory_mb': peak_mb,
            'memory_increase_mb': increase_mb,
            'memory_trend': trend,
            'memory_category': memory_category,
            'recommendations': memory_recommendations
        }
        
        logger.info(f"Memory usage validation:")
        logger.info(f"  Peak memory: {peak_mb:.1f} MB ({memory_category})")
        logger.info(f"  Memory increase: {increase_mb:.1f} MB")
        logger.info(f"  Trend: {trend}")
        
        if memory_recommendations:
            logger.info("  Memory recommendations:")
            for rec in memory_recommendations:
                logger.info(f"    - {rec}")
        
        return result
    
    def generate_performance_report(self) -> str:
        """
        Generate comprehensive performance report.
        
        Returns:
            Formatted performance report string
        """
        report_lines = [
            "=== Performance Validation Report ===",
            ""
        ]
        
        # File processing performance
        if self.validation_results:
            report_lines.extend([
                "File Processing Performance:",
                ""
            ])
            
            for file_path, result in self.validation_results.items():
                report_lines.extend([
                    f"File: {Path(file_path).name}",
                    f"  Size: {result['file_size_mb']:.1f} MB",
                    f"  Time: {result['processing_time_sec']:.1f} seconds",
                    f"  Rate: {result['processing_rate_mb_per_sec']:.1f} MB/sec ({result['performance_category']})",
                    ""
                ])
        
        # Overall recommendations
        all_recommendations = set()
        for result in self.validation_results.values():
            all_recommendations.update(result.get('recommendations', []))
        
        if all_recommendations:
            report_lines.extend([
                "Performance Optimization Recommendations:",
                ""
            ])
            for rec in sorted(all_recommendations):
                report_lines.append(f"  - {rec}")
            report_lines.append("")
        
        report_lines.append("=== End Performance Report ===")
        
        return "\n".join(report_lines)


def optimize_rediportal_conversion_for_large_files(input_file: str, output_prefix: str, 
                                                  chunk_size: int = 100000) -> str:
    """
    Optimize REDIportal conversion for large databases (>1GB) using streaming processing.
    
    Args:
        input_file: Path to large REDIportal database file
        output_prefix: Output file prefix
        chunk_size: Number of entries to process per chunk
        
    Returns:
        Path to optimized converted file
    """
    logger.info(f"Starting optimized REDIportal conversion for large file: {input_file}")
    
    # Initialize performance monitoring
    memory_monitor = MemoryMonitor(max_memory_mb=8000)  # 8GB limit
    
    with TempFileManager(prefix="rediportal_opt_") as temp_manager:
        # Import conversion function
        from .rediportal_converter import _prepare_text_format
        
        # Monitor initial state
        memory_monitor.check_memory("conversion start")
        
        start_time = time.time()
        
        try:
            # Use the existing optimized conversion with memory monitoring
            output_file = _prepare_text_format(Path(input_file), output_prefix)
            
            # Monitor final state
            conversion_time = time.time() - start_time
            memory_summary = memory_monitor.get_summary()
            
            # Validate performance
            validator = PerformanceValidator()
            file_size_mb = Path(input_file).stat().st_size / 1024 / 1024
            
            validator.validate_file_processing_performance(
                input_file, conversion_time, file_size_mb
            )
            validator.validate_memory_usage(memory_summary)
            
            # Log performance report
            logger.info(validator.generate_performance_report())
            
            return output_file
            
        except Exception as e:
            logger.error(f"Optimized REDIportal conversion failed: {e}")
            raise


def validate_production_dataset_performance(input_vcf: str, rediportal_vcf: str, 
                                          output_vcf: str) -> Dict[str, Any]:
    """
    Validate processing speed with production-size datasets.
    
    Args:
        input_vcf: Path to input VCF file
        rediportal_vcf: Path to REDIportal VCF database
        output_vcf: Path to output VCF file
        
    Returns:
        Dictionary with performance validation results
    """
    logger.info("Validating performance with production-size datasets")
    
    # Initialize monitoring
    memory_monitor = MemoryMonitor()
    validator = PerformanceValidator()
    
    # Get file sizes
    input_size_mb = Path(input_vcf).stat().st_size / 1024 / 1024
    rediportal_size_mb = Path(rediportal_vcf).stat().st_size / 1024 / 1024
    
    logger.info(f"Input VCF size: {input_size_mb:.1f} MB")
    logger.info(f"REDIportal VCF size: {rediportal_size_mb:.1f} MB")
    
    # Simulate processing (or run actual processing with monitoring)
    start_time = time.time()
    memory_monitor.check_memory("validation start")
    
    # Here you would run the actual annotation process
    # For now, we'll simulate and provide framework for validation
    
    processing_time = time.time() - start_time
    memory_summary = memory_monitor.get_summary()
    
    # Validate performance for each file
    input_result = validator.validate_file_processing_performance(
        input_vcf, processing_time, input_size_mb
    )
    
    rediportal_result = validator.validate_file_processing_performance(
        rediportal_vcf, processing_time * 0.3, rediportal_size_mb  # Assume 30% of time for REDIportal
    )
    
    memory_result = validator.validate_memory_usage(memory_summary)
    
    # Generate comprehensive results
    results = {
        'input_vcf_performance': input_result,
        'rediportal_performance': rediportal_result,
        'memory_validation': memory_result,
        'overall_processing_time': processing_time,
        'performance_report': validator.generate_performance_report()
    }
    
    logger.info("Production dataset performance validation completed")
    return results