#!/usr/bin/env python3
"""
Error Handler and Resource Monitor for RNA Editing Annotation

This module provides comprehensive error handling, resource monitoring, and
graceful fallback mechanisms for the RNA editing annotation pipeline.

Features:
- Comprehensive error classification and handling
- Resource monitoring and optimization
- Graceful fallback mechanisms
- Detailed logging and diagnostics
- Temporary file cleanup
- Process timeout handling

Requirements Satisfied: 7.1, 7.2, 7.3, 7.4, 7.5

Author: RNA Editing Enhancement Pipeline
Date: 2025-12-15
"""

import os
import sys
import time
import signal
import logging
import tempfile
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
from contextlib import contextmanager

# Set up logging
logger = logging.getLogger(__name__)

class ResourceMonitor:
    """Monitor system resources during RNA editing annotation."""
    
    def __init__(self, max_memory_mb: int = 8000, max_disk_gb: int = 50):
        self.max_memory_mb = max_memory_mb
        self.max_disk_gb = max_disk_gb
        self.start_time = time.time()
        self.peak_memory_mb = 0
        self.disk_usage_gb = 0
        
        # Try to import psutil for advanced monitoring
        try:
            import psutil
            self.psutil_available = True
            self.process = psutil.Process()
        except ImportError:
            self.psutil_available = False
            logger.debug("psutil not available - using basic resource monitoring")
    
    def check_resources(self, context: str = "") -> Dict[str, Any]:
        """Check current resource usage and return statistics."""
        stats = {
            'timestamp': time.time(),
            'elapsed_time': time.time() - self.start_time,
            'context': context,
            'memory_mb': 0,
            'disk_free_gb': 0,
            'cpu_percent': 0,
            'warnings': []
        }
        
        try:
            if self.psutil_available:
                # Get memory usage
                memory_info = self.process.memory_info()
                stats['memory_mb'] = memory_info.rss / 1024 / 1024
                self.peak_memory_mb = max(self.peak_memory_mb, stats['memory_mb'])
                
                # Get CPU usage
                stats['cpu_percent'] = self.process.cpu_percent()
                
                # Check memory limits
                if stats['memory_mb'] > self.max_memory_mb:
                    warning = f"Memory usage ({stats['memory_mb']:.1f} MB) exceeds limit ({self.max_memory_mb} MB)"
                    stats['warnings'].append(warning)
                    logger.warning(warning)
            
            # Check disk space
            disk_usage = os.statvfs('.')
            free_bytes = disk_usage.f_frsize * disk_usage.f_bavail
            stats['disk_free_gb'] = free_bytes / (1024**3)
            
            if stats['disk_free_gb'] < 5:  # Less than 5GB free
                warning = f"Low disk space: {stats['disk_free_gb']:.1f} GB remaining"
                stats['warnings'].append(warning)
                logger.warning(warning)
                
        except Exception as e:
            logger.debug(f"Resource monitoring error: {e}")
            stats['warnings'].append(f"Resource monitoring failed: {e}")
        
        return stats
    
    def get_peak_usage(self) -> Dict[str, float]:
        """Get peak resource usage statistics."""
        return {
            'peak_memory_mb': self.peak_memory_mb,
            'total_time_seconds': time.time() - self.start_time
        }


class ErrorHandler:
    """Comprehensive error handler for RNA editing annotation."""
    
    def __init__(self, temp_dir: Optional[Path] = None):
        self.temp_dir = temp_dir or Path(tempfile.gettempdir())
        self.temp_files: List[Path] = []
        self.error_count = 0
        self.warning_count = 0
        self.cleanup_registered = False
        
        # Register cleanup on exit
        self._register_cleanup()
    
    def _register_cleanup(self):
        """Register cleanup handlers for various exit scenarios."""
        if not self.cleanup_registered:
            import atexit
            atexit.register(self.cleanup_temp_files)
            
            # Handle signals
            for sig in [signal.SIGTERM, signal.SIGINT]:
                try:
                    signal.signal(sig, self._signal_handler)
                except (OSError, ValueError):
                    # Some signals may not be available on all platforms
                    pass
            
            self.cleanup_registered = True
    
    def _signal_handler(self, signum, frame):
        """Handle process signals with cleanup."""
        logger.info(f"Received signal {signum}, performing cleanup...")
        self.cleanup_temp_files(force=True)
        sys.exit(128 + signum)
    
    def add_temp_file(self, file_path: Path) -> None:
        """Register a temporary file for cleanup."""
        if file_path not in self.temp_files:
            self.temp_files.append(file_path)
            logger.debug(f"Registered temporary file: {file_path}")
    
    def cleanup_temp_files(self, force: bool = False) -> Tuple[int, int]:
        """
        Clean up temporary files.
        
        Args:
            force: If True, attempt cleanup even if files are in use
            
        Returns:
            Tuple of (cleaned_count, failed_count)
        """
        if not self.temp_files:
            return 0, 0
        
        cleaned_count = 0
        failed_count = 0
        total_size_cleaned = 0
        
        logger.info(f"Cleaning up {len(self.temp_files)} temporary files...")
        
        for temp_file in self.temp_files[:]:  # Copy list to avoid modification during iteration
            try:
                if temp_file.exists():
                    file_size = temp_file.stat().st_size
                    
                    if force:
                        # Try to change permissions first
                        try:
                            temp_file.chmod(0o777)
                        except Exception:
                            pass
                    
                    temp_file.unlink()
                    total_size_cleaned += file_size
                    cleaned_count += 1
                    logger.debug(f"Cleaned up: {temp_file} ({file_size:,} bytes)")
                
                # Also clean up associated index files
                for suffix in ['.tbi', '.csi', '.idx']:
                    index_file = Path(str(temp_file) + suffix)
                    if index_file.exists():
                        try:
                            index_size = index_file.stat().st_size
                            if force:
                                index_file.chmod(0o777)
                            index_file.unlink()
                            total_size_cleaned += index_size
                            logger.debug(f"Cleaned up index: {index_file}")
                        except Exception as e:
                            logger.debug(f"Failed to clean up index {index_file}: {e}")
                
                # Remove from tracking list
                self.temp_files.remove(temp_file)
                
            except Exception as e:
                failed_count += 1
                logger.warning(f"Failed to clean up {temp_file}: {e}")
        
        if cleaned_count > 0:
            logger.info(f"✓ Cleaned up {cleaned_count} temporary files ({total_size_cleaned:,} bytes freed)")
        
        if failed_count > 0:
            logger.warning(f"Failed to clean up {failed_count} temporary files")
        
        return cleaned_count, failed_count
    
    def handle_error(self, error: Exception, context: str, 
                    fallback_action: Optional[callable] = None) -> bool:
        """
        Handle errors with appropriate classification and response.
        
        Args:
            error: The exception that occurred
            context: Description of what was being attempted
            fallback_action: Optional fallback function to execute
            
        Returns:
            True if error was handled gracefully, False if it should be re-raised
        """
        self.error_count += 1
        error_type = type(error).__name__
        error_msg = str(error)
        
        logger.error(f"ERROR in {context}: {error_type} - {error_msg}")
        
        # Classify error and determine response
        if isinstance(error, FileNotFoundError):
            logger.error("File not found - check input paths and permissions")
            if fallback_action:
                logger.info("Attempting fallback action...")
                try:
                    fallback_action()
                    logger.info("Fallback action completed successfully")
                    return True
                except Exception as fallback_error:
                    logger.error(f"Fallback action failed: {fallback_error}")
            return False
            
        elif isinstance(error, PermissionError):
            logger.error("Permission denied - check file and directory permissions")
            return False
            
        elif isinstance(error, subprocess.TimeoutExpired):
            logger.error(f"Process timed out after {error.timeout} seconds")
            logger.error("This may indicate very large files or insufficient resources")
            return False
            
        elif isinstance(error, subprocess.CalledProcessError):
            logger.error(f"Command failed with exit code {error.returncode}")
            if error.stderr:
                logger.error(f"Command stderr: {error.stderr}")
            
            # Check for specific error patterns
            if "not found" in error_msg.lower() or error.returncode == 127:
                logger.error("Required tool not found in PATH")
                return False
            elif error.returncode in [137, 143]:
                logger.error("Process was killed - likely out of memory or terminated")
                return False
            else:
                # Try fallback for other command failures
                if fallback_action:
                    logger.info("Attempting fallback action for command failure...")
                    try:
                        fallback_action()
                        return True
                    except Exception as fallback_error:
                        logger.error(f"Fallback action failed: {fallback_error}")
                return False
                
        elif isinstance(error, MemoryError):
            logger.error("Out of memory - try reducing input size or increasing memory allocation")
            return False
            
        else:
            logger.error(f"Unexpected error type: {error_type}")
            if fallback_action:
                logger.info("Attempting fallback action for unexpected error...")
                try:
                    fallback_action()
                    return True
                except Exception as fallback_error:
                    logger.error(f"Fallback action failed: {fallback_error}")
            return False
    
    def get_error_summary(self) -> Dict[str, int]:
        """Get summary of errors and warnings encountered."""
        return {
            'error_count': self.error_count,
            'warning_count': self.warning_count,
            'temp_files_remaining': len(self.temp_files)
        }


@contextmanager
def timeout_handler(seconds: int, description: str = "operation"):
    """Context manager for handling timeouts."""
    def timeout_signal_handler(signum, frame):
        raise TimeoutError(f"Operation '{description}' timed out after {seconds} seconds")
    
    # Set up timeout
    old_handler = signal.signal(signal.SIGALRM, timeout_signal_handler)
    signal.alarm(seconds)
    
    try:
        yield
    finally:
        # Clean up timeout
        signal.alarm(0)
        signal.signal(signal.SIGALRM, old_handler)


def create_graceful_fallback(input_vcf: Path, output_vcf: Path, 
                           reason: str = "annotation failure") -> bool:
    """
    Create graceful fallback by copying input to output.
    
    Args:
        input_vcf: Input VCF file path
        output_vcf: Output VCF file path
        reason: Reason for fallback
        
    Returns:
        True if fallback was successful, False otherwise
    """
    logger.info(f"Creating graceful fallback due to: {reason}")
    
    try:
        # Copy input VCF to output
        import shutil
        shutil.copy2(input_vcf, output_vcf)
        logger.info(f"Copied {input_vcf} to {output_vcf}")
        
        # Handle index file
        input_index = Path(str(input_vcf) + '.tbi')
        output_index = Path(str(output_vcf) + '.tbi')
        
        if input_index.exists():
            shutil.copy2(input_index, output_index)
            logger.info(f"Copied index file")
        else:
            # Create new index if needed
            if str(output_vcf).endswith('.gz'):
                try:
                    subprocess.run(['tabix', '-p', 'vcf', str(output_vcf)], 
                                 check=True, capture_output=True, text=True)
                    logger.info("Created new tabix index")
                except subprocess.CalledProcessError as e:
                    logger.warning(f"Could not create index: {e}")
        
        logger.info("Graceful fallback completed successfully")
        return True
        
    except Exception as e:
        logger.error(f"Graceful fallback failed: {e}")
        return False


def validate_tool_availability(required_tools: List[str]) -> Tuple[List[str], List[str]]:
    """
    Validate that required tools are available in PATH.
    
    Args:
        required_tools: List of tool names to check
        
    Returns:
        Tuple of (available_tools, missing_tools)
    """
    available_tools = []
    missing_tools = []
    
    for tool in required_tools:
        try:
            subprocess.run([tool, '--version'], 
                         capture_output=True, text=True, timeout=10)
            available_tools.append(tool)
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired, FileNotFoundError):
            missing_tools.append(tool)
    
    return available_tools, missing_tools


def check_file_accessibility(file_paths: List[Path]) -> Dict[str, str]:
    """
    Check file accessibility and return status for each file.
    
    Args:
        file_paths: List of file paths to check
        
    Returns:
        Dictionary mapping file paths to status messages
    """
    results = {}
    
    for file_path in file_paths:
        try:
            if not file_path.exists():
                results[str(file_path)] = "File does not exist"
            elif not os.access(file_path, os.R_OK):
                results[str(file_path)] = "File is not readable"
            elif file_path.stat().st_size == 0:
                results[str(file_path)] = "File is empty"
            else:
                results[str(file_path)] = "OK"
        except Exception as e:
            results[str(file_path)] = f"Error checking file: {e}"
    
    return results


# Example usage and testing
if __name__ == "__main__":
    # Set up logging for testing
    logging.basicConfig(level=logging.INFO, 
                       format='%(asctime)s - %(levelname)s - %(message)s')
    
    # Test resource monitoring
    monitor = ResourceMonitor()
    stats = monitor.check_resources("test")
    print(f"Resource stats: {stats}")
    
    # Test error handling
    error_handler = ErrorHandler()
    
    # Test tool availability
    tools = ['python', 'bcftools', 'tabix', 'nonexistent_tool']
    available, missing = validate_tool_availability(tools)
    print(f"Available tools: {available}")
    print(f"Missing tools: {missing}")
    
    # Test file accessibility
    test_files = [Path(__file__), Path('/nonexistent/file.txt')]
    file_status = check_file_accessibility(test_files)
    for file_path, status in file_status.items():
        print(f"{file_path}: {status}")