#!/usr/bin/env python3
"""
bcftools VCF-to-VCF Annotation Engine

This module provides a dedicated bcftools annotation engine for VCF-to-VCF annotation
without requiring separate header files or column specifications. It implements exact
coordinate and allele matching between input VCF and annotation VCF database.

FEATURES:
- Direct VCF-to-VCF annotation using bcftools annotate
- Automatic header conflict detection and resolution
- Exact coordinate and allele matching (CHROM, POS, REF, ALT)
- Compressed VCF input/output handling with proper indexing
- Comprehensive error handling and logging

Requirements Satisfied: 3.1, 3.2, 3.3, 3.5

Author: RNA Editing Enhancement Pipeline
Date: 2025-12-13
"""

import logging
import shutil
import subprocess
import os
import time
from pathlib import Path
from typing import Dict, Any, List, Optional, Tuple

logger = logging.getLogger(__name__)

# Enhanced logging configuration for bcftools operations
class BcftoolsOperationLogger:
    """Enhanced logging for bcftools annotation operations."""
    
    def __init__(self):
        self.stats = {
            'start_time': time.time(),
            'command_executions': 0,
            'successful_commands': 0,
            'failed_commands': 0,
            'total_stdout_lines': 0,
            'total_stderr_lines': 0,
            'command_history': [],
            'error_details': [],
            'performance_metrics': {}
        }
    
    def log_command_start(self, cmd: List[str], description: str = ""):
        """Log the start of a bcftools command execution."""
        cmd_info = {
            'command': ' '.join(str(x) for x in cmd),
            'description': description,
            'start_time': time.time(),
            'pid': None
        }
        
        self.stats['command_executions'] += 1
        self.stats['command_history'].append(cmd_info)
        
        logger.info(f"Executing bcftools command: {description}")
        logger.debug(f"Command: {cmd_info['command']}")
        
        return len(self.stats['command_history']) - 1  # Return index for tracking
    
    def log_command_completion(self, cmd_index: int, returncode: int, stdout: str, stderr: str, 
                             execution_time: float):
        """Log the completion of a bcftools command execution."""
        if cmd_index < len(self.stats['command_history']):
            cmd_info = self.stats['command_history'][cmd_index]
            cmd_info['returncode'] = returncode
            cmd_info['execution_time'] = execution_time
            cmd_info['stdout_lines'] = len(stdout.split('\n')) if stdout else 0
            cmd_info['stderr_lines'] = len(stderr.split('\n')) if stderr else 0
            
            self.stats['total_stdout_lines'] += cmd_info['stdout_lines']
            self.stats['total_stderr_lines'] += cmd_info['stderr_lines']
            
            if returncode == 0:
                self.stats['successful_commands'] += 1
                logger.info(f"✓ Command completed successfully in {execution_time:.2f}s")
                
                # Log informational stderr (bcftools often writes progress to stderr)
                if stderr:
                    self._log_bcftools_stderr(stderr, is_error=False)
                
                # Log stdout if it contains useful information
                if stdout and len(stdout.strip()) > 0:
                    logger.debug(f"Command stdout ({cmd_info['stdout_lines']} lines):")
                    for line in stdout.strip().split('\n')[:10]:  # Log first 10 lines
                        logger.debug(f"  {line}")
                    if cmd_info['stdout_lines'] > 10:
                        logger.debug(f"  ... ({cmd_info['stdout_lines'] - 10} more lines)")
            else:
                self.stats['failed_commands'] += 1
                error_info = {
                    'command': cmd_info['command'],
                    'returncode': returncode,
                    'stderr': stderr,
                    'stdout': stdout,
                    'execution_time': execution_time,
                    'timestamp': time.time()
                }
                self.stats['error_details'].append(error_info)
                
                logger.error(f"✗ Command failed with return code {returncode} after {execution_time:.2f}s")
                logger.error(f"Command: {cmd_info['command']}")
                
                # Log detailed error information
                if stderr:
                    self._log_bcftools_stderr(stderr, is_error=True)
                if stdout:
                    logger.error(f"Command stdout: {stdout}")
    
    def _log_bcftools_stderr(self, stderr: str, is_error: bool = False):
        """Log bcftools stderr output with appropriate level."""
        if not stderr:
            return
        
        stderr_lines = stderr.strip().split('\n')
        
        for line in stderr_lines:
            line = line.strip()
            if not line:
                continue
            
            # Classify stderr messages
            if any(keyword in line.lower() for keyword in ['error', 'failed', 'cannot', 'invalid']):
                logger.error(f"bcftools error: {line}")
            elif any(keyword in line.lower() for keyword in ['warning', 'warn']):
                logger.warning(f"bcftools warning: {line}")
            elif any(keyword in line.lower() for keyword in ['processed', 'written', 'lines', 'records']):
                logger.info(f"bcftools progress: {line}")
            elif is_error:
                logger.error(f"bcftools stderr: {line}")
            else:
                logger.debug(f"bcftools info: {line}")
    
    def log_file_operation(self, operation: str, file_path: Path, success: bool, 
                          details: str = "", execution_time: float = 0):
        """Log file operations (compression, indexing, etc.)."""
        if success:
            logger.info(f"✓ {operation} completed for {file_path.name}")
            if execution_time > 0:
                logger.info(f"  Time: {execution_time:.2f}s")
            if details:
                logger.info(f"  Details: {details}")
        else:
            logger.error(f"✗ {operation} failed for {file_path.name}")
            if details:
                logger.error(f"  Error: {details}")
    
    def log_performance_metrics(self):
        """Log performance metrics for bcftools operations."""
        total_time = time.time() - self.stats['start_time']
        
        logger.info("=== bcftools Performance Metrics ===")
        logger.info(f"Total execution time: {total_time:.2f} seconds")
        logger.info(f"Commands executed: {self.stats['command_executions']}")
        logger.info(f"Successful commands: {self.stats['successful_commands']}")
        logger.info(f"Failed commands: {self.stats['failed_commands']}")
        
        if self.stats['command_executions'] > 0:
            success_rate = (self.stats['successful_commands'] / self.stats['command_executions']) * 100
            logger.info(f"Success rate: {success_rate:.1f}%")
        
        # Log command timing breakdown
        if self.stats['command_history']:
            logger.info("Command execution times:")
            for i, cmd_info in enumerate(self.stats['command_history']):
                if 'execution_time' in cmd_info:
                    status = "✓" if cmd_info.get('returncode', -1) == 0 else "✗"
                    logger.info(f"  {status} {cmd_info.get('description', f'Command {i+1}')}: {cmd_info['execution_time']:.2f}s")
        
        # Log output statistics
        logger.info(f"Total stdout lines captured: {self.stats['total_stdout_lines']:,}")
        logger.info(f"Total stderr lines captured: {self.stats['total_stderr_lines']:,}")
        
        # Log error summary
        if self.stats['error_details']:
            logger.warning(f"Errors encountered: {len(self.stats['error_details'])}")
            for error in self.stats['error_details']:
                logger.warning(f"  Command failed: {error['command'][:100]}...")
                logger.warning(f"    Return code: {error['returncode']}")
                logger.warning(f"    Time: {error['execution_time']:.2f}s")
    
    def get_statistics(self) -> Dict:
        """Get operation statistics dictionary."""
        return self.stats.copy()


class BcftoolsAnnotator:
    """
    bcftools VCF-to-VCF annotation engine with exact matching and conflict resolution.
    
    This class provides streamlined VCF annotation using bcftools annotate for exact
    coordinate and allele matching between input VCF and annotation VCF database.
    It handles header conflicts gracefully and supports compressed VCF files.
    """
    
    def __init__(self, input_vcf: str, annotation_vcf: str, output_vcf: str):
        """
        Initialize bcftools annotator.
        
        Args:
            input_vcf: Path to input VCF file to be annotated
            annotation_vcf: Path to annotation VCF database
            output_vcf: Path to output annotated VCF file
        """
        self.input_vcf = Path(input_vcf)
        self.annotation_vcf = Path(annotation_vcf)
        self.output_vcf = Path(output_vcf)
        
        # Tool paths discovered via system PATH
        self.tool_paths: Dict[str, Optional[str]] = {}
        
        # Enhanced logging for bcftools operations
        self.operation_logger = BcftoolsOperationLogger()
        
        # Statistics and timing
        self.stats = {
            'start_time': time.time(),
            'processing_steps': [],
            'errors': [],
            'warnings': []
        }
        
        logger.info("=== bcftools VCF-to-VCF Annotation Engine Initialized ===")
        logger.info(f"Input VCF: {self.input_vcf}")
        logger.info(f"Annotation VCF: {self.annotation_vcf}")
        logger.info(f"Output VCF: {self.output_vcf}")
        
        # Log file sizes for diagnostics
        try:
            if self.input_vcf.exists():
                input_size = self.input_vcf.stat().st_size
                logger.info(f"Input VCF size: {input_size:,} bytes ({input_size/1024/1024:.1f} MB)")
            
            if self.annotation_vcf.exists():
                annotation_size = self.annotation_vcf.stat().st_size
                logger.info(f"Annotation VCF size: {annotation_size:,} bytes ({annotation_size/1024/1024:.1f} MB)")
        except Exception as e:
            logger.debug(f"Could not get file sizes: {e}")
        
        # Validate inputs and tools
        self.validate_inputs()
        self.validate_tools()
    
    def validate_inputs(self) -> None:
        """Validate input files exist and are accessible."""
        step_start = time.time()
        logger.info("Validating input files...")
        
        # Check input VCF existence and accessibility
        if not self.input_vcf.exists():
            error_msg = f"Input VCF not found: {self.input_vcf}"
            logger.error(error_msg)
            raise FileNotFoundError(error_msg)
        
        if not self.annotation_vcf.exists():
            error_msg = f"Annotation VCF not found: {self.annotation_vcf}"
            logger.error(error_msg)
            raise FileNotFoundError(error_msg)
        
        # Check file permissions
        if not os.access(self.input_vcf, os.R_OK):
            error_msg = f"Cannot read input VCF: {self.input_vcf}"
            logger.error(error_msg)
            raise PermissionError(error_msg)
            
        if not os.access(self.annotation_vcf, os.R_OK):
            error_msg = f"Cannot read annotation VCF: {self.annotation_vcf}"
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
        annotation_size = self.annotation_vcf.stat().st_size
        logger.info(f"Input VCF size: {input_size:,} bytes ({input_size/1024/1024:.1f} MB)")
        logger.info(f"Annotation VCF size: {annotation_size:,} bytes ({annotation_size/1024/1024:.1f} MB)")
        
        step_time = time.time() - step_start
        self.stats['processing_steps'].append(('validate_inputs', step_time))
        logger.info(f"✓ Input files validated successfully ({step_time:.2f}s)")
    
    def validate_tools(self) -> None:
        """Validate required tools are available in system PATH."""
        step_start = time.time()
        logger.info("Validating required tools in system PATH...")
        
        required_tools = ['bcftools', 'tabix', 'bgzip']
        missing_tools = []
        tool_versions = {}
        
        for tool in required_tools:
            tool_path = shutil.which(tool)
            if tool_path:
                self.tool_paths[tool] = tool_path
                
                # Get tool version for diagnostics
                try:
                    version_result = subprocess.run(
                        [tool_path, '--version'],
                        capture_output=True,
                        text=True,
                        timeout=10
                    )
                    if version_result.returncode == 0:
                        version_line = version_result.stdout.split('\n')[0]
                        tool_versions[tool] = version_line
                        logger.info(f"✓ Found {tool}: {tool_path} ({version_line})")
                    else:
                        logger.warning(f"✓ Found {tool}: {tool_path} (version check failed)")
                        tool_versions[tool] = "version unknown"
                except Exception as e:
                    logger.warning(f"✓ Found {tool}: {tool_path} (version check error: {e})")
                    tool_versions[tool] = f"version error: {e}"
            else:
                self.tool_paths[tool] = None
                missing_tools.append(tool)
                logger.error(f"✗ {tool} not found in system PATH")
        
        if missing_tools:
            error_msg = f"Required tools not found in system PATH: {missing_tools}"
            logger.error(error_msg)
            logger.error("Installation instructions:")
            logger.error("  Ubuntu/Debian: apt install bcftools tabix")
            logger.error("  CentOS/RHEL: yum install bcftools htslib")
            logger.error("  Conda: conda install bcftools htslib")
            logger.error("  Homebrew: brew install bcftools htslib")
            raise RuntimeError(error_msg)
        
        # Store tool versions in stats
        self.stats['tool_versions'] = tool_versions
        
        step_time = time.time() - step_start
        self.stats['processing_steps'].append(('validate_tools', step_time))
        logger.info(f"✓ All required tools found in system PATH ({step_time:.2f}s)")
    
    def check_header_compatibility(self) -> Dict[str, Any]:
        """
        Check VCF header for existing fields and determine annotation strategy.
        
        Returns:
            Dict containing header analysis results and annotation strategy
        """
        step_start = time.time()
        logger.info("Checking VCF header compatibility...")
        
        try:
            # Read VCF header to check for existing INFO fields
            cmd = [str(self.tool_paths['bcftools']), 'view', '-h', str(self.input_vcf)]
            logger.debug(f"Running command: {' '.join(cmd)}")
            
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True,
                timeout=60
            )
            
            header_lines = result.stdout.strip().split('\n')
            existing_info_fields = set()
            header_stats = {
                'total_lines': len(header_lines),
                'info_lines': 0,
                'format_lines': 0,
                'contig_lines': 0,
                'sample_count': 0
            }
            
            for line in header_lines:
                if line.startswith('##INFO=<ID='):
                    header_stats['info_lines'] += 1
                    # Extract field ID
                    try:
                        start = line.find('ID=') + 3
                        end = line.find(',', start)
                        if end == -1:
                            end = line.find('>', start)
                        field_id = line[start:end]
                        existing_info_fields.add(field_id)
                    except Exception as e:
                        logger.warning(f"Failed to parse INFO field from line: {line[:100]}... Error: {e}")
                elif line.startswith('##FORMAT='):
                    header_stats['format_lines'] += 1
                elif line.startswith('##contig='):
                    header_stats['contig_lines'] += 1
                elif line.startswith('#CHROM'):
                    # Count samples
                    columns = line.split('\t')
                    if len(columns) > 9:
                        header_stats['sample_count'] = len(columns) - 9
            
            logger.info(f"VCF header analysis: {header_stats['total_lines']} lines, "
                       f"{header_stats['info_lines']} INFO fields, "
                       f"{header_stats['sample_count']} samples")
            
            # Check for potential conflicts with annotation fields
            annotation_fields = self._get_annotation_fields()
            conflicts = existing_info_fields.intersection(annotation_fields)
            
            # Determine annotation strategy based on conflicts
            if conflicts:
                logger.warning(f"Found existing INFO fields that may conflict: {conflicts}")
                strategy = 'handle_conflicts'
                logger.info("Will use conflict-aware annotation strategy")
                self.stats['warnings'].append(f"Header conflicts detected: {conflicts}")
            else:
                logger.info("✓ No header conflicts detected")
                strategy = 'standard'
            
            analysis_result = {
                'existing_fields': existing_info_fields,
                'conflicting_fields': conflicts,
                'strategy': strategy,
                'requires_preprocessing': len(conflicts) > 0,
                'header_stats': header_stats
            }
            
            step_time = time.time() - step_start
            self.stats['processing_steps'].append(('check_header_compatibility', step_time))
            logger.info(f"✓ Header compatibility check completed ({step_time:.2f}s)")
            
            return analysis_result
                
        except subprocess.CalledProcessError as e:
            error_msg = f"Failed to check VCF header: {e}"
            logger.error(error_msg)
            logger.error(f"Command: {' '.join(cmd)}")
            logger.error(f"Return code: {e.returncode}")
            logger.error(f"stdout: {e.stdout}")
            logger.error(f"stderr: {e.stderr}")
            raise RuntimeError(error_msg)
        except subprocess.TimeoutExpired:
            error_msg = "VCF header check timed out after 60 seconds"
            logger.error(error_msg)
            raise RuntimeError(error_msg)
        except Exception as e:
            logger.error(f"Unexpected error during header compatibility check: {e}")
            raise
    
    def _get_annotation_fields(self) -> set:
        """
        Get annotation fields from annotation VCF header.
        
        Returns:
            Set of INFO field IDs from annotation VCF
        """
        try:
            cmd = [str(self.tool_paths['bcftools']), 'view', '-h', str(self.annotation_vcf)]
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True,
                timeout=60
            )
            
            annotation_fields = set()
            for line in result.stdout.split('\n'):
                if line.startswith('##INFO=<ID='):
                    try:
                        start = line.find('ID=') + 3
                        end = line.find(',', start)
                        if end == -1:
                            end = line.find('>', start)
                        field_id = line[start:end]
                        annotation_fields.add(field_id)
                    except Exception as e:
                        logger.warning(f"Failed to parse annotation field from line: {line[:100]}... Error: {e}")
            
            logger.debug(f"Found annotation fields: {annotation_fields}")
            return annotation_fields
            
        except Exception as e:
            logger.warning(f"Could not determine annotation fields: {e}")
            # Return common REDIportal fields as fallback
            return {'REDI_ACCESSION', 'REDI_DB', 'REDI_TYPE', 'REDI_REPEAT', 'REDI_FUNC', 'REDI_STRAND', 'DB'}
    
    def build_annotation_command(self) -> Tuple[List[str], Optional[Path]]:
        """
        Build bcftools annotate command using VCF-to-VCF annotation approach.
        
        Returns:
            Tuple of (command arguments, temporary header file path if created)
        """
        logger.info("Building bcftools annotate command for VCF-to-VCF annotation...")
        
        # Check header compatibility first
        header_analysis = self.check_header_compatibility()
        
        if header_analysis['strategy'] == 'handle_conflicts':
            return self._build_conflict_aware_command(header_analysis)
        else:
            return self._build_standard_command()
    
    def _build_standard_command(self) -> Tuple[List[str], Optional[Path]]:
        """Build standard bcftools annotate command for VCF-to-VCF annotation."""
        logger.info("Building standard VCF-to-VCF annotation command...")
        
        # Use VCF-to-VCF annotation without requiring separate header files
        # bcftools will automatically match INFO fields from annotation VCF
        cmd = [
            str(self.tool_paths['bcftools']), 'annotate',
            '-a', str(self.annotation_vcf),
            '-c', 'CHROM,POS,REF,ALT,INFO',  # Copy all INFO fields from annotation VCF
            '-o', str(self.output_vcf),
            str(self.input_vcf)
        ]
        
        # Add compression if output is .gz
        if str(self.output_vcf).endswith('.gz'):
            cmd.extend(['-O', 'z'])
        
        logger.info(f"Standard command: {' '.join(cmd)}")
        return cmd, None
    
    def _build_conflict_aware_command(self, header_analysis: Dict[str, Any]) -> Tuple[List[str], Optional[Path]]:
        """
        Build bcftools annotate command that handles header conflicts.
        
        Args:
            header_analysis: Results from header compatibility check
            
        Returns:
            Tuple of (command arguments, temporary header file path if created)
        """
        conflicting_fields = header_analysis['conflicting_fields']
        logger.info(f"Building conflict-aware command for fields: {conflicting_fields}")
        
        # Get available annotation fields
        annotation_fields = self._get_annotation_fields()
        safe_fields = annotation_fields - conflicting_fields
        
        if safe_fields:
            # Annotate only non-conflicting fields
            columns = 'CHROM,POS,REF,ALT,' + ','.join(f'INFO/{field}' for field in safe_fields)
            
            cmd = [
                str(self.tool_paths['bcftools']), 'annotate',
                '-a', str(self.annotation_vcf),
                '-c', columns,
                '-o', str(self.output_vcf),
                str(self.input_vcf)
            ]
            
            logger.info(f"Selective annotation for fields: {safe_fields}")
        else:
            # If all fields conflict, just mark sites without adding INFO fields
            cmd = [
                str(self.tool_paths['bcftools']), 'annotate',
                '-a', str(self.annotation_vcf),
                '--mark-sites', '+ANNOTATION_MATCH',
                '-o', str(self.output_vcf),
                str(self.input_vcf)
            ]
            
            logger.warning("All annotation fields conflict - using site marking only")
        
        # Add compression if output is .gz
        if str(self.output_vcf).endswith('.gz'):
            cmd.extend(['-O', 'z'])
        
        logger.info(f"Conflict-aware command: {' '.join(cmd)}")
        return cmd, None
    
    def execute_annotation(self) -> None:
        """
        Execute bcftools annotate command with comprehensive error handling.
        """
        step_start = time.time()
        logger.info("Executing bcftools VCF-to-VCF annotation...")
        
        # Build annotation command
        cmd, temp_header_file = self.build_annotation_command()
        
        # Log command execution start
        cmd_index = self.operation_logger.log_command_start(cmd, "VCF-to-VCF annotation")
        
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
            
            # Log command completion
            self.operation_logger.log_command_completion(
                cmd_index, result.returncode, result.stdout, result.stderr, execution_time
            )
            
            # Parse and log bcftools output statistics
            self._parse_bcftools_output_statistics(result.stdout, result.stderr)
            
            # Clean up temporary header file if created
            if temp_header_file and temp_header_file.exists():
                temp_header_file.unlink()
                logger.debug(f"Cleaned up temporary header file: {temp_header_file}")
            
            self.stats['processing_steps'].append(('execute_annotation', execution_time))
            
            logger.info(f"✓ bcftools annotation completed successfully ({execution_time:.2f}s)")
            
        except subprocess.CalledProcessError as e:
            execution_time = time.time() - step_start
            
            # Log command failure
            self.operation_logger.log_command_completion(
                cmd_index, e.returncode, e.stdout, e.stderr, execution_time
            )
            
            error_context = f"bcftools annotate failed with exit code {e.returncode}"
            logger.error(error_context)
            logger.error(f"Command: {' '.join(str(x) for x in cmd)}")
            logger.error(f"Execution time: {execution_time:.2f}s")
            
            # Log detailed output
            if e.stdout:
                logger.error(f"stdout ({len(e.stdout.split())} lines):")
                for i, line in enumerate(e.stdout.split('\n')[:20]):  # First 20 lines
                    logger.error(f"  {i+1}: {line}")
                if len(e.stdout.split('\n')) > 20:
                    logger.error(f"  ... ({len(e.stdout.split('\n')) - 20} more lines)")
            
            if e.stderr:
                logger.error(f"stderr ({len(e.stderr.split())} lines):")
                for i, line in enumerate(e.stderr.split('\n')[:20]):  # First 20 lines
                    logger.error(f"  {i+1}: {line}")
                if len(e.stderr.split('\n')) > 20:
                    logger.error(f"  ... ({len(e.stderr.split('\n')) - 20} more lines)")
            
            # Provide specific error guidance based on stderr content
            self._diagnose_bcftools_error(e.stderr)
            
            # Clean up temporary files on error
            if temp_header_file and temp_header_file.exists():
                temp_header_file.unlink()
                logger.debug("Cleaned up temporary header file after error")
            
            raise RuntimeError(f"bcftools annotate failed: {e}")
            
        except subprocess.TimeoutExpired:
            execution_time = time.time() - step_start
            error_msg = f"bcftools annotate timed out after {execution_time:.0f} seconds (limit: 3600s)"
            
            logger.error(error_msg)
            logger.error("This may indicate:")
            logger.error("  - Very large input files requiring more processing time")
            logger.error("  - System performance issues or resource constraints")
            logger.error("  - Possible infinite loop or deadlock in bcftools")
            logger.error("  - Insufficient memory causing excessive swapping")
            
            # Log system resource information if available
            try:
                import psutil
                memory = psutil.virtual_memory()
                disk = psutil.disk_usage('/')
                logger.error(f"System memory: {memory.percent}% used ({memory.available/1024/1024/1024:.1f} GB available)")
                logger.error(f"Disk space: {disk.percent}% used ({disk.free/1024/1024/1024:.1f} GB available)")
            except ImportError:
                logger.debug("psutil not available for system resource monitoring")
            except Exception as e:
                logger.debug(f"Could not get system resource info: {e}")
            
            # Clean up temporary files on timeout
            if temp_header_file and temp_header_file.exists():
                temp_header_file.unlink()
                logger.debug("Cleaned up temporary header file after timeout")
            
            raise RuntimeError(error_msg)
            
        except Exception as e:
            execution_time = time.time() - step_start
            logger.error(f"Unexpected error during bcftools execution after {execution_time:.2f}s: {e}")
            logger.error(f"Error type: {type(e).__name__}")
            
            # Clean up temporary files on error
            if temp_header_file and temp_header_file.exists():
                temp_header_file.unlink()
                logger.debug("Cleaned up temporary header file after unexpected error")
            
            raise
    
    def _parse_bcftools_output_statistics(self, stdout: str, stderr: str):
        """Parse and log statistics from bcftools output."""
        try:
            # Look for common bcftools statistics patterns
            combined_output = (stdout or "") + "\n" + (stderr or "")
            
            # Parse lines processed
            import re
            lines_match = re.search(r'(\d+)\s+lines?\s+processed', combined_output, re.IGNORECASE)
            if lines_match:
                lines_processed = int(lines_match.group(1))
                logger.info(f"bcftools processed {lines_processed:,} lines")
            
            # Parse records written
            records_match = re.search(r'(\d+)\s+records?\s+written', combined_output, re.IGNORECASE)
            if records_match:
                records_written = int(records_match.group(1))
                logger.info(f"bcftools wrote {records_written:,} records")
            
            # Parse annotations added
            annotations_match = re.search(r'(\d+)\s+annotations?\s+added', combined_output, re.IGNORECASE)
            if annotations_match:
                annotations_added = int(annotations_match.group(1))
                logger.info(f"bcftools added {annotations_added:,} annotations")
            
        except Exception as e:
            logger.debug(f"Could not parse bcftools output statistics: {e}")
    
    def _diagnose_bcftools_error(self, stderr: str) -> None:
        """
        Provide specific error guidance based on bcftools stderr content.
        
        Args:
            stderr: Standard error output from bcftools command
        """
        if not stderr:
            return
        
        stderr_lower = stderr.lower()
        
        if "header already exists" in stderr_lower:
            logger.error("DIAGNOSIS: VCF header conflict detected")
            logger.error("SOLUTION: The input VCF already contains annotation fields that conflict")
            logger.error("  - Check if the VCF was previously annotated")
            logger.error("  - Consider removing existing annotation fields first")
            logger.error("  - Use bcftools annotate -x to remove conflicting fields")
            
        elif "no such file" in stderr_lower or "not found" in stderr_lower:
            logger.error("DIAGNOSIS: File not found error")
            logger.error("SOLUTION: Check that all input files exist and are accessible")
            logger.error(f"  - Input VCF: {self.input_vcf} (exists: {self.input_vcf.exists()})")
            logger.error(f"  - Annotation VCF: {self.annotation_vcf} (exists: {self.annotation_vcf.exists()})")
            
        elif "permission denied" in stderr_lower:
            logger.error("DIAGNOSIS: Permission denied error")
            logger.error("SOLUTION: Check file permissions and output directory access")
            logger.error(f"  - Output directory: {self.output_vcf.parent}")
            logger.error(f"  - Directory writable: {os.access(self.output_vcf.parent, os.W_OK)}")
            
        elif "not compressed with bgzip" in stderr_lower or "not indexed" in stderr_lower:
            logger.error("DIAGNOSIS: VCF file format issue")
            logger.error("SOLUTION: Files should be bgzip compressed and tabix indexed for bcftools")
            logger.error("  - Use 'bgzip file.vcf' to compress")
            logger.error("  - Use 'tabix -p vcf file.vcf.gz' to index")
            
        elif "malformed" in stderr_lower or "invalid" in stderr_lower:
            logger.error("DIAGNOSIS: VCF format validation error")
            logger.error("SOLUTION: Input VCF may have format issues")
            logger.error("  - Use 'bcftools view -h' to check header format")
            logger.error("  - Use 'bcftools stats' to validate VCF structure")
            
        elif "out of memory" in stderr_lower or "memory" in stderr_lower:
            logger.error("DIAGNOSIS: Memory allocation error")
            logger.error("SOLUTION: Insufficient memory for processing")
            logger.error("  - Try processing smaller VCF chunks")
            logger.error("  - Increase available system memory")
            
        else:
            logger.error("DIAGNOSIS: Unknown bcftools error")
            logger.error("SOLUTION: Check bcftools documentation and input file formats")
    
    def validate_output(self) -> None:
        """Validate that output file was created successfully and is accessible."""
        step_start = time.time()
        logger.info("Validating output file creation and accessibility...")
        
        try:
            # Check file existence
            if not self.output_vcf.exists():
                raise RuntimeError(f"Output file was not created: {self.output_vcf}")
            
            # Check file size
            file_size = self.output_vcf.stat().st_size
            if file_size == 0:
                raise RuntimeError(f"Output file is empty: {self.output_vcf}")
            
            # Check file permissions
            if not os.access(self.output_vcf, os.R_OK):
                raise RuntimeError(f"Output file is not readable: {self.output_vcf}")
            
            # Validate VCF format
            try:
                cmd = [str(self.tool_paths['bcftools']), 'view', '-h', str(self.output_vcf)]
                result = subprocess.run(cmd, capture_output=True, text=True, check=True, timeout=60)
                
                # Check for VCF header
                if not result.stdout.startswith('##fileformat=VCF'):
                    raise RuntimeError("Output file does not have valid VCF header")
                
                logger.debug("✓ Output VCF format validation passed")
                    
            except subprocess.CalledProcessError as e:
                logger.warning(f"VCF format validation failed: {e}")
                logger.warning("Output file may have format issues")
            except Exception as e:
                logger.warning(f"Could not validate VCF format: {e}")
            
            # Create index if output is compressed
            if str(self.output_vcf).endswith('.gz'):
                self._create_output_index()
            
            step_time = time.time() - step_start
            self.stats['processing_steps'].append(('validate_output', step_time))
            
            logger.info(f"✓ Output file validation completed ({file_size:,} bytes, {step_time:.2f}s)")
            
        except Exception as e:
            logger.error(f"Output file validation failed: {e}")
            raise
    
    def _create_output_index(self) -> None:
        """Create tabix index for compressed output VCF."""
        logger.info("Creating tabix index for compressed output...")
        
        try:
            index_file = Path(str(self.output_vcf) + '.tbi')
            
            # Check if index already exists and is valid
            if index_file.exists():
                vcf_mtime = self.output_vcf.stat().st_mtime
                index_mtime = index_file.stat().st_mtime
                
                if index_mtime >= vcf_mtime:
                    logger.info("Valid tabix index already exists")
                    return
                else:
                    logger.info("Existing index is outdated, recreating...")
                    index_file.unlink()
            
            cmd = [str(self.tool_paths['tabix']), '-p', 'vcf', str(self.output_vcf)]
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
            logger.info(f"✓ Tabix index created successfully ({index_size:,} bytes)")
            
        except subprocess.CalledProcessError as e:
            error_msg = f"Failed to create tabix index: {e}"
            logger.error(error_msg)
            logger.error(f"Command stderr: {e.stderr}")
            raise RuntimeError(error_msg)
            
        except subprocess.TimeoutExpired:
            error_msg = "Tabix indexing timed out after 5 minutes"
            logger.error(error_msg)
            raise RuntimeError(error_msg)
            
        except Exception as e:
            logger.error(f"Unexpected error during indexing: {e}")
            raise
    
    def run_annotation(self) -> None:
        """
        Run the complete bcftools VCF-to-VCF annotation pipeline.
        """
        pipeline_start = time.time()
        logger.info("Starting bcftools VCF-to-VCF annotation pipeline...")
        
        temp_files = []
        
        try:
            # Prepare input files (compress if needed)
            logger.info("Step 1: Preparing input files...")
            prepared_input = self._prepare_input_vcf()
            if prepared_input != self.input_vcf:
                temp_files.append(prepared_input)
            
            logger.info("Step 2: Preparing annotation VCF...")
            prepared_annotation = self._prepare_annotation_vcf()
            if prepared_annotation != self.annotation_vcf:
                temp_files.append(prepared_annotation)
            
            # Update paths for annotation
            original_input = self.input_vcf
            original_annotation = self.annotation_vcf
            self.input_vcf = prepared_input
            self.annotation_vcf = prepared_annotation
            
            # Execute annotation
            logger.info("Step 3: Executing annotation...")
            self.execute_annotation()
            
            # Restore original paths
            self.input_vcf = original_input
            self.annotation_vcf = original_annotation
            
            # Validate output
            logger.info("Step 4: Validating output...")
            self.validate_output()
            
            # Clean up temporary files
            logger.info("Step 5: Cleaning up temporary files...")
            self._cleanup_temp_files(temp_files)
            
            # Generate final statistics
            pipeline_time = time.time() - pipeline_start
            self.stats['processing_steps'].append(('total_pipeline', pipeline_time))
            
            logger.info("✓ bcftools VCF-to-VCF annotation completed successfully!")
            self._generate_summary_statistics()
            
        except Exception as e:
            logger.error(f"bcftools annotation pipeline failed: {e}")
            
            # Clean up any partial output files
            if self.output_vcf.exists():
                try:
                    self.output_vcf.unlink()
                    logger.info("Cleaned up partial output file")
                except Exception as cleanup_error:
                    logger.warning(f"Failed to clean up partial output: {cleanup_error}")
            
            # Clean up temporary files
            try:
                self._cleanup_temp_files(temp_files)
            except Exception as cleanup_error:
                logger.warning(f"Failed to clean up temporary files: {cleanup_error}")
            
            # Generate error summary
            self._generate_summary_statistics()
            raise
    
    def _prepare_input_vcf(self) -> Path:
        """
        Prepare input VCF file for bcftools (compress if needed).
        
        Returns:
            Path to prepared (possibly compressed) input VCF
        """
        step_start = time.time()
        
        if str(self.input_vcf).endswith('.gz'):
            logger.info("Input VCF is already compressed")
            # Check if index exists
            index_file = Path(str(self.input_vcf) + '.tbi')
            if not index_file.exists():
                logger.info("Creating missing tabix index for input VCF...")
                try:
                    subprocess.run(
                        [str(self.tool_paths['tabix']), '-p', 'vcf', str(self.input_vcf)],
                        capture_output=True,
                        text=True,
                        check=True,
                        timeout=300
                    )
                    logger.info("✓ Tabix index created")
                except Exception as e:
                    logger.error(f"Failed to create tabix index for input VCF: {e}")
                    raise
            return self.input_vcf
        
        logger.info("Compressing input VCF for bcftools compatibility...")
        compressed_input = Path(str(self.input_vcf) + '.gz')
        
        try:
            # Compress with bgzip
            logger.info(f"Compressing {self.input_vcf} to {compressed_input}...")
            cmd = [str(self.tool_paths['bgzip']), '-c', str(self.input_vcf)]
            
            with open(compressed_input, 'wb') as f:
                subprocess.run(
                    cmd, 
                    stdout=f, 
                    stderr=subprocess.PIPE, 
                    text=True,
                    check=True,
                    timeout=1800  # 30 minute timeout for compression
                )
            
            # Verify compression worked
            if not compressed_input.exists() or compressed_input.stat().st_size == 0:
                raise RuntimeError("Compression failed - output file is missing or empty")
            
            logger.info(f"✓ Compression completed, size: {compressed_input.stat().st_size:,} bytes")
            
            # Index with tabix
            logger.info("Creating tabix index...")
            subprocess.run(
                [str(self.tool_paths['tabix']), '-p', 'vcf', str(compressed_input)],
                capture_output=True,
                text=True,
                check=True,
                timeout=300
            )
            
            # Verify index was created
            index_file = Path(str(compressed_input) + '.tbi')
            if not index_file.exists():
                raise RuntimeError("Tabix indexing failed - index file not created")
            
            logger.info(f"✓ Input VCF compressed and indexed: {compressed_input}")
            
            step_time = time.time() - step_start
            self.stats['processing_steps'].append(('prepare_input_vcf', step_time))
            return compressed_input
            
        except subprocess.CalledProcessError as e:
            error_msg = f"Failed to prepare input VCF: {e}"
            logger.error(error_msg)
            logger.error(f"Command stderr: {e.stderr}")
            
            # Clean up partial files
            if compressed_input.exists():
                try:
                    compressed_input.unlink()
                    logger.info("Cleaned up partial compressed file")
                except Exception:
                    pass
            
            raise RuntimeError(error_msg)
        except Exception as e:
            logger.error(f"Unexpected error during input VCF preparation: {e}")
            raise
    
    def _prepare_annotation_vcf(self) -> Path:
        """
        Prepare annotation VCF file for bcftools (compress if needed).
        
        Returns:
            Path to prepared (possibly compressed) annotation VCF
        """
        step_start = time.time()
        
        if str(self.annotation_vcf).endswith('.gz'):
            logger.info("Annotation VCF is already compressed")
            # Check if index exists
            index_file = Path(str(self.annotation_vcf) + '.tbi')
            if not index_file.exists():
                logger.info("Creating missing tabix index for annotation VCF...")
                try:
                    subprocess.run(
                        [str(self.tool_paths['tabix']), '-p', 'vcf', str(self.annotation_vcf)],
                        capture_output=True,
                        text=True,
                        check=True,
                        timeout=300
                    )
                    logger.info("✓ Tabix index created")
                except Exception as e:
                    logger.error(f"Failed to create tabix index for annotation VCF: {e}")
                    raise
            return self.annotation_vcf
        
        logger.info("Compressing annotation VCF for bcftools compatibility...")
        compressed_annotation = Path(str(self.annotation_vcf) + '.gz')
        
        try:
            # Compress with bgzip
            logger.info(f"Compressing {self.annotation_vcf} to {compressed_annotation}...")
            cmd = [str(self.tool_paths['bgzip']), '-c', str(self.annotation_vcf)]
            
            with open(compressed_annotation, 'wb') as f:
                subprocess.run(
                    cmd, 
                    stdout=f, 
                    stderr=subprocess.PIPE, 
                    text=True,
                    check=True,
                    timeout=1800  # 30 minute timeout for compression
                )
            
            # Verify compression worked
            if not compressed_annotation.exists() or compressed_annotation.stat().st_size == 0:
                raise RuntimeError("Compression failed - output file is missing or empty")
            
            logger.info(f"✓ Compression completed, size: {compressed_annotation.stat().st_size:,} bytes")
            
            # Index with tabix
            logger.info("Creating tabix index...")
            subprocess.run(
                [str(self.tool_paths['tabix']), '-p', 'vcf', str(compressed_annotation)],
                capture_output=True,
                text=True,
                check=True,
                timeout=300
            )
            
            # Verify index was created
            index_file = Path(str(compressed_annotation) + '.tbi')
            if not index_file.exists():
                raise RuntimeError("Tabix indexing failed - index file not created")
            
            logger.info(f"✓ Annotation VCF compressed and indexed: {compressed_annotation}")
            
            step_time = time.time() - step_start
            self.stats['processing_steps'].append(('prepare_annotation_vcf', step_time))
            return compressed_annotation
            
        except subprocess.CalledProcessError as e:
            error_msg = f"Failed to prepare annotation VCF: {e}"
            logger.error(error_msg)
            logger.error(f"Command stderr: {e.stderr}")
            
            # Clean up partial files
            if compressed_annotation.exists():
                try:
                    compressed_annotation.unlink()
                    logger.info("Cleaned up partial compressed file")
                except Exception:
                    pass
            
            raise RuntimeError(error_msg)
        except Exception as e:
            logger.error(f"Unexpected error during annotation VCF preparation: {e}")
            raise
    
    def _cleanup_temp_files(self, temp_files: List[Path]) -> None:
        """
        Clean up temporary files created during processing.
        
        Args:
            temp_files: List of temporary files to clean up
        """
        if not temp_files:
            logger.debug("No temporary files to clean up")
            return
        
        logger.info(f"Cleaning up {len(temp_files)} temporary files...")
        cleaned_count = 0
        failed_count = 0
        total_size_cleaned = 0
        
        for temp_file in temp_files:
            try:
                if temp_file.exists():
                    file_size = temp_file.stat().st_size
                    total_size_cleaned += file_size
                    temp_file.unlink()
                    logger.debug(f"Cleaned up temporary file: {temp_file} ({file_size:,} bytes)")
                    cleaned_count += 1
                
                # Also clean up associated index files
                for suffix in ['.tbi', '.csi']:
                    index_file = Path(str(temp_file) + suffix)
                    if index_file.exists():
                        index_size = index_file.stat().st_size
                        total_size_cleaned += index_size
                        index_file.unlink()
                        logger.debug(f"Cleaned up index file: {index_file} ({index_size:,} bytes)")
                        
            except Exception as e:
                logger.warning(f"Failed to clean up temporary file {temp_file}: {e}")
                failed_count += 1
        
        if cleaned_count > 0:
            logger.info(f"✓ Cleaned up {cleaned_count} temporary files ({total_size_cleaned:,} bytes freed)")
        if failed_count > 0:
            logger.warning(f"Failed to clean up {failed_count} temporary files")

    def _generate_summary_statistics(self) -> None:
        """Generate and log comprehensive summary statistics."""
        total_time = time.time() - self.stats['start_time']
        
        logger.info("=== bcftools VCF-to-VCF Annotation Summary ===")
        logger.info(f"Total processing time: {total_time:.2f} seconds")
        logger.info(f"Input VCF: {self.input_vcf}")
        logger.info(f"Annotation VCF: {self.annotation_vcf}")
        logger.info(f"Output VCF: {self.output_vcf}")
        
        # Log file size information
        try:
            if self.input_vcf.exists():
                input_size = self.input_vcf.stat().st_size
                logger.info(f"Input VCF size: {input_size:,} bytes ({input_size/1024/1024:.1f} MB)")
            
            if self.annotation_vcf.exists():
                annotation_size = self.annotation_vcf.stat().st_size
                logger.info(f"Annotation VCF size: {annotation_size:,} bytes ({annotation_size/1024/1024:.1f} MB)")
            
            if self.output_vcf.exists():
                output_size = self.output_vcf.stat().st_size
                logger.info(f"Output VCF size: {output_size:,} bytes ({output_size/1024/1024:.1f} MB)")
        except Exception as e:
            logger.debug(f"Could not get file size information: {e}")
        
        # Log processing step timings
        if self.stats['processing_steps']:
            logger.info("Processing step timings:")
            total_step_time = 0
            for step_name, step_time in self.stats['processing_steps']:
                percentage = (step_time / total_time * 100) if total_time > 0 else 0
                logger.info(f"  {step_name}: {step_time:.2f}s ({percentage:.1f}%)")
                total_step_time += step_time
        
        # Log bcftools operation performance metrics
        self.operation_logger.log_performance_metrics()
        
        # Log tool versions used
        if 'tool_versions' in self.stats:
            logger.info("Tool versions used:")
            for tool, version in self.stats['tool_versions'].items():
                logger.info(f"  {tool}: {version}")
        
        # Log warnings and errors with categorization
        if self.stats['warnings']:
            logger.warning(f"Warnings encountered: {len(self.stats['warnings'])}")
            for warning in self.stats['warnings']:
                logger.warning(f"  {warning}")
        
        if self.stats['errors']:
            logger.error(f"Errors encountered: {len(self.stats['errors'])}")
            for error_info in self.stats['errors']:
                logger.error(f"  {error_info}")
        
        # Log operation statistics from enhanced logger
        operation_stats = self.operation_logger.get_statistics()
        if operation_stats['command_executions'] > 0:
            logger.info("bcftools operation summary:")
            logger.info(f"  Commands executed: {operation_stats['command_executions']}")
            logger.info(f"  Successful commands: {operation_stats['successful_commands']}")
            logger.info(f"  Failed commands: {operation_stats['failed_commands']}")
            logger.info(f"  Total stdout lines: {operation_stats['total_stdout_lines']:,}")
            logger.info(f"  Total stderr lines: {operation_stats['total_stderr_lines']:,}")
        
        logger.info("=== Summary Complete ===")


def main():
    """Main entry point for standalone usage."""
    import argparse
    import sys
    
    parser = argparse.ArgumentParser(
        description="bcftools VCF-to-VCF annotation engine",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
bcftools VCF-to-VCF Annotation Engine:
  This module provides a dedicated bcftools annotation engine for VCF-to-VCF
  annotation without requiring separate header files or column specifications.
  It implements exact coordinate and allele matching between input VCF and
  annotation VCF database.

Key Features:
  - Direct VCF-to-VCF annotation using bcftools annotate
  - Automatic header conflict detection and resolution
  - Exact coordinate and allele matching (CHROM, POS, REF, ALT)
  - Compressed VCF input/output handling with proper indexing
  - Comprehensive error handling and logging

Tool Requirements:
  - bcftools (available in system PATH)
  - tabix (available in system PATH)
  - bgzip (available in system PATH)

Examples:
  # Basic VCF-to-VCF annotation
  python bcftools_annotator.py -i input.vcf.gz -a annotation.vcf.gz -o output.vcf.gz
  
  # Verbose logging
  python bcftools_annotator.py -i input.vcf.gz -a annotation.vcf.gz -o output.vcf.gz --verbose
        """
    )
    
    parser.add_argument(
        '-i', '--input',
        required=True,
        metavar='VCF_FILE',
        help='Input VCF file to be annotated (can be compressed with .gz)'
    )
    
    parser.add_argument(
        '-a', '--annotation',
        required=True,
        metavar='VCF_FILE',
        help='Annotation VCF database file (can be compressed with .gz)'
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
    
    # Run annotation
    try:
        annotator = BcftoolsAnnotator(
            input_vcf=args.input,
            annotation_vcf=args.annotation,
            output_vcf=args.output
        )
        
        annotator.run_annotation()
        
        logger.info("✓ bcftools VCF-to-VCF annotation completed successfully!")
        sys.exit(0)
        
    except KeyboardInterrupt:
        logger.error("Annotation interrupted by user")
        sys.exit(130)
    except Exception as e:
        logger.error(f"Annotation failed: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()