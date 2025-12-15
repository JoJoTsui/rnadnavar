#!/usr/bin/env python3
"""
Logging Configuration for RNA Editing Annotation Pipeline

This module provides comprehensive logging configuration for operational monitoring
and integration with Nextflow's logging and reporting system.

Features:
- Configurable log levels for different operational environments
- Structured logging for metrics and monitoring
- Integration with Nextflow reporting
- Performance metrics collection
- Operational status monitoring

Requirements Satisfied: 8.1, 8.2, 8.3, 8.4, 8.5
"""

import logging
import os
import sys
import time
from pathlib import Path
from typing import Dict, Optional, Union
import json

class OperationalLogger:
    """
    Operational logger for RNA editing annotation with Nextflow integration.
    
    Provides structured logging, performance metrics, and operational monitoring
    capabilities designed for production use and integration with monitoring systems.
    """
    
    def __init__(self, 
                 name: str = "rna_editing_annotation",
                 log_level: str = "INFO",
                 output_dir: Optional[Path] = None,
                 enable_metrics: bool = True,
                 enable_performance_monitoring: bool = True):
        """
        Initialize operational logger.
        
        Args:
            name: Logger name
            log_level: Logging level (DEBUG, INFO, WARNING, ERROR)
            output_dir: Directory for log files
            enable_metrics: Enable metrics collection
            enable_performance_monitoring: Enable performance monitoring
        """
        self.name = name
        self.log_level = getattr(logging, log_level.upper(), logging.INFO)
        self.output_dir = output_dir or Path.cwd()
        self.enable_metrics = enable_metrics
        self.enable_performance_monitoring = enable_performance_monitoring
        
        # Initialize logging components
        self.logger = logging.getLogger(name)
        self.metrics = {}
        self.performance_data = {}
        self.start_time = time.time()
        
        # Set up logging configuration
        self._setup_logging()
        
        # Initialize metrics collection
        if self.enable_metrics:
            self._initialize_metrics()
    
    def _setup_logging(self) -> None:
        """Set up comprehensive logging configuration."""
        # Clear any existing handlers
        self.logger.handlers.clear()
        self.logger.setLevel(self.log_level)
        
        # Console handler with operational format
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(self.log_level)
        
        # Operational console format for Nextflow integration
        console_format = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        console_handler.setFormatter(console_format)
        self.logger.addHandler(console_handler)
        
        # File handler for detailed logging
        if self.output_dir:
            try:
                self.output_dir.mkdir(parents=True, exist_ok=True)
                
                # Detailed log file
                log_file = self.output_dir / f"{self.name}_{int(time.time())}.log"
                file_handler = logging.FileHandler(log_file)
                file_handler.setLevel(logging.DEBUG)
                
                # Detailed file format
                file_format = logging.Formatter(
                    '%(asctime)s - %(name)s - %(levelname)s - [%(filename)s:%(lineno)d] - %(funcName)s - %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S'
                )
                file_handler.setFormatter(file_format)
                self.logger.addHandler(file_handler)
                
                self.logger.info(f"Detailed logging to: {log_file}")
                
            except Exception as e:
                self.logger.warning(f"Could not set up file logging: {e}")
    
    def _initialize_metrics(self) -> None:
        """Initialize metrics collection system."""
        self.metrics = {
            'pipeline_start_time': self.start_time,
            'pipeline_start_date': time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(self.start_time)),
            'log_level': logging.getLevelName(self.log_level),
            'metrics_enabled': self.enable_metrics,
            'performance_monitoring_enabled': self.enable_performance_monitoring,
            'counters': {},
            'timers': {},
            'gauges': {},
            'events': []
        }
        
        self.logger.debug("Metrics collection initialized")
    
    def log_metric(self, metric_name: str, value: Union[int, float], 
                   metric_type: str = "gauge", unit: str = "") -> None:
        """
        Log a metric for operational monitoring.
        
        Args:
            metric_name: Name of the metric
            value: Metric value
            metric_type: Type of metric (counter, gauge, timer)
            unit: Unit of measurement
        """
        if not self.enable_metrics:
            return
        
        timestamp = time.time()
        metric_data = {
            'value': value,
            'unit': unit,
            'timestamp': timestamp,
            'type': metric_type
        }
        
        # Store in appropriate metrics category
        if metric_type == "counter":
            self.metrics['counters'][metric_name] = metric_data
        elif metric_type == "timer":
            self.metrics['timers'][metric_name] = metric_data
        else:  # gauge
            self.metrics['gauges'][metric_name] = metric_data
        
        # Log for immediate visibility
        self.logger.info(f"METRIC: {metric_name}={value}{unit} [{metric_type}]")
    
    def log_event(self, event_name: str, event_data: Dict = None, level: str = "INFO") -> None:
        """
        Log an operational event.
        
        Args:
            event_name: Name of the event
            event_data: Additional event data
            level: Log level for the event
        """
        timestamp = time.time()
        event_record = {
            'event': event_name,
            'timestamp': timestamp,
            'data': event_data or {},
            'level': level
        }
        
        if self.enable_metrics:
            self.metrics['events'].append(event_record)
        
        # Log the event
        log_method = getattr(self.logger, level.lower(), self.logger.info)
        log_method(f"EVENT: {event_name} - {event_data or {}}")
    
    def log_performance(self, operation: str, duration: float, 
                       items_processed: int = 0, additional_data: Dict = None) -> None:
        """
        Log performance metrics for an operation.
        
        Args:
            operation: Name of the operation
            duration: Duration in seconds
            items_processed: Number of items processed
            additional_data: Additional performance data
        """
        if not self.enable_performance_monitoring:
            return
        
        rate = items_processed / duration if duration > 0 and items_processed > 0 else 0
        
        performance_record = {
            'operation': operation,
            'duration_seconds': duration,
            'items_processed': items_processed,
            'processing_rate': rate,
            'timestamp': time.time(),
            'additional_data': additional_data or {}
        }
        
        self.performance_data[operation] = performance_record
        
        # Log performance metrics
        self.logger.info(f"PERFORMANCE: {operation} - duration: {duration:.2f}s, "
                        f"items: {items_processed}, rate: {rate:.0f}/s")
        
        # Store as metrics
        self.log_metric(f"{operation}_duration", duration, "timer", "s")
        if items_processed > 0:
            self.log_metric(f"{operation}_items", items_processed, "counter")
            self.log_metric(f"{operation}_rate", rate, "gauge", "/s")
    
    def log_stage(self, stage_name: str, status: str = "START", 
                  additional_data: Dict = None) -> None:
        """
        Log processing stage for pipeline monitoring.
        
        Args:
            stage_name: Name of the processing stage
            status: Stage status (START, COMPLETE, FAILED)
            additional_data: Additional stage data
        """
        timestamp = time.time()
        elapsed = timestamp - self.start_time
        
        stage_data = {
            'stage': stage_name,
            'status': status,
            'elapsed_time': elapsed,
            'timestamp': timestamp,
            'additional_data': additional_data or {}
        }
        
        self.log_event(f"STAGE_{status}", stage_data)
        self.logger.info(f"STAGE: {stage_name} - {status} (elapsed: {elapsed:.2f}s)")
    
    def log_resource_usage(self, context: str) -> Dict:
        """
        Log current resource usage.
        
        Args:
            context: Context for resource usage logging
            
        Returns:
            Dictionary with resource usage data
        """
        resource_data = {'context': context, 'timestamp': time.time()}
        
        try:
            import psutil
            process = psutil.Process()
            
            # Memory usage
            memory_info = process.memory_info()
            memory_mb = memory_info.rss / 1024 / 1024
            resource_data['memory_mb'] = memory_mb
            
            # CPU usage
            cpu_percent = process.cpu_percent()
            resource_data['cpu_percent'] = cpu_percent
            
            # System memory
            system_memory = psutil.virtual_memory()
            resource_data['system_memory_percent'] = system_memory.percent
            resource_data['system_memory_available_mb'] = system_memory.available / 1024 / 1024
            
            # Disk usage
            disk_usage = psutil.disk_usage('.')
            resource_data['disk_free_gb'] = disk_usage.free / 1024 / 1024 / 1024
            resource_data['disk_used_percent'] = (disk_usage.used / disk_usage.total) * 100
            
            self.logger.info(f"RESOURCES: {context} - memory: {memory_mb:.1f}MB, "
                           f"cpu: {cpu_percent:.1f}%, disk_free: {resource_data['disk_free_gb']:.1f}GB")
            
            # Store as metrics
            self.log_metric(f"memory_usage_mb_{context}", memory_mb, "gauge", "MB")
            self.log_metric(f"cpu_usage_percent_{context}", cpu_percent, "gauge", "%")
            
        except ImportError:
            resource_data['error'] = 'psutil not available'
            self.logger.debug(f"RESOURCES: {context} - psutil not available for resource monitoring")
        except Exception as e:
            resource_data['error'] = str(e)
            self.logger.warning(f"RESOURCES: {context} - failed to get resource usage: {e}")
        
        return resource_data
    
    def get_metrics_summary(self) -> Dict:
        """Get comprehensive metrics summary for reporting."""
        if not self.enable_metrics:
            return {'metrics_disabled': True}
        
        total_elapsed = time.time() - self.start_time
        
        summary = {
            'pipeline_info': {
                'name': self.name,
                'start_time': self.start_time,
                'total_elapsed_seconds': total_elapsed,
                'total_elapsed_minutes': total_elapsed / 60,
                'log_level': logging.getLevelName(self.log_level)
            },
            'metrics': self.metrics,
            'performance_data': self.performance_data,
            'summary_generated_at': time.time()
        }
        
        return summary
    
    def save_metrics_report(self, output_file: Optional[Path] = None) -> Path:
        """
        Save comprehensive metrics report to JSON file.
        
        Args:
            output_file: Output file path (optional)
            
        Returns:
            Path to saved metrics file
        """
        if not output_file:
            output_file = self.output_dir / f"{self.name}_metrics_{int(time.time())}.json"
        
        try:
            metrics_summary = self.get_metrics_summary()
            
            with open(output_file, 'w') as f:
                json.dump(metrics_summary, f, indent=2, default=str)
            
            self.logger.info(f"Metrics report saved to: {output_file}")
            return output_file
            
        except Exception as e:
            self.logger.error(f"Failed to save metrics report: {e}")
            raise
    
    def configure_for_nextflow(self) -> None:
        """Configure logger for optimal Nextflow integration."""
        # Ensure INFO level for Nextflow visibility
        if self.log_level > logging.INFO:
            self.log_level = logging.INFO
            self.logger.setLevel(logging.INFO)
        
        # Log configuration for Nextflow
        self.logger.info("Logger configured for Nextflow integration")
        self.logger.info(f"Log level: {logging.getLevelName(self.log_level)}")
        self.logger.info(f"Metrics enabled: {self.enable_metrics}")
        self.logger.info(f"Performance monitoring enabled: {self.enable_performance_monitoring}")


def get_operational_logger(name: str = "rna_editing_annotation",
                          log_level: Optional[str] = None,
                          output_dir: Optional[Path] = None) -> OperationalLogger:
    """
    Get configured operational logger instance.
    
    Args:
        name: Logger name
        log_level: Log level (from environment or default to INFO)
        output_dir: Output directory for log files
        
    Returns:
        Configured OperationalLogger instance
    """
    # Get log level from environment or use default
    if not log_level:
        log_level = os.environ.get('RNA_EDITING_LOG_LEVEL', 'INFO')
    
    # Get output directory from environment or use provided/default
    if not output_dir:
        output_dir_str = os.environ.get('RNA_EDITING_LOG_DIR')
        if output_dir_str:
            output_dir = Path(output_dir_str)
        else:
            output_dir = Path.cwd()
    
    # Check for metrics and performance monitoring settings
    enable_metrics = os.environ.get('RNA_EDITING_ENABLE_METRICS', 'true').lower() == 'true'
    enable_performance = os.environ.get('RNA_EDITING_ENABLE_PERFORMANCE', 'true').lower() == 'true'
    
    logger = OperationalLogger(
        name=name,
        log_level=log_level,
        output_dir=output_dir,
        enable_metrics=enable_metrics,
        enable_performance_monitoring=enable_performance
    )
    
    # Configure for Nextflow if running in Nextflow environment
    if os.environ.get('NXF_TASK_WORKDIR') or os.environ.get('NEXTFLOW_TASK_WORKDIR'):
        logger.configure_for_nextflow()
    
    return logger


# Environment variable configuration guide
ENVIRONMENT_VARIABLES = {
    'RNA_EDITING_LOG_LEVEL': 'Set logging level (DEBUG, INFO, WARNING, ERROR)',
    'RNA_EDITING_LOG_DIR': 'Directory for log files',
    'RNA_EDITING_ENABLE_METRICS': 'Enable metrics collection (true/false)',
    'RNA_EDITING_ENABLE_PERFORMANCE': 'Enable performance monitoring (true/false)'
}

if __name__ == '__main__':
    # Example usage and testing
    logger = get_operational_logger()
    
    logger.log_stage("TEST_PIPELINE", "START")
    logger.log_metric("test_metric", 42, "gauge", "units")
    logger.log_performance("test_operation", 1.5, 100)
    logger.log_resource_usage("test_context")
    logger.log_event("test_event", {"key": "value"})
    logger.log_stage("TEST_PIPELINE", "COMPLETE")
    
    # Save metrics report
    metrics_file = logger.save_metrics_report()
    print(f"Test metrics saved to: {metrics_file}")