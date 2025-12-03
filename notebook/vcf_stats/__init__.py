#!/usr/bin/env python3
"""
VCF Statistics Core Module

Core functionality for VCF statistics extraction and analysis.
Extracted from notebook/vcf_statistics.ipynb for better organization.
"""

import warnings
warnings.filterwarnings("ignore")

# Standard library imports
import os
import sys
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Tuple, Optional, Any

# Constants
TOOLS = ["strelka", "deepsomatic", "mutect2"]
MODALITIES = ["DNA_TUMOR_vs_DNA_NORMAL", "RNA_TUMOR_vs_DNA_NORMAL"]
CATEGORY_ORDER = ["Somatic", "Germline", "Reference", "Artifact"]

# Import classes from modules
from .file_discovery import VCFFileDiscovery
from .vcf_processor import VCFStatisticsExtractor, process_all_vcfs
from .classifiers import (
    classify_strelka_variant,
    classify_deepsomatic_variant,
    classify_mutect2_variant,
    classify_variant,
    get_sample_indices
)
from .visualizer import VCFVisualizer
from .statistics_aggregator import StatisticsAggregator
from .bam_validator import BAMValidator
from .rescue_analyzer import analyze_rescue_vcf, export_rescue_analysis, create_resuce_transition_matrix

# Export all public classes
__all__ = [
    'VCFFileDiscovery',
    'VCFStatisticsExtractor',
    'process_all_vcfs',
    'VCFVisualizer',
    'StatisticsAggregator',
    'BAMValidator',
    'classify_strelka_variant',
    'classify_deepsomatic_variant',
    'classify_mutect2_variant',
    'classify_variant',
    'get_sample_indices',
    'analyze_rescue_vcf',
    'export_rescue_analysis',
    'create_resuce_transition_matrix',
    'TOOLS',
    'MODALITIES',
    'CATEGORY_ORDER'
]

print("✓ VCF statistics core module initialized")