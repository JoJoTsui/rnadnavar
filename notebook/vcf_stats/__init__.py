#!/usr/bin/env python3
"""
VCF Statistics Analysis Package (Refactored v2.0)

Comprehensive analysis of VCF files through the annotation pipeline with support for:

PROCESSING STAGES (in order):
  1. normalized       - VCFs from standalone variant callers (Strelka, DeepSomatic, Mutect2)
  2. consensus        - DNA and RNA consensus VCFs
  3. rescue           - Combined DNA + RNA rescue VCFs
  4. cosmic_gnomad    - Cosmic/gnomAD annotated rescue VCFs
  5. rna_editing      - RNA editing annotated rescue VCFs
  6. filtered_rescue  - Final filtered rescue VCFs (post-RaVeX filtering)

VARIANT CATEGORIES:
  - Somatic       - Somatic variants (biological classification)
  - Germline      - Germline variants
  - Reference     - Reference calls
  - Artifact      - Filtering artifacts
  - RNA_Edit      - RNA editing events (NEW)
  - NoConsensus   - Non-consensus variants (NEW)

KEY FEATURES:
  ✓ Automated VCF discovery with support for sample-pair naming (DT/DN/RT)
  ✓ Multi-caller support (Strelka, DeepSomatic, Mutect2)
  ✓ Count distribution analysis replacing pass/filtered metrics
  ✓ Unified color legend across all visualizations (legend duplication fixed)
  ✓ Annotation progression tracking through pipeline stages
  ✓ BAM validation using final filtered rescue VCF
  ✓ Category-based tiering and IGV-reports visualization
  ✓ Comprehensive aggregation and export (Excel/CSV)

ORGANIZATION:
  - file_discovery.py: Automated VCF and BAM file discovery
  - vcf_processor.py: Statistics extraction with new categories
  - classifiers.py: Variant classification (biological + annotation-based)
  - visualizer.py: 5 main plots + annotation progression visualization
  - statistics_aggregator.py: Count distribution and stage progression tracking
  - bam_validator.py: BAM-based variant validation
  - tiering.py: Fine-grained tier assignment
  - tier_visualizer.py: Tier-based visualizations
  - igv_reports_wrapper.py: IGV-reports integration

REFACTORING CHANGES (v2.0):
  1. Removed raw variant_calling VCF collection (use normalized only)
  2. Added 3 new annotation stages (cosmic_gnomad, rna_editing, filtered_rescue)
  3. Added 2 new variant categories (RNA_Edit, NoConsensus)
  4. Replaced pass/filtered metrics with category count distribution
  5. Fixed color legend duplication in Plot 2
  6. Added Plot 5: Annotation progression visualization
  7. Updated BAM validation to use final filtered rescue VCF
  8. Enhanced aggregation with stage progression tracking
  9. Improved docstrings and inline documentation

For detailed usage, see notebook/vcf_statistics_P2374372.ipynb
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

# Variant categories for classification and aggregation
CATEGORY_ORDER = ["Somatic", "Germline", "Reference", "Artifact", "RNA_Edit", "NoConsensus"]

# VCF processing stages in order (for annotation pipeline)
VCF_STAGE_ORDER = ["rescue", "cosmic_gnomad", "rna_editing", "filtered_rescue"]

# Color scheme for variant categories and annotations
CATEGORY_COLORS = {
    "Somatic": "#636EFA",        # Blue
    "Germline": "#00CC96",       # Green
    "Reference": "#FFA15A",      # Orange
    "Artifact": "#EF553B",       # Red
    "RNA_Edit": "#AB63FA",       # Purple
    "NoConsensus": "#8A8A8A",    # Gray
    "PASS": "#636EFA",           # Same as Somatic
    "LowQual": "#EF553B",        # Same as Artifact
    "StrandBias": "#AB63FA",     # Purple
    "Clustered": "#FFA500",      # Orange
}

# Import classes from modules
from .file_discovery import VCFFileDiscovery
from .vcf_processor import VCFStatisticsExtractor, process_all_vcfs
from .classifiers import (
    classify_strelka_variant,
    classify_deepsomatic_variant,
    classify_mutect2_variant,
    classify_variant,
    detect_rna_edit_variant,
    detect_cosmic_gnomad_annotation,
    detect_no_consensus_variant,
    classify_annotated_variant,
    get_sample_indices
)
from .visualizer import VCFVisualizer
from .statistics_aggregator import StatisticsAggregator
from .bam_validator import BAMValidator
from .rescue_analyzer import analyze_rescue_vcf, export_rescue_analysis, create_resuce_transition_matrix
from .tiering import tier_rescue_variants, sample_tier_representatives
from .igv_like import IGVLikePlotter, visualize_variants_igv_like
from .igv_reports_wrapper import (
    create_subset_vcf,
    generate_igv_report_subprocess,
    create_navigation_page,
    organize_by_category_tier,
    visualize_rescue_variants_with_igvreports,
    check_igv_reports_available,
    get_alignment_index_path,
)
from .tier_visualizer import TierVisualizer, create_tier_visualization_report

# Export all public classes and constants
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
    'detect_rna_edit_variant',
    'detect_cosmic_gnomad_annotation',
    'detect_no_consensus_variant',
    'classify_annotated_variant',
    'get_sample_indices',
    'analyze_rescue_vcf',
    'export_rescue_analysis',
    'create_resuce_transition_matrix',
    'tier_rescue_variants',
    'sample_tier_representatives',
    'IGVLikePlotter',
    'visualize_variants_igv_like',
    'create_subset_vcf',
    'generate_igv_report_subprocess',
    'create_navigation_page',
    'organize_by_category_tier',
    'visualize_rescue_variants_with_igvreports',
    'check_igv_reports_available',
    'get_alignment_index_path',
    'TierVisualizer',
    'create_tier_visualization_report',
    'TOOLS',
    'MODALITIES',
    'CATEGORY_ORDER',
    'VCF_STAGE_ORDER',
    'CATEGORY_COLORS'
]

print("✓ VCF statistics core module initialized")