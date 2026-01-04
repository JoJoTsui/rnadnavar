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
  10. Integrated with shared bin/common/vcf_config for consistency across pipeline
  11. Fixed plot_filter_status to include RNA_Edit in secondary view
  12. Fixed plot_tier_distribution to properly display all categories in stacked bars
  13. Fixed IGV reports to eliminate duplicate sample names
  14. Ensured all categories have proper tier organization in IGV reports

INTEGRATION WITH VCF_UTILS:
  - Shares CATEGORY_ORDER, VCF_STAGE_ORDER, and CATEGORY_COLORS with bin/common/vcf_config
  - Falls back to local definitions if shared module unavailable
  - Ready for seamless integration with vcf_utils Nextflow pipeline

For detailed usage, see notebook/vcf_statistics_P2374372.ipynb
"""

import warnings

warnings.filterwarnings("ignore")

# Standard library imports
import os
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

# Constants
TOOLS = ["strelka", "deepsomatic", "mutect2"]
MODALITIES = ["DNA_TUMOR_vs_DNA_NORMAL", "RNA_TUMOR_vs_DNA_NORMAL"]

# Try to import from shared vcf_config module (for integration with vcf_utils)
# Fall back to local definitions if shared module not available
try:
    # Attempt to import from bin/common/vcf_config
    import sys
    from pathlib import Path as _Path

    # Add parent directories to path for shared config discovery
    _current_dir = _Path(__file__).parent.parent.parent
    _bin_common = _current_dir / "bin" / "common"
    if _bin_common.exists() and str(_bin_common) not in sys.path:
        sys.path.insert(0, str(_bin_common))

    from vcf_config import (
        CATEGORY_COLORS as _SHARED_CATEGORY_COLORS,
    )
    from vcf_config import (
        CATEGORY_ORDER as _SHARED_CATEGORY_ORDER,
    )
    from vcf_config import (
        STAGE_DISPLAY_NAMES as _SHARED_STAGE_DISPLAY_NAMES,
    )
    from vcf_config import (
        VCF_STAGE_ORDER as _SHARED_VCF_STAGE_ORDER,
    )

    # Use shared definitions
    CATEGORY_ORDER = _SHARED_CATEGORY_ORDER
    VCF_STAGE_ORDER = _SHARED_VCF_STAGE_ORDER
    STAGE_DISPLAY_NAMES = _SHARED_STAGE_DISPLAY_NAMES
    CATEGORY_COLORS = _SHARED_CATEGORY_COLORS
    print("✓ Using shared vcf_config module from bin/common/")
except (ImportError, ModuleNotFoundError):
    # Fall back to local definitions
    # Variant categories for classification and aggregation
    CATEGORY_ORDER = [
        "Somatic",
        "Germline",
        "Reference",
        "Artifact",
        "RNA_Edit",
        "NoConsensus",
    ]

    # VCF processing stages in order (complete pipeline from normalization to final filtering)
    VCF_STAGE_ORDER = [
        "normalized",  # Individual caller VCFs (Strelka, DeepSomatic, Mutect2)
        "dna_consensus",  # DNA consensus VCF (combined DNA callers)
        "rna_consensus",  # RNA consensus VCF (combined RNA callers)
        "rescue",  # Rescue VCF (DNA + RNA combined)
        "cosmic_gnomad",  # Cosmic/gnomAD annotated rescue VCF
        "rna_editing",  # RNA editing annotated rescue VCF
        "filtered_rescue",  # Final filtered rescue VCF
    ]

    # Human-readable display names for stages
    STAGE_DISPLAY_NAMES = {
        "normalized": "Normalized",
        "dna_consensus": "DNA Consensus",
        "rna_consensus": "RNA Consensus",
        "rescue": "Rescued",
        "cosmic_gnomad": "COSMIC/GnomAD",
        "rna_editing": "RNA Editing",
        "filtered_rescue": "Filtered",
        # Alias for consensus detection
        "consensus": "Consensus",
    }

    # Color scheme for variant categories and annotations
    CATEGORY_COLORS = {
        "Somatic": "#636EFA",  # Blue
        "Germline": "#00CC96",  # Green
        "Reference": "#FFA15A",  # Orange
        "Artifact": "#EF553B",  # Red
        "RNA_Edit": "#AB63FA",  # Purple
        "NoConsensus": "#8A8A8A",  # Gray
        "PASS": "#636EFA",  # Same as Somatic
        "LowQual": "#EF553B",  # Same as Artifact
        "StrandBias": "#AB63FA",  # Purple
        "Clustered": "#FFA500",  # Orange
    }

# Import utilities first (no circular dependency)
from .bam_validator import BAMValidator, RealignmentBAMValidator
from .classifiers import (
    classify_annotated_variant,
    classify_by_filter,
    get_sample_indices,
)

# Import workflow comparison
from .comparison import WorkflowComparator

# Import classes from modules
from .file_discovery import VCFFileDiscovery
from .igv_like import IGVLikePlotter, visualize_variants_igv_like
from .igv_reports_wrapper import (
    check_igv_reports_available,
    create_navigation_page,
    create_subset_vcf,
    generate_igv_report_subprocess,
    get_alignment_index_path,
    organize_by_category_tier,
    visualize_rescue_variants_with_igvreports,
)
from .rescue_analyzer import (
    analyze_rescue_vcf,
    create_resuce_transition_matrix,
    export_rescue_analysis,
)
from .statistics_aggregator import StatisticsAggregator
from .tier_visualizer import TierVisualizer, create_tier_visualization_report
from .tiering import sample_tier_representatives, tier_rescue_variants
from .utils import (
    collect_stage_statistics,
    create_dual_view_plots,
    create_legend_config,
    ensure_all_categories_in_legend,
    filter_no_consensus,
    get_stage_display_name,
    get_stage_order_key,
    should_show_legend,
    sort_stages,
)
from .vcf_processor import VCFStatisticsExtractor, process_all_vcfs
from .visualizer import VCFVisualizer

# Import workflow abstraction
from .workflow import (
    ALL_STAGES,
    STAGE_CONSENSUS,
    STAGE_COSMIC_GNOMAD,
    STAGE_FILTERED_RESCUE,
    STAGE_NORMALIZED,
    STAGE_RESCUE,
    STAGE_RNA_EDITING,
    WORKFLOW_TYPES,
    WorkflowConfig,
    WorkflowManager,
    WorkflowType,
)

# Export all public classes and constants
__all__ = [
    # Utility functions
    "get_stage_order_key",
    "sort_stages",
    "should_show_legend",
    "create_legend_config",
    "filter_no_consensus",
    "get_stage_display_name",
    "create_dual_view_plots",
    "ensure_all_categories_in_legend",
    "collect_stage_statistics",
    # Workflow abstraction
    "WorkflowType",
    "WorkflowConfig",
    "WorkflowManager",
    "WorkflowComparator",
    "WORKFLOW_TYPES",
    "STAGE_NORMALIZED",
    "STAGE_CONSENSUS",
    "STAGE_RESCUE",
    "STAGE_COSMIC_GNOMAD",
    "STAGE_RNA_EDITING",
    "STAGE_FILTERED_RESCUE",
    "ALL_STAGES",
    # Core classes
    "VCFFileDiscovery",
    "VCFStatisticsExtractor",
    "process_all_vcfs",
    "VCFVisualizer",
    "StatisticsAggregator",
    "BAMValidator",
    "RealignmentBAMValidator",
    # Classification functions
    "classify_by_filter",
    "classify_annotated_variant",
    "get_sample_indices",
    # Analysis functions
    "analyze_rescue_vcf",
    "export_rescue_analysis",
    "create_resuce_transition_matrix",
    # Tiering functions
    "tier_rescue_variants",
    "sample_tier_representatives",
    # IGV functions
    "IGVLikePlotter",
    "visualize_variants_igv_like",
    "create_subset_vcf",
    "generate_igv_report_subprocess",
    "create_navigation_page",
    "organize_by_category_tier",
    "visualize_rescue_variants_with_igvreports",
    "check_igv_reports_available",
    "get_alignment_index_path",
    # Tier visualization
    "TierVisualizer",
    "create_tier_visualization_report",
    # Constants
    "TOOLS",
    "MODALITIES",
    "CATEGORY_ORDER",
    "VCF_STAGE_ORDER",
    "STAGE_DISPLAY_NAMES",
    "CATEGORY_COLORS",
]

print("✓ VCF statistics core module initialized")
