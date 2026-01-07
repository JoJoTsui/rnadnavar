"""
Caller Comparison Module for VCF Statistics Analysis

This module provides functionality for comparing variant calling results across
different callers (Strelka, DeepSomatic, Mutect2) using:
- Variant set extraction with FILTER category classification
- 3-way Venn diagrams showing SNP/INDEL/shared overlap
- Interactive plotly visualizations with HTML export
- Comprehensive comparison statistics

Author: VCF Statistics Pipeline
Date: 2026-01-07
"""

from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

import pandas as pd
import plotly.graph_objects as go
from cyvcf2 import VCF

# Import shared modules
from vcf_config import VARIANT_TYPE_ORDER, get_variant_type

from .classifiers import classify_by_filter, get_sample_indices

# Type aliases for clarity
VariantKey = Tuple[str, int, str, str, str]  # (CHROM, POS, REF, ALT, FILTER_CATEGORY)
VariantSet = Set[VariantKey]

# Caller color scheme - more distinguishable colors
CALLER_COLORS = {
    "strelka": "#E74C3C",  # Red
    "deepsomatic": "#3498DB",  # Blue
    "mutect2": "#2ECC71",  # Green
}


def extract_classified_variant_sets(
    vcf_path: Path, stage: str = "normalized", caller_name: Optional[str] = None
) -> Dict[str, VariantSet]:
    """
    Extract variant sets from a VCF file with FILTER category classification.

    For normalized VCFs, variants are first classified by FILTER using classify_by_filter()
    to assign them to categories (Somatic, Germline, Reference, Artifact, etc.).
    Variants are then keyed by (CHROM, POS, REF, ALT, FILTER_CATEGORY).

    Args:
        vcf_path: Path to VCF file (must be normalized for caller comparison)
        stage: VCF processing stage ("normalized", "consensus", "rescue")
        caller_name: Name of the caller (strelka, deepsomaticdna, mutect2, etc.)
                    If None, will attempt to infer from filename

    Returns:
        Dictionary mapping variant type to sets of variants:
        {
            "SNP": {(chr1, 100, A, G, Somatic), ...},
            "INDEL": {(chr2, 200, AT, A, Germline), ...},
            "OTHER": {(chr3, 300, A, G,AT, Artifact), ...}
        }
    """
    if not vcf_path.exists():
        raise FileNotFoundError(f"VCF file not found: {vcf_path}")

    # Infer caller name from filename if not provided
    if caller_name is None:
        filename_lower = vcf_path.name.lower()
        if "strelka" in filename_lower:
            caller_name = "strelka"
        elif "deepsomatic" in filename_lower:
            caller_name = "deepsomatic"
        elif "mutect2" in filename_lower:
            caller_name = "mutect2"
        else:
            caller_name = "unknown"

    # Initialize variant sets by type
    variant_sets = {"SNP": set(), "INDEL": set(), "OTHER": set()}

    # Open VCF and get sample indices
    vcf = VCF(str(vcf_path))
    sample_indices = get_sample_indices(vcf, caller_name)

    # Process each variant
    for variant in vcf:
        # Classify variant by FILTER (done once per variant position)
        filter_category = classify_by_filter(
            variant,
            stage=stage,
            caller_name=caller_name,
            sample_indices=sample_indices,
        )

        # Use get_variant_type for consistent type detection
        var_type = get_variant_type(variant)

        # For multi-allelic sites, create separate keys for each ALT
        if len(variant.ALT) > 1:
            for alt_allele in variant.ALT:
                var_key = (
                    variant.CHROM,
                    variant.POS,
                    variant.REF,
                    alt_allele,
                    filter_category,
                )
                variant_sets[var_type].add(var_key)
        else:
            # Single ALT allele
            var_key = (
                variant.CHROM,
                variant.POS,
                variant.REF,
                variant.ALT[0],
                filter_category,
            )
            variant_sets[var_type].add(var_key)

    vcf.close()

    return variant_sets


def compute_venn_overlaps(
    set1: VariantSet, set2: VariantSet, set3: VariantSet, labels: Tuple[str, str, str]
) -> Dict[str, Any]:
    """
    Compute overlaps for 3-way Venn diagram.

    Args:
        set1, set2, set3: Variant sets from three callers
        labels: Tuple of caller names (e.g., ("Strelka", "DeepSomatic", "Mutect2"))

    Returns:
        Dictionary with overlap counts and sets:
        {
            "label1": str, "label2": str, "label3": str,
            "only_1": int, "only_2": int, "only_3": int,
            "overlap_1_2": int, "overlap_1_3": int, "overlap_2_3": int,
            "overlap_all": int,
            "set1": int, "set2": int, "set3": int,
            "variants_only_1": set, "variants_only_2": set, ...
        }
    """
    # Compute intersections
    overlap_1_2 = set1 & set2 - set3
    overlap_1_3 = set1 & set3 - set2
    overlap_2_3 = set2 & set3 - set1
    overlap_all = set1 & set2 & set3

    only_1 = set1 - set2 - set3
    only_2 = set2 - set1 - set3
    only_3 = set3 - set1 - set2

    return {
        # Labels
        "label1": labels[0],
        "label2": labels[1],
        "label3": labels[2],
        # Counts
        "only_1": len(only_1),
        "only_2": len(only_2),
        "only_3": len(only_3),
        "overlap_1_2": len(overlap_1_2),
        "overlap_1_3": len(overlap_1_3),
        "overlap_2_3": len(overlap_2_3),
        "overlap_all": len(overlap_all),
        "set1": len(set1),
        "set2": len(set2),
        "set3": len(set3),
        # Actual variant sets (for detailed inspection)
        "variants_only_1": only_1,
        "variants_only_2": only_2,
        "variants_only_3": only_3,
        "variants_overlap_1_2": overlap_1_2,
        "variants_overlap_1_3": overlap_1_3,
        "variants_overlap_2_3": overlap_2_3,
        "variants_overlap_all": overlap_all,
    }


def compute_3way_venn_plotly(
    overlaps: Dict[str, Any],
    title: str = "3-Way Venn Diagram",
    colors: Tuple[str, str, str] = ("#FF6B6B", "#4ECDC4", "#45B7D1"),
) -> go.Figure:
    """
    Create interactive 3-way Venn diagram using plotly.

    Uses three overlapping circles with hover information showing variant counts.
    Supports HTML export for embedding in notebooks.

    Args:
        overlaps: Dictionary from compute_venn_overlaps()
        title: Plot title
        colors: Tuple of RGB/hex colors for three circles

    Returns:
        Plotly Figure object with interactive Venn diagram
    """
    # Circle parameters (positions and radius optimized for overlap visualization)
    # Circle 1 (left): centered at (0.5, 0.8)
    # Circle 2 (right): centered at (1.5, 0.8)
    # Circle 3 (bottom): centered at (1.0, 0.2)
    radius = 0.6

    circle_params = [
        {"x": 0.5, "y": 0.8, "color": colors[0], "label": overlaps["label1"]},
        {"x": 1.5, "y": 0.8, "color": colors[1], "label": overlaps["label2"]},
        {"x": 1.0, "y": 0.2, "color": colors[2], "label": overlaps["label3"]},
    ]

    # Create figure
    fig = go.Figure()

    # Add three circles as shapes
    for params in circle_params:
        # Add circle as a shape
        fig.add_shape(
            type="circle",
            x0=params["x"] - radius,
            y0=params["y"] - radius,
            x1=params["x"] + radius,
            y1=params["y"] + radius,
            line=dict(color=params["color"], width=2),
            fillcolor=params["color"],
            opacity=0.3,
        )

    # Add text annotations for region counts
    # Position annotations at approximate centers of each region
    annotations = []

    # Only set 1 (left circle, upper left)
    annotations.append(
        dict(
            x=0.2,
            y=1.0,
            text=f"<b>{overlaps['only_1']}</b>",
            showarrow=False,
            font=dict(size=14, color="black"),
            xanchor="center",
        )
    )

    # Only set 2 (right circle, upper right)
    annotations.append(
        dict(
            x=1.8,
            y=1.0,
            text=f"<b>{overlaps['only_2']}</b>",
            showarrow=False,
            font=dict(size=14, color="black"),
            xanchor="center",
        )
    )

    # Only set 3 (bottom circle, lower center)
    annotations.append(
        dict(
            x=1.0,
            y=-0.1,
            text=f"<b>{overlaps['only_3']}</b>",
            showarrow=False,
            font=dict(size=14, color="black"),
            xanchor="center",
        )
    )

    # Overlap 1-2 (top center)
    annotations.append(
        dict(
            x=1.0,
            y=1.0,
            text=f"<b>{overlaps['overlap_1_2']}</b>",
            showarrow=False,
            font=dict(size=14, color="black"),
            xanchor="center",
        )
    )

    # Overlap 1-3 (left-bottom)
    annotations.append(
        dict(
            x=0.6,
            y=0.5,
            text=f"<b>{overlaps['overlap_1_3']}</b>",
            showarrow=False,
            font=dict(size=14, color="black"),
            xanchor="center",
        )
    )

    # Overlap 2-3 (right-bottom)
    annotations.append(
        dict(
            x=1.4,
            y=0.5,
            text=f"<b>{overlaps['overlap_2_3']}</b>",
            showarrow=False,
            font=dict(size=14, color="black"),
            xanchor="center",
        )
    )

    # Overlap all three (center)
    annotations.append(
        dict(
            x=1.0,
            y=0.65,
            text=f"<b>{overlaps['overlap_all']}</b>",
            showarrow=False,
            font=dict(size=16, color="darkblue"),
            xanchor="center",
        )
    )

    # Add circle labels
    annotations.append(
        dict(
            x=0.5,
            y=1.6,
            text=f"<b>{overlaps['label1']}</b><br>({overlaps['set1']})",
            showarrow=False,
            font=dict(size=12, color=colors[0]),
            xanchor="center",
        )
    )

    annotations.append(
        dict(
            x=1.5,
            y=1.6,
            text=f"<b>{overlaps['label2']}</b><br>({overlaps['set2']})",
            showarrow=False,
            font=dict(size=12, color=colors[1]),
            xanchor="center",
        )
    )

    annotations.append(
        dict(
            x=1.0,
            y=-0.6,
            text=f"<b>{overlaps['label3']}</b><br>({overlaps['set3']})",
            showarrow=False,
            font=dict(size=12, color=colors[2]),
            xanchor="center",
        )
    )

    # Update layout
    fig.update_layout(
        title=dict(text=title, x=0.5, xanchor="center", font=dict(size=16)),
        xaxis=dict(
            showgrid=False, showticklabels=False, zeroline=False, range=[-0.5, 2.5]
        ),
        yaxis=dict(
            showgrid=False,
            showticklabels=False,
            zeroline=False,
            range=[-1.0, 2.0],
            scaleanchor="x",
            scaleratio=1,
        ),
        annotations=annotations,
        plot_bgcolor="white",
        width=600,
        height=600,
        margin=dict(l=50, r=50, t=80, b=50),
    )

    return fig


class CallerComparator:
    """
    Comprehensive caller comparison analysis for VCF statistics.

    Compares variant calling results across Strelka, DeepSomatic, and Mutect2
    for both DNA and RNA modalities using:
    - Classified variant sets (with FILTER category)
    - 3-way Venn diagrams for SNPs, INDELs, and shared variants
    - Summary statistics tables
    - Interactive HTML visualizations

    Attributes:
        all_vcfs: Dictionary of all discovered VCF files organized by stage
        variant_sets: Cached variant sets by modality and caller
        comparison_results: Cached comparison results
    """

    def __init__(self, all_vcfs: Dict[str, Dict[str, Dict[str, Path]]]):
        """
        Initialize caller comparator with discovered VCF files.

        Args:
            all_vcfs: Nested dict from VCFFileDiscovery.discover_all_workflows()
                     Format: {workflow: {stage: {file_id: path}}}
        """
        self.all_vcfs = all_vcfs
        self.variant_sets = {}
        self.comparison_results = {}

    def _extract_caller_from_filename(self, filename: str) -> Optional[str]:
        """Extract caller name from VCF filename."""
        filename_lower = filename.lower()
        if "strelka" in filename_lower:
            return "strelka"
        elif "deepsomatic" in filename_lower:
            return "deepsomatic"
        elif "mutect2" in filename_lower:
            return "mutect2"
        return None

    def _extract_modality_from_filename(self, filename: str) -> Optional[str]:
        """Extract modality (DNA/RNA) from VCF filename."""
        filename_upper = filename.upper()
        if "DT" in filename_upper or "DNA_TUMOR" in filename_upper:
            return "DNA"
        elif "RT" in filename_upper or "RNA_TUMOR" in filename_upper:
            return "RNA"
        return None

    def _generate_focused_caller_venns(
        self,
        variant_sets_dict: Dict[str, Dict[str, VariantSet]],
        labels: List[str],
        comparison_title: str,
    ) -> Tuple[List[Dict], Dict[str, go.Figure]]:
        """
        Generate focused Venn diagrams for caller comparison:
        - Somatic SNP and Somatic INDEL (per-category)
        - All SNP and All INDEL (aggregated across all categories)

        Args:
            variant_sets_dict: {label: {var_type: variant_set}} for each comparison group
            labels: List of labels (e.g., ["strelka", "deepsomatic", "mutect2"])
            comparison_title: Title prefix for Venn diagrams

        Returns:
            Tuple of (summary_rows, figures_dict)
        """
        figures = {}
        summary_rows = []

        # Get caller colors
        colors = tuple(
            [CALLER_COLORS.get(label.lower(), "#999999") for label in labels]
        )

        # 1. Generate Somatic-only Venn diagrams (SNP and INDEL)
        print("  Generating Somatic SNP and INDEL Venn diagrams...")
        for var_type in ["SNP", "INDEL"]:
            # Filter variant sets by Somatic category
            somatic_sets = []
            for label in labels:
                somatic_set = {
                    vk
                    for vk in variant_sets_dict[label][var_type]
                    if vk[4] == "Somatic"
                }
                somatic_sets.append(somatic_set)

            # Skip if all sets are empty
            if all(len(s) == 0 for s in somatic_sets):
                print(f"    Skipping Somatic {var_type} (all callers have 0 variants)")
                continue

            # Compute overlaps
            overlaps = compute_venn_overlaps(
                somatic_sets[0],
                somatic_sets[1],
                somatic_sets[2],
                labels=tuple(labels),
            )

            # Create Venn diagram
            fig_key = f"{comparison_title}_Somatic_{var_type}"
            fig = compute_3way_venn_plotly(
                overlaps,
                title=f"{comparison_title}\nSomatic {var_type}",
                colors=colors,
            )
            figures[fig_key] = fig

            # Add to summary
            summary_rows.append(
                {
                    "Comparison": comparison_title,
                    "Category": "Somatic",
                    "Variant_Type": var_type,
                    f"{labels[0]}_count": overlaps["set1"],
                    f"{labels[1]}_count": overlaps["set2"],
                    f"{labels[2]}_count": overlaps["set3"],
                    "All_Three": overlaps["overlap_all"],
                    "Any_Two_or_More": overlaps["overlap_all"]
                    + overlaps["overlap_1_2"]
                    + overlaps["overlap_1_3"]
                    + overlaps["overlap_2_3"],
                }
            )

        # 2. Generate All-categories Venn diagrams (SNP and INDEL)
        print("  Generating All-categories SNP and INDEL Venn diagrams...")
        for var_type in ["SNP", "INDEL"]:
            # Aggregate all categories for each caller
            # Keep the category in the key so variants must match position AND category
            all_sets = []
            for label in labels:
                # Keep full key including category: (CHROM, POS, REF, ALT, CATEGORY)
                all_set = set(variant_sets_dict[label][var_type])
                all_sets.append(all_set)

            # Skip if all sets are empty
            if all(len(s) == 0 for s in all_sets):
                print(f"    Skipping All {var_type} (all callers have 0 variants)")
                continue

            # Compute overlaps
            overlaps = compute_venn_overlaps(
                all_sets[0],
                all_sets[1],
                all_sets[2],
                labels=tuple(labels),
            )

            # Create Venn diagram
            fig_key = f"{comparison_title}_All_{var_type}"
            fig = compute_3way_venn_plotly(
                overlaps,
                title=f"{comparison_title}\nAll Categories - {var_type}",
                colors=colors,
            )
            figures[fig_key] = fig

            # Add to summary
            summary_rows.append(
                {
                    "Comparison": comparison_title,
                    "Category": "All",
                    "Variant_Type": var_type,
                    f"{labels[0]}_count": overlaps["set1"],
                    f"{labels[1]}_count": overlaps["set2"],
                    f"{labels[2]}_count": overlaps["set3"],
                    "All_Three": overlaps["overlap_all"],
                    "Any_Two_or_More": overlaps["overlap_all"]
                    + overlaps["overlap_1_2"]
                    + overlaps["overlap_1_3"]
                    + overlaps["overlap_2_3"],
                }
            )

        return summary_rows, figures

    def _generate_per_category_venns(
        self,
        variant_sets_dict: Dict[str, Dict[str, VariantSet]],
        labels: List[str],
        comparison_title: str,
        colors: Tuple[str, str, str] = ("#E74C3C", "#3498DB", "#2ECC71"),
    ) -> Tuple[List[Dict], Dict[str, go.Figure]]:
        """
        Generate per-category Venn diagrams for consensus/workflow comparison.
        Keeps all categories (used for consensus comparison, not caller comparison).

        Args:
            variant_sets_dict: {label: {var_type: variant_set}} for each comparison group
            labels: List of labels (e.g., ["DNA_Standard", "RNA_Standard", "RNA_Realignment"])
            comparison_title: Title prefix for Venn diagrams
            colors: Tuple of colors for the three sets

        Returns:
            Tuple of (summary_rows, figures_dict)
        """
        figures = {}
        summary_rows = []

        # Get all variant sets for determining categories
        all_categories = set()
        for var_type in VARIANT_TYPE_ORDER:
            for label in labels:
                for var_key in variant_sets_dict[label][var_type]:
                    all_categories.add(var_key[4])  # FILTER_CATEGORY is index 4

        # Generate Venn diagram for each category and variant type
        for category in sorted(all_categories):
            for var_type in ["SNP", "INDEL"]:  # Only SNP and INDEL for consensus
                # Filter variant sets by category
                category_sets = []
                for label in labels:
                    cat_set = {
                        vk
                        for vk in variant_sets_dict[label][var_type]
                        if vk[4] == category
                    }
                    category_sets.append(cat_set)

                # Skip if all sets are empty
                if all(len(s) == 0 for s in category_sets):
                    continue

                # Compute overlaps
                overlaps = compute_venn_overlaps(
                    category_sets[0],
                    category_sets[1],
                    category_sets[2],
                    labels=tuple(labels),
                )

                # Create Venn diagram
                fig_key = f"{comparison_title}_{category}_{var_type}"
                fig = compute_3way_venn_plotly(
                    overlaps,
                    title=f"{comparison_title} - {category} {var_type}",
                    colors=colors,
                )
                figures[fig_key] = fig

                # Add to summary
                summary_rows.append(
                    {
                        "Comparison": comparison_title,
                        "Category": category,
                        "Variant_Type": var_type,
                        f"{labels[0]}_count": overlaps["set1"],
                        f"{labels[1]}_count": overlaps["set2"],
                        f"{labels[2]}_count": overlaps["set3"],
                        "All_Three": overlaps["overlap_all"],
                        "Any_Two_or_More": overlaps["overlap_all"]
                        + overlaps["overlap_1_2"]
                        + overlaps["overlap_1_3"]
                        + overlaps["overlap_2_3"],
                    }
                )

        return summary_rows, figures

    def compare_dna_callers(
        self, workflow: str = "standard", stage: str = "normalized"
    ) -> Tuple[pd.DataFrame, Dict[str, go.Figure]]:
        """
        Compare DNA variant callers (Strelka, DeepSomatic, Mutect2) with per-category Venn diagrams.

        Args:
            workflow: Workflow name ("standard" or "realignment")
            stage: Processing stage ("normalized" recommended for caller comparison)

        Returns:
            Tuple of (summary_df, figures_dict) where:
            - summary_df: DataFrame with variant counts and overlaps per category
            - figures_dict: Dictionary of plotly figures keyed by "{category}_{vartype}"
        """
        # Extract DNA normalized VCFs
        if workflow not in self.all_vcfs:
            raise ValueError(f"Workflow '{workflow}' not found in VCF files")

        if stage not in self.all_vcfs[workflow]:
            raise ValueError(f"Stage '{stage}' not found in workflow '{workflow}'")

        stage_vcfs = self.all_vcfs[workflow][stage]

        # Group by caller
        dna_vcfs = {}
        for file_id, vcf_info in stage_vcfs.items():
            # Extract path from metadata dict or use directly
            if isinstance(vcf_info, dict):
                vcf_path = vcf_info.get("path", vcf_info)
            else:
                vcf_path = vcf_info

            # Ensure we have a Path object
            if not isinstance(vcf_path, Path):
                vcf_path = Path(vcf_path)

            caller = self._extract_caller_from_filename(vcf_path.name)
            modality = self._extract_modality_from_filename(vcf_path.name)

            if caller and modality == "DNA":
                dna_vcfs[caller] = vcf_path

        if len(dna_vcfs) < 2:
            raise ValueError(
                f"Need at least 2 DNA callers, found {len(dna_vcfs)}: {list(dna_vcfs.keys())}"
            )

        print(f"Found DNA callers: {list(dna_vcfs.keys())}")

        # Extract variant sets for each caller
        caller_variant_sets = {}
        for caller, vcf_path in dna_vcfs.items():
            print(f"  Extracting variants from {caller}...")
            caller_variant_sets[caller] = extract_classified_variant_sets(
                vcf_path, stage=stage, caller_name=caller
            )

        # Prepare for 3-way comparison (if we have exactly 3 callers)
        callers = list(caller_variant_sets.keys())
        if len(callers) == 3:
            print("  Generating focused Venn diagrams (Somatic + All categories)...")
            summary_rows, figures = self._generate_focused_caller_venns(
                caller_variant_sets,
                callers,
                f"DNA_{workflow}_{stage}",
            )
            summary_df = pd.DataFrame(summary_rows)

        else:
            # For 2 callers, create pairwise comparison
            summary_rows = []
            figures = {}

            for var_type in VARIANT_TYPE_ORDER:
                set1 = caller_variant_sets[callers[0]][var_type]
                set2 = caller_variant_sets[callers[1]][var_type]

                overlap = set1 & set2
                only_1 = set1 - set2
                only_2 = set2 - set1

                summary_rows.append(
                    {
                        "Variant_Type": var_type,
                        "Caller_1": callers[0],
                        "Caller_1_Count": len(set1),
                        "Caller_2": callers[1],
                        "Caller_2_Count": len(set2),
                        "Overlap": len(overlap),
                        "Only_Caller_1": len(only_1),
                        "Only_Caller_2": len(only_2),
                    }
                )

            summary_df = pd.DataFrame(summary_rows)

        return summary_df, figures

    def compare_rna_callers(
        self, workflow: str = "standard", stage: str = "normalized"
    ) -> Tuple[pd.DataFrame, Dict[str, go.Figure]]:
        """
        Compare RNA variant callers with per-category Venn diagrams.

        Similar to compare_dna_callers() but for RNA modality.

        Args:
            workflow: Workflow name ("standard" or "realignment")
            stage: Processing stage ("normalized" recommended)

        Returns:
            Tuple of (summary_df, figures_dict) where:
            - summary_df: DataFrame with variant counts and overlaps per category
            - figures_dict: Dictionary of plotly figures keyed by "{category}_{vartype}"
        """
        # Extract RNA normalized VCFs
        if workflow not in self.all_vcfs:
            raise ValueError(f"Workflow '{workflow}' not found in VCF files")

        if stage not in self.all_vcfs[workflow]:
            raise ValueError(f"Stage '{stage}' not found in workflow '{workflow}'")

        stage_vcfs = self.all_vcfs[workflow][stage]

        # Group by caller
        rna_vcfs = {}
        for file_id, vcf_info in stage_vcfs.items():
            # Extract path from metadata dict or use directly
            if isinstance(vcf_info, dict):
                vcf_path = vcf_info.get("path", vcf_info)
            else:
                vcf_path = vcf_info

            # Ensure we have a Path object
            if not isinstance(vcf_path, Path):
                vcf_path = Path(vcf_path)

            caller = self._extract_caller_from_filename(vcf_path.name)
            modality = self._extract_modality_from_filename(vcf_path.name)

            if caller and modality == "RNA":
                rna_vcfs[caller] = vcf_path

        if len(rna_vcfs) < 2:
            raise ValueError(
                f"Need at least 2 RNA callers, found {len(rna_vcfs)}: {list(rna_vcfs.keys())}"
            )

        print(f"Found RNA callers: {list(rna_vcfs.keys())}")

        # Extract variant sets for each caller
        caller_variant_sets = {}
        for caller, vcf_path in rna_vcfs.items():
            print(f"  Extracting variants from {caller}...")
            caller_variant_sets[caller] = extract_classified_variant_sets(
                vcf_path, stage=stage, caller_name=caller
            )

        # Prepare for 3-way comparison
        callers = list(caller_variant_sets.keys())
        if len(callers) == 3:
            print("  Generating focused Venn diagrams (Somatic + All categories)...")
            summary_rows, figures = self._generate_focused_caller_venns(
                caller_variant_sets,
                callers,
                f"RNA_{workflow}_{stage}",
            )
            summary_df = pd.DataFrame(summary_rows)

        else:
            # Pairwise comparison
            summary_rows = []
            figures = {}

            for var_type in VARIANT_TYPE_ORDER:
                set1 = caller_variant_sets[callers[0]][var_type]
                set2 = caller_variant_sets[callers[1]][var_type]

                overlap = set1 & set2
                only_1 = set1 - set2
                only_2 = set2 - set1

                summary_rows.append(
                    {
                        "Variant_Type": var_type,
                        "Caller_1": callers[0],
                        "Caller_1_Count": len(set1),
                        "Caller_2": callers[1],
                        "Caller_2_Count": len(set2),
                        "Overlap": len(overlap),
                        "Only_Caller_1": len(only_1),
                        "Only_Caller_2": len(only_2),
                    }
                )

            summary_df = pd.DataFrame(summary_rows)

        return summary_df, figures

    def compare_consensus_vcfs(
        self, include_realignment: bool = True
    ) -> Tuple[pd.DataFrame, Dict[str, go.Figure]]:
        """
        Create 3-way Venn diagrams comparing DNA, RNA standard, and RNA realignment consensus VCFs.

        Compares all three workflows together with per-category Venn diagrams:
        - DNA standard consensus
        - RNA standard consensus
        - RNA realignment consensus (if available)

        Args:
            include_realignment: Whether to include realignment consensus comparison

        Returns:
            Tuple of (summary_df, figures_dict) where figures are keyed by "{category}_{vartype}"
        """
        stage = "consensus"

        # Get standard consensus VCFs
        if "standard" not in self.all_vcfs:
            raise ValueError("Standard workflow not found")
        if stage not in self.all_vcfs["standard"]:
            raise ValueError(f"Stage '{stage}' not found in standard workflow")

        standard_vcfs = self.all_vcfs["standard"][stage]
        dna_consensus = None
        rna_standard_consensus = None

        for file_id, vcf_info in standard_vcfs.items():
            if isinstance(vcf_info, dict):
                vcf_path = vcf_info.get("path", vcf_info)
            else:
                vcf_path = vcf_info
            if not isinstance(vcf_path, Path):
                vcf_path = Path(vcf_path)

            modality = self._extract_modality_from_filename(vcf_path.name)
            if modality == "DNA":
                dna_consensus = vcf_path
            elif modality == "RNA":
                rna_standard_consensus = vcf_path

        # Get realignment consensus VCF if available
        rna_realign_consensus = None
        if include_realignment and "realignment" in self.all_vcfs:
            if stage in self.all_vcfs["realignment"]:
                realign_vcfs = self.all_vcfs["realignment"][stage]
                for file_id, vcf_info in realign_vcfs.items():
                    if isinstance(vcf_info, dict):
                        vcf_path = vcf_info.get("path", vcf_info)
                    else:
                        vcf_path = vcf_info
                    if not isinstance(vcf_path, Path):
                        vcf_path = Path(vcf_path)

                    modality = self._extract_modality_from_filename(vcf_path.name)
                    if modality == "RNA":
                        rna_realign_consensus = vcf_path
                        break

        # Check if we have all three consensus VCFs for 3-way comparison
        if dna_consensus and rna_standard_consensus and rna_realign_consensus:
            print(
                "\n3-Way Consensus Comparison: DNA vs RNA Standard vs RNA Realignment"
            )
            print(f"  DNA: {dna_consensus.name}")
            print(f"  RNA Standard: {rna_standard_consensus.name}")
            print(f"  RNA Realignment: {rna_realign_consensus.name}")

            # Extract variant sets
            consensus_sets = {
                "DNA_Standard": extract_classified_variant_sets(
                    dna_consensus, stage=stage, caller_name="consensus"
                ),
                "RNA_Standard": extract_classified_variant_sets(
                    rna_standard_consensus, stage=stage, caller_name="consensus"
                ),
                "RNA_Realignment": extract_classified_variant_sets(
                    rna_realign_consensus, stage=stage, caller_name="consensus"
                ),
            }

            # Generate 3-way per-category Venn diagrams
            print("  Generating per-category 3-way Venn diagrams...")
            summary_rows, figures = self._generate_per_category_venns(
                consensus_sets,
                ["DNA_Standard", "RNA_Standard", "RNA_Realignment"],
                "Consensus 3way",
                colors=("#FF6B6B", "#4ECDC4", "#45B7D1"),
            )
            summary_df = pd.DataFrame(summary_rows)

        elif dna_consensus and rna_standard_consensus:
            # Fallback to 2-way comparison if no realignment
            print("\n2-Way Consensus Comparison: DNA vs RNA Standard")
            print(f"  DNA: {dna_consensus.name}")
            print(f"  RNA Standard: {rna_standard_consensus.name}")

            dna_sets = extract_classified_variant_sets(
                dna_consensus, stage=stage, caller_name="consensus"
            )
            rna_std_sets = extract_classified_variant_sets(
                rna_standard_consensus, stage=stage, caller_name="consensus"
            )

            # Get all unique categories
            all_categories = set()
            for var_type in VARIANT_TYPE_ORDER:
                for var_key in dna_sets[var_type] | rna_std_sets[var_type]:
                    all_categories.add(var_key[4])

            figures = {}
            all_summary_rows = []

            # Create per-category 2-way Venn diagrams
            for category in sorted(all_categories):
                for var_type in VARIANT_TYPE_ORDER:
                    dna_cat = {vk for vk in dna_sets[var_type] if vk[4] == category}
                    rna_cat = {vk for vk in rna_std_sets[var_type] if vk[4] == category}

                    # Skip if both empty
                    if len(dna_cat) == 0 and len(rna_cat) == 0:
                        continue

                    # Skip OTHER if all zeros
                    if var_type == "OTHER" and len(dna_cat) == 0 and len(rna_cat) == 0:
                        print(
                            f"  Skipping {category} OTHER Venn (all workflows have 0 variants)"
                        )
                        continue

                    overlap = dna_cat & rna_cat
                    only_dna = dna_cat - rna_cat
                    only_rna = rna_cat - dna_cat

                    all_summary_rows.append(
                        {
                            "Comparison": "DNA_vs_RNA_Standard",
                            "Category": category,
                            "Variant_Type": var_type,
                            "DNA_Consensus": len(dna_cat),
                            "RNA_Consensus": len(rna_cat),
                            "Overlap": len(overlap),
                            "Only_DNA": len(only_dna),
                            "Only_RNA": len(only_rna),
                        }
                    )

                    # Create 2-way Venn
                    if len(dna_cat) > 0 or len(rna_cat) > 0:
                        fig_key = f"Consensus_2way_{category}_{var_type}"
                        figures[fig_key] = self._create_2way_venn(
                            dna_cat,
                            rna_cat,
                            "DNA Standard",
                            "RNA Standard",
                            f"Consensus: DNA vs RNA Standard - {category} {var_type}",
                        )

            summary_df = pd.DataFrame(all_summary_rows)
        else:
            raise ValueError("Need at least DNA and RNA standard consensus VCFs")

        return summary_df, figures

    def _create_2way_venn(
        self, set1: VariantSet, set2: VariantSet, label1: str, label2: str, title: str
    ) -> go.Figure:
        """Create a 2-way Venn diagram."""
        overlap = set1 & set2
        only_1 = set1 - set2
        only_2 = set2 - set1

        fig = go.Figure()

        # Two overlapping circles
        radius = 0.5
        fig.add_shape(
            type="circle",
            x0=0.3 - radius,
            y0=0.5 - radius,
            x1=0.3 + radius,
            y1=0.5 + radius,
            line=dict(color="#FF6B6B", width=2),
            fillcolor="#FF6B6B",
            opacity=0.3,
        )
        fig.add_shape(
            type="circle",
            x0=0.7 - radius,
            y0=0.5 - radius,
            x1=0.7 + radius,
            y1=0.5 + radius,
            line=dict(color="#4ECDC4", width=2),
            fillcolor="#4ECDC4",
            opacity=0.3,
        )

        # Add annotations
        annotations = [
            dict(
                x=0.15,
                y=0.5,
                text=f"<b>{len(only_1)}</b>",
                showarrow=False,
                font=dict(size=14),
            ),
            dict(
                x=0.85,
                y=0.5,
                text=f"<b>{len(only_2)}</b>",
                showarrow=False,
                font=dict(size=14),
            ),
            dict(
                x=0.5,
                y=0.5,
                text=f"<b>{len(overlap)}</b>",
                showarrow=False,
                font=dict(size=16, color="darkblue"),
            ),
            dict(
                x=0.3,
                y=1.2,
                text=f"<b>{label1}</b><br>({len(set1)})",
                showarrow=False,
                font=dict(size=12, color="#FF6B6B"),
            ),
            dict(
                x=0.7,
                y=1.2,
                text=f"<b>{label2}</b><br>({len(set2)})",
                showarrow=False,
                font=dict(size=12, color="#4ECDC4"),
            ),
        ]

        fig.update_layout(
            title=dict(text=title, x=0.5, xanchor="center", font=dict(size=14)),
            xaxis=dict(
                showgrid=False, showticklabels=False, zeroline=False, range=[0, 1]
            ),
            yaxis=dict(
                showgrid=False,
                showticklabels=False,
                zeroline=False,
                range=[0, 1.5],
                scaleanchor="x",
                scaleratio=1,
            ),
            annotations=annotations,
            plot_bgcolor="white",
            width=500,
            height=500,
            margin=dict(l=50, r=50, t=80, b=50),
        )

        return fig
