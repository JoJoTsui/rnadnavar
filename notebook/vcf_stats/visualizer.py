#!/usr/bin/env python3
"""
VCF Statistics Visualizer Module

Create interactive visualizations for VCF statistics including variant counts,
quality distributions, category distributions, and workflow comparisons.

Key Features:
    - Interactive Plotly visualizations
    - Multi-tool and multi-modality comparisons
    - Stage progression tracking
    - Workflow comparison plots (standard vs realignment)
    - RNA-focused comparison visualizations
    - Integrative view plots (DNA + RNA standard + RNA realignment)
    - Customizable color schemes and layouts

Enhanced Features (v2.0):
    - RNA workflow comparison plots
    - RNA stage progression comparison
    - RNA annotation impact heatmaps
    - Integrative view visualizations
    - Workflow-aware plot generation

Usage Example:
    >>> from vcf_stats.visualizer import VCFVisualizer
    >>>
    >>> # Create visualizer for standard workflow
    >>> visualizer = VCFVisualizer(
    ...     all_stats=standard_workflow_stats,
    ...     workflow_type="standard"
    ... )
    >>>
    >>> # Generate standard plots
    >>> fig1 = visualizer.plot_variant_counts_by_tool()
    >>> fig1.show()
    >>>
    >>> fig2 = visualizer.plot_category_distribution()
    >>> fig2.show()
    >>>
    >>> # Generate workflow comparison plots (if realignment available)
    >>> fig3 = visualizer.plot_rna_workflow_comparison(
    ...     rna_standard_stats=rna_standard_stats,
    ...     rna_realignment_stats=rna_realignment_stats
    ... )
    >>> fig3.show()
    >>>
    >>> # Generate integrative view
    >>> fig4 = visualizer.plot_integrative_view(
    ...     dna_stats=dna_stats,
    ...     rna_standard_stats=rna_standard_stats,
    ...     rna_realignment_stats=rna_realignment_stats
    ... )
    >>> fig4.show()

Visualization Types:
    Standard Plots:
        - Variant counts by tool: Compare variant callers
        - Category distribution: Somatic, Germline, Artifact, etc.
        - Stage progression: Track variants through pipeline stages
        - Quality distributions: Box plots of quality scores

    Workflow Comparison Plots (NEW):
        - RNA workflow comparison: Side-by-side RNA standard vs realignment
        - RNA stage progression: Line plots showing stage-to-stage changes
        - RNA annotation impact: Heatmap of annotation stage differences
        - Integrative view: Comprehensive DNA + RNA standard + RNA realignment

Color Scheme:
    - Somatic: Blue (#1f77b4)
    - Germline: Orange (#ff7f0e)
    - Reference: Green (#2ca02c)
    - Artifact: Red (#d62728)
    - RNA_Edit: Purple (#9467bd)
    - NoConsensus: Gray (#7f7f7f)

Design Principles:
    - Interactive: All plots use Plotly for interactivity
    - Consistent: Unified color scheme across all plots
    - Informative: Clear labels, legends, and annotations
    - Workflow-aware: Support for multiple workflow types
"""

from typing import Any, Dict, Optional

# Import visualization dependencies with error handling
try:
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import plotly.express as px
    import plotly.graph_objects as go
    import seaborn as sns
    from plotly.subplots import make_subplots

    VISUALIZATION_AVAILABLE = True
except ImportError as e:
    print(f"Warning: Visualization dependencies not available: {e}")
    VISUALIZATION_AVAILABLE = False

# Import constants and utilities from main module
from . import CATEGORY_COLORS, CATEGORY_ORDER, STAGE_DISPLAY_NAMES, VCF_STAGE_ORDER
from .plot_utils import (
    build_legend_tracker,
    heatmap_matrix,
    legend_config,
    percentile_cap,
    should_add_to_legend,
)
from .utils import sort_stages


class VCFVisualizer:
    """
    Create interactive visualizations for VCF statistics.

    This class provides comprehensive visualization capabilities for VCF statistics
    across all processing stages and workflow types. It supports both standard and
    realignment workflows, enabling detailed comparison and analysis.

    Key Capabilities:
        - Variant count visualizations by tool and modality
        - Category distribution plots
        - Stage progression tracking
        - Quality score distributions
        - Workflow comparison plots (standard vs realignment)
        - RNA-focused comparison visualizations
        - Integrative view plots

    Workflow Support:
        - Standard workflow: DNA + RNA samples
        - Realignment workflow: RNA samples only
        - Comparison: RNA standard vs RNA realignment

    Usage Example:
        >>> # Create visualizer
        >>> visualizer = VCFVisualizer(
        ...     all_stats=workflow_stats,
        ...     workflow_type="standard"
        ... )
        >>>
        >>> # Generate standard plots
        >>> fig1 = visualizer.plot_variant_counts_by_tool()
        >>> fig1.show()
        >>>
        >>> fig2 = visualizer.plot_category_distribution()
        >>> fig2.show()
        >>>
        >>> # Generate workflow comparison (if realignment available)
        >>> fig3 = visualizer.plot_rna_workflow_comparison(
        ...     rna_standard_stats, rna_realignment_stats
        ... )
        >>> fig3.show()

    Attributes:
        all_stats: Dictionary containing all VCF statistics
        workflow_type: Type of workflow ("standard" or "realignment")
        CATEGORY_COLORS: Color scheme for variant categories
    """

    def __init__(self, all_stats: Dict[str, Any], workflow_type: str = "standard"):
        """
        Initialize visualizer with statistics data.

        Args:
            all_stats: Dictionary containing all VCF statistics
                Format: {stage: {sample_name: {stats_dict}}}
            workflow_type: Type of workflow ("standard" or "realignment")
                Used for labeling plots and organizing visualizations

        Example:
            >>> stats = {
            ...     'filtered_rescue': {
            ...         'DNA_TUMOR_vs_DNA_NORMAL': {
            ...             'stats': {
            ...                 'basic': {'total_variants': 1234, ...}
            ...             }
            ...         }
            ...     }
            ... }
            >>> visualizer = VCFVisualizer(stats, workflow_type="standard")
        """
        self.all_stats = all_stats
        self.workflow_type = workflow_type

        # Reuse unified color scheme exported by the package
        self.CATEGORY_COLORS = CATEGORY_COLORS.copy()

    def plot_variant_counts_by_tool(self):
        """
        Bar plot comparing variant counts with FILTER categories across tools and modalities.

        Note: Does not create NoConsensus-free view since normalized caller VCFs
        do not have consensus/rescue classification.
        """
        if not VISUALIZATION_AVAILABLE:
            print("Visualization libraries not available. Skipping plot.")
            return None

        # Collect data - use normalized caller VCFs (variant_calling removed)
        data = []

        if "normalized" in self.all_stats:
            for name, vcf_data in self.all_stats["normalized"].items():
                if "stats" in vcf_data and "basic" in vcf_data["stats"]:
                    basic = vcf_data["stats"]["basic"]
                    classification = basic.get("classification", {})
                    parts = name.split("_")
                    tool = parts[0] if parts else name
                    modality = (
                        "DNA" if "DT_vs_DN" in name or "DNA_TUMOR" in name else "RNA"
                    )

                    # Add data for each FILTER category
                    for filter_cat in CATEGORY_ORDER:
                        count = classification.get(filter_cat, 0)
                        if count > 0:
                            data.append(
                                {
                                    "Tool": tool,
                                    "Modality": modality,
                                    "Category": filter_cat,
                                    "Count": count,
                                }
                            )

        if not data:
            print("No data available for plotting")
            return None

        df = pd.DataFrame(data)

        # Create subplots for DNA and RNA
        fig = make_subplots(
            rows=1,
            cols=2,
            subplot_titles=("DNA Modality", "RNA Modality"),
            horizontal_spacing=0.12,
            shared_yaxes=True,
        )

        # Plot DNA modality
        df_dna = df[df["Modality"] == "DNA"]
        if not df_dna.empty:
            tools = sorted(df_dna["Tool"].unique())
            for filter_cat in CATEGORY_ORDER:
                df_cat = df_dna[df_dna["Category"] == filter_cat]
                counts = [
                    df_cat[df_cat["Tool"] == t]["Count"].sum()
                    if not df_cat[df_cat["Tool"] == t].empty
                    else 0
                    for t in tools
                ]
                if sum(counts) > 0:
                    fig.add_trace(
                        go.Bar(
                            name=filter_cat,
                            x=tools,
                            y=counts,
                            marker_color=self.CATEGORY_COLORS.get(
                                filter_cat, "#8A8A8A"
                            ),
                            text=counts,
                            textposition="inside",
                            showlegend=True,
                            legendgroup=filter_cat,
                        ),
                        row=1,
                        col=1,
                    )

        # Plot RNA modality
        df_rna = df[df["Modality"] == "RNA"]
        if not df_rna.empty:
            tools = sorted(df_rna["Tool"].unique())
            for filter_cat in CATEGORY_ORDER:
                df_cat = df_rna[df_rna["Category"] == filter_cat]
                counts = [
                    df_cat[df_cat["Tool"] == t]["Count"].sum()
                    if not df_cat[df_cat["Tool"] == t].empty
                    else 0
                    for t in tools
                ]
                if sum(counts) > 0:
                    fig.add_trace(
                        go.Bar(
                            name=filter_cat,
                            x=tools,
                            y=counts,
                            marker_color=self.CATEGORY_COLORS.get(
                                filter_cat, "#8A8A8A"
                            ),
                            text=counts,
                            textposition="inside",
                            showlegend=False,
                            legendgroup=filter_cat,
                        ),
                        row=1,
                        col=2,
                    )

        fig.update_xaxes(title_text="Tool", row=1, col=1)
        fig.update_xaxes(title_text="Tool", row=1, col=2)
        fig.update_yaxes(title_text="Number of Variants", row=1, col=1)
        fig.update_yaxes(title_text="Number of Variants", row=1, col=2)

        fig.update_layout(
            title="Variant Counts by Tool and FILTER Category",
            template="plotly_white",
            barmode="stack",
            height=500,
            showlegend=True,
        )

        return fig

    def plot_variant_type_distribution(
        self, exclude_no_consensus_view: bool = True, debug_rna_edit: bool = False
    ):
        """
        Stacked bar charts showing SNP vs INDEL distribution with FILTER categories.

        Shows variant type breakdown across all pipeline stages including:
        - DNA Consensus, RNA Consensus, Rescue
        - COSMIC/GnomAD, RNA Editing, Filtered

        Args:
            exclude_no_consensus_view: If True, creates second plot excluding NoConsensus
            debug_rna_edit: If True, prints debug info about RNA_Edit variant type distribution
        """
        if not VISUALIZATION_AVAILABLE:
            print("Visualization libraries not available. Skipping plot.")
            return None

        data = []
        debug_info = []  # Collect debug information

        # Stage mapping for consistent naming
        stage_mapping = {
            "consensus_DNA_TUMOR": "DNA Consensus",
            "consensus_RNA_TUMOR": "RNA Consensus",
            "rescue": "Rescue",
            "cosmic_gnomad": "COSMIC/GnomAD",
            "rna_editing": "RNA Editing",
            "filtered_rescue": "Filtered",
        }

        # Collect data from all annotation pipeline stages
        for vcf_type in [
            "consensus",
            "rescue",
            "cosmic_gnomad",
            "rna_editing",
            "filtered_rescue",
        ]:
            if vcf_type in self.all_stats:
                for name, vcf_data in self.all_stats[vcf_type].items():
                    if "stats" in vcf_data and "basic" in vcf_data["stats"]:
                        basic = vcf_data["stats"]["basic"]
                        classification = basic.get("classification", {})

                        # Determine stage name
                        if vcf_type == "consensus":
                            if "DNA_TUMOR" in name:
                                stage_name = "DNA Consensus"
                            elif "RNA_TUMOR" in name:
                                stage_name = "RNA Consensus"
                            else:
                                continue
                        else:
                            stage_name = STAGE_DISPLAY_NAMES.get(
                                vcf_type, vcf_type.title()
                            )

                        # Get variant types
                        variant_types = basic.get("variant_types", {})
                        snps = variant_types.get("SNP", 0)
                        indels = variant_types.get("DEL", 0) + variant_types.get(
                            "INS", 0
                        )

                        # For each FILTER category, calculate proportional SNP/INDEL split
                        total_vars = basic.get("total_variants", 1)
                        for filter_cat in CATEGORY_ORDER:
                            count = classification.get(filter_cat, 0)
                            if count > 0:
                                # Proportionally split into SNPs and INDELs
                                snp_count = int(count * (snps / total_vars))
                                indel_count = count - snp_count

                                if snp_count > 0:
                                    data.append(
                                        {
                                            "Stage": stage_name,
                                            "Type": "SNP",
                                            "Category": filter_cat,
                                            "Count": snp_count,
                                        }
                                    )
                                if indel_count > 0:
                                    data.append(
                                        {
                                            "Stage": stage_name,
                                            "Type": "INDEL",
                                            "Category": filter_cat,
                                            "Count": indel_count,
                                        }
                                    )

                                # Collect debug info for RNA_Edit
                                if debug_rna_edit and filter_cat == "RNA_Edit":
                                    debug_info.append(
                                        {
                                            "stage": stage_name,
                                            "total_variants": total_vars,
                                            "total_snps": snps,
                                            "total_indels": indels,
                                            "rna_edit_count": count,
                                            "rna_edit_snps": snp_count,
                                            "rna_edit_indels": indel_count,
                                        }
                                    )

        # Debug output for RNA_Edit INDEL handling
        if debug_rna_edit and debug_info:
            print("\n" + "=" * 80)
            print("DEBUG: RNA_Edit Variant Type Distribution")
            print("=" * 80)
            for info in debug_info:
                print(f"\n{info['stage']}:")
                print(f"  Total variants: {info['total_variants']:,}")
                print(
                    f"  SNPs (total): {info['total_snps']:,} ({info['total_snps'] / info['total_variants'] * 100:.1f}%)"
                )
                print(
                    f"  INDELs (total): {info['total_indels']:,} ({info['total_indels'] / info['total_variants'] * 100:.1f}%)"
                )
                print(f"  RNA_Edit count: {info['rna_edit_count']:,}")
                print(
                    f"    - SNPs: {info['rna_edit_snps']:,} ({info['rna_edit_snps'] / info['rna_edit_count'] * 100:.1f}%)"
                )
                print(
                    f"    - INDELs: {info['rna_edit_indels']:,} ({info['rna_edit_indels'] / info['rna_edit_count'] * 100:.1f}%)"
                )
            print("=" * 80 + "\n")

        if not data:
            print("No data available for plotting")
            return None

        df = pd.DataFrame(data)

        def _create_variant_type_plot(plot_df, title_suffix=""):
            """Helper to create variant type distribution plot."""
            # Determine stages present in data and sort by pipeline order
            stage_order = [
                "DNA Consensus",
                "RNA Consensus",
                "Rescue",
                "COSMIC/GnomAD",
                "RNA Editing",
                "Filtered",
            ]
            available_stages = [s for s in stage_order if s in plot_df["Stage"].values]
            n_stages = len(available_stages)

            if n_stages == 0:
                return None

            fig = make_subplots(
                rows=1,
                cols=n_stages,
                subplot_titles=available_stages,
                horizontal_spacing=0.08,
                shared_yaxes=True,
            )

            # Track legend to avoid duplication
            categories_in_legend = set()

            for i, stage in enumerate(available_stages, 1):
                df_stage = plot_df[plot_df["Stage"] == stage]

                # Group by Type (SNP/INDEL) and stack by Category
                for var_type in ["SNP", "INDEL"]:
                    df_type = df_stage[df_stage["Type"] == var_type]
                    if not df_type.empty:
                        for filter_cat in CATEGORY_ORDER:
                            df_cat = df_type[df_type["Category"] == filter_cat]
                            if not df_cat.empty:
                                count = df_cat["Count"].sum()

                                # Show in legend only once per category
                                show_legend = filter_cat not in categories_in_legend
                                if show_legend:
                                    categories_in_legend.add(filter_cat)

                                fig.add_trace(
                                    go.Bar(
                                        name=filter_cat,
                                        x=[var_type],
                                        y=[count],
                                        marker_color=self.CATEGORY_COLORS.get(
                                            filter_cat, "#8A8A8A"
                                        ),
                                        text=[count],
                                        textposition="inside",
                                        showlegend=show_legend,
                                        legendgroup=filter_cat,
                                        hovertemplate=f"<b>{var_type}</b><br>{filter_cat}: %{{y}}<extra></extra>",
                                    ),
                                    row=1,
                                    col=i,
                                )

                fig.update_xaxes(title_text="Type", row=1, col=i)
                if i == 1:
                    fig.update_yaxes(title_text="Number of Variants", row=1, col=i)

            fig.update_layout(
                title_text=f"Variant Type Distribution (SNP vs INDEL) Across Pipeline Stages{title_suffix}",
                height=500,
                barmode="stack",
                template="plotly_white",
                showlegend=True,
            )
            return fig

        # Create main plot with all categories
        full_fig = _create_variant_type_plot(df)

        # Create NoConsensus-free plot if requested
        if exclude_no_consensus_view:
            df_filtered = df[df["Category"] != "NoConsensus"].copy()
            if not df_filtered.empty:
                filtered_fig = _create_variant_type_plot(
                    df_filtered, " (excluding NoConsensus)"
                )

    def plot_consensus_comparison(self):
        """
        Compare consensus variants to individual tools with FILTER category breakdown.
        """
        if not VISUALIZATION_AVAILABLE:
            print("Visualization libraries not available. Skipping plot.")
            return None

        data = []

        # Get tool-level classification data from normalized VCFs
        if "normalized" in self.all_stats:
            for name, vcf_data in self.all_stats["normalized"].items():
                if "stats" in vcf_data and "basic" in vcf_data["stats"]:
                    basic = vcf_data["stats"]["basic"]
                    classification = basic.get("classification", {})
                    parts = name.split("_")
                    tool = parts[0] if parts else name
                    modality = (
                        "DNA" if "DT_vs_DN" in name or "DNA_TUMOR" in name else "RNA"
                    )

                    for filter_cat in CATEGORY_ORDER:
                        count = classification.get(filter_cat, 0)
                        if count > 0:
                            data.append(
                                {
                                    "Category": tool,
                                    "Modality": modality,
                                    "FilterCat": filter_cat,
                                    "Count": count,
                                }
                            )

        # Get consensus classification data
        if "consensus" in self.all_stats:
            for modality_key, vcf_data in self.all_stats["consensus"].items():
                if "stats" in vcf_data and "basic" in vcf_data["stats"]:
                    basic = vcf_data["stats"]["basic"]
                    classification = basic.get("classification", {})

                    # Explicitly detect modality using multiple patterns
                    modality_key_upper = modality_key.upper()
                    if any(
                        pattern in modality_key_upper
                        for pattern in ["RT_VS_DN", "RNA_TUMOR", "RNA_VS", "_RT_"]
                    ):
                        modality = "RNA"
                    elif any(
                        pattern in modality_key_upper
                        for pattern in ["DT_VS_DN", "DNA_TUMOR", "DNA_VS", "_DT_"]
                    ):
                        modality = "DNA"
                    else:
                        # Fallback: if contains RNA, it's RNA, otherwise DNA
                        modality = "RNA" if "RNA" in modality_key_upper else "DNA"

                    for filter_cat in CATEGORY_ORDER:
                        count = classification.get(filter_cat, 0)
                        if count > 0:
                            data.append(
                                {
                                    "Category": "consensus",
                                    "Modality": modality,
                                    "FilterCat": filter_cat,
                                    "Count": count,
                                }
                            )

        if not data:
            print("No comparison data available")
            return None

        df = pd.DataFrame(data)

        # Create subplots for DNA and RNA modalities
        fig = make_subplots(
            rows=1,
            cols=2,
            subplot_titles=("DNA Modality", "RNA Modality"),
            horizontal_spacing=0.12,
            shared_yaxes=True,
        )

        # Get all unique categories and sort them (tools alphabetically, consensus last)
        categories = sorted([c for c in df["Category"].unique() if c != "consensus"])
        if "consensus" in df["Category"].unique():
            categories.append("consensus")

        # DNA modality
        df_dna = df[df["Modality"] == "DNA"]
        if not df_dna.empty:
            for filter_cat in CATEGORY_ORDER:
                df_filter = df_dna[df_dna["FilterCat"] == filter_cat]
                if not df_filter.empty:
                    counts = [
                        df_filter[df_filter["Category"] == cat]["Count"].sum()
                        if not df_filter[df_filter["Category"] == cat].empty
                        else 0
                        for cat in categories
                    ]
                    fig.add_trace(
                        go.Bar(
                            name=filter_cat,
                            x=categories,
                            y=counts,
                            marker_color=self.CATEGORY_COLORS.get(
                                filter_cat, "#8A8A8A"
                            ),
                            text=counts,
                            textposition="inside",
                            showlegend=True,
                            legendgroup=filter_cat,
                        ),
                        row=1,
                        col=1,
                    )

        # RNA modality
        df_rna = df[df["Modality"] == "RNA"]
        if not df_rna.empty:
            for filter_cat in CATEGORY_ORDER:
                df_filter = df_rna[df_rna["FilterCat"] == filter_cat]
                if not df_filter.empty:
                    counts = [
                        df_filter[df_filter["Category"] == cat]["Count"].sum()
                        if not df_filter[df_filter["Category"] == cat].empty
                        else 0
                        for cat in categories
                    ]
                    fig.add_trace(
                        go.Bar(
                            name=filter_cat,
                            x=categories,
                            y=counts,
                            marker_color=self.CATEGORY_COLORS.get(
                                filter_cat, "#8A8A8A"
                            ),
                            text=counts,
                            textposition="inside",
                            showlegend=False,
                            legendgroup=filter_cat,
                        ),
                        row=1,
                        col=2,
                    )

        fig.update_xaxes(title_text="Category", row=1, col=1)
        fig.update_xaxes(title_text="Category", row=1, col=2)
        fig.update_yaxes(title_text="Number of Variants", row=1, col=1)

        fig.update_layout(
            title_text="Variant Counts by Category with FILTER Classification",
            template="plotly_white",
            height=500,
            barmode="stack",
        )

        return fig

    def plot_filter_status(self, dual_view_no_consensus: bool = True):
        """
        Stacked bar chart showing all 6 FILTER categories:
        Somatic, Germline, Reference, Artifact, RNA_Edit, NoConsensus.

        Ensures all categories appear in legend with their assigned colors.
        """
        if not VISUALIZATION_AVAILABLE:
            print("Visualization libraries not available. Skipping plot.")
            return None

        # Collect data from all VCF stats
        data = []

        for category, files in self.all_stats.items():
            for name, vcf_data in files.items():
                if "stats" not in vcf_data or "basic" not in vcf_data["stats"]:
                    continue

                basic = vcf_data["stats"]["basic"]
                classification = basic.get("classification", {})

                if not classification:
                    continue

                # Determine subplot category and tool name
                parts = name.split("_")

                if category == "normalized":
                    tool = parts[0] if parts else "unknown"
                    if "DT_vs_DN" in name or "DNA_TUMOR" in name:
                        subplot_cat = "DNA"
                    elif "RT_vs_DN" in name or "RNA_TUMOR" in name:
                        subplot_cat = "RNA"
                    else:
                        continue

                elif category == "consensus":
                    tool = "consensus"
                    if "DNA_TUMOR" in name:
                        subplot_cat = "DNA"
                    elif "RNA_TUMOR" in name:
                        subplot_cat = "RNA"
                    else:
                        continue

                elif category in {
                    "rescue",
                    "cosmic_gnomad",
                    "rna_editing",
                    "filtered_rescue",
                }:
                    tool = category
                    # Show these together in a single subplot
                    subplot_cat = "Annotation Stages"

                else:
                    continue

                # Add counts for each classification category
                for filter_cat in CATEGORY_ORDER:
                    count = classification.get(filter_cat, 0)
                    if count > 0:  # Only add if there are variants
                        data.append(
                            {
                                "SubplotCategory": subplot_cat,
                                "Tool": tool,
                                "FilterCategory": filter_cat,
                                "Count": count,
                                "Name": name,
                            }
                        )

        if not data:
            print("No classification data available for plotting")
            return None

        df = pd.DataFrame(data)

        # Determine available subplot categories
        subplot_categories = ["DNA", "RNA", "Annotation Stages"]
        available_subplots = [
            cat for cat in subplot_categories if cat in df["SubplotCategory"].values
        ]
        n_subplots = len(available_subplots)

        if n_subplots == 0:
            print("No data to plot")
            return None

        # Create subplots with unified y-axis
        fig = make_subplots(
            rows=1,
            cols=n_subplots,
            subplot_titles=available_subplots,
            horizontal_spacing=0.10,
            shared_yaxes=True,
        )

        # Find global y-axis max for consistent scaling
        max_y = 0

        # Track categories added to legend globally
        categories_seen = build_legend_tracker()

        # Plot each subplot
        for i, subplot_cat in enumerate(available_subplots, 1):
            df_subplot = df[df["SubplotCategory"] == subplot_cat]

            # Get unique tools in this subplot and sort them
            tools = df_subplot["Tool"].unique()

            # Sort annotation stages by pipeline order, others alphabetically
            if subplot_cat == "Annotation Stages":
                tools = sort_stages(list(tools))
            else:
                tools = sorted(tools)

            # For each filter category, add a stacked bar
            for filter_cat in CATEGORY_ORDER:
                df_filter = df_subplot[df_subplot["FilterCategory"] == filter_cat]

                # Create counts array aligned with tools
                counts = []
                for tool in tools:
                    tool_data = df_filter[df_filter["Tool"] == tool]
                    count = tool_data["Count"].sum() if not tool_data.empty else 0
                    counts.append(count)

                # Only add trace if there are non-zero counts
                if sum(counts) > 0:
                    fig.add_trace(
                        go.Bar(
                            name=filter_cat,
                            x=tools,
                            y=counts,
                            marker_color=self.CATEGORY_COLORS.get(
                                filter_cat, "#8A8A8A"
                            ),
                            text=counts,
                            textposition="inside",
                            textfont=dict(color="white", size=10),
                            showlegend=should_add_to_legend(
                                categories_seen, filter_cat
                            ),
                            legendgroup=filter_cat,
                            hovertemplate=f"<b>%{{x}}</b><br>{filter_cat}: %{{y}}<extra></extra>",
                        ),
                        row=1,
                        col=i,
                    )

            # Calculate total height for this subplot
            for tool in tools:
                tool_data = df_subplot[df_subplot["Tool"] == tool]
                total = tool_data["Count"].sum()
                max_y = max(max_y, total)

            # Update axes labels
            fig.update_xaxes(title_text="Caller", row=1, col=i)
            if i == 1:
                fig.update_yaxes(title_text="Number of Variants", row=1, col=i)

        # Update layout with unified settings
        fig.update_layout(
            title_text="Variant Classification by FILTER Category (Stacked)",
            template="plotly_white",
            height=500,
            barmode="stack",
            showlegend=True,
            legend=legend_config("right"),
        )

        # Set unified y-axis range with some padding
        for i in range(1, n_subplots + 1):
            fig.update_yaxes(range=[0, max_y * 1.1], row=1, col=i)

        # Optional secondary view excluding NoConsensus to improve minor category visibility
        if dual_view_no_consensus and data:
            df_small = df[df["FilterCategory"] != "NoConsensus"].copy()
            if not df_small.empty:
                fig_small = make_subplots(
                    rows=1,
                    cols=len(available_subplots),
                    subplot_titles=[
                        f"{t} (No NoConsensus)" for t in available_subplots
                    ],
                    horizontal_spacing=0.10,
                    shared_yaxes=True,
                )

                max_y2 = 0
                # Track categories added to legend for second plot
                categories_seen_small = build_legend_tracker()

                for i, subplot_cat in enumerate(available_subplots, 1):
                    df_subplot = df_small[df_small["SubplotCategory"] == subplot_cat]
                    tools = df_subplot["Tool"].unique()

                    # Sort annotation stages by pipeline order
                    if subplot_cat == "Annotation Stages":
                        tools = sort_stages(list(tools))
                    else:
                        tools = sorted(tools)

                    for filter_cat in [c for c in CATEGORY_ORDER if c != "NoConsensus"]:
                        df_filter = df_subplot[
                            df_subplot["FilterCategory"] == filter_cat
                        ]
                        counts = []
                        for tool in tools:
                            tool_data = df_filter[df_filter["Tool"] == tool]
                            count = (
                                tool_data["Count"].sum() if not tool_data.empty else 0
                            )
                            counts.append(count)
                        if sum(counts) > 0:
                            fig_small.add_trace(
                                go.Bar(
                                    name=filter_cat,
                                    x=tools,
                                    y=counts,
                                    marker_color=self.CATEGORY_COLORS.get(
                                        filter_cat, "#8A8A8A"
                                    ),
                                    text=counts,
                                    textposition="inside",
                                    textfont=dict(color="white", size=10),
                                    showlegend=should_add_to_legend(
                                        categories_seen_small, filter_cat
                                    ),
                                    legendgroup=filter_cat,
                                ),
                                row=1,
                                col=i,
                            )
                    for tool in tools:
                        tool_data = df_subplot[df_subplot["Tool"] == tool]
                        total = tool_data["Count"].sum()
                        max_y2 = max(max_y2, total)

                    fig_small.update_xaxes(title_text="Caller/Stage", row=1, col=i)
                    if i == 1:
                        fig_small.update_yaxes(
                            title_text="Number of Variants", row=1, col=i
                        )

                fig_small.update_layout(
                    title_text="Variant Classification (excluding NoConsensus)",
                    template="plotly_white",
                    height=500,
                    barmode="stack",
                    showlegend=True,
                    legend=legend_config("right"),
                )
                for i in range(1, len(available_subplots) + 1):
                    fig_small.update_yaxes(range=[0, max_y2 * 1.1], row=1, col=i)
                fig_small.show()

    def plot_annotation_progression(self, dual_view_no_consensus: bool = True):
        """
        Plot variant count progression through annotation pipeline stages:
        rescue → cosmic_gnomad → rna_editing → filtered_rescue

        Note: Excludes 'normalized' and consensus stages (not part of annotation pipeline).
        Includes an optional secondary view that hides NoConsensus to
        highlight smaller categories when NoConsensus dominates.
        """
        if not VISUALIZATION_AVAILABLE:
            print("Visualization libraries not available. Skipping plot.")
            return None

        data = []

        # Only include annotation pipeline stages, not normalized or consensus
        annotation_stages = [
            "rescue",
            "cosmic_gnomad",
            "rna_editing",
            "filtered_rescue",
        ]

        # Collect data from each annotation stage in order
        for stage in annotation_stages:
            if stage not in self.all_stats:
                continue

            for name, vcf_data in self.all_stats[stage].items():
                if "stats" not in vcf_data or "basic" not in vcf_data["stats"]:
                    continue

                basic = vcf_data["stats"]["basic"]
                classification = basic.get("classification", {})
                total_variants = basic.get("total_variants", 0)

                # For each category, record count at this stage
                for category in CATEGORY_ORDER:
                    count = classification.get(category, 0)
                    if count > 0 or stage == "rescue":  # Always include rescue stage
                        data.append(
                            {
                                "Stage": stage,
                                "Stage_Order": annotation_stages.index(stage),
                                "Category": category,
                                "Count": count,
                                "Total": total_variants,
                            }
                        )

        if not data:
            print("No annotation stage data available for progression plot")
            return None

        df = pd.DataFrame(data)

        def _build_progression(fig_df, title_suffix=""):
            fig = go.Figure()

            stages_ordered = sorted(
                fig_df["Stage"].unique(),
                key=lambda x: annotation_stages.index(x)
                if x in annotation_stages
                else 999,
            )

            categories_shown = set()

            for category in CATEGORY_ORDER:
                df_cat = fig_df[fig_df["Category"] == category]
                if df_cat.empty:
                    continue

                counts_by_stage = []
                for stage in stages_ordered:
                    stage_data = df_cat[df_cat["Stage"] == stage]
                    counts_by_stage.append(
                        stage_data["Count"].sum() if not stage_data.empty else 0
                    )

                if sum(counts_by_stage) > 0:
                    show_legend = category not in categories_shown
                    categories_shown.add(category)

                    fig.add_trace(
                        go.Bar(
                            name=category,
                            x=stages_ordered,
                            y=counts_by_stage,
                            marker_color=self.CATEGORY_COLORS.get(category, "#8A8A8A"),
                            text=counts_by_stage,
                            textposition="inside",
                            textfont=dict(color="white", size=10),
                            showlegend=show_legend,
                            legendgroup=category,
                            hovertemplate=f"<b>%{{x}}</b><br>{category}: %{{y}}<extra></extra>",
                        )
                    )

            fig.update_layout(
                title_text=f"Variant Count Progression Through Annotation Stages{title_suffix}",
                xaxis_title="Processing Stage",
                yaxis_title="Number of Variants",
                barmode="stack",
                template="plotly_white",
                height=500,
                showlegend=True,
                legend=dict(
                    orientation="v",
                    yanchor="top",
                    y=1.0,
                    xanchor="left",
                    x=1.02,
                    title="Category",
                ),
            )
            return fig

        # Full view including NoConsensus
        full_fig = _build_progression(df)

        # Optional view excluding NoConsensus to improve readability
        if dual_view_no_consensus:
            df_small = df[df["Category"] != "NoConsensus"].copy()
            if not df_small.empty:
                small_fig = _build_progression(df_small, " (excluding NoConsensus)")

    def plot_category_heatmap(
        self, stages: Optional[list] = None, dual_view_without_no_consensus: bool = True
    ):
        """
        Heatmap showing variant counts across stages and categories.

        Creates dual heatmaps: one with all categories (including NoConsensus),
        one without NoConsensus to improve visibility of other categories.
        Uses per-stage relative coloring to handle NoConsensus domination.

        Args:
            stages: List of stages to include (defaults to all annotation pipeline stages)
            dual_view_without_noconse nsus: If True, creates secondary view without NoConsensus
        """
        if not VISUALIZATION_AVAILABLE:
            print("Visualization libraries not available. Skipping plot.")
            return None

        if stages is None:
            stages = ["rescue", "cosmic_gnomad", "rna_editing", "filtered_rescue"]

        # Collect data
        data = []
        for stage in stages:
            if stage in self.all_stats:
                for name, vcf_data in self.all_stats[stage].items():
                    if "stats" in vcf_data and "basic" in vcf_data["stats"]:
                        basic = vcf_data["stats"]["basic"]
                        classification = basic.get("classification", {})

                        stage_name = STAGE_DISPLAY_NAMES.get(stage, stage.title())

                        for category in CATEGORY_ORDER:
                            count = classification.get(category, 0)
                            data.append(
                                {
                                    "Stage": stage_name,
                                    "Category": category,
                                    "Count": count,
                                }
                            )

        if not data:
            print("No data available for heatmap")
            return None

        df = pd.DataFrame(data)

        def _create_heatmap(plot_df, title_suffix="", normalize: Optional[str] = None):
            """Helper to create heatmap visualization (counts or percent per stage)."""
            # Build matrix and sort axes
            mat = heatmap_matrix(plot_df, normalize=normalize)
            # Sort columns by stage display order
            stage_names = [
                STAGE_DISPLAY_NAMES.get(s, s.title())
                for s in stages
                if STAGE_DISPLAY_NAMES.get(s, s.title()) in mat.columns
            ]
            mat = mat[stage_names]
            # Sort rows by CATEGORY_ORDER
            mat = mat.reindex([c for c in CATEGORY_ORDER if c in mat.index])

            zvals = mat.values
            zmax = None if normalize == "stage" else percentile_cap(zvals, 95.0)

            fig = go.Figure(
                data=go.Heatmap(
                    z=zvals,
                    x=mat.columns,
                    y=mat.index,
                    colorscale="YlGnBu",
                    zmax=zmax,
                    text=np.round(zvals, 1) if normalize == "stage" else zvals,
                    texttemplate="%{text}" if normalize == "stage" else "%{text:,}",
                    textfont={"size": 10},
                    colorbar=dict(
                        title="Percent" if normalize == "stage" else "Variant Count"
                    ),
                    hovertemplate=(
                        "<b>%{y}</b><br>%{x}<br>"
                        + (
                            "Percent: %{z:.1f}%"
                            if normalize == "stage"
                            else "Count: %{z:,}"
                        )
                        + "<extra></extra>"
                    ),
                    zauto=(normalize == "stage"),
                )
            )

            fig.update_layout(
                title=f"Variant Category Distribution Across Pipeline Stages{title_suffix}",
                xaxis_title="Processing Stage",
                yaxis_title="Category",
                height=500 if "NoConsensus" in plot_df["Category"].values else 400,
                template="plotly_white",
            )
            return fig

        # Full view (counts) with all categories
        fig_full = _create_heatmap(df)
        fig_full.show()

        # Percent-per-stage normalized view
        fig_pct = _create_heatmap(df, " (percent per stage)", normalize="stage")
        fig_pct.show()

        # Secondary views without NoConsensus to improve visibility
        if dual_view_without_no_consensus:
            df_no_nc = df[df["Category"] != "NoConsensus"].copy()
            if not df_no_nc.empty:
                fig_no_nc = _create_heatmap(df_no_nc, " (excluding NoConsensus)")
                fig_no_nc.show()
                fig_no_nc_pct = _create_heatmap(
                    df_no_nc,
                    " (percent per stage, excluding NoConsensus)",
                    normalize="stage",
                )
                fig_no_nc_pct.show()

    def plot_stage_transition_sankey(self, categories_to_show: Optional[list] = None):
        """
        Sankey diagram showing variant flow through pipeline stages.

        Visualizes how variants transition between stages, with flow thickness
        representing variant counts. Uses category colors for consistent styling.

        Args:
            categories_to_show: List of categories to include (defaults to [Somatic, Germline, RNA_Edit, NoConsensus])
        """
        if not VISUALIZATION_AVAILABLE:
            print("Visualization libraries not available. Skipping plot.")
            return None

        # Default to showing major categories with visual distinction
        if categories_to_show is None:
            categories_to_show = ["Somatic", "Germline", "RNA_Edit", "NoConsensus"]

        stages = ["rescue", "cosmic_gnomad", "rna_editing", "filtered_rescue"]

        # Collect stage statistics
        from .utils import collect_stage_statistics

        stage_stats = collect_stage_statistics(self.all_stats, stages)

        # Build nodes and links for Sankey
        nodes = []
        node_map = {}
        node_idx = 0

        # Create nodes for each stage-category combination
        for stage in stages:
            stage_name = STAGE_DISPLAY_NAMES.get(stage, stage.title())
            for category in categories_to_show:
                node_label = f"{stage_name}\n{category}"
                nodes.append(node_label)
                node_map[(stage, category)] = node_idx
                node_idx += 1

        # Create links between consecutive stages
        links = []
        for i in range(len(stages) - 1):
            stage_from = stages[i]
            stage_to = stages[i + 1]

            for category in categories_to_show:
                count_from = stage_stats[stage_from]["classification"].get(category, 0)
                count_to = stage_stats[stage_to]["classification"].get(category, 0)

                # Use minimum count for flow (represents variants that continue)
                flow_value = min(count_from, count_to)

                if flow_value > 0:
                    links.append(
                        {
                            "source": node_map[(stage_from, category)],
                            "target": node_map[(stage_to, category)],
                            "value": flow_value,
                            "category": category,
                            "label": f"{category}: {flow_value:,}",
                        }
                    )

        if not links:
            print("No transition data available for Sankey diagram")
            return None

        # Create node colors based on category colors
        node_colors = []
        for node_label in nodes:
            # Extract category from node label (format: "StageName\nCategoryName")
            parts = node_label.split("\n")
            if len(parts) > 1:
                category = parts[1]
                node_colors.append(self.CATEGORY_COLORS.get(category, "#8A8A8A"))
            else:
                node_colors.append("#D3D3D3")

        # Create link colors (semi-transparent version of source node color)
        link_colors = []
        for link in links:
            source_idx = link["source"]
            source_color = (
                node_colors[source_idx] if source_idx < len(node_colors) else "#8A8A8A"
            )
            # Convert hex to RGBA with transparency
            hex_color = source_color.lstrip("#")
            try:
                rgb = tuple(int(hex_color[i : i + 2], 16) for i in (0, 2, 4))
                link_colors.append(f"rgba({rgb[0]}, {rgb[1]}, {rgb[2]}, 0.5)")
            except:
                link_colors.append("rgba(128, 128, 128, 0.5)")

        # Create Sankey diagram with styled colors and improved layout
        fig = go.Figure(
            data=[
                go.Sankey(
                    node=dict(
                        pad=15,
                        thickness=20,
                        line=dict(color="darkgray", width=1),
                        label=nodes,
                        color=node_colors,
                    ),
                    link=dict(
                        source=[l["source"] for l in links],
                        target=[l["target"] for l in links],
                        value=[l["value"] for l in links],
                        color=link_colors,
                        label=[l["label"] for l in links],
                    ),
                )
            ]
        )

        fig.update_layout(
            title=f"Variant Flow Through Pipeline Stages ({', '.join(categories_to_show)})",
            height=650,
            font=dict(size=11, family="Arial"),
            template="plotly_white",
        )

    def plot_tier_distribution(
        self, rescue_vcf_path: Optional[str] = None, dual_view_no_consensus: bool = True
    ):
        """
        Bar chart showing tier distribution (C1D0-C7D1) per category for rescue VCF.

        NEW HYBRID TIERING SYSTEM:
        - Caller tiers (C1-C7): Based on DNA/RNA caller counts
          C1: ≥2 DNA + ≥2 RNA | C2: ≥2 DNA + (0-1) RNA | C3: ≥2 RNA + (0-1) DNA
          C4: 1 DNA + 1 RNA | C5: 1 DNA only | C6: 1 RNA only | C7: 0 DNA + 0 RNA

        - Database tiers (D0-D1):
          D1: Has database support | D0: No database support

        - Final tiers: CxDy format (14 total tiers from C1D1 to C7D0)

        Args:
            rescue_vcf_path: Path to rescue VCF (auto-detected if None)
            dual_view_no_consensus: If True, creates secondary plot excluding NoConsensus
        """
        if not VISUALIZATION_AVAILABLE:
            print("Visualization libraries not available. Skipping plot.")
            return None

        try:
            import sys
            from pathlib import Path as PathLib

            from .tiering import tier_rescue_variants
            from .tiering_engine import TieringEngine

            # Add bin/common to path
            bin_common = PathLib(__file__).parent.parent.parent / "bin" / "common"
            if str(bin_common) not in sys.path:
                sys.path.insert(0, str(bin_common))

            from tier_config import TIER_COLORS, TIER_ORDER, get_tier_display_name
        except ImportError as e:
            print(f"Required modules not available: {e}")
            return None

        # Auto-detect rescue VCF if not provided
        if rescue_vcf_path is None:
            if "rescue" in self.all_stats and self.all_stats["rescue"]:
                first_rescue = next(iter(self.all_stats["rescue"].values()))
                rescue_vcf_path = first_rescue.get("path")

        if not rescue_vcf_path:
            print("No rescue VCF found for tier distribution")
            return None

        # Get tiered data
        from pathlib import Path

        tiered_df = tier_rescue_variants(Path(rescue_vcf_path))

        if tiered_df.empty:
            print("No tiered data available")
            return None

        def _create_tier_plot(plot_df, title_suffix=""):
            """Helper function to create tier distribution plot with CxDy tiers."""
            # Create stacked bar chart data by tier
            data = []
            for category in CATEGORY_ORDER:
                cat_data = plot_df[plot_df["filter_category"] == category]
                if not cat_data.empty:
                    tier_counts = cat_data["tier"].value_counts().sort_index()
                    for tier, count in tier_counts.items():
                        data.append(
                            {"Category": category, "Tier": tier, "Count": count}
                        )

            if not data:
                return None

            df = pd.DataFrame(data)

            # Create stacked bar chart (stacked by category, x-axis = tier)
            fig = go.Figure()

            # Add traces for each category (will be stacked)
            for category in CATEGORY_ORDER:
                cat_df = df[df["Category"] == category]
                if not cat_df.empty:
                    # Ensure we have all tiers represented (fill with 0 if missing)
                    # Use TIER_ORDER from tier_config (C1D1, C1D0, C2D1, ...)
                    tier_counts = {tier: 0 for tier in TIER_ORDER}
                    for _, row in cat_df.iterrows():
                        if row["Tier"] in tier_counts:
                            tier_counts[row["Tier"]] = row["Count"]

                    fig.add_trace(
                        go.Bar(
                            name=category,
                            x=TIER_ORDER,
                            y=[tier_counts[t] for t in TIER_ORDER],
                            marker_color=self.CATEGORY_COLORS.get(category, "#8A8A8A"),
                            text=[
                                tier_counts[t] if tier_counts[t] > 0 else ""
                                for t in TIER_ORDER
                            ],
                            textposition="inside",
                            showlegend=True,
                            legendgroup=category,
                            hovertemplate="<b>%{x}</b><br>"
                            + category
                            + ": %{y}<extra></extra>",
                        )
                    )

            # Filter to show only tiers that have data
            existing_tiers = [t for t in TIER_ORDER if t in df["Tier"].values]

            fig.update_layout(
                title=f"Rescue VCF Tier Distribution by Category{title_suffix}<br><sub>CxDy format: Caller tier (C1-C7) + Database tier (D0-D1)</sub>",
                xaxis=dict(
                    title="Tier (Caller × Database)",
                    categoryorder="array",
                    categoryarray=existing_tiers if existing_tiers else TIER_ORDER,
                    tickangle=-45,
                ),
                yaxis_title="Number of Variants",
                barmode="stack",
                template="plotly_white",
                height=600,
                showlegend=True,
                legend=dict(
                    orientation="v",
                    yanchor="top",
                    y=1.0,
                    xanchor="left",
                    x=1.02,
                    title="Category",
                ),
            )
            return fig

        # Full view with all categories
        fig_full = _create_tier_plot(tiered_df)
        if fig_full:
            fig_full.show()

        # Secondary view excluding NoConsensus
        if dual_view_no_consensus:
            tiered_df_no_nc = tiered_df[
                tiered_df["filter_category"] != "NoConsensus"
            ].copy()
            if not tiered_df_no_nc.empty:
                fig_no_nc = _create_tier_plot(
                    tiered_df_no_nc, " (excluding NoConsensus)"
                )
                if fig_no_nc:
                    fig_no_nc.show()

    def plot_rna_workflow_comparison(
        self, rna_standard_stats: Dict[str, Any], rna_realignment_stats: Dict[str, Any]
    ):
        """
        Side-by-side comparison of RNA standard vs realignment workflows.

        Creates subplots showing:
        - RNA variant counts by stage
        - RNA category distribution
        - RNA SNP/INDEL ratios

        Note: Only compares RNA modality as realignment only applies to RNA.

        Args:
            rna_standard_stats: Statistics from standard RNA workflow
            rna_realignment_stats: Statistics from realignment RNA workflow
        """
        if not VISUALIZATION_AVAILABLE:
            print("Visualization libraries not available. Skipping plot.")
            return None

        # Import comparison module
        from .comparison import WorkflowComparator

        # Create comparator
        comparator = WorkflowComparator(rna_standard_stats, rna_realignment_stats)

        # Get comparison data
        variant_counts_df = comparator.compare_rna_variant_counts()
        category_dist_df = comparator.compare_rna_category_distribution()

        # Create figure with 3 subplots
        fig = make_subplots(
            rows=1,
            cols=3,
            subplot_titles=(
                "RNA Variant Counts by Stage",
                "RNA Category Distribution",
                "RNA SNP/INDEL Ratios",
            ),
            horizontal_spacing=0.12,
            specs=[[{"type": "bar"}, {"type": "bar"}, {"type": "bar"}]],
        )

        # Subplot 1: Variant counts by stage
        stages = variant_counts_df["Stage"].unique()
        standard_counts = variant_counts_df["Standard_RNA_Count"].values
        realignment_counts = variant_counts_df["Realignment_RNA_Count"].values

        fig.add_trace(
            go.Bar(
                name="Standard RNA",
                x=stages,
                y=standard_counts,
                marker_color="#4472C4",
                text=standard_counts,
                textposition="inside",
                showlegend=True,
            ),
            row=1,
            col=1,
        )

        fig.add_trace(
            go.Bar(
                name="Realignment RNA",
                x=stages,
                y=realignment_counts,
                marker_color="#ED7D31",
                text=realignment_counts,
                textposition="inside",
                showlegend=True,
            ),
            row=1,
            col=1,
        )

        # Subplot 2: Category distribution (stacked bars for final stage)
        final_stage = "filtered_rescue"
        final_stage_data = category_dist_df[category_dist_df["Stage"] == final_stage]

        categories_in_legend = set()

        for category in CATEGORY_ORDER:
            cat_data = final_stage_data[final_stage_data["Category"] == category]
            if not cat_data.empty:
                standard_count = cat_data["Standard_RNA_Count"].values[0]
                realignment_count = cat_data["Realignment_RNA_Count"].values[0]

                show_legend = category not in categories_in_legend
                categories_in_legend.add(category)

                fig.add_trace(
                    go.Bar(
                        name=category,
                        x=["Standard", "Realignment"],
                        y=[standard_count, realignment_count],
                        marker_color=self.CATEGORY_COLORS.get(category, "#8A8A8A"),
                        text=[standard_count, realignment_count],
                        textposition="inside",
                        showlegend=show_legend,
                        legendgroup=category,
                    ),
                    row=1,
                    col=2,
                )

        # Subplot 3: SNP/INDEL ratios
        # Calculate SNP/INDEL ratios from final stage
        def _get_snp_indel_ratio(stats, stage):
            """Helper to calculate SNP/INDEL ratio for a stage."""
            if stage not in stats:
                return 0, 0

            total_snps = 0
            total_indels = 0

            for name, data in stats[stage].items():
                if "stats" in data and "basic" in data["stats"]:
                    basic = data["stats"]["basic"]
                    variant_types = basic.get("variant_types", {})
                    total_snps += variant_types.get("SNP", 0)
                    total_indels += variant_types.get("DEL", 0) + variant_types.get(
                        "INS", 0
                    )

            return total_snps, total_indels

        std_snps, std_indels = _get_snp_indel_ratio(rna_standard_stats, final_stage)
        real_snps, real_indels = _get_snp_indel_ratio(
            rna_realignment_stats, final_stage
        )

        fig.add_trace(
            go.Bar(
                name="SNPs",
                x=["Standard", "Realignment"],
                y=[std_snps, real_snps],
                marker_color="#70AD47",
                text=[std_snps, real_snps],
                textposition="inside",
                showlegend=True,
            ),
            row=1,
            col=3,
        )

        fig.add_trace(
            go.Bar(
                name="INDELs",
                x=["Standard", "Realignment"],
                y=[std_indels, real_indels],
                marker_color="#FFC000",
                text=[std_indels, real_indels],
                textposition="inside",
                showlegend=True,
            ),
            row=1,
            col=3,
        )

        # Update axes
        fig.update_xaxes(title_text="Stage", row=1, col=1, tickangle=-45)
        fig.update_xaxes(title_text="Workflow", row=1, col=2)
        fig.update_xaxes(title_text="Workflow", row=1, col=3)

        fig.update_yaxes(title_text="Variant Count", row=1, col=1)
        fig.update_yaxes(title_text="Variant Count", row=1, col=2)
        fig.update_yaxes(title_text="Variant Count", row=1, col=3)

        # Update layout
        fig.update_layout(
            title_text="RNA Workflow Comparison: Standard vs Realignment",
            template="plotly_white",
            height=600,
            barmode="stack",
            showlegend=True,
        )

        return fig

    def plot_rna_stage_progression_comparison(
        self, rna_standard_stats: Dict[str, Any], rna_realignment_stats: Dict[str, Any]
    ):
        """
        Line plot showing RNA variant count changes through stages.

        Two lines (standard RNA and realignment RNA) showing how counts
        change from normalized → filtered_rescue.

        Args:
            rna_standard_stats: Statistics from standard RNA workflow (already filtered to RNA-only)
            rna_realignment_stats: Statistics from realignment RNA workflow (RNA-only)
        """
        if not VISUALIZATION_AVAILABLE:
            print("Visualization libraries not available. Skipping plot.")
            return None

        # Aggregate counts directly from the pre-filtered RNA stats
        def aggregate_total_variants(stats: Dict[str, Any], stage: str) -> int:
            """Aggregate total variants for a stage."""
            total = 0
            if stage in stats:
                stage_data = stats[stage]
                # Check if this is the expected nested structure
                if isinstance(stage_data, dict):
                    for name, data in stage_data.items():
                        if isinstance(data, dict):
                            if "stats" in data and "basic" in data["stats"]:
                                total += data["stats"]["basic"].get("total_variants", 0)
                            # Handle case where data has basic stats directly
                            elif "basic" in data:
                                total += data["basic"].get("total_variants", 0)
            return total

        # Build comparison data - include ALL stages present in either workflow
        all_stages = set()
        for stage in VCF_STAGE_ORDER:
            if stage in rna_standard_stats or stage in rna_realignment_stats:
                all_stages.add(stage)

        # Also check for any stages that might not be in VCF_STAGE_ORDER
        all_stages.update(rna_standard_stats.keys())
        all_stages.update(rna_realignment_stats.keys())

        # Sort stages according to VCF_STAGE_ORDER
        sorted_stages = [s for s in VCF_STAGE_ORDER if s in all_stages]
        # Add any remaining stages not in VCF_STAGE_ORDER
        for stage in all_stages:
            if stage not in sorted_stages:
                sorted_stages.append(stage)

        rows = []
        for stage in sorted_stages:
            standard_count = aggregate_total_variants(rna_standard_stats, stage)
            realignment_count = aggregate_total_variants(rna_realignment_stats, stage)

            # Only include stages with data
            if standard_count > 0 or realignment_count > 0:
                rows.append(
                    {
                        "Stage": stage,
                        "Standard_RNA_Count": standard_count,
                        "Realignment_RNA_Count": realignment_count,
                    }
                )

        if not rows:
            print("No data available for RNA stage progression comparison")
            return None

        import pandas as pd

        variant_counts_df = pd.DataFrame(rows)

        # Create line plot
        fig = go.Figure()

        stages = variant_counts_df["Stage"].values
        standard_counts = variant_counts_df["Standard_RNA_Count"].values
        realignment_counts = variant_counts_df["Realignment_RNA_Count"].values

        # Add standard workflow line
        fig.add_trace(
            go.Scatter(
                name="Standard RNA",
                x=stages,
                y=standard_counts,
                mode="lines+markers",
                line=dict(color="#4472C4", width=3),
                marker=dict(size=10, symbol="circle"),
                text=[f"{c:,}" for c in standard_counts],
                textposition="top center",
                hovertemplate="<b>%{x}</b><br>Standard: %{y:,}<extra></extra>",
            )
        )

        # Add realignment workflow line
        fig.add_trace(
            go.Scatter(
                name="Realignment RNA",
                x=stages,
                y=realignment_counts,
                mode="lines+markers",
                line=dict(color="#ED7D31", width=3),
                marker=dict(size=10, symbol="diamond"),
                text=[f"{c:,}" for c in realignment_counts],
                textposition="bottom center",
                hovertemplate="<b>%{x}</b><br>Realignment: %{y:,}<extra></extra>",
            )
        )

        # Update layout
        fig.update_layout(
            title_text="RNA Variant Count Progression Through Pipeline Stages",
            xaxis_title="Processing Stage",
            yaxis_title="Number of Variants",
            template="plotly_white",
            height=600,
            showlegend=True,
            legend=dict(
                orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1
            ),
            hovermode="x unified",
        )

        # Rotate x-axis labels for readability
        fig.update_xaxes(tickangle=-45)

        return fig

    def plot_rna_annotation_impact_comparison(
        self, rna_standard_stats: Dict[str, Any], rna_realignment_stats: Dict[str, Any]
    ):
        """
        Heatmap showing annotation stage differences for RNA samples.

        Rows: Stages (cosmic_gnomad, rna_editing, filtered_rescue)
        Columns: Categories (Somatic, Germline, etc.)
        Values: Difference in counts (realignment - standard)

        Args:
            rna_standard_stats: Statistics from standard RNA workflow
            rna_realignment_stats: Statistics from realignment RNA workflow
        """
        if not VISUALIZATION_AVAILABLE:
            print("Visualization libraries not available. Skipping plot.")
            return None

        # Import comparison module
        from .comparison import WorkflowComparator

        # Create comparator
        comparator = WorkflowComparator(rna_standard_stats, rna_realignment_stats)

        # Get annotation stage comparison data
        annotation_df = comparator.compare_rna_annotation_stages()

        if annotation_df.empty:
            print("No annotation stage data available for comparison")
            return None

        # Check if we have the required columns
        if (
            "Stage" not in annotation_df.columns
            or "Category" not in annotation_df.columns
        ):
            print(
                f"Missing required columns in annotation data. Available columns: {annotation_df.columns.tolist()}"
            )
            return None

        # Pivot data for heatmap
        # Rows = Stages, Columns = Categories, Values = Difference
        try:
            heatmap_data = annotation_df.pivot(
                index="Stage", columns="Category", values="Difference"
            )
        except Exception as e:
            print(f"Error creating pivot table: {e}")
            print(f"DataFrame shape: {annotation_df.shape}")
            print(f"DataFrame columns: {annotation_df.columns.tolist()}")
            print(f"First few rows:\n{annotation_df.head()}")
            return None

        # Reorder columns by CATEGORY_ORDER
        available_categories = [c for c in CATEGORY_ORDER if c in heatmap_data.columns]
        heatmap_data = heatmap_data[available_categories]

        # Create heatmap
        fig = go.Figure(
            data=go.Heatmap(
                z=heatmap_data.values,
                x=heatmap_data.columns,
                y=heatmap_data.index,
                colorscale="RdBu_r",  # Red for negative (fewer), Blue for positive (more)
                zmid=0,  # Center colorscale at 0
                text=heatmap_data.values,
                texttemplate="%{text:,}",
                textfont={"size": 12},
                colorbar=dict(title="Difference<br>(Realignment - Standard)"),
                hovertemplate="<b>%{y}</b><br>%{x}<br>Difference: %{z:,}<extra></extra>",
            )
        )

        # Update layout
        fig.update_layout(
            title_text="RNA Annotation Stage Impact: Realignment vs Standard<br><sub>Positive values = more variants in realignment, Negative = fewer variants</sub>",
            xaxis_title="Category",
            yaxis_title="Processing Stage",
            template="plotly_white",
            height=500,
            width=1000,
        )

        return fig

    def plot_integrative_view(
        self,
        dna_stats: Dict[str, Any],
        rna_standard_stats: Dict[str, Any],
        rna_realignment_stats: Dict[str, Any],
    ):
        """
        Comprehensive visualization of all three modalities.

        Creates grouped bar charts showing:
        - DNA-tumor (standard only)
        - RNA-tumor (standard)
        - RNA-tumor (realignment)

        This provides a complete picture of variant calling across
        all samples and workflows, highlighting the realignment impact
        in the context of DNA results.

        Args:
            dna_stats: Statistics from DNA standard workflow
            rna_standard_stats: Statistics from RNA standard workflow
            rna_realignment_stats: Statistics from RNA realignment workflow
        """
        if not VISUALIZATION_AVAILABLE:
            print("Visualization libraries not available. Skipping plot.")
            return None

        # Import comparison module
        from .comparison import WorkflowComparator

        # Create comparator (using dummy standard stats since we're providing all three separately)
        comparator = WorkflowComparator({}, {})

        # Get integrative view data
        integrative_df = comparator.create_integrative_view(
            dna_stats=dna_stats,
            rna_standard_stats=rna_standard_stats,
            rna_realignment_stats=rna_realignment_stats,
        )

        # Create figure with subplots for each stage
        annotation_stages = [
            "rescue",
            "cosmic_gnomad",
            "rna_editing",
            "filtered_rescue",
        ]
        available_stages = [
            s for s in annotation_stages if s in integrative_df["Stage"].values
        ]
        n_stages = len(available_stages)

        if n_stages == 0:
            print("No data available for integrative view")
            return None

        fig = make_subplots(
            rows=1,
            cols=n_stages,
            subplot_titles=available_stages,
            horizontal_spacing=0.10,
            shared_yaxes=True,
        )

        # Track categories for legend
        categories_in_legend = set()

        # For each stage, create grouped bars by category
        for i, stage in enumerate(available_stages, 1):
            stage_data = integrative_df[integrative_df["Stage"] == stage]

            # For each category, add grouped bars
            for category in CATEGORY_ORDER:
                cat_data = stage_data[stage_data["Category"] == category]
                if not cat_data.empty:
                    dna_count = cat_data["DNA_Standard"].values[0]
                    rna_std_count = cat_data["RNA_Standard"].values[0]
                    rna_real_count = cat_data["RNA_Realignment"].values[0]

                    # Only show in legend once
                    show_legend = category not in categories_in_legend
                    if show_legend:
                        categories_in_legend.add(category)

                    # Add three bars for this category (DNA, RNA Standard, RNA Realignment)
                    fig.add_trace(
                        go.Bar(
                            name=category,
                            x=["DNA", "RNA-Std", "RNA-Real"],
                            y=[dna_count, rna_std_count, rna_real_count],
                            marker_color=self.CATEGORY_COLORS.get(category, "#8A8A8A"),
                            text=[dna_count, rna_std_count, rna_real_count],
                            textposition="inside",
                            showlegend=show_legend,
                            legendgroup=category,
                            hovertemplate=f"<b>%{{x}}</b><br>{category}: %{{y}}<extra></extra>",
                        ),
                        row=1,
                        col=i,
                    )

            # Update axes
            fig.update_xaxes(title_text="Sample Type", row=1, col=i)
            if i == 1:
                fig.update_yaxes(title_text="Number of Variants", row=1, col=i)

        # Update layout
        fig.update_layout(
            title_text="Integrative View: DNA + RNA Standard + RNA Realignment<br><sub>Comparing all modalities across annotation stages</sub>",
            template="plotly_white",
            height=600,
            barmode="stack",
            showlegend=True,
            legend=dict(
                orientation="v",
                yanchor="top",
                y=1.0,
                xanchor="left",
                x=1.02,
                title="Category",
            ),
        )

        return fig


print("✓ VCF Visualizer module loaded successfully")
