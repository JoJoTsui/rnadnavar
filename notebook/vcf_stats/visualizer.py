#!/usr/bin/env python3
"""
VCF Statistics Visualizer Module (Fixed Version)

Create visualizations for VCF statistics including variant counts,
quality distributions, and rescue analysis plots.
Modified to work with the actual data structure returned by vcf_processor.py.
"""

from typing import Dict, Any, Optional

# Import visualization dependencies with error handling
try:
    import matplotlib.pyplot as plt
    import seaborn as sns
    import plotly.express as px
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    import pandas as pd
    VISUALIZATION_AVAILABLE = True
except ImportError as e:
    print(f"Warning: Visualization dependencies not available: {e}")
    VISUALIZATION_AVAILABLE = False

# Import constants and utilities from main module
from . import CATEGORY_ORDER, CATEGORY_COLORS, VCF_STAGE_ORDER, STAGE_DISPLAY_NAMES
from .utils import sort_stages, should_show_legend, create_legend_config


class VCFVisualizer:
    """Create visualizations for VCF statistics"""

    def __init__(self, all_stats: Dict[str, Any]):
        """
        Initialize visualizer with statistics data.

        Args:
            all_stats: Dictionary containing all VCF statistics
        """
        self.all_stats = all_stats

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
                    modality = "DNA" if "DT_vs_DN" in name or "DNA_TUMOR" in name else "RNA"

                    # Add data for each FILTER category
                    for filter_cat in CATEGORY_ORDER:
                        count = classification.get(filter_cat, 0)
                        if count > 0:
                            data.append({
                                "Tool": tool,
                                "Modality": modality,
                                "Category": filter_cat,
                                "Count": count
                            })

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
            shared_yaxes=True
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
                            marker_color=self.CATEGORY_COLORS.get(filter_cat, "#8A8A8A"),
                            text=counts,
                            textposition="inside",
                            showlegend=True,
                            legendgroup=filter_cat
                        ),
                        row=1,
                        col=1
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
                            marker_color=self.CATEGORY_COLORS.get(filter_cat, "#8A8A8A"),
                            text=counts,
                            textposition="inside",
                            showlegend=False,
                            legendgroup=filter_cat
                        ),
                        row=1,
                        col=2
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
            showlegend=True
        )

        fig.show()

    def plot_variant_type_distribution(self, exclude_no_consensus_view: bool = True, debug_rna_edit: bool = False):
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
            "filtered_rescue": "Filtered"
        }

        # Collect data from all annotation pipeline stages
        for vcf_type in ["consensus", "rescue", "cosmic_gnomad", "rna_editing", "filtered_rescue"]:
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
                            stage_name = STAGE_DISPLAY_NAMES.get(vcf_type, vcf_type.title())

                        # Get variant types
                        variant_types = basic.get("variant_types", {})
                        snps = variant_types.get("SNP", 0)
                        indels = variant_types.get("DEL", 0) + variant_types.get("INS", 0)

                        # For each FILTER category, calculate proportional SNP/INDEL split
                        total_vars = basic.get("total_variants", 1)
                        for filter_cat in CATEGORY_ORDER:
                            count = classification.get(filter_cat, 0)
                            if count > 0:
                                # Proportionally split into SNPs and INDELs
                                snp_count = int(count * (snps / total_vars))
                                indel_count = count - snp_count

                                if snp_count > 0:
                                    data.append({
                                        "Stage": stage_name,
                                        "Type": "SNP",
                                        "Category": filter_cat,
                                        "Count": snp_count
                                    })
                                if indel_count > 0:
                                    data.append({
                                        "Stage": stage_name,
                                        "Type": "INDEL",
                                        "Category": filter_cat,
                                        "Count": indel_count
                                    })
                                
                                # Collect debug info for RNA_Edit
                                if debug_rna_edit and filter_cat == "RNA_Edit":
                                    debug_info.append({
                                        "stage": stage_name,
                                        "total_variants": total_vars,
                                        "total_snps": snps,
                                        "total_indels": indels,
                                        "rna_edit_count": count,
                                        "rna_edit_snps": snp_count,
                                        "rna_edit_indels": indel_count
                                    })

        # Debug output for RNA_Edit INDEL handling
        if debug_rna_edit and debug_info:
            print("\n" + "="*80)
            print("DEBUG: RNA_Edit Variant Type Distribution")
            print("="*80)
            for info in debug_info:
                print(f"\n{info['stage']}:")
                print(f"  Total variants: {info['total_variants']:,}")
                print(f"  SNPs (total): {info['total_snps']:,} ({info['total_snps']/info['total_variants']*100:.1f}%)")
                print(f"  INDELs (total): {info['total_indels']:,} ({info['total_indels']/info['total_variants']*100:.1f}%)")
                print(f"  RNA_Edit count: {info['rna_edit_count']:,}")
                print(f"    - SNPs: {info['rna_edit_snps']:,} ({info['rna_edit_snps']/info['rna_edit_count']*100:.1f}%)")
                print(f"    - INDELs: {info['rna_edit_indels']:,} ({info['rna_edit_indels']/info['rna_edit_count']*100:.1f}%)")
            print("="*80 + "\n")

        if not data:
            print("No data available for plotting")
            return None

        df = pd.DataFrame(data)

        def _create_variant_type_plot(plot_df, title_suffix=""):
            """Helper to create variant type distribution plot."""
            # Determine stages present in data and sort by pipeline order
            stage_order = ["DNA Consensus", "RNA Consensus", "Rescue", 
                          "COSMIC/GnomAD", "RNA Editing", "Filtered"]
            available_stages = [s for s in stage_order if s in plot_df["Stage"].values]
            n_stages = len(available_stages)

            if n_stages == 0:
                return None

            fig = make_subplots(
                rows=1,
                cols=n_stages,
                subplot_titles=available_stages,
                horizontal_spacing=0.08,
                shared_yaxes=True
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
                                        marker_color=self.CATEGORY_COLORS.get(filter_cat, "#8A8A8A"),
                                        text=[count],
                                        textposition="inside",
                                        showlegend=show_legend,
                                        legendgroup=filter_cat,
                                        hovertemplate=f"<b>{var_type}</b><br>{filter_cat}: %{{y}}<extra></extra>"
                                    ),
                                    row=1,
                                    col=i
                                )

                fig.update_xaxes(title_text="Type", row=1, col=i)
                if i == 1:
                    fig.update_yaxes(title_text="Number of Variants", row=1, col=i)

            fig.update_layout(
                title_text=f"Variant Type Distribution (SNP vs INDEL) Across Pipeline Stages{title_suffix}",
                height=500,
                barmode="stack",
                template="plotly_white",
                showlegend=True
            )
            return fig

        # Create main plot with all categories
        full_fig = _create_variant_type_plot(df)
        if full_fig:
            full_fig.show()

        # Create NoConsensus-free plot if requested
        if exclude_no_consensus_view:
            df_filtered = df[df["Category"] != "NoConsensus"].copy()
            if not df_filtered.empty:
                filtered_fig = _create_variant_type_plot(df_filtered, " (excluding NoConsensus)")
                if filtered_fig:
                    filtered_fig.show()

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
                    modality = "DNA" if "DT_vs_DN" in name or "DNA_TUMOR" in name else "RNA"

                    for filter_cat in CATEGORY_ORDER:
                        count = classification.get(filter_cat, 0)
                        if count > 0:
                            data.append({
                                "Category": tool,
                                "Modality": modality,
                                "FilterCat": filter_cat,
                                "Count": count
                            })

        # Get consensus classification data
        if "consensus" in self.all_stats:
            for modality_key, vcf_data in self.all_stats["consensus"].items():
                if "stats" in vcf_data and "basic" in vcf_data["stats"]:
                    basic = vcf_data["stats"]["basic"]
                    classification = basic.get("classification", {})
                    # Explicitly detect modality using suffix keys
                    if "RT_vs_DN" in modality_key or "RNA_TUMOR" in modality_key:
                        modality = "RNA"
                    else:
                        modality = "DNA"

                    for filter_cat in CATEGORY_ORDER:
                        count = classification.get(filter_cat, 0)
                        if count > 0:
                            data.append({
                                "Category": "consensus",
                                "Modality": modality,
                                "FilterCat": filter_cat,
                                "Count": count
                            })

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
            shared_yaxes=True
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
                            marker_color=self.CATEGORY_COLORS.get(filter_cat, "#8A8A8A"),
                            text=counts,
                            textposition="inside",
                            showlegend=True,
                            legendgroup=filter_cat
                        ),
                        row=1,
                        col=1
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
                            marker_color=self.CATEGORY_COLORS.get(filter_cat, "#8A8A8A"),
                            text=counts,
                            textposition="inside",
                            showlegend=False,
                            legendgroup=filter_cat
                        ),
                        row=1,
                        col=2
                    )

        fig.update_xaxes(title_text="Category", row=1, col=1)
        fig.update_xaxes(title_text="Category", row=1, col=2)
        fig.update_yaxes(title_text="Number of Variants", row=1, col=1)

        fig.update_layout(
            title_text="Variant Counts by Category with FILTER Classification",
            template="plotly_white",
            height=500,
            barmode="stack"
        )

        fig.show()

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

                elif category in {"rescue", "cosmic_gnomad", "rna_editing", "filtered_rescue"}:
                    tool = category
                    # Show these together in a single subplot
                    subplot_cat = "Annotation Stages"

                else:
                    continue

                # Add counts for each classification category
                for filter_cat in CATEGORY_ORDER:
                    count = classification.get(filter_cat, 0)
                    if count > 0:  # Only add if there are variants
                        data.append({
                            "SubplotCategory": subplot_cat,
                            "Tool": tool,
                            "FilterCategory": filter_cat,
                            "Count": count,
                            "Name": name
                        })

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
            shared_yaxes=True
        )

        # Find global y-axis max for consistent scaling
        max_y = 0

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
                            marker_color=self.CATEGORY_COLORS.get(filter_cat, "#8A8A8A"),
                            text=counts,
                            textposition="inside",
                            textfont=dict(color="white", size=10),
                            showlegend=(i == 1),
                            legendgroup=filter_cat,
                            hovertemplate=f"<b>%{{x}}</b><br>{filter_cat}: %{{y}}<extra></extra>"
                        ),
                        row=1,
                        col=i
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
            legend=dict(
                orientation="v",
                yanchor="top",
                y=1.0,
                xanchor="left",
                x=1.02,
                title="FILTER Category"
            )
        )

        # Set unified y-axis range with some padding
        for i in range(1, n_subplots + 1):
            fig.update_yaxes(range=[0, max_y * 1.1], row=1, col=i)

        fig.show()

        # Optional secondary view excluding NoConsensus to improve minor category visibility
        if dual_view_no_consensus and data:
            df_small = df[df["FilterCategory"] != "NoConsensus"].copy()
            if not df_small.empty:
                fig_small = make_subplots(
                    rows=1,
                    cols=len(available_subplots),
                    subplot_titles=[f"{t} (No NoConsensus)" for t in available_subplots],
                    horizontal_spacing=0.10,
                    shared_yaxes=True
                )

                max_y2 = 0
                for i, subplot_cat in enumerate(available_subplots, 1):
                    df_subplot = df_small[df_small["SubplotCategory"] == subplot_cat]
                    tools = df_subplot["Tool"].unique()
                    
                    # Sort annotation stages by pipeline order
                    if subplot_cat == "Annotation Stages":
                        tools = sort_stages(list(tools))
                    else:
                        tools = sorted(tools)
                    
                    for filter_cat in [c for c in CATEGORY_ORDER if c != "NoConsensus"]:
                        df_filter = df_subplot[df_subplot["FilterCategory"] == filter_cat]
                        counts = []
                        for tool in tools:
                            tool_data = df_filter[df_filter["Tool"] == tool]
                            count = tool_data["Count"].sum() if not tool_data.empty else 0
                            counts.append(count)
                        if sum(counts) > 0:
                            fig_small.add_trace(
                                go.Bar(
                                    name=filter_cat,
                                    x=tools,
                                    y=counts,
                                    marker_color=self.CATEGORY_COLORS.get(filter_cat, "#8A8A8A"),
                                    text=counts,
                                    textposition="inside",
                                    textfont=dict(color="white", size=10),
                                    showlegend=(i == 1),
                                    legendgroup=filter_cat,
                                ),
                                row=1,
                                col=i
                            )
                    for tool in tools:
                        tool_data = df_subplot[df_subplot["Tool"] == tool]
                        total = tool_data["Count"].sum()
                        max_y2 = max(max_y2, total)

                    fig_small.update_xaxes(title_text="Caller/Stage", row=1, col=i)
                    if i == 1:
                        fig_small.update_yaxes(title_text="Number of Variants", row=1, col=i)

                fig_small.update_layout(
                    title_text="Variant Classification (excluding NoConsensus)",
                    template="plotly_white",
                    height=500,
                    barmode="stack",
                    showlegend=True
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

        from . import VCF_STAGE_ORDER
        
        data = []
        
        # Only include annotation pipeline stages, not normalized or consensus
        annotation_stages = ["rescue", "cosmic_gnomad", "rna_editing", "filtered_rescue"]

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
                        data.append({
                            "Stage": stage,
                            "Stage_Order": annotation_stages.index(stage),
                            "Category": category,
                            "Count": count,
                            "Total": total_variants
                        })

        if not data:
            print("No annotation stage data available for progression plot")
            return None

        df = pd.DataFrame(data)

        def _build_progression(fig_df, title_suffix=""):
            fig = go.Figure()

            stages_ordered = sorted(
                fig_df["Stage"].unique(),
                key=lambda x: annotation_stages.index(x) if x in annotation_stages else 999,
            )

            categories_shown = set()

            for category in CATEGORY_ORDER:
                df_cat = fig_df[fig_df["Category"] == category]
                if df_cat.empty:
                    continue

                counts_by_stage = []
                for stage in stages_ordered:
                    stage_data = df_cat[df_cat["Stage"] == stage]
                    counts_by_stage.append(stage_data["Count"].sum() if not stage_data.empty else 0)

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
        full_fig.show()

        # Optional view excluding NoConsensus to improve readability
        if dual_view_no_consensus:
            df_small = df[df["Category"] != "NoConsensus"].copy()
            if not df_small.empty:
                small_fig = _build_progression(df_small, " (excluding NoConsensus)")
                small_fig.show()

    def plot_category_heatmap(self, stages: Optional[list] = None, dual_view_without_no_consensus: bool = True):
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
                            data.append({
                                "Stage": stage_name,
                                "Category": category,
                                "Count": count
                            })

        if not data:
            print("No data available for heatmap")
            return None

        df = pd.DataFrame(data)
        
        def _create_heatmap(plot_df, title_suffix=""):
            """Helper to create heatmap visualization."""
            # Pivot for heatmap
            heatmap_data = plot_df.pivot_table(
                index="Category", 
                columns="Stage", 
                values="Count", 
                fill_value=0
            )
            
            # Sort stages by pipeline order
            stage_names = [STAGE_DISPLAY_NAMES.get(s, s.title()) for s in stages if STAGE_DISPLAY_NAMES.get(s, s.title()) in heatmap_data.columns]
            heatmap_data = heatmap_data[stage_names]
            
            # Sort categories by CATEGORY_ORDER
            heatmap_data = heatmap_data.reindex([c for c in CATEGORY_ORDER if c in heatmap_data.index])

            # Create heatmap with automatic scaling
            fig = go.Figure(data=go.Heatmap(
                z=heatmap_data.values,
                x=heatmap_data.columns,
                y=heatmap_data.index,
                colorscale="Viridis",
                text=heatmap_data.values,
                texttemplate="%{text:,}",
                textfont={"size": 10},
                colorbar=dict(title="Variant Count"),
                hovertemplate="<b>%{y}</b><br>%{x}<br>Count: %{z:,}<extra></extra>",
                zauto=True  # Auto-scale based on data range
            ))

            fig.update_layout(
                title=f"Variant Category Distribution Across Pipeline Stages{title_suffix}",
                xaxis_title="Processing Stage",
                yaxis_title="Category",
                height=500 if "NoConsensus" in plot_df["Category"].values else 400,
                template="plotly_white"
            )
            return fig

        # Full view with all categories
        fig_full = _create_heatmap(df)
        fig_full.show()
        
        # Secondary view without NoConsensus to improve visibility
        if dual_view_without_no_consensus:
            df_no_nc = df[df["Category"] != "NoConsensus"].copy()
            if not df_no_nc.empty:
                fig_no_nc = _create_heatmap(df_no_nc, " (excluding NoConsensus)")
                fig_no_nc.show()

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
                    links.append({
                        "source": node_map[(stage_from, category)],
                        "target": node_map[(stage_to, category)],
                        "value": flow_value,
                        "category": category,
                        "label": f"{category}: {flow_value:,}"
                    })

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
            source_color = node_colors[source_idx] if source_idx < len(node_colors) else "#8A8A8A"
            # Convert hex to RGBA with transparency
            hex_color = source_color.lstrip('#')
            try:
                rgb = tuple(int(hex_color[i:i+2], 16) for i in (0, 2, 4))
                link_colors.append(f"rgba({rgb[0]}, {rgb[1]}, {rgb[2]}, 0.5)")
            except:
                link_colors.append("rgba(128, 128, 128, 0.5)")
        
        # Create Sankey diagram with styled colors and improved layout
        fig = go.Figure(data=[go.Sankey(
            node=dict(
                pad=15,
                thickness=20,
                line=dict(color="darkgray", width=1),
                label=nodes,
                color=node_colors
            ),
            link=dict(
                source=[l["source"] for l in links],
                target=[l["target"] for l in links],
                value=[l["value"] for l in links],
                color=link_colors,
                label=[l["label"] for l in links]
            )
        )])

        fig.update_layout(
            title=f"Variant Flow Through Pipeline Stages ({', '.join(categories_to_show)})",
            height=650,
            font=dict(size=11, family="Arial"),
            template="plotly_white"
        )

        fig.show()

    def plot_tier_distribution(self, rescue_vcf_path: Optional[str] = None, dual_view_no_consensus: bool = True):
        """
        Bar chart showing tier distribution (T1-T7) per category for rescue VCF.
        
        Tier definitions:
            T1: 2+ DNA + 2+ RNA (strongest consensus)
            T2: 2+ DNA + 1 RNA
            T3: 2+ DNA only
            T4: 1 DNA + 1+ RNA
            T5: 1 DNA only
            T6: 2+ RNA only
            T7: 1 RNA (weakest)
        
        Args:
            rescue_vcf_path: Path to rescue VCF (auto-detected if None)
            dual_view_no_consensus: If True, creates secondary plot excluding NoConsensus
        """
        if not VISUALIZATION_AVAILABLE:
            print("Visualization libraries not available. Skipping plot.")
            return None

        try:
            from .tiering import tier_rescue_variants
        except ImportError:
            print("Tiering module not available")
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
            """Helper function to create tier distribution plot."""
            # Create grouped bar chart data
            data = []
            for category in CATEGORY_ORDER:
                cat_data = plot_df[plot_df["filter_category"] == category]
                if not cat_data.empty:
                    tier_counts = cat_data["tier"].value_counts().sort_index()
                    for tier, count in tier_counts.items():
                        data.append({
                            "Category": category,
                            "Tier": tier,
                            "Count": count
                        })

            if not data:
                return None

            df = pd.DataFrame(data)

            # Create grouped bar chart
            fig = go.Figure()

            for category in CATEGORY_ORDER:
                cat_df = df[df["Category"] == category]
                if not cat_df.empty:
                    fig.add_trace(go.Bar(
                        name=category,
                        x=cat_df["Tier"],
                        y=cat_df["Count"],
                        marker_color=self.CATEGORY_COLORS.get(category, "#8A8A8A"),
                        text=cat_df["Count"],
                        textposition="inside",
                        showlegend=(category in CATEGORY_ORDER),
                        legendgroup=category,
                        hovertemplate="<b>%{x}</b><br>%{fullData.name}: %{y}<extra></extra>"
                    ))

            # Ensure tier order (T1-T7)
            tier_order = ["T1", "T2", "T3", "T4", "T5", "T6", "T7"]
            existing_tiers = [t for t in tier_order if t in df["Tier"].values]

            fig.update_layout(
                title=f"Rescue VCF Tier Distribution by Category{title_suffix}",
                xaxis=dict(
                    title="Tier",
                    categoryorder="array",
                    categoryarray=existing_tiers
                ),
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
                    title="Category"
                )
            )
            return fig

        # Full view with all categories
        fig_full = _create_tier_plot(tiered_df)
        if fig_full:
            fig_full.show()
        
        # Secondary view excluding NoConsensus
        if dual_view_no_consensus:
            tiered_df_no_nc = tiered_df[tiered_df["filter_category"] != "NoConsensus"].copy()
            if not tiered_df_no_nc.empty:
                fig_no_nc = _create_tier_plot(tiered_df_no_nc, " (excluding NoConsensus)")
                if fig_no_nc:
                    fig_no_nc.show()