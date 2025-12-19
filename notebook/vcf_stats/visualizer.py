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

# Import constants from main module
from . import CATEGORY_ORDER


class VCFVisualizer:
    """Create visualizations for VCF statistics"""

    def __init__(self, all_stats: Dict[str, Any]):
        """
        Initialize visualizer with statistics data.

        Args:
            all_stats: Dictionary containing all VCF statistics
        """
        self.all_stats = all_stats

        # Define unified color scheme (matching original notebook)
        self.CATEGORY_COLORS = {
            "Somatic": "#636EFA",
            "Germline": "#00CC96",
            "Reference": "#FFA15A",
            "Artifact": "#EF553B",
            "RNA_Edit": "#AB63FA",
            "NoConsensus": "#8A8A8A",
            "PASS": "#636EFA",
            "LowQual": "#EF553B",
            "StrandBias": "#AB63FA",
            "Clustered": "#FFA500",
            "Other": "#8A8A8A"
        }

    def plot_variant_counts_by_tool(self, exclude_no_consensus_view: bool = True):
        """
        Bar plot comparing variant counts with FILTER categories across tools and modalities.
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

        # Optional: a second view excluding NoConsensus to highlight minor categories
        if exclude_no_consensus_view:
            df2 = df[df["Category"] != "NoConsensus"].copy()
            if not df2.empty:
                fig2 = make_subplots(
                    rows=1,
                    cols=2,
                    subplot_titles=("DNA Modality (No NoConsensus)", "RNA Modality (No NoConsensus)"),
                    horizontal_spacing=0.12,
                    shared_yaxes=True
                )

                # DNA
                df_dna2 = df2[df2["Modality"] == "DNA"]
                if not df_dna2.empty:
                    tools = sorted(df_dna2["Tool"].unique())
                    for filter_cat in [c for c in CATEGORY_ORDER if c != "NoConsensus"]:
                        df_cat = df_dna2[df_dna2["Category"] == filter_cat]
                        counts = [
                            df_cat[df_cat["Tool"] == t]["Count"].sum()
                            if not df_cat[df_cat["Tool"] == t].empty
                            else 0
                            for t in tools
                        ]
                        if sum(counts) > 0:
                            fig2.add_trace(
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

                # RNA
                df_rna2 = df2[df2["Modality"] == "RNA"]
                if not df_rna2.empty:
                    tools = sorted(df_rna2["Tool"].unique())
                    for filter_cat in [c for c in CATEGORY_ORDER if c != "NoConsensus"]:
                        df_cat = df_rna2[df_rna2["Category"] == filter_cat]
                        counts = [
                            df_cat[df_cat["Tool"] == t]["Count"].sum()
                            if not df_cat[df_cat["Tool"] == t].empty
                            else 0
                            for t in tools
                        ]
                        if sum(counts) > 0:
                            fig2.add_trace(
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

                fig2.update_layout(
                    title="Variant Counts by Tool (excluding NoConsensus)",
                    template="plotly_white",
                    barmode="stack",
                    height=500,
                    showlegend=True
                )
                fig2.show()

    def plot_variant_type_distribution(self):
        """
        Stacked bar charts showing SNP vs INDEL distribution with FILTER categories.
        
        Fixed: Unified color legend (no duplication) using legendgroup.
        """
        if not VISUALIZATION_AVAILABLE:
            print("Visualization libraries not available. Skipping plot.")
            return None

        data = []

        # Collect data from consensus VCFs for cleaner view
        for vcf_type in ["consensus", "rescue"]:
            if vcf_type in self.all_stats:
                for name, vcf_data in self.all_stats[vcf_type].items():
                    if "stats" in vcf_data and "basic" in vcf_data["stats"]:
                        basic = vcf_data["stats"]["basic"]
                        classification = basic.get("classification", {})

                        if vcf_type == "consensus":
                            if "DNA_TUMOR" in name:
                                modality = "DNA Consensus"
                            elif "RNA_TUMOR" in name:
                                modality = "RNA Consensus"
                            else:
                                continue
                        else:
                            modality = "Rescued"

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
                                        "Modality": modality,
                                        "Type": "SNP",
                                        "Category": filter_cat,
                                        "Count": snp_count
                                    })
                                if indel_count > 0:
                                    data.append({
                                        "Modality": modality,
                                        "Type": "INDEL",
                                        "Category": filter_cat,
                                        "Count": indel_count
                                    })

        if not data:
            print("No data available for plotting")
            return None

        df = pd.DataFrame(data)

        # Create subplots
        modalities = ["DNA Consensus", "RNA Consensus", "Rescued"]
        available_mods = [m for m in modalities if m in df["Modality"].values]
        n_mods = len(available_mods)

        fig = make_subplots(
            rows=1,
            cols=n_mods,
            subplot_titles=available_mods,
            horizontal_spacing=0.10,
            shared_yaxes=True
        )

        # Track which categories have been shown in legend (to avoid duplication)
        categories_in_legend = set()

        for i, modality in enumerate(available_mods, 1):
            df_mod = df[df["Modality"] == modality]

            # Group by Type (SNP/INDEL) and stack by Category
            for var_type in ["SNP", "INDEL"]:
                df_type = df_mod[df_mod["Type"] == var_type]
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

            fig.update_xaxes(title_text="Variant Type", row=1, col=i)
            if i == 1:
                fig.update_yaxes(title_text="Number of Variants", row=1, col=i)

        fig.update_layout(
            title_text="Variant Type Distribution (SNP vs INDEL) by FILTER Category",
            height=500,
            barmode="stack",
            template="plotly_white",
            showlegend=True
        )

        fig.show()

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
        Stacked bar chart showing unified FILTER categories (Somatic, Germline, Reference, Artifact).
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
            tools = sorted(df_subplot["Tool"].unique())

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
                    tools = sorted(df_subplot["Tool"].unique())
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

    def plot_annotation_progression(self):
        """
        Plot variant count progression through annotation stages:
        rescue → cosmic_gnomad → rna_editing → filtered_rescue
        
        Shows how variant counts and categories change at each stage,
        allowing visualization of filtering/annotation effects.
        """
        if not VISUALIZATION_AVAILABLE:
            print("Visualization libraries not available. Skipping plot.")
            return None

        from . import VCF_STAGE_ORDER
        
        data = []

        # Collect data from each annotation stage in order
        for stage in VCF_STAGE_ORDER:
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
                            "Stage_Order": VCF_STAGE_ORDER.index(stage),
                            "Category": category,
                            "Count": count,
                            "Total": total_variants
                        })

        if not data:
            print("No annotation stage data available for progression plot")
            return None

        df = pd.DataFrame(data)

        # Create figure with stacked bars showing progression
        fig = go.Figure()

        # Get unique stages in order
        stages_ordered = sorted(df["Stage"].unique(), 
                               key=lambda x: VCF_STAGE_ORDER.index(x) if x in VCF_STAGE_ORDER else 999)

        # Track which categories to show in legend
        categories_shown = set()

        # Add bars for each category
        for category in CATEGORY_ORDER:
            df_cat = df[df["Category"] == category]
            if df_cat.empty:
                continue

            # Get counts for each stage
            counts_by_stage = []
            for stage in stages_ordered:
                stage_data = df_cat[df_cat["Stage"] == stage]
                if not stage_data.empty:
                    counts_by_stage.append(stage_data["Count"].sum())
                else:
                    counts_by_stage.append(0)

            # Only add if there are non-zero values
            if sum(counts_by_stage) > 0:
                show_legend = category not in categories_shown
                categories_shown.add(category)
                
                fig.add_trace(go.Bar(
                    name=category,
                    x=stages_ordered,
                    y=counts_by_stage,
                    marker_color=self.CATEGORY_COLORS.get(category, "#8A8A8A"),
                    text=counts_by_stage,
                    textposition="inside",
                    textfont=dict(color="white", size=10),
                    showlegend=show_legend,
                    legendgroup=category,
                    hovertemplate=f"<b>%{{x}}</b><br>{category}: %{{y}}<extra></extra>"
                ))

        fig.update_layout(
            title_text="Variant Count Progression Through Annotation Stages",
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
                title="Category"
            )
        )

        fig.show()