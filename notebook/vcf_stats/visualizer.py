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

        # Define unified color scheme
        self.CATEGORY_COLORS = {
            "Somatic": "#636EFA",
            "Germline": "#00CC96",
            "Reference": "#EF553B",
            "Artifact": "#AB63FA",
            "PASS": "#636EFA",
            "LowQual": "#EF553B",
            "StrandBias": "#AB63FA",
            "Clustered": "#FFA500",
            "Other": "#8A8A8A"
        }

    def plot_variant_counts_by_tool(self):
        """
        Bar plot comparing variant counts with FILTER categories across tools and modalities.
        """
        if not VISUALIZATION_AVAILABLE:
            print("Visualization libraries not available. Skipping plot.")
            return None

        # Collect data
        data = []
        print(f"DEBUG: Starting plot_variant_counts_by_tool with {len(self.all_stats)} categories")

        for category, files in self.all_stats.items():
            for name, file_data in files.items():
                # The correct data structure has stats nested under the file_data
                if "stats" not in file_data:
                    print(f"DEBUG: No 'stats' for {name}")
                    continue
                    
                stats = file_data["stats"]
                if "basic" not in stats:
                    print(f"DEBUG: No 'basic' stats for {name}")
                    continue

                basic = stats["basic"]
                classification = basic.get("classification", {})
                
                # Only process files that have classification data
                if not classification:
                    continue

                parts = name.split("_")
                tool = parts[0] if parts else name
                modality = "DNA" if "DNA_TUMOR" in name else "RNA"

                print(f"DEBUG: {name}: classification = {classification}")

                # Add data for each classification category
                for filter_cat in CATEGORY_ORDER:
                    count = classification.get(filter_cat, 0)
                    if count > 0:
                        data.append({
                            "Tool": tool,
                            "Modality": modality,
                            "Category": filter_cat,
                            "Count": count,
                            "File": name
                        })

        print(f"DEBUG: Collected {len(data)} data entries for plotting")
        
        if not data:
            print("No data available for plotting.")
            return None

        # Create DataFrame
        df = pd.DataFrame(data)

        # Create subplots for DNA and RNA
        fig = make_subplots(
            rows=1, cols=2,
            subplot_titles=("DNA Tumor", "RNA Tumor"),
            specs=[[{"secondary_y": False}, {"secondary_y": False}]]
        )

        # Plot for each modality
        for i, modality in enumerate(["DNA", "RNA"]):
            col_idx = i + 1
            modality_data = df[df["Modality"] == modality]

            if modality_data.empty:
                continue

            # Get unique tools for this modality
            tools = sorted(modality_data["Tool"].unique())
            
            # Create stacked bar chart
            for filter_cat in CATEGORY_ORDER:
                y_values = []
                for tool in tools:
                    tool_data = modality_data[modality_data["Tool"] == tool]
                    if not tool_data.empty:
                        count = tool_data.iloc[0]["Count"] if tool_data["Category"].iloc[0] == filter_cat else 0
                        y_values.append(count)
                
                # Only add trace if there are values
                if any(y_values):
                    fig.add_trace(
                        go.Bar(
                            name=filter_cat,
                            x=tools,
                            y=y_values,
                            marker_color=self.CATEGORY_COLORS.get(filter_cat, "#8A8A8A"),
                            showlegend=(i == 0),  # Only show legend for first subplot
                            hovertemplate=f"<b>{filter_cat}</b><br>Tool: %{{x}}<br>Count: %{{y}}"
                        ),
                        row=1, col=col_idx
                    )

            # Update subplot layout
            fig.update_xaxes(title_text="Variant Calling Tool", row=1, col=col_idx)
            fig.update_yaxes(title_text="Number of Variants", row=1, col=col_idx)

        # Update overall layout
        fig.update_layout(
            title="Variant Counts by Tool and Modality",
            barmode="stack",
            height=600,
            width=1200,
            legend=dict(
                orientation="h",
                yanchor="bottom",
                y=1.02,
                xanchor="right",
                x=1
            )
        )

        fig.show()
        return fig

    def plot_variant_type_distribution(self):
        """
        Plot variant type distribution (SNP vs INDEL) by FILTER category.
        """
        if not VISUALIZATION_AVAILABLE:
            print("Visualization libraries not available. Skipping plot.")
            return None

        data = []

        # Collect variant type data from consensus and rescue categories
        for category in ["consensus", "rescue"]:
            if category not in self.all_stats:
                continue

            for name, file_data in self.all_stats[category].items():
                if "stats" not in file_data or "basic" not in file_data["stats"]:
                    continue

                basic = file_data["stats"]["basic"]
                classification = basic.get("classification", {})
                snps = basic.get("snps", 0)
                indels = basic.get("indels", 0)

                # Determine modality
                if category == "consensus":
                    modality = "DNA Consensus" if "DNA" in name else "RNA Consensus"
                elif category == "rescue":
                    modality = "Rescued"
                else:
                    continue

                # For each classification category, split by SNP/INDEL
                # Note: We approximate the split based on overall SNP/INDEL ratio
                total = snps + indels
                if total == 0:
                    continue

                snp_ratio = snps / total
                indel_ratio = indels / total

                for filter_cat in CATEGORY_ORDER:
                    count = classification.get(filter_cat, 0)
                    if count > 0:
                        snp_count = int(count * snp_ratio)
                        indel_count = int(count * indel_ratio)

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
            print("No variant type data available for plotting")
            return None

        df = pd.DataFrame(data)

        # Create subplots
        modalities = ["DNA Consensus", "RNA Consensus", "Rescued"]
        available_mods = [m for m in modalities if m in df["Modality"].values]
        n_mods = len(available_mods)

        if n_mods == 0:
            print("No modality data available")
            return None

        fig = make_subplots(
            rows=1,
            cols=n_mods,
            subplot_titles=available_mods,
            horizontal_spacing=0.10
        )

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
                            fig.add_trace(
                                go.Bar(
                                    name=filter_cat,
                                    x=[var_type],
                                    y=[count],
                                    marker_color=self.CATEGORY_COLORS.get(filter_cat, "#8A8A8A"),
                                    text=[count],
                                    textposition="inside",
                                    showlegend=(i == 1),
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
        return fig

    def plot_consensus_comparison(self):
        """
        Compare consensus variants to individual tools with FILTER category breakdown.
        """
        if not VISUALIZATION_AVAILABLE:
            print("Visualization libraries not available. Skipping plot.")
            return None

        data = []

        # Get tool-level classification data
        if "variant_calling" in self.all_stats:
            for name, vcf_data in self.all_stats["variant_calling"].items():
                if "stats" in vcf_data and "basic" in vcf_data["stats"]:
                    basic = vcf_data["stats"]["basic"]
                    classification = basic.get("classification", {})
                    parts = name.split("_")
                    tool = parts[0] if parts else name
                    modality = "DNA" if "DNA_TUMOR" in name else "RNA"

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
                    modality = "DNA" if "DNA" in modality_key else "RNA"

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
        return fig

    def plot_filter_status(self):
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

                if category == "variant_calling":
                    tool = parts[0] if parts else "unknown"
                    if "DNA_TUMOR" in name:
                        subplot_cat = "DNA"
                    elif "RNA_TUMOR" in name:
                        subplot_cat = "RNA"
                    else:
                        continue  # Skip if not DNA or RNA

                elif category == "consensus":
                    tool = "consensus"
                    if "DNA" in name:
                        subplot_cat = "DNA"
                    elif "RNA" in name:
                        subplot_cat = "RNA"
                    else:
                        continue

                elif category == "rescue":
                    tool = "rescue"
                    subplot_cat = "Rescue"

                else:
                    continue  # Skip other categories

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
        subplot_categories = ["DNA", "RNA", "Rescue"]
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
        return fig