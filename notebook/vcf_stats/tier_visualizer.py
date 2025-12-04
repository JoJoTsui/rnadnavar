#!/usr/bin/env python3
"""
Tiered Variant Visualization Module

Generate HTML/Plotly visualizations for tiered variants from rescue VCF analysis.
Includes tier distribution plots, category-tier heatmaps, and statistical summaries.
"""

from typing import Dict, Optional, Tuple, List, Any
from pathlib import Path
import pandas as pd

# Import visualization dependencies with error handling
try:
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    import plotly.express as px
    VISUALIZATION_AVAILABLE = True
except ImportError as e:
    print(f"Warning: Visualization dependencies not available: {e}")
    VISUALIZATION_AVAILABLE = False


class TierVisualizer:
    """Create visualizations for tiered rescue variant statistics"""

    def __init__(self, tiered_df: pd.DataFrame):
        """
        Initialize tier visualizer with tiered variant data.

        Args:
            tiered_df: DataFrame with columns: chrom, pos, ref, alt, filter_category,
                      tier, dna_callers, rna_callers, etc.
        """
        self.tiered_df = tiered_df

        # Define unified color scheme for tiers (progression from strong to weak support)
        self.TIER_COLORS = {
            "T1": "#006400",  # Dark green - Strongest consensus
            "T2": "#228B22",  # Forest green
            "T3": "#32CD32",  # Lime green
            "T4": "#90EE90",  # Light green
            "T5": "#FFD700",  # Gold - Single caller
            "T6": "#FF8C00",  # Dark orange
            "T7": "#FF6347",  # Tomato
            "T8": "#DC143C",  # Crimson - Weakest
        }

        # Category colors (biological classification)
        self.CATEGORY_COLORS = {
            "PASS": "#636EFA",
            "Somatic": "#636EFA",
            "Germline": "#00CC96",
            "Reference": "#FFA15A",
            "Artifact": "#EF553B",
            "LowQual": "#AB63FA",
            "StrandBias": "#FFA500",
            "Clustered": "#8B4789",
            "Other": "#8A8A8A"
        }

    def plot_tier_distribution(self) -> Optional[go.Figure]:
        """
        Stacked bar chart showing variant counts for each tier.
        Each bar represents a tier, colored by biological category.
        """
        if not VISUALIZATION_AVAILABLE:
            print("Visualization libraries not available. Skipping plot.")
            return None

        if self.tiered_df.empty:
            print("No tiered variant data available")
            return None

        # Group by tier and filter_category
        data = []
        for tier in [f"T{i}" for i in range(1, 9)]:
            tier_data = self.tiered_df[self.tiered_df["tier"] == tier]
            if not tier_data.empty:
                for category in sorted(tier_data["filter_category"].unique()):
                    cat_count = len(tier_data[tier_data["filter_category"] == category])
                    data.append({
                        "Tier": tier,
                        "Category": category,
                        "Count": cat_count
                    })

        if not data:
            print("No tier distribution data available")
            return None

        df = pd.DataFrame(data)

        # Create stacked bar chart
        fig = go.Figure()

        # Get unique categories in sorted order
        categories = sorted(df["Category"].unique())

        for category in categories:
            df_cat = df[df["Category"] == category]
            tiers = [f"T{i}" for i in range(1, 9)]
            counts = [
                df_cat[df_cat["Tier"] == t]["Count"].sum()
                if not df_cat[df_cat["Tier"] == t].empty else 0
                for t in tiers
            ]

            fig.add_trace(go.Bar(
                name=category,
                x=tiers,
                y=counts,
                marker_color=self.CATEGORY_COLORS.get(category, "#8A8A8A"),
                text=counts,
                textposition="inside",
                textfont=dict(color="white", size=11),
                hovertemplate=f"<b>%{{x}}</b><br>{category}: %{{y}} variants<extra></extra>"
            ))

        fig.update_layout(
            title={
                "text": "Rescue Variant Distribution by Tier and Biological Category",
                "font": {"size": 16, "color": "#333"}
            },
            xaxis_title="Tier (Caller Support Strength)",
            yaxis_title="Number of Variants",
            template="plotly_white",
            barmode="stack",
            height=500,
            hovermode="x unified",
            legend=dict(
                orientation="v",
                yanchor="top",
                y=0.99,
                xanchor="right",
                x=0.99,
                bgcolor="rgba(255, 255, 255, 0.8)",
                bordercolor="rgba(0, 0, 0, 0.2)",
                borderwidth=1
            )
        )

        return fig

    def plot_category_tier_heatmap(self) -> Optional[go.Figure]:
        """
        Heatmap showing variant counts for each Category-Tier combination.
        Rows are biological categories, columns are tiers.
        """
        if not VISUALIZATION_AVAILABLE:
            print("Visualization libraries not available. Skipping plot.")
            return None

        if self.tiered_df.empty:
            print("No tiered variant data available")
            return None

        # Create pivot table: Categories x Tiers
        categories = sorted(self.tiered_df["filter_category"].unique())
        tiers = [f"T{i}" for i in range(1, 9)]

        data_matrix = []
        for category in categories:
            row = []
            cat_df = self.tiered_df[self.tiered_df["filter_category"] == category]
            for tier in tiers:
                count = len(cat_df[cat_df["tier"] == tier])
                row.append(count)
            data_matrix.append(row)

        # Create heatmap
        fig = go.Figure(data=go.Heatmap(
            z=data_matrix,
            x=tiers,
            y=categories,
            colorscale="YlOrRd",
            text=data_matrix,
            texttemplate="%{text}",
            textfont={"size": 12},
            colorbar=dict(
                title="Variant<br>Count",
                tickfont={"size": 11}
            ),
            hovertemplate="<b>Category:</b> %{y}<br><b>Tier:</b> %{x}<br><b>Count:</b> %{z}<extra></extra>"
        ))

        fig.update_layout(
            title={
                "text": "Variant Count Heatmap: Biological Category × Tier",
                "font": {"size": 16, "color": "#333"}
            },
            xaxis_title="Tier (Caller Support)",
            yaxis_title="Biological Category",
            template="plotly_white",
            height=400,
            width=900,
            xaxis={"side": "bottom"},
            yaxis={"side": "left"}
        )

        return fig

    def plot_tier_statistics_summary(self) -> Optional[go.Figure]:
        """
        Multi-panel visualization showing:
        1. Tier distribution pie chart
        2. DNA caller distribution
        3. RNA caller distribution
        4. Category distribution
        """
        if not VISUALIZATION_AVAILABLE:
            print("Visualization libraries not available. Skipping plot.")
            return None

        if self.tiered_df.empty:
            print("No tiered variant data available")
            return None

        # Create subplots: 2x2 grid
        fig = make_subplots(
            rows=2, cols=2,
            specs=[[{"type": "pie"}, {"type": "bar"}],
                   [{"type": "bar"}, {"type": "bar"}]],
            subplot_titles=(
                "Tier Distribution",
                "DNA Caller Support",
                "RNA Caller Support",
                "Biological Category"
            ),
            vertical_spacing=0.15,
            horizontal_spacing=0.12
        )

        # Panel 1: Tier distribution pie chart
        tier_counts = self.tiered_df["tier"].value_counts().sort_index()
        tier_labels = [f"{t}: {tier_counts.get(t, 0)} vars" for t in [f"T{i}" for i in range(1, 9)] if t in tier_counts.index]
        tier_values = [tier_counts[t] for t in sorted(tier_counts.index)]

        fig.add_trace(
            go.Pie(
                labels=sorted(self.tiered_df["tier"].unique()),
                values=tier_values,
                marker=dict(colors=[self.TIER_COLORS.get(t, "#8A8A8A") for t in sorted(self.tiered_df["tier"].unique())]),
                textposition="inside",
                textinfo="label+percent",
                hovertemplate="<b>%{label}</b><br>Count: %{value}<br>Percent: %{percent}<extra></extra>"
            ),
            row=1, col=1
        )

        # Panel 2: DNA caller distribution bar chart
        dna_dist = self.tiered_df["dna_callers"].value_counts().sort_index()
        fig.add_trace(
            go.Bar(
                x=[f"{int(x)} callers" for x in dna_dist.index],
                y=dna_dist.values,
                marker_color="#636EFA",
                text=dna_dist.values,
                textposition="outside",
                name="DNA Callers",
                hovertemplate="<b>%{x}</b><br>Count: %{y}<extra></extra>"
            ),
            row=1, col=2
        )

        # Panel 3: RNA caller distribution bar chart
        rna_dist = self.tiered_df["rna_callers"].value_counts().sort_index()
        fig.add_trace(
            go.Bar(
                x=[f"{int(x)} callers" for x in rna_dist.index],
                y=rna_dist.values,
                marker_color="#00CC96",
                text=rna_dist.values,
                textposition="outside",
                name="RNA Callers",
                hovertemplate="<b>%{x}</b><br>Count: %{y}<extra></extra>"
            ),
            row=2, col=1
        )

        # Panel 4: Biological category distribution bar chart
        cat_dist = self.tiered_df["filter_category"].value_counts().sort_values(ascending=False)
        fig.add_trace(
            go.Bar(
                x=cat_dist.index,
                y=cat_dist.values,
                marker_color=[self.CATEGORY_COLORS.get(c, "#8A8A8A") for c in cat_dist.index],
                text=cat_dist.values,
                textposition="outside",
                name="Category",
                hovertemplate="<b>%{x}</b><br>Count: %{y}<extra></extra>"
            ),
            row=2, col=2
        )

        # Update axes labels
        fig.update_xaxes(title_text="", row=1, col=2)
        fig.update_yaxes(title_text="Count", row=1, col=2)
        fig.update_xaxes(title_text="", row=2, col=1)
        fig.update_yaxes(title_text="Count", row=2, col=1)
        fig.update_xaxes(title_text="", row=2, col=2)
        fig.update_yaxes(title_text="Count", row=2, col=2)

        fig.update_layout(
            title={
                "text": "Rescue Variant Tiering Statistics Summary",
                "font": {"size": 16, "color": "#333"}
            },
            template="plotly_white",
            height=800,
            showlegend=False,
            hovermode="closest"
        )

        return fig

    def plot_caller_support_distribution(self) -> Optional[go.Figure]:
        """
        Scatter plot showing relationship between DNA and RNA caller support.
        Each point represents variants at a specific (dna_callers, rna_callers) coordinate.
        """
        if not VISUALIZATION_AVAILABLE:
            print("Visualization libraries not available. Skipping plot.")
            return None

        if self.tiered_df.empty:
            print("No tiered variant data available")
            return None

        # Count variants at each (dna_callers, rna_callers) position
        support_counts = self.tiered_df.groupby(
            ["dna_callers", "rna_callers", "tier"]
        ).size().reset_index(name="count")

        fig = go.Figure()

        # Add traces for each tier
        for tier in sorted(support_counts["tier"].unique()):
            tier_data = support_counts[support_counts["tier"] == tier]
            fig.add_trace(go.Scatter(
                x=tier_data["dna_callers"],
                y=tier_data["rna_callers"],
                mode="markers",
                name=tier,
                marker=dict(
                    size=tier_data["count"] ** 0.5 * 5,  # Size proportional to count
                    color=self.TIER_COLORS.get(tier, "#8A8A8A"),
                    opacity=0.7,
                    line=dict(width=1, color="white")
                ),
                text=[f"{tier}<br>DNA: {d}, RNA: {r}<br>Variants: {c}" 
                      for d, r, c in zip(tier_data["dna_callers"], tier_data["rna_callers"], tier_data["count"])],
                hovertemplate="<b>%{text}</b><extra></extra>"
            ))

        fig.update_layout(
            title={
                "text": "Variant Distribution by DNA and RNA Caller Support",
                "font": {"size": 16, "color": "#333"}
            },
            xaxis_title="DNA Callers",
            yaxis_title="RNA Callers",
            template="plotly_white",
            height=500,
            width=800,
            hovermode="closest",
            legend=dict(
                title="Tier",
                yanchor="top",
                y=0.99,
                xanchor="right",
                x=0.99
            ),
            xaxis=dict(dtick=1, range=[-0.5, self.tiered_df["dna_callers"].max() + 0.5]),
            yaxis=dict(dtick=1, range=[-0.5, self.tiered_df["rna_callers"].max() + 0.5])
        )

        return fig

    def plot_tier_composition_by_category(self) -> Optional[go.Figure]:
        """
        Grouped bar chart showing tier composition within each biological category.
        X-axis: biological categories, grouped by tier colors.
        """
        if not VISUALIZATION_AVAILABLE:
            print("Visualization libraries not available. Skipping plot.")
            return None

        if self.tiered_df.empty:
            print("No tiered variant data available")
            return None

        # Group by category and tier
        data = []
        for category in sorted(self.tiered_df["filter_category"].unique()):
            cat_df = self.tiered_df[self.tiered_df["filter_category"] == category]
            for tier in [f"T{i}" for i in range(1, 9)]:
                count = len(cat_df[cat_df["tier"] == tier])
                if count > 0:
                    data.append({
                        "Category": category,
                        "Tier": tier,
                        "Count": count
                    })

        if not data:
            print("No tier composition data available")
            return None

        df = pd.DataFrame(data)

        fig = go.Figure()

        # Add bars for each tier
        tiers = sorted(df["Tier"].unique())
        for tier in tiers:
            tier_data = df[df["Tier"] == tier]
            categories = tier_data["Category"].values
            counts = tier_data["Count"].values

            fig.add_trace(go.Bar(
                name=tier,
                x=categories,
                y=counts,
                marker_color=self.TIER_COLORS.get(tier, "#8A8A8A"),
                text=counts,
                textposition="outside",
                hovertemplate="<b>%{x}</b><br>" + tier + ": %{y} variants<extra></extra>"
            ))

        fig.update_layout(
            title={
                "text": "Tier Composition within Each Biological Category",
                "font": {"size": 16, "color": "#333"}
            },
            xaxis_title="Biological Category",
            yaxis_title="Number of Variants",
            template="plotly_white",
            barmode="group",
            height=500,
            hovermode="x unified",
            legend=dict(
                title="Tier",
                yanchor="top",
                y=0.99,
                xanchor="right",
                x=0.99
            )
        )

        return fig


def create_tier_visualization_report(
    tiered_df: pd.DataFrame,
    output_dir: Optional[Path] = None
) -> Dict[str, go.Figure]:
    """
    Create comprehensive visualization report for tiered variants.
    
    Args:
        tiered_df: DataFrame with tiered variant data
        output_dir: Optional directory to save HTML reports
        
    Returns:
        Dictionary mapping plot names to figure objects
    """
    visualizer = TierVisualizer(tiered_df)

    report = {
        "tier_distribution": visualizer.plot_tier_distribution(),
        "category_tier_heatmap": visualizer.plot_category_tier_heatmap(),
        "statistics_summary": visualizer.plot_tier_statistics_summary(),
        "caller_support": visualizer.plot_caller_support_distribution(),
        "tier_composition": visualizer.plot_tier_composition_by_category(),
    }

    # Save HTML reports if output directory specified
    if output_dir:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        for plot_name, fig in report.items():
            if fig is not None:
                try:
                    html_path = output_dir / f"{plot_name}.html"
                    fig.write_html(str(html_path))
                    print(f"✓ Saved: {html_path}")
                except Exception as e:
                    print(f"✗ Error saving {plot_name}: {e}")

    return report
