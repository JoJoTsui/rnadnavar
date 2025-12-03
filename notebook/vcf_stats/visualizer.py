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