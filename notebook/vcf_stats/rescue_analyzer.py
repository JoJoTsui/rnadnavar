#!/usr/bin/env python3
"""
Rescue Analyzer Module

Analyze and visualize rescue VCF statistics with FILTER category tracking
and transition analysis between DNA consensus, RNA consensus, and rescued variants.
"""

from typing import Dict, List, Any, Optional, Tuple
import pandas as pd

# Import visualization dependencies with error handling
try:
    import plotly.express as px
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    VISUALIZATION_AVAILABLE = True
except ImportError:
    VISUALIZATION_AVAILABLE = False

# Import constants and utilities from main module
from . import CATEGORY_ORDER, CATEGORY_COLORS, STAGE_DISPLAY_NAMES
from .utils import should_show_legend


def analyze_rescue_vcf(all_vcf_stats: Dict[str, Any], show_plot: bool = True) -> Dict[str, Any]:
    """
    Analyze and visualize rescue VCF statistics with FILTER category tracking and transition analysis.

    This function:
    1. Shows FILTER category breakdown for DNA consensus, RNA consensus, and rescued variants
    2. Analyzes transition patterns between stages
    3. Creates visualizations of the rescue process

    Args:
        all_vcf_stats: Dictionary containing all VCF statistics from VCFStatisticsExtractor
        show_plot: Whether to display interactive plots

    Returns:
        Dictionary containing rescue analysis results
    """
    if not VISUALIZATION_AVAILABLE and show_plot:
        print("Visualization libraries not available. Skipping plots.")
        show_plot = False

    # Initialize result dictionary
    rescue_analysis = {
        "dna_consensus": {},
        "rna_consensus": {},
        "rescue": {},
        "cosmic_gnomad": {},
        "rna_editing": {},
        "filtered_rescue": {},
        "transitions": {},
        "summary": {}
    }

    # Collect classification data for each stage
    dna_classification = {}
    rna_classification = {}
    rescue_classification = {}

    dna_total = 0
    rna_total = 0
    rescue_total = 0

    # Get DNA consensus classification
    if "consensus" in all_vcf_stats:
        for name, data in all_vcf_stats["consensus"].items():
            if "DNA_TUMOR" in name:
                basic = data.get("stats", {}).get("basic", {})
                dna_classification = basic.get("classification", {})
                dna_total = basic.get("total_variants", 0)
                rescue_analysis["dna_consensus"] = basic

    # Get RNA consensus classification
    if "consensus" in all_vcf_stats:
        for name, data in all_vcf_stats["consensus"].items():
            if "RNA_TUMOR" in name:
                basic = data.get("stats", {}).get("basic", {})
                rna_classification = basic.get("classification", {})
                rna_total = basic.get("total_variants", 0)
                rescue_analysis["rna_consensus"] = basic

    # Get rescued classification
    if "rescue" in all_vcf_stats:
        for name, data in all_vcf_stats["rescue"].items():
            basic = data.get("stats", {}).get("basic", {})
            rescue_classification = basic.get("classification", {})
            rescue_total = basic.get("total_variants", 0)
            rescue_analysis["rescue"] = basic

    # Store classifications in result
    rescue_analysis["dna_consensus"]["classification"] = dna_classification
    rescue_analysis["rna_consensus"]["classification"] = rna_classification
    rescue_analysis["rescue"]["classification"] = rescue_classification

    # Collect downstream annotation stages if available
    for stage in ["cosmic_gnomad", "rna_editing", "filtered_rescue"]:
        stage_class = {}
        stage_total = 0
        if stage in all_vcf_stats:
            for _, data in all_vcf_stats[stage].items():
                basic = data.get("stats", {}).get("basic", {})
                cls = basic.get("classification", {})
                tot = basic.get("total_variants", 0)
                # Merge counts across files of the same stage
                stage_total += tot
                for k, v in cls.items():
                    stage_class[k] = stage_class.get(k, 0) + v
        rescue_analysis[stage] = {
            "classification": stage_class,
            "total_variants": stage_total,
        }

    # Print summary statistics
    print(
        f"\n{'Category':<15} {'DNA Consensus':<15} {'RNA Consensus':<15} {'Rescued':<15} {'COSMIC':<15} {'RNA_Edit':<15} {'Filtered':<15}"
    )
    print("-" * 60)
    for filter_cat in CATEGORY_ORDER + ["PASS", "LowQual", "StrandBias", "Clustered", "Other"]:
        dna_count = dna_classification.get(filter_cat, 0)
        rna_count = rna_classification.get(filter_cat, 0)
        rescue_count = rescue_classification.get(filter_cat, 0)
        cg_count = rescue_analysis["cosmic_gnomad"]["classification"].get(filter_cat, 0)
        re_count = rescue_analysis["rna_editing"]["classification"].get(filter_cat, 0)
        flt_count = rescue_analysis["filtered_rescue"]["classification"].get(filter_cat, 0)
        print(f"{filter_cat:<15} {dna_count:<15,} {rna_count:<15,} {rescue_count:<15,} {cg_count:<15,} {re_count:<15,} {flt_count:<15,}")

    print("-" * 60)
    print(f"{'TOTAL':<15} {dna_total:<15,} {rna_total:<15,} {rescue_total:<15,}")
    
    # Add more informative output if all values are 0
    if dna_total == 0 and rna_total == 0 and rescue_total == 0:
        print("\nNote: No consensus or rescue data found. This dataset may not contain these files.")
        print("The rescue analysis requires consensus or rescue VCF files.")
        print("Check if the 'consensus' and 'rescue' directories exist in your data.")
        if "consensus" not in all_vcf_stats or len(all_vcf_stats.get("consensus", {})) == 0:
            print("  - No consensus data found")
        if "rescue" not in all_vcf_stats or len(all_vcf_stats.get("rescue", {})) == 0:
            print("  - No rescue data found")

    # Detailed breakdown for non-zero categories
    print("\nDetailed breakdown:")
    for filter_cat in CATEGORY_ORDER:
        dna_count = dna_classification.get(filter_cat, 0)
        rna_count = rna_classification.get(filter_cat, 0)
        rescue_count = rescue_classification.get(filter_cat, 0)
        cg_count = rescue_analysis["cosmic_gnomad"]["classification"].get(filter_cat, 0)
        re_count = rescue_analysis["rna_editing"]["classification"].get(filter_cat, 0)
        flt_count = rescue_analysis["filtered_rescue"]["classification"].get(filter_cat, 0)

        if any([dna_count, rna_count, rescue_count, cg_count, re_count, flt_count]):
            print(f"\n{filter_cat} Category:")
            print(f"  DNA Consensus: {dna_count:,}")
            print(f"  RNA Consensus: {rna_count:,}")
            print(f"  Rescued: {rescue_count:,}")
            print(f"  COSMIC/GnomAD: {cg_count:,}")
            print(f"  RNA Editing: {re_count:,}")
            print(f"  Filtered: {flt_count:,}")

            # Calculate rescue effectiveness
            max_dna_rna = max(dna_count, rna_count)
            if max_dna_rna > 0:
                rescue_gain = rescue_count - max_dna_rna
                rescue_percent = (rescue_gain / max_dna_rna) * 100
                print(f"  Rescue Gain: {rescue_gain:,} ({rescue_percent:.1f}%)")

    # Calculate transition statistics
    transitions = {}
    for filter_cat in CATEGORY_ORDER:
        dna_count = dna_classification.get(filter_cat, 0)
        rna_count = rna_classification.get(filter_cat, 0)
        rescue_count = rescue_classification.get(filter_cat, 0)

        max_dna_rna = max(dna_count, rna_count)
        if max_dna_rna > 0:
            rescue_gain = rescue_count - max_dna_rna
            transitions[filter_cat] = {
                "dna": dna_count,
                "rna": rna_count,
                "max_input": max_dna_rna,
                "rescue": rescue_count,
                "gain": rescue_gain,
                "gain_percent": (rescue_gain / max_dna_rna) * 100
            }

    rescue_analysis["transitions"] = transitions

    # Create summary statistics
    summary = {
        "dna_total": dna_total,
        "rna_total": rna_total,
        "rescue_total": rescue_total,
        "cosmic_total": rescue_analysis["cosmic_gnomad"].get("total_variants", 0),
        "rna_edit_total": rescue_analysis["rna_editing"].get("total_variants", 0),
        "filtered_total": rescue_analysis["filtered_rescue"].get("total_variants", 0),
        "max_input_total": max(dna_total, rna_total),
        "total_gain": rescue_total - max(dna_total, rna_total),
        "total_gain_percent": ((rescue_total - max(dna_total, rna_total)) / max(dna_total, rna_total)) * 100 if max(dna_total, rna_total) > 0 else 0
    }
    rescue_analysis["summary"] = summary

    print(f"\nOverall Rescue Summary:")
    print(f"  DNA Consensus: {dna_total:,}")
    print(f"  RNA Consensus: {rna_total:,}")
    print(f"  Rescued Total: {rescue_total:,}")
    print(f"  Overall Rescue Gain: {summary['total_gain']:,} ({summary['total_gain_percent']:.1f}%)")

    # Print category color legend for visualization interpretation
    print("\n" + "="*60)
    print("Category Color Legend (for plots):")
    print("="*60)
    for category in CATEGORY_ORDER:
        color = CATEGORY_COLORS.get(category, "#8A8A8A")
        print(f"  {category:<15} : {color}")
    print("="*60)

    # Create visualization if requested
    if show_plot and VISUALIZATION_AVAILABLE:
        try:
            # Create stacked bar chart
            fig = go.Figure()

            # Unified color palette (matches visualizer)
            category_colors = CATEGORY_COLORS.copy()

            # Use stage display names from constants
            stages = [
                (STAGE_DISPLAY_NAMES.get("dna_consensus", "DNA Consensus"), dna_classification, 1),
                (STAGE_DISPLAY_NAMES.get("rna_consensus", "RNA Consensus"), rna_classification, 2),
                (STAGE_DISPLAY_NAMES.get("rescue", "Rescued"), rescue_classification, 3),
                (STAGE_DISPLAY_NAMES.get("cosmic_gnomad", "COSMIC/GnomAD"), rescue_analysis["cosmic_gnomad"]["classification"], 4),
                (STAGE_DISPLAY_NAMES.get("rna_editing", "RNA Editing"), rescue_analysis["rna_editing"]["classification"], 5),
                (STAGE_DISPLAY_NAMES.get("filtered_rescue", "Filtered"), rescue_analysis["filtered_rescue"]["classification"], 6),
            ]

            def _create_rescue_plot(include_no_consensus=True, title_suffix=""):
                """Helper to create rescue analysis plot."""
                fig_rescue = go.Figure()
                
                for stage_name, classification, col_idx in stages:
                    for filter_cat in CATEGORY_ORDER:
                        # Skip NoConsensus if not included
                        if not include_no_consensus and filter_cat == "NoConsensus":
                            continue
                            
                        count = classification.get(filter_cat, 0)
                        if count > 0:
                            fig_rescue.add_trace(
                                go.Bar(
                                    name=filter_cat,
                                    x=[stage_name],
                                    y=[count],
                                    marker_color=category_colors.get(filter_cat, "#8A8A8A"),
                                    showlegend=should_show_legend(col_idx, 1),
                                    legendgroup=filter_cat,
                                )
                            )

                fig_rescue.update_layout(
                    title=f"Rescue Analysis: Pipeline Progression{title_suffix}",
                    xaxis_title="Processing Stage",
                    yaxis_title="Number of Variants",
                    barmode="stack",
                    height=600,
                    width=1000,
                    legend=dict(
                        orientation="h",
                        yanchor="bottom",
                        y=1.02,
                        xanchor="right",
                        x=1
                    )
                )
                return fig_rescue

            # Create main plot with all categories
            fig = _create_rescue_plot(include_no_consensus=True)
            fig.show()
            
            # Create NoConsensus-free plot for better visibility of other categories
            fig_no_nc = _create_rescue_plot(include_no_consensus=False, title_suffix=" (excluding NoConsensus)")
            fig_no_nc.show()

        except Exception as e:
            print(f"Error creating rescue analysis plot: {e}")

    return rescue_analysis


def create_resuce_transition_matrix(rescue_analysis: Dict[str, Any]) -> pd.DataFrame:
    """
    Create a transition matrix showing how variants move through the rescue process.

    Args:
        rescue_analysis: Output from analyze_rescue_vcf function

    Returns:
        DataFrame with transition matrix
    """
    transitions = rescue_analysis.get("transitions", {})

    rows = []
    for category, transition_data in transitions.items():
        rows.append({
            "Category": category,
            "DNA_Count": transition_data["dna"],
            "RNA_Count": transition_data["rna"],
            "Max_Input": transition_data["max_input"],
            "Rescue_Count": transition_data["rescue"],
            "Rescue_Gain": transition_data["gain"],
            "Rescue_Gain_Percent": transition_data["gain_percent"]
        })

    df = pd.DataFrame(rows)
    return df.sort_values("Rescue_Gain", ascending=False)


def compare_rescue_strategies(all_vcf_stats: Dict[str, Any]) -> Dict[str, Any]:
    """
    Compare different rescue strategies if multiple are present.

    Args:
        all_vcf_stats: Dictionary containing all VCF statistics

    Returns:
        Dictionary with strategy comparison results
    """
    comparison = {
        "strategies": {},
        "summary": {}
    }

    # Identify different rescue strategies
    rescue_categories = {}
    for category in all_vcf_stats.keys():
        if "rescue" in category.lower():
            rescue_categories[category] = all_vcf_stats[category]

    if len(rescue_categories) <= 1:
        # Only one rescue strategy, or none
        return {"message": "Only one or no rescue strategies found for comparison"}

    # Analyze each strategy
    for strategy_name, strategy_data in rescue_categories.items():
        strategy_analysis = {
            "total_variants": 0,
            "classification": {},
            "files": list(strategy_data.keys())
        }

        for file_name, file_data in strategy_data.items():
            if "basic" in file_data:
                basic = file_data["basic"]
                strategy_analysis["total_variants"] += basic.get("total_variants", 0)

                # Merge classifications
                file_classification = basic.get("classification", {})
                for class_name, count in file_classification.items():
                    strategy_analysis["classification"][class_name] = (
                        strategy_analysis["classification"].get(class_name, 0) + count
                    )

        comparison["strategies"][strategy_name] = strategy_analysis

    # Create summary comparison
    summary_rows = []
    for strategy_name, strategy_data in comparison["strategies"].items():
        total = strategy_data["total_variants"]
        classification = strategy_data["classification"]

        row = {
            "Strategy": strategy_name,
            "Total_Variants": total,
            "Files": len(strategy_data["files"]),
        }

        # Add classification counts and percentages
        for class_name in CATEGORY_ORDER:
            count = classification.get(class_name, 0)
            row[f"{class_name}_Count"] = count
            row[f"{class_name}_Percent"] = (count / total * 100) if total > 0 else 0

        summary_rows.append(row)

    comparison["summary"] = pd.DataFrame(summary_rows)

    return comparison


def export_rescue_analysis(
    rescue_analysis: Dict[str, Any],
    output_path: str,
    format: str = "excel"
):
    """
    Export rescue analysis results to file.

    Args:
        rescue_analysis: Output from analyze_rescue_vcf function
        output_path: Directory to save the results
        format: Export format ('excel', 'csv', 'both')
    """
    from pathlib import Path
    output_dir = Path(output_path)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Create transition matrix
    transition_df = create_resuce_transition_matrix(rescue_analysis)

    # Create summary data
    summary = rescue_analysis.get("summary", {})
    summary_df = pd.DataFrame([summary])

    if format in ["excel", "both"]:
        excel_path = output_dir / "rescue_analysis.xlsx"
        with pd.ExcelWriter(excel_path, engine='openpyxl') as writer:
            transition_df.to_excel(writer, sheet_name="Transition_Matrix", index=False)
            summary_df.to_excel(writer, sheet_name="Summary", index=False)

            # Add detailed classification data
            for stage in ["dna_consensus", "rna_consensus", "rescue"]:
                if stage in rescue_analysis:
                    stage_data = rescue_analysis[stage]
                    if "classification" in stage_data:
                        classification_df = pd.DataFrame([
                            {"Category": cat, "Count": count}
                            for cat, count in stage_data["classification"].items()
                        ])
                        classification_df.to_excel(writer, sheet_name=f"{stage.title()}_Classification", index=False)

        print(f"✓ Rescue analysis exported to Excel: {excel_path}")

    if format in ["csv", "both"]:
        csv_dir = output_dir / "csv_reports"
        csv_dir.mkdir(exist_ok=True)

        transition_df.to_csv(csv_dir / "transition_matrix.csv", index=False)
        summary_df.to_csv(csv_dir / "summary.csv", index=False)

        print(f"✓ Rescue analysis exported to CSV files in: {csv_dir}")

    return rescue_analysis