#!/usr/bin/env python3
"""
VCF Statistics Aggregator Module

Aggregate and summarize VCF statistics across all processing stages:
normalized → consensus → rescue → cosmic_gnomad → rna_editing → filtered_rescue.

Features:
- Category count distribution (replaces pass/filtered metrics)
- Stage-to-stage progression tracking
- Support for new categories (RNA_Edit, NoConsensus)
"""

from typing import Dict, Any
import pandas as pd
from pathlib import Path

# Import constants from main module
from . import CATEGORY_ORDER, VCF_STAGE_ORDER


class StatisticsAggregator:
    """Aggregate and summarize VCF statistics"""

    def __init__(self, all_stats: Dict[str, Any]):
        """
        Initialize aggregator with statistics data.

        Args:
            all_stats: Dictionary containing all VCF statistics
        """
        self.all_stats = all_stats

    def create_variant_count_summary(self) -> pd.DataFrame:
        """
        Create summary table of variant counts across all VCFs.
        
        Includes:
        - Total variant counts
        - Variant type breakdown (SNPs, Indels)
        - Category distribution (counts and percentages)
        - Removed: pass/filtered counts
        
        Returns:
            DataFrame with variant counts and category distributions for each VCF
        """
        rows = []

        for category, files in self.all_stats.items():
            for name, data in files.items():
                if "stats" not in data:
                    continue

                basic = data["stats"].get("basic", {})

                # Parse tool and modality from name
                parts = name.split("_")
                if len(parts) >= 2:
                    tool = parts[0]
                    modality = "_".join(parts[1:])
                else:
                    tool = name
                    modality = "Unknown"

                row = {
                    "Category": category,
                    "Tool": tool,
                    "Modality": modality,
                    "File": name,
                    "Total_Variants": basic.get("total_variants", 0),
                    "SNPs": basic.get("snps", 0),
                    "Indels": basic.get("indels", 0),
                }

                # Add classification counts for each category
                if "classification" in basic:
                    for class_name, count in basic["classification"].items():
                        row[class_name] = count
                        
                        # Calculate percentage for this category
                        total = basic.get("total_variants", 1)
                        pct_col = f"{class_name}_pct"
                        row[pct_col] = (count / total * 100) if total > 0 else 0

                rows.append(row)

        # Create DataFrame with all expected columns
        expected_cols = [
            "Category", "Tool", "Modality", "File",
            "Total_Variants", "SNPs", "Indels"
        ]
        
        # Add category columns and their percentages
        for class_name in CATEGORY_ORDER:
            if class_name not in expected_cols:
                expected_cols.append(class_name)
            # Add percentage column
            pct_col = f"{class_name}_pct"
            if pct_col not in expected_cols:
                expected_cols.append(pct_col)
                
        if rows:
            df = pd.DataFrame(rows)
            
            # Add any missing columns with default values
            for col in expected_cols:
                if col not in df.columns:
                    df[col] = 0
                    
            return df[expected_cols]
        else:
            # Return empty DataFrame with expected columns
            return pd.DataFrame(columns=expected_cols)

    def create_category_distribution_summary(self) -> pd.DataFrame:
        """
        Create detailed category distribution summary.
        
        Shows percentage breakdown of each category across all VCFs.
        Useful for understanding variant composition.
        
        Returns:
            DataFrame with category distribution analysis
        """
        rows = []
        
        for category, files in self.all_stats.items():
            for name, data in files.items():
                if "stats" not in data:
                    continue
                    
                basic = data["stats"].get("basic", {})
                total = basic.get("total_variants", 0)
                
                if total == 0:
                    continue
                
                row = {
                    "VCF_Category": category,
                    "VCF_Name": name,
                    "Total_Variants": total
                }
                
                # Add percentage for each category
                classification = basic.get("classification", {})
                for cat in CATEGORY_ORDER:
                    count = classification.get(cat, 0)
                    pct = (count / total * 100) if total > 0 else 0
                    row[f"{cat}_pct"] = pct
                    row[f"{cat}_count"] = count
                
                rows.append(row)
        
        if rows:
            return pd.DataFrame(rows)
        else:
            return pd.DataFrame()

    def create_stage_progression_summary(self) -> pd.DataFrame:
        """
        Create annotation stage progression summary.
        
        Shows how variant counts change through each stage:
        rescue → cosmic_gnomad → rna_editing → filtered_rescue
        
        Returns:
            DataFrame showing variant retention/filtering at each stage
        """
        stage_data = []
        
        for stage in VCF_STAGE_ORDER:
            if stage not in self.all_stats:
                continue
                
            for name, data in self.all_stats[stage].items():
                if "stats" not in data:
                    continue
                    
                basic = data["stats"].get("basic", {})
                classification = basic.get("classification", {})
                
                row = {
                    "Stage": stage,
                    "VCF_Name": name,
                    "Total_Variants": basic.get("total_variants", 0),
                    "SNPs": basic.get("snps", 0),
                    "Indels": basic.get("indels", 0),
                }
                
                # Add category counts
                for cat in CATEGORY_ORDER:
                    row[cat] = classification.get(cat, 0)
                
                stage_data.append(row)
        
        if stage_data:
            df = pd.DataFrame(stage_data)
            
            # Add progression metrics (if multiple samples)
            if len(df) > 0 and "Stage" in df.columns:
                # Sort by stage order
                stage_order_map = {stage: i for i, stage in enumerate(VCF_STAGE_ORDER)}
                df["Stage_Order"] = df["Stage"].map(stage_order_map)
                df = df.sort_values("Stage_Order").drop("Stage_Order", axis=1)
            
            return df
        else:
            return pd.DataFrame()

    def create_summary_report(self) -> Dict[str, pd.DataFrame]:
        """
        Create a comprehensive summary report with all key analyses.

        Returns:
            Dictionary containing all summary DataFrames:
            - variant_count_summary: Count and category distribution
            - category_distribution: Percentage breakdown by category
            - stage_progression: Variant counts through annotation stages (if available)
        """
        report = {
            "variant_count_summary": self.create_variant_count_summary(),
            "category_distribution": self.create_category_distribution_summary(),
        }
        
        # Add stage progression if annotation stages are available
        stage_prog = self.create_stage_progression_summary()
        if not stage_prog.empty:
            report["stage_progression"] = stage_prog
        
        return report

    def export_report(self, output_dir: str, format: str = "excel"):
        """
        Export summary report to files.

        Args:
            output_dir: Directory to save report
            format: Export format ('excel', 'csv', 'both')
        """
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        report = self.create_summary_report()

        if format in ["excel", "both"]:
            # Export to Excel with multiple sheets
            excel_path = output_path / "vcf_statistics_report.xlsx"
            with pd.ExcelWriter(excel_path, engine='openpyxl') as writer:
                for sheet_name, df in report.items():
                    if not df.empty:
                        df.to_excel(writer, sheet_name=sheet_name, index=False)
            print(f"✓ Report exported to Excel: {excel_path}")

        if format in ["csv", "both"]:
            # Export to separate CSV files
            csv_dir = output_path / "csv_reports"
            csv_dir.mkdir(exist_ok=True)

            for name, df in report.items():
                if not df.empty:
                    csv_path = csv_dir / f"{name}.csv"
                    df.to_csv(csv_path, index=False)
            print(f"✓ Report exported to CSV files in: {csv_dir}")

        return report


print("✓ Enhanced Statistics Aggregator imported successfully")