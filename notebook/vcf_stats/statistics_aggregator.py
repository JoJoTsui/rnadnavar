#!/usr/bin/env python3
"""
VCF Statistics Aggregator Module

Clean version of statistics aggregator based on notebook implementation.
"""

from typing import Dict, Any
import pandas as pd
from pathlib import Path

# Import constants from main module
from . import CATEGORY_ORDER


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

        Returns:
            DataFrame with variant counts and classifications for each VCF
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

                # Add classification counts if available
                if "classification" in basic:
                    for class_name, count in basic["classification"].items():
                        row[class_name] = count

                rows.append(row)

        # Create DataFrame with all expected columns
        expected_cols = [
            "Category", "Tool", "Modality", "File",
            "Total_Variants", "SNPs", "Indels"
        ]
        
        # Add classification columns that might be present
        for class_name in CATEGORY_ORDER:
            if class_name not in expected_cols:
                expected_cols.append(class_name)
                
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

    def create_summary_report(self) -> Dict[str, pd.DataFrame]:
        """
        Create a comprehensive summary report with all key analyses.

        Returns:
            Dictionary containing all summary DataFrames
        """
        report = {
            "variant_count_summary": self.create_variant_count_summary()
        }
        
        return report

    def export_report(self, output_dir: str, format: str = "excel"):
        """
        Export summary report to files.

        Args:
            output_dir: Directory to save report
            format: Export format ('excel', 'csv', 'both')
        """
        from pathlib import Path
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


print("✓ Clean Statistics Aggregator imported successfully")