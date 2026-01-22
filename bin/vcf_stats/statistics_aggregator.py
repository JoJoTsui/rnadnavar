#!/usr/bin/env python3
"""
VCF Statistics Aggregator Module

Aggregate and summarize VCF statistics across all processing stages:
normalized → consensus → rescue → cosmic_gnomad → rna_editing → filtered_rescue.

Key Features:
    - Category count distribution (Somatic, Germline, Reference, Artifact, RNA_Edit, NoConsensus)
    - Stage-to-stage progression tracking
    - Workflow comparison support (standard vs realignment)
    - Comprehensive summary tables and reports
    - Export capabilities for downstream analysis

Enhanced Features (v2.0):
    - Workflow-aware aggregation (standard and realignment workflows)
    - RNA-focused comparison summaries
    - Integrative view support (DNA + RNA standard + RNA realignment)
    - Realignment impact metrics

Usage Example:
    >>> from vcf_stats.statistics_aggregator import StatisticsAggregator
    >>>
    >>> # Create aggregator for standard workflow
    >>> aggregator = StatisticsAggregator(
    ...     all_stats=standard_workflow_stats,
    ...     workflow_type="standard"
    ... )
    >>>
    >>> # Generate variant count summary
    >>> summary = aggregator.create_variant_count_summary()
    >>> print(summary[['Category', 'Tool', 'Modality', 'Total_Variants']])
    >>>
    >>> # Create workflow comparison (if realignment available)
    >>> comparison = aggregator.create_workflow_comparison_summary(
    ...     standard_stats=standard_workflow_stats,
    ...     realignment_stats=realignment_workflow_stats
    ... )
    >>> print(comparison[['Stage', 'Category', 'Standard_Count',
    ...                    'Realignment_Count', 'Difference']])
    >>>
    >>> # Export results
    >>> aggregator.export_report(output_dir)

Processing Stages:
    1. normalized: Individual variant caller outputs
    2. consensus: Within-modality consensus (DNA, RNA)
    3. rescue: Cross-modality rescue (DNA + RNA)
    4. cosmic_gnomad: COSMIC/gnomAD annotation
    5. rna_editing: RNA editing annotation
    6. filtered_rescue: Final filtered variants

Variant Categories:
    - Somatic: Tumor-specific variants (FILTER=PASS)
    - Germline: Germline variants (FILTER=germline)
    - Reference: Reference calls (FILTER=reference)
    - Artifact: Likely artifacts (other FILTER values)
    - RNA_Edit: RNA editing sites (FILTER=RNA_Edit)
    - NoConsensus: No consensus variants (FILTER=NoConsensus)

Design Principles:
    - Comprehensive: Aggregate all relevant statistics
    - Workflow-aware: Support multiple workflow types
    - Exportable: Generate reports for downstream analysis
    - Flexible: Customizable aggregation and filtering
"""

from pathlib import Path
from typing import Any, Dict, Optional

import pandas as pd

# Import constants from main module
from . import CATEGORY_ORDER, VCF_STAGE_ORDER


class StatisticsAggregator:
    """
    Aggregate and summarize VCF statistics.

    This class provides comprehensive aggregation and summarization of VCF statistics
    across all processing stages and workflow types. It supports both standard and
    realignment workflows, enabling detailed comparison and analysis.

    Key Capabilities:
        - Variant count summaries across all stages
        - Category distribution analysis
        - Stage-to-stage progression tracking
        - Workflow comparison (standard vs realignment)
        - RNA-focused comparison summaries
        - Export to CSV/Excel formats

    Workflow Support:
        - Standard workflow: DNA + RNA samples
        - Realignment workflow: RNA samples only
        - Comparison: RNA standard vs RNA realignment

    Usage Example:
        >>> # Create aggregator
        >>> aggregator = StatisticsAggregator(
        ...     all_stats=workflow_stats,
        ...     workflow_type="standard"
        ... )
        >>>
        >>> # Generate summaries
        >>> variant_summary = aggregator.create_variant_count_summary()
        >>> stage_progression = aggregator.create_stage_progression_summary()
        >>>
        >>> # Workflow comparison (if realignment available)
        >>> comparison = aggregator.create_workflow_comparison_summary(
        ...     standard_stats, realignment_stats
        ... )
        >>>
        >>> # Export results
        >>> aggregator.export_report(output_dir)

    Attributes:
        all_stats: Dictionary containing all VCF statistics
        workflow_type: Type of workflow ("standard" or "realignment")
    """

    def __init__(self, all_stats: Dict[str, Any], workflow_type: str = "standard"):
        """
        Initialize aggregator with statistics data.

        Args:
            all_stats: Dictionary containing all VCF statistics
                Format: {stage: {sample_name: {stats_dict}}}
            workflow_type: Type of workflow ("standard" or "realignment")
                Used for labeling and organizing output

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
            >>> aggregator = StatisticsAggregator(stats, workflow_type="standard")
        """
        self.all_stats = all_stats
        self.workflow_type = workflow_type

    def create_variant_count_summary(self) -> pd.DataFrame:
        """
        Create summary table of variant counts across all VCFs.

        Includes:
        - Total variant counts
        - Variant type breakdown (SNPs, INDELs, OTHER)
        - Category distribution (counts and percentages)
        - Validation: RNAedit must have zero INDEL + OTHER counts

        Returns:
            DataFrame with variant counts and category distributions for each VCF
        """
        # Import variant type constants
        try:
            from vcf_config import VARIANT_TYPE_ORDER
        except ImportError:
            pass

        rows = []

        for category, files in self.all_stats.items():
            for name, data in files.items():
                if "stats" not in data:
                    continue

                basic = data["stats"].get("basic", {})

                # Parse Tool/Modality with stage-aware logic
                if category in {
                    "rescue",
                    "cosmic_gnomad",
                    "rna_editing",
                    "filtered_rescue",
                }:
                    tool = category
                    modality = name
                elif category == "consensus":
                    # For consensus files, tool is "consensus" and modality is the full name
                    tool = "consensus"
                    modality = name
                else:
                    parts = name.split("_")
                    if len(parts) >= 2:
                        tool = parts[0]
                        modality = "_".join(parts[1:])
                    else:
                        tool = name
                        modality = "Unknown"

                # Get 3-category variant type counts
                category_variant_types = basic.get("category_variant_types", {})

                # Calculate total SNP/INDEL/OTHER across all categories
                total_snp = 0
                total_indel = 0
                total_other = 0
                for cat_name, type_counts in category_variant_types.items():
                    total_snp += type_counts.get("SNP", 0)
                    total_indel += type_counts.get("INDEL", 0)
                    total_other += type_counts.get("OTHER", 0)

                row = {
                    "Category": category,
                    "Tool": tool,
                    "Modality": modality,
                    "File": name,
                    "Total_Variants": basic.get("total_variants", 0),
                    "SNP": total_snp,
                    "INDEL": total_indel,
                    "OTHER": total_other,
                }

                # Validation: Check RNAedit has zero INDEL + OTHER
                if "classification" in basic and "RNAedit" in basic["classification"]:
                    rna_edit_types = category_variant_types.get("RNAedit", {})
                    rna_edit_indel = rna_edit_types.get("INDEL", 0)
                    rna_edit_other = rna_edit_types.get("OTHER", 0)
                    if rna_edit_indel > 0 or rna_edit_other > 0:
                        row["RNAedit_VALIDATION_ERROR"] = (
                            f"INDEL={rna_edit_indel}, OTHER={rna_edit_other}"
                        )

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
            "Category",
            "Tool",
            "Modality",
            "File",
            "Total_Variants",
            "SNP",
            "INDEL",
            "OTHER",
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

            # Ensure numeric columns are numeric (avoid NaN for missing cats)
            numeric_cols = [
                c
                for c in expected_cols
                if c not in {"Category", "Tool", "Modality", "File"}
            ]
            for c in numeric_cols:
                df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0)

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
                    "Total_Variants": total,
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

    def create_workflow_comparison_summary(
        self, standard_stats: Dict[str, Any], realignment_stats: Dict[str, Any]
    ) -> pd.DataFrame:
        """
        Compare variant counts between standard and realignment workflows.

        This method creates a comprehensive comparison table showing how variant
        counts differ between the two workflows across all stages and categories.

        Args:
            standard_stats: Statistics from standard workflow
            realignment_stats: Statistics from realignment workflow

        Returns:
            DataFrame with columns:
            - Stage
            - Category
            - Standard_Count
            - Realignment_Count
            - Difference
            - Percent_Change
        """
        from .comparison import WorkflowComparator

        # Create comparator instance
        comparator = WorkflowComparator(standard_stats, realignment_stats)

        # Get RNA category distribution comparison (most comprehensive)
        comparison_df = comparator.compare_rna_category_distribution()

        return comparison_df

    def create_rna_stage_comparison_summary(
        self, standard_stats: Dict[str, Any], realignment_stats: Dict[str, Any]
    ) -> pd.DataFrame:
        """
        Compare RNA annotation stages between standard and realignment workflows.

        Focuses on annotation stages (cosmic_gnomad, rna_editing, filtered_rescue)
        to show how realignment affects variant annotation and filtering.

        Args:
            standard_stats: Statistics from standard workflow
            realignment_stats: Statistics from realignment workflow

        Returns:
            DataFrame with detailed annotation stage comparison
        """
        from .comparison import WorkflowComparator

        # Create comparator instance
        comparator = WorkflowComparator(standard_stats, realignment_stats)

        # Get RNA annotation stage comparison
        annotation_comparison = comparator.compare_rna_annotation_stages()

        return annotation_comparison

    def export_report(
        self,
        output_dir: str,
        format: str = "excel",
        comparison_data: Optional[Dict[str, pd.DataFrame]] = None,
    ):
        """
        Export summary report to files.

        Args:
            output_dir: Directory to save report
            format: Export format ('excel', 'csv', 'both')
            comparison_data: Optional dictionary of comparison DataFrames to include
                           (e.g., from create_workflow_comparison_summary)
        """
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        report = self.create_summary_report()

        # Add comparison data if provided
        if comparison_data:
            report.update(comparison_data)

        if format in ["excel", "both"]:
            # Export to Excel with multiple sheets
            excel_filename = f"vcf_statistics_report_{self.workflow_type}.xlsx"
            excel_path = output_path / excel_filename
            with pd.ExcelWriter(excel_path, engine="openpyxl") as writer:
                for sheet_name, df in report.items():
                    if not df.empty:
                        df.to_excel(writer, sheet_name=sheet_name, index=False)
            print(f"✓ Report exported to Excel: {excel_path}")

        if format in ["csv", "both"]:
            # Export to separate CSV files
            csv_dir = output_path / f"csv_reports_{self.workflow_type}"
            csv_dir.mkdir(exist_ok=True)

            for name, df in report.items():
                if not df.empty:
                    csv_path = csv_dir / f"{name}.csv"
                    df.to_csv(csv_path, index=False)
            print(f"✓ Report exported to CSV files in: {csv_dir}")

        return report

    def aggregate_dna_only_stats(
        self, vcf_files: Dict[str, Dict[str, Any]]
    ) -> pd.DataFrame:
        """
        Aggregate statistics for DNA-only workflow.

        This method handles the two-stage DNA-only pipeline (variant_calling → consensus)
        without expecting rescue or annotation stages. It computes per-caller statistics
        at the variant_calling stage and consensus statistics.

        The DNA-only workflow has only two stages:
        - variant_calling: Individual caller outputs (Strelka, DeepSomatic, Mutect2)
        - consensus: Combined consensus VCF

        Args:
            vcf_files: Dictionary mapping file_id to VCF metadata
                Format: {file_id: {path, stage, tool, sample, file_id}}
                Example:
                    {
                        'mutect2_WES_T_vs_N': {
                            'path': Path('/path/to/vcf.gz'),
                            'stage': 'variant_calling',
                            'tool': 'mutect2',
                            'sample': 'WES_T_vs_N',
                            'file_id': 'mutect2_WES_T_vs_N'
                        },
                        'consensus_WES_T_vs_N': {
                            'path': Path('/path/to/consensus.vcf.gz'),
                            'stage': 'consensus',
                            'tool': 'consensus',
                            'sample': 'WES_T_vs_N',
                            'file_id': 'consensus_WES_T_vs_N'
                        }
                    }

        Returns:
            DataFrame with columns:
            - Stage: Processing stage (variant_calling, consensus)
            - Tool: Variant caller name (strelka, deepsomatic, mutect2, consensus)
            - Sample: Sample pair identifier
            - Total_Variants: Total variant count
            - SNP: SNP count
            - INDEL: INDEL count
            - OTHER: Other variant type count
            - Somatic: Somatic variant count
            - Germline: Germline variant count
            - Reference: Reference call count
            - Artifact: Artifact count
            - Somatic_pct: Somatic percentage
            - Germline_pct: Germline percentage
            - Reference_pct: Reference percentage
            - Artifact_pct: Artifact percentage

        Example:
            >>> aggregator = StatisticsAggregator({}, workflow_type="dna_only")
            >>> vcf_files = {
            ...     'mutect2_sample': {'path': Path('...'), 'stage': 'variant_calling', ...},
            ...     'consensus_sample': {'path': Path('...'), 'stage': 'consensus', ...}
            ... }
            >>> df = aggregator.aggregate_dna_only_stats(vcf_files)
            >>> print(df[['Stage', 'Tool', 'Total_Variants', 'Somatic']])

        Requirements:
            - 5.1: Per-caller statistics at variant_calling stage
            - 5.2: Consensus statistics
            - 5.3: Handle two-stage pipeline without rescue/annotation stages
            - 5.4: Return DataFrame compatible with existing aggregation format
        """
        from .vcf_processor import VCFStatisticsExtractor

        rows = []

        # DNA-only stages in order
        dna_only_stages = ["variant_calling", "consensus"]

        # Group VCF files by stage
        files_by_stage: Dict[str, Dict[str, Any]] = {}
        for file_id, metadata in vcf_files.items():
            stage = metadata.get("stage", "unknown")
            if stage not in files_by_stage:
                files_by_stage[stage] = {}
            files_by_stage[stage][file_id] = metadata

        # Process each stage in order
        for stage in dna_only_stages:
            if stage not in files_by_stage:
                continue

            stage_files = files_by_stage[stage]

            for file_id, metadata in stage_files.items():
                vcf_path = metadata.get("path")
                tool = metadata.get("tool", "unknown")
                sample = metadata.get("sample", "unknown")

                if vcf_path is None:
                    continue

                # Extract statistics from VCF file
                try:
                    extractor = VCFStatisticsExtractor(
                        vcf_path, caller_name=tool, stage=stage
                    )
                    stats = extractor.extract_basic_stats()

                    if stats is None:
                        continue

                    # Get variant type counts from category_variant_types
                    category_variant_types = stats.get("category_variant_types", {})

                    # Calculate total SNP/INDEL/OTHER across all categories
                    total_snp = 0
                    total_indel = 0
                    total_other = 0
                    for cat_name, type_counts in category_variant_types.items():
                        total_snp += type_counts.get("SNP", 0)
                        total_indel += type_counts.get("INDEL", 0)
                        total_other += type_counts.get("OTHER", 0)

                    # Build row with basic info
                    row = {
                        "Stage": stage,
                        "Tool": tool,
                        "Sample": sample,
                        "File_ID": file_id,
                        "Total_Variants": stats.get("total_variants", 0),
                        "SNP": total_snp,
                        "INDEL": total_indel,
                        "OTHER": total_other,
                    }

                    # Add classification counts and percentages
                    classification = stats.get("classification", {})
                    total_variants = stats.get("total_variants", 1)

                    # DNA-only mode uses these categories
                    for cat in CATEGORY_ORDER:
                        count = classification.get(cat, 0)
                        row[cat] = count
                        # Calculate percentage
                        pct = (count / total_variants * 100) if total_variants > 0 else 0
                        row[f"{cat}_pct"] = pct

                    rows.append(row)

                except Exception as e:
                    print(f"Warning: Failed to process {vcf_path}: {e}")
                    continue

        # Define expected columns for consistent output
        expected_cols = [
            "Stage",
            "Tool",
            "Sample",
            "File_ID",
            "Total_Variants",
            "SNP",
            "INDEL",
            "OTHER",
        ]

        # Add category columns and their percentages
        for cat in CATEGORY_ORDER:
            expected_cols.append(cat)
            expected_cols.append(f"{cat}_pct")

        if rows:
            df = pd.DataFrame(rows)

            # Add any missing columns with default values
            for col in expected_cols:
                if col not in df.columns:
                    df[col] = 0

            # Ensure numeric columns are numeric
            numeric_cols = [
                c
                for c in expected_cols
                if c not in {"Stage", "Tool", "Sample", "File_ID"}
            ]
            for c in numeric_cols:
                df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0)

            # Sort by stage order, then by tool
            stage_order_map = {stage: i for i, stage in enumerate(dna_only_stages)}
            df["_stage_order"] = df["Stage"].map(
                lambda x: stage_order_map.get(x, 999)
            )
            df = df.sort_values(["_stage_order", "Tool"]).drop("_stage_order", axis=1)

            return df[expected_cols]
        else:
            # Return empty DataFrame with expected columns
            return pd.DataFrame(columns=expected_cols)


print("✓ Enhanced Statistics Aggregator imported successfully")
