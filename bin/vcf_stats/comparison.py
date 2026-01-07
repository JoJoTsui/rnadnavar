#!/usr/bin/env python3
"""
Workflow Comparison Module

Provides comprehensive comparison capabilities between standard and realignment workflows.
Focus on RNA-modality comparisons since realignment only applies to RNA samples.

Key Features:
    - RNA-focused variant count comparisons
    - RNA category distribution analysis
    - RNA annotation stage comparisons
    - Integrative view (DNA + RNA standard + RNA realignment)
    - Realignment impact calculation
    - Comprehensive comparison report export

Important Note:
    Realignment workflow only applies to RNA-tumor samples, not DNA samples.
    Therefore, all comparisons focus on RNA modality:
    - Primary comparison: RNA standard vs RNA realignment
    - Integrative view: DNA (standard) + RNA (standard) + RNA (realignment)

Usage Example:
    >>> from vcf_stats.comparison import WorkflowComparator
    >>>
    >>> # Create comparator with statistics from both workflows
    >>> comparator = WorkflowComparator(
    ...     standard_stats=standard_workflow_stats,
    ...     realignment_stats=realignment_workflow_stats
    ... )
    >>>
    >>> # Compare RNA variant counts
    >>> rna_counts = comparator.compare_rna_variant_counts()
    >>> print(rna_counts)
    >>>
    >>> # Compare RNA category distribution
    >>> rna_categories = comparator.compare_rna_category_distribution()
    >>> print(rna_categories)
    >>>
    >>> # Create integrative view (DNA + RNA standard + RNA realignment)
    >>> integrative = comparator.create_integrative_view(
    ...     dna_stats=dna_standard_stats,
    ...     rna_standard_stats=rna_standard_stats,
    ...     rna_realignment_stats=realignment_stats
    ... )
    >>> print(integrative)
    >>>
    >>> # Calculate realignment impact
    >>> impact = comparator.calculate_realignment_impact()
    >>> print(f"Variants added: {impact['rna_variants_added']}")
    >>> print(f"Improvement: {impact['realignment_improvement']:.2f}%")
    >>>
    >>> # Export comprehensive report
    >>> comparator.export_comparison_report(output_dir)

Comparison Metrics:
    - Variant counts: Total variants at each processing stage
    - Category distribution: Somatic, Germline, Reference, Artifact, RNA_Edit, NoConsensus
    - Annotation stage impacts: cosmic_gnomad, rna_editing, filtered_rescue
    - Realignment effectiveness: Net improvement in variant calling
    - VAF improvements: Variant allele frequency changes

Design Principles:
    - RNA-focused: All comparisons focus on RNA modality (realignment scope)
    - Comprehensive: Multiple comparison perspectives (counts, categories, stages)
    - Quantitative: Calculate concrete impact metrics
    - Exportable: Generate reports for downstream analysis
"""

from pathlib import Path
from typing import Any, Dict, Optional

import pandas as pd

# Import constants from main module
from . import CATEGORY_ORDER, VCF_STAGE_ORDER


class WorkflowComparator:
    """
    Compare statistics between standard and realignment workflows.

    This class provides comprehensive comparison capabilities for analyzing the impact
    of realignment on RNA variant calling. Since realignment only applies to RNA samples,
    all comparisons focus on RNA modality.

    Key Comparison Types:
        1. RNA variant counts: Compare total variants at each stage
        2. RNA category distribution: Compare variant categories (Somatic, Germline, etc.)
        3. RNA annotation stages: Compare annotation stage results
        4. Integrative view: Combine DNA + RNA standard + RNA realignment
        5. Realignment impact: Quantify effectiveness of realignment

    Workflow Context:
        - Standard workflow: DNA + RNA samples with standard alignment
        - Realignment workflow: RNA samples only with focused realignment
        - Comparison focus: RNA standard vs RNA realignment

    Usage Example:
        >>> comparator = WorkflowComparator(
        ...     standard_stats=standard_workflow_stats,
        ...     realignment_stats=realignment_workflow_stats
        ... )
        >>>
        >>> # Compare RNA variant counts
        >>> rna_counts = comparator.compare_rna_variant_counts()
        >>> print(rna_counts[['Stage', 'Category', 'Standard_RNA_Count',
        ...                    'Realignment_RNA_Count', 'Difference']])
        >>>
        >>> # Calculate realignment impact
        >>> impact = comparator.calculate_realignment_impact()
        >>> print(f"Net improvement: {impact['rna_variants_added'] - impact['rna_variants_removed']}")

    Attributes:
        standard_stats: Statistics from standard workflow (DNA + RNA)
        realignment_stats: Statistics from realignment workflow (RNA only)
    """

    def __init__(
        self, standard_stats: Dict[str, Any], realignment_stats: Dict[str, Any]
    ):
        """
        Initialize with statistics from both workflows.

        Note: Realignment stats should only contain RNA modality data,
        as realignment only applies to RNA samples.

        Args:
            standard_stats: Statistics from standard workflow (DNA + RNA)
                Format: {stage: {sample_name: {stats_dict}}}
            realignment_stats: Statistics from realignment workflow (RNA only)
                Format: {stage: {sample_name: {stats_dict}}}

        Example:
            >>> standard_stats = {
            ...     'filtered_rescue': {
            ...         'DNA_TUMOR_vs_DNA_NORMAL': {...},
            ...         'RNA_TUMOR_vs_DNA_NORMAL': {...}
            ...     }
            ... }
            >>> realignment_stats = {
            ...     'filtered_rescue': {
            ...         'RNA_TUMOR_realign_vs_DNA_NORMAL': {...}
            ...     }
            ... }
            >>> comparator = WorkflowComparator(standard_stats, realignment_stats)
        """
        self.standard_stats = standard_stats
        self.realignment_stats = realignment_stats

    def _extract_rna_stats(self, stats: Dict[str, Any]) -> Dict[str, Any]:
        """
        Extract RNA-specific statistics from a workflow stats dictionary.

        Args:
            stats: Complete workflow statistics

        Returns:
            Dictionary containing only RNA sample statistics
        """
        rna_stats = {}

        for stage, files in stats.items():
            rna_files = {}
            for name, data in files.items():
                # Check if this is an RNA sample
                # RNA samples can have "RNA" or "RT" in the name
                name_upper = name.upper()
                if "RNA" in name_upper or "_RT" in name or name.startswith("RT"):
                    rna_files[name] = data

            if rna_files:
                rna_stats[stage] = rna_files

        return rna_stats

    def _extract_dna_stats(self, stats: Dict[str, Any]) -> Dict[str, Any]:
        """
        Extract DNA-specific statistics from a workflow stats dictionary.

        Args:
            stats: Complete workflow statistics

        Returns:
            Dictionary containing only DNA sample statistics
        """
        dna_stats = {}

        for stage, files in stats.items():
            dna_files = {}
            for name, data in files.items():
                # Check if this is a DNA sample
                # DNA samples have "DNA" or "DT" in the name, but not RNA/RT
                name_upper = name.upper()
                has_dna = "DNA" in name_upper or "_DT" in name or name.startswith("DT")
                has_rna = (
                    "RNA" in name_upper or "_RT" in name or "RESCUED_RT" in name_upper
                )

                if has_dna and not has_rna:
                    dna_files[name] = data

            if dna_files:
                dna_stats[stage] = dna_files

        return dna_stats

    def _aggregate_stage_counts(
        self, stats: Dict[str, Any], stage: str
    ) -> Dict[str, int]:
        """
        Aggregate variant counts for a specific stage across all samples.

        Args:
            stats: Workflow statistics
            stage: Stage name

        Returns:
            Dictionary with aggregated counts by category
        """
        aggregated = {"total_variants": 0, "snps": 0, "indels": 0, "classification": {}}

        if stage not in stats:
            return aggregated

        for name, data in stats[stage].items():
            if "stats" not in data or "basic" not in data["stats"]:
                continue

            basic = data["stats"]["basic"]

            # Aggregate basic counts
            aggregated["total_variants"] += basic.get("total_variants", 0)
            aggregated["snps"] += basic.get("snps", 0)
            aggregated["indels"] += basic.get("indels", 0)

            # Aggregate classification counts
            classification = basic.get("classification", {})
            for category, count in classification.items():
                if category not in aggregated["classification"]:
                    aggregated["classification"][category] = 0
                aggregated["classification"][category] += count

        return aggregated

    def compare_rna_variant_counts(self) -> pd.DataFrame:
        """
        Compare RNA variant counts between standard and realignment workflows.

        Focuses on RNA_TUMOR samples only, as realignment only applies to RNA.
        Compares counts across all stages.

        Returns:
            DataFrame with columns:
            - Stage
            - Standard_RNA_Count
            - Realignment_RNA_Count
            - Difference (Realignment - Standard)
            - Percent_Change
        """
        # Extract RNA-specific stats
        rna_standard = self._extract_rna_stats(self.standard_stats)
        rna_realignment = self._extract_rna_stats(self.realignment_stats)

        rows = []

        # Iterate through all stages
        for stage in VCF_STAGE_ORDER:
            # Get aggregated counts for this stage
            standard_counts = self._aggregate_stage_counts(rna_standard, stage)
            realignment_counts = self._aggregate_stage_counts(rna_realignment, stage)

            standard_total = standard_counts["total_variants"]
            realignment_total = realignment_counts["total_variants"]

            # Calculate difference and percent change
            difference = realignment_total - standard_total
            if standard_total > 0:
                percent_change = (difference / standard_total) * 100
            else:
                percent_change = 0.0 if realignment_total == 0 else float("inf")

            rows.append(
                {
                    "Stage": stage,
                    "Standard_RNA_Count": standard_total,
                    "Realignment_RNA_Count": realignment_total,
                    "Difference": difference,
                    "Percent_Change": percent_change,
                }
            )

        return pd.DataFrame(rows)

    def compare_rna_category_distribution(self) -> pd.DataFrame:
        """
        Compare RNA category percentages between workflows.

        Returns DataFrame with columns:
        - Stage
        - Category
        - Standard_RNA_Count
        - Realignment_RNA_Count
        - Difference
        - Percent_Change
        """
        # Extract RNA-specific stats
        rna_standard = self._extract_rna_stats(self.standard_stats)
        rna_realignment = self._extract_rna_stats(self.realignment_stats)

        rows = []

        # Iterate through all stages
        for stage in VCF_STAGE_ORDER:
            # Get aggregated counts for this stage
            standard_counts = self._aggregate_stage_counts(rna_standard, stage)
            realignment_counts = self._aggregate_stage_counts(rna_realignment, stage)

            standard_classification = standard_counts["classification"]
            realignment_classification = realignment_counts["classification"]

            # Get all categories present in either workflow
            all_categories = set(standard_classification.keys()) | set(
                realignment_classification.keys()
            )

            # Add rows for each category
            for category in CATEGORY_ORDER:
                if category not in all_categories:
                    continue

                standard_count = standard_classification.get(category, 0)
                realignment_count = realignment_classification.get(category, 0)

                # Calculate difference and percent change
                difference = realignment_count - standard_count
                if standard_count > 0:
                    percent_change = (difference / standard_count) * 100
                else:
                    percent_change = 0.0 if realignment_count == 0 else float("inf")

                rows.append(
                    {
                        "Stage": stage,
                        "Category": category,
                        "Standard_RNA_Count": standard_count,
                        "Realignment_RNA_Count": realignment_count,
                        "Difference": difference,
                        "Percent_Change": percent_change,
                    }
                )

        return pd.DataFrame(rows)

    def compare_rna_annotation_stages(self) -> pd.DataFrame:
        """
        Compare RNA annotation stage results.

        Focus on cosmic_gnomad, rna_editing, filtered_rescue stages
        for RNA samples only.

        Returns:
            DataFrame with detailed annotation stage comparison
        """
        # Focus on annotation stages
        annotation_stages = ["cosmic_gnomad", "rna_editing", "filtered_rescue"]

        # Extract RNA-specific stats
        rna_standard = self._extract_rna_stats(self.standard_stats)
        rna_realignment = self._extract_rna_stats(self.realignment_stats)

        rows = []

        for stage in annotation_stages:
            # Get aggregated counts for this stage
            standard_counts = self._aggregate_stage_counts(rna_standard, stage)
            realignment_counts = self._aggregate_stage_counts(rna_realignment, stage)

            standard_classification = standard_counts["classification"]
            realignment_classification = realignment_counts["classification"]

            # Get all categories
            all_categories = set(standard_classification.keys()) | set(
                realignment_classification.keys()
            )

            for category in CATEGORY_ORDER:
                if category not in all_categories:
                    continue

                standard_count = standard_classification.get(category, 0)
                realignment_count = realignment_classification.get(category, 0)
                difference = realignment_count - standard_count

                if standard_count > 0:
                    percent_change = (difference / standard_count) * 100
                else:
                    percent_change = 0.0 if realignment_count == 0 else float("inf")

                rows.append(
                    {
                        "Stage": stage,
                        "Category": category,
                        "Standard_RNA_Count": standard_count,
                        "Realignment_RNA_Count": realignment_count,
                        "Difference": difference,
                        "Percent_Change": percent_change,
                    }
                )

        return pd.DataFrame(rows)

    def create_integrative_view(
        self,
        dna_stats: Optional[Dict[str, Any]] = None,
        rna_standard_stats: Optional[Dict[str, Any]] = None,
        rna_realignment_stats: Optional[Dict[str, Any]] = None,
    ) -> pd.DataFrame:
        """
        Create integrative comparison of all three modalities:
        - DNA-tumor (standard only)
        - RNA-tumor (standard)
        - RNA-tumor (realignment)

        This provides a comprehensive view of variant calling across
        all samples and workflows.

        Args:
            dna_stats: DNA-specific stats (if None, extracted from standard_stats)
            rna_standard_stats: RNA standard stats (if None, extracted from standard_stats)
            rna_realignment_stats: RNA realignment stats (if None, uses realignment_stats)

        Returns:
            DataFrame with columns:
            - Stage
            - Category
            - DNA_Standard
            - RNA_Standard
            - RNA_Realignment
            - RNA_Difference (Realignment - Standard)
        """
        # Extract stats if not provided
        if dna_stats is None:
            dna_stats = self._extract_dna_stats(self.standard_stats)
        if rna_standard_stats is None:
            rna_standard_stats = self._extract_rna_stats(self.standard_stats)
        if rna_realignment_stats is None:
            rna_realignment_stats = self._extract_rna_stats(self.realignment_stats)

        rows = []

        # Iterate through all stages
        for stage in VCF_STAGE_ORDER:
            # Get aggregated counts for this stage
            dna_counts = self._aggregate_stage_counts(dna_stats, stage)
            rna_std_counts = self._aggregate_stage_counts(rna_standard_stats, stage)
            rna_real_counts = self._aggregate_stage_counts(rna_realignment_stats, stage)

            dna_classification = dna_counts["classification"]
            rna_std_classification = rna_std_counts["classification"]
            rna_real_classification = rna_real_counts["classification"]

            # Get all categories present in any workflow
            all_categories = (
                set(dna_classification.keys())
                | set(rna_std_classification.keys())
                | set(rna_real_classification.keys())
            )

            # Add rows for each category
            for category in CATEGORY_ORDER:
                if category not in all_categories:
                    continue

                dna_count = dna_classification.get(category, 0)
                rna_std_count = rna_std_classification.get(category, 0)
                rna_real_count = rna_real_classification.get(category, 0)

                # Calculate RNA difference
                rna_difference = rna_real_count - rna_std_count

                rows.append(
                    {
                        "Stage": stage,
                        "Category": category,
                        "DNA_Standard": dna_count,
                        "RNA_Standard": rna_std_count,
                        "RNA_Realignment": rna_real_count,
                        "RNA_Difference": rna_difference,
                    }
                )

        return pd.DataFrame(rows)

    def calculate_realignment_impact(self) -> Dict[str, Any]:
        """
        Calculate the impact of realignment on RNA samples.

        Returns:
            {
                "rna_variants_added": int,
                "rna_variants_removed": int,
                "category_changes": {category: delta},
                "stage_impacts": {stage: impact_metrics},
                "realignment_improvement": float  # percentage
            }
        """
        # Extract RNA-specific stats
        rna_standard = self._extract_rna_stats(self.standard_stats)
        rna_realignment = self._extract_rna_stats(self.realignment_stats)

        # Calculate overall impact using final filtered stage
        final_stage = "filtered_rescue"

        standard_final = self._aggregate_stage_counts(rna_standard, final_stage)
        realignment_final = self._aggregate_stage_counts(rna_realignment, final_stage)

        standard_total = standard_final["total_variants"]
        realignment_total = realignment_final["total_variants"]

        # Calculate variants added/removed
        variants_added = max(0, realignment_total - standard_total)
        variants_removed = max(0, standard_total - realignment_total)

        # Calculate category changes
        category_changes = {}
        standard_classification = standard_final["classification"]
        realignment_classification = realignment_final["classification"]

        all_categories = set(standard_classification.keys()) | set(
            realignment_classification.keys()
        )

        for category in all_categories:
            std_count = standard_classification.get(category, 0)
            real_count = realignment_classification.get(category, 0)
            category_changes[category] = real_count - std_count

        # Calculate stage-by-stage impacts
        stage_impacts = {}

        for stage in VCF_STAGE_ORDER:
            std_counts = self._aggregate_stage_counts(rna_standard, stage)
            real_counts = self._aggregate_stage_counts(rna_realignment, stage)

            std_total = std_counts["total_variants"]
            real_total = real_counts["total_variants"]

            difference = real_total - std_total
            if std_total > 0:
                percent_change = (difference / std_total) * 100
            else:
                percent_change = 0.0 if real_total == 0 else float("inf")

            stage_impacts[stage] = {
                "standard_count": std_total,
                "realignment_count": real_total,
                "difference": difference,
                "percent_change": percent_change,
            }

        # Calculate overall improvement metric
        # Positive if realignment found more variants, negative if fewer
        if standard_total > 0:
            realignment_improvement = (
                (realignment_total - standard_total) / standard_total
            ) * 100
        else:
            realignment_improvement = 0.0 if realignment_total == 0 else 100.0

        return {
            "rna_variants_added": variants_added,
            "rna_variants_removed": variants_removed,
            "category_changes": category_changes,
            "stage_impacts": stage_impacts,
            "realignment_improvement": realignment_improvement,
        }

    def export_comparison_report(self, output_dir: Path):
        """
        Export comprehensive comparison report to CSV/Excel.

        Creates separate sheets/files for:
        - RNA standard vs realignment comparison
        - Integrative view (DNA + RNA standard + RNA realignment)
        - Realignment impact summary

        Args:
            output_dir: Directory to save comparison reports
        """
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        # Generate all comparison data
        rna_variant_counts = self.compare_rna_variant_counts()
        rna_category_dist = self.compare_rna_category_distribution()
        rna_annotation = self.compare_rna_annotation_stages()
        integrative_view = self.create_integrative_view()
        realignment_impact = self.calculate_realignment_impact()

        # Convert realignment impact to DataFrame for export
        impact_rows = []

        # Overall metrics
        impact_rows.append(
            {
                "Metric": "RNA Variants Added",
                "Value": realignment_impact["rna_variants_added"],
            }
        )
        impact_rows.append(
            {
                "Metric": "RNA Variants Removed",
                "Value": realignment_impact["rna_variants_removed"],
            }
        )
        impact_rows.append(
            {
                "Metric": "Realignment Improvement (%)",
                "Value": realignment_impact["realignment_improvement"],
            }
        )

        # Category changes
        for category, delta in realignment_impact["category_changes"].items():
            impact_rows.append(
                {"Metric": f"Category Change: {category}", "Value": delta}
            )

        impact_df = pd.DataFrame(impact_rows)

        # Stage impacts DataFrame
        stage_impact_rows = []
        for stage, metrics in realignment_impact["stage_impacts"].items():
            stage_impact_rows.append(
                {
                    "Stage": stage,
                    "Standard_Count": metrics["standard_count"],
                    "Realignment_Count": metrics["realignment_count"],
                    "Difference": metrics["difference"],
                    "Percent_Change": metrics["percent_change"],
                }
            )
        stage_impact_df = pd.DataFrame(stage_impact_rows)

        # Export to Excel with multiple sheets
        excel_path = output_path / "workflow_comparison_report.xlsx"
        with pd.ExcelWriter(excel_path, engine="openpyxl") as writer:
            rna_variant_counts.to_excel(
                writer, sheet_name="RNA_Variant_Counts", index=False
            )
            rna_category_dist.to_excel(
                writer, sheet_name="RNA_Category_Distribution", index=False
            )
            rna_annotation.to_excel(
                writer, sheet_name="RNA_Annotation_Stages", index=False
            )
            integrative_view.to_excel(
                writer, sheet_name="Integrative_View", index=False
            )
            impact_df.to_excel(writer, sheet_name="Realignment_Impact", index=False)
            stage_impact_df.to_excel(writer, sheet_name="Stage_Impacts", index=False)

        print(f"✓ Workflow comparison report exported to Excel: {excel_path}")

        # Also export to separate CSV files
        csv_dir = output_path / "comparison_csv"
        csv_dir.mkdir(exist_ok=True)

        rna_variant_counts.to_csv(csv_dir / "rna_variant_counts.csv", index=False)
        rna_category_dist.to_csv(csv_dir / "rna_category_distribution.csv", index=False)
        rna_annotation.to_csv(csv_dir / "rna_annotation_stages.csv", index=False)
        integrative_view.to_csv(csv_dir / "integrative_view.csv", index=False)
        impact_df.to_csv(csv_dir / "realignment_impact.csv", index=False)
        stage_impact_df.to_csv(csv_dir / "stage_impacts.csv", index=False)

        print(f"✓ Workflow comparison CSV files exported to: {csv_dir}")

        return {
            "rna_variant_counts": rna_variant_counts,
            "rna_category_distribution": rna_category_dist,
            "rna_annotation_stages": rna_annotation,
            "integrative_view": integrative_view,
            "realignment_impact": impact_df,
            "stage_impacts": stage_impact_df,
        }


print("✓ Workflow Comparison module loaded successfully")
