#!/usr/bin/env python3
"""
BAM Validator Module

Validate variants using BAM/CRAM alignment files to check read support
and provide independent validation of VCF results.
"""

from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
import pandas as pd

# Try to import pysam for BAM handling
try:
    import pysam
    BAM_SUPPORT = True
except ImportError:
    BAM_SUPPORT = False
    print("Warning: pysam not available. BAM validation will be disabled.")


class BAMValidator:
    """Validate variants using BAM/CRAM alignment files"""

    def __init__(self, reference_fasta: Optional[str] = None):
        """
        Initialize BAM validator.

        Args:
            reference_fasta: Optional path to reference FASTA file
        """
        self.reference_fasta = reference_fasta

    def validate_variants(
        self,
        vcf_path: Path,
        bam_paths: Dict[str, Path],
        max_variants: int = 100
    ) -> List[Dict[str, Any]]:
        """
        Validate variants by checking read support in BAM files.

        Args:
            vcf_path: Path to VCF file with variants to validate
            bam_paths: Dictionary mapping sample names to BAM paths
            max_variants: Maximum number of variants to validate (for performance)

        Returns:
            List of validation results for each variant
        """
        if not BAM_SUPPORT:
            print("BAM validation not available (pysam not installed)")
            return []

        validation_results = []

        try:
            # Open VCF file
            vcf = pysam.VariantFile(str(vcf_path))

            variant_count = 0
            for variant in vcf:
                if variant_count >= max_variants:
                    break

                variant_count += 1
                validation_result = self._validate_single_variant(variant, bam_paths)
                validation_results.append(validation_result)

            vcf.close()

        except Exception as e:
            print(f"Error validating variants in {vcf_path}: {e}")

        return validation_results

    def _validate_single_variant(
        self,
        variant: Any,
        bam_paths: Dict[str, Path]
    ) -> Dict[str, Any]:
        """
        Validate a single variant against BAM files.

        Args:
            variant: pysam.VariantRecord object
            bam_paths: Dictionary mapping sample names to BAM paths

        Returns:
            Dictionary with validation results
        """
        result = {
            "chrom": variant.chrom,
            "pos": variant.pos,
            "ref": str(variant.ref) if variant.ref else "",
            "alt": str(variant.alts[0]) if variant.alts else "",
            "qual": variant.qual,
            "filter": variant.filter.keys() if variant.filter else [],
            "sample_results": {}
        }

        # Validate against each BAM file
        for sample_name, bam_path in bam_paths.items():
            sample_result = self._validate_variant_in_bam(
                variant, bam_path, sample_name
            )
            result["sample_results"][sample_name] = sample_result

        return result

    def _validate_variant_in_bam(
        self,
        variant: Any,
        bam_path: Path,
        sample_name: str
    ) -> Dict[str, Any]:
        """
        Validate a variant in a single BAM file.

        Args:
            variant: pysam.VariantRecord object
            bam_path: Path to BAM/CRAM file
            sample_name: Sample identifier

        Returns:
            Dictionary with validation results for this sample
        """
        result = {
            "bam_file": bam_path.name,
            "total_depth": 0,
            "ref_depth": 0,
            "alt_depth": 0,
            "alt_fraction": 0.0,
            "support": "unknown",
            "quality": 0.0,
            "mapping_quality": 0.0,
            "read_positions": []
        }

        try:
            # Open BAM file
            with pysam.AlignmentFile(str(bam_path), "rb") as bam:
                # Get reference name
                chrom = variant.chrom
                if chrom.startswith('chr') and any(ref.startswith('chr') for ref in bam.references):
                    pass  # Keep chr prefix
                elif not chrom.startswith('chr') and any(not ref.startswith('chr') for ref in bam.references):
                    pass  # Keep no chr prefix
                elif chrom.startswith('chr') and any(not ref.startswith('chr') for ref in bam.references):
                    chrom = chrom[3:]  # Remove chr prefix
                elif not chrom.startswith('chr') and any(ref.startswith('chr') for ref in bam.references):
                    chrom = f"chr{chrom}"  # Add chr prefix

                # Fetch reads at variant position
                reads = list(bam.fetch(chrom, variant.pos - 1, variant.pos + len(variant.ref) - 1))

                # Process reads
                total_depth = len(reads)
                ref_depth = 0
                alt_depth = 0
                qualities = []
                mapping_qualities = []
                read_positions = []

                for read in reads:
                    # Skip low quality or unmapped reads
                    if read.mapping_quality < 20:
                        continue
                    if read.is_unmapped:
                        continue

                    qualities.append(read.query_qualification or 0)
                    mapping_qualities.append(read.mapping_quality)

                    # Get read position relative to variant
                    read_pos = read.get_reference_positions()
                    if variant.pos - 1 in read_pos:
                        idx = read_pos.index(variant.pos - 1)
                        base = read.query_sequence[idx]
                        read_positions.append(base)

                        # Count reference vs alternate
                        if idx < len(variant.ref) and base == variant.ref:
                            ref_depth += 1
                        elif base == str(variant.alts[0]) if variant.alts else False:
                            alt_depth += 1

                # Calculate metrics
                result["total_depth"] = total_depth
                result["ref_depth"] = ref_depth
                result["alt_depth"] = alt_depth
                result["alt_fraction"] = alt_depth / max(total_depth, 1)
                result["quality"] = sum(qualities) / max(len(qualities), 1)
                result["mapping_quality"] = sum(mapping_qualities) / max(len(mapping_qualities), 1)
                result["read_positions"] = read_positions

                # Determine support
                if total_depth >= 10 and alt_depth >= 3:
                    result["support"] = "supported"
                elif total_depth >= 5 and alt_depth >= 2:
                    result["support"] = "weak_support"
                elif total_depth >= 1 and alt_depth >= 1:
                    result["support"] = "minimal_support"
                else:
                    result["support"] = "unsupported"

        except Exception as e:
            result["error"] = str(e)
            result["support"] = "error"

        return result

    def summarize_validation(self, validation_results: List[Dict]) -> pd.DataFrame:
        """
        Convert validation results to DataFrame for analysis.

        Args:
            validation_results: List of validation result dictionaries

        Returns:
            DataFrame with validation summary
        """
        if not validation_results:
            return pd.DataFrame()

        rows = []

        for result in validation_results:
            base_row = {
                "chrom": result["chrom"],
                "pos": result["pos"],
                "ref": result["ref"],
                "alt": result["alt"],
                "qual": result["qual"],
                "filter": ";".join(result["filter"]) if result["filter"] else "PASS"
            }

            # Add sample-specific results
            for sample_name, sample_result in result["sample_results"].items():
                row = base_row.copy()
                row.update({
                    "sample": sample_name,
                    "bam_file": sample_result["bam_file"],
                    "total_depth": sample_result["total_depth"],
                    "ref_depth": sample_result["ref_depth"],
                    "alt_depth": sample_result["alt_depth"],
                    "alt_fraction": sample_result["alt_fraction"],
                    "support": sample_result["support"],
                    "quality": sample_result["quality"],
                    "mapping_quality": sample_result["mapping_quality"],
                })

                if "error" in sample_result:
                    row["error"] = sample_result["error"]

                rows.append(row)

        df = pd.DataFrame(rows)
        return df

    def analyze_validation_quality(
        self,
        validation_df: pd.DataFrame
    ) -> Dict[str, Any]:
        """
        Analyze validation quality metrics.

        Args:
            validation_df: DataFrame from summarize_validation

        Returns:
            Dictionary with quality analysis
        """
        if validation_df.empty:
            return {}

        analysis = {
            "total_variants": len(validation_df),
            "variants_per_sample": {},
            "support_rates": {},
            "depth_statistics": {},
            "quality_statistics": {}
        }

        # Analysis per sample
        for sample in validation_df["sample"].unique():
            sample_data = validation_df[validation_df["sample"] == sample]

            analysis["variants_per_sample"][sample] = len(sample_data)

            # Support rates
            support_counts = sample_data["support"].value_counts()
            total_sample = len(sample_data)

            analysis["support_rates"][sample] = {
                "supported": support_counts.get("supported", 0) / total_sample * 100,
                "weak_support": support_counts.get("weak_support", 0) / total_sample * 100,
                "minimal_support": support_counts.get("minimal_support", 0) / total_sample * 100,
                "unsupported": support_counts.get("unsupported", 0) / total_sample * 100,
            }

            # Depth statistics
            analysis["depth_statistics"][sample] = {
                "mean_total_depth": sample_data["total_depth"].mean(),
                "median_total_depth": sample_data["total_depth"].median(),
                "mean_alt_depth": sample_data["alt_depth"].mean(),
                "median_alt_depth": sample_data["alt_depth"].median(),
                "mean_alt_fraction": sample_data["alt_fraction"].mean(),
                "median_alt_fraction": sample_data["alt_fraction"].median(),
            }

            # Quality statistics
            analysis["quality_statistics"][sample] = {
                "mean_quality": sample_data["quality"].mean(),
                "median_quality": sample_data["quality"].median(),
                "mean_mapping_quality": sample_data["mapping_quality"].mean(),
                "median_mapping_quality": sample_data["mapping_quality"].median(),
            }

        return analysis

    def export_validation_results(
        self,
        validation_results: List[Dict],
        output_path: str,
        format: str = "excel"
    ):
        """
        Export validation results to file.

        Args:
            validation_results: List of validation result dictionaries
            output_path: Path to save the results
            format: Export format ('excel', 'csv', 'json')
        """
        from pathlib import Path
        output_dir = Path(output_path)
        output_dir.mkdir(parents=True, exist_ok=True)

        # Create DataFrame
        df = self.summarize_validation(validation_results)

        if df.empty:
            print("No validation results to export")
            return

        # Export in requested format
        if format in ["excel", "both"]:
            excel_path = output_dir / "bam_validation_results.xlsx"
            with pd.ExcelWriter(excel_path, engine='openpyxl') as writer:
                df.to_excel(writer, sheet_name="Validation_Results", index=False)

                # Add quality analysis
                analysis = self.analyze_validation_quality(df)
                if analysis:
                    analysis_df = pd.DataFrame([
                        {
                            "Metric": metric_type,
                            "Value": str(metric_value)
                        }
                        for metric_type, metric_data in analysis.items()
                        if isinstance(metric_data, dict)
                        for key, value in metric_data.items()
                    ])
                    analysis_df.to_excel(writer, sheet_name="Quality_Analysis", index=False)

            print(f"✓ Validation results exported to Excel: {excel_path}")

        if format in ["csv", "both"]:
            csv_path = output_dir / "bam_validation_results.csv"
            df.to_csv(csv_path, index=False)
            print(f"✓ Validation results exported to CSV: {csv_path}")

        if format == "json":
            import json
            json_path = output_dir / "bam_validation_results.json"
            with open(json_path, 'w') as f:
                json.dump(validation_results, f, indent=2, default=str)
            print(f"✓ Validation results exported to JSON: {json_path}")

        return df