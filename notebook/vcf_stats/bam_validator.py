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
        max_variants: int = 100,
        stage_hint: Optional[str] = None,
        selected_modalities: Optional[List[str]] = None,
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

        # Apply modality gating based on stage hint or explicit selection
        filtered_bams: Dict[str, Path] = {}
        if selected_modalities:
            allow = set(selected_modalities)
            filtered_bams = {k: v for k, v in bam_paths.items() if k in allow}
        else:
            if stage_hint == "rna_editing":
                allow = {"RNA_TUMOR"}
            elif stage_hint in {"rescue", "filtered_rescue"}:
                allow = {"DNA_TUMOR", "RNA_TUMOR"}
            else:
                allow = set(bam_paths.keys())
            filtered_bams = {k: v for k, v in bam_paths.items() if k in allow}

        try:
            # Open VCF file
            vcf = pysam.VariantFile(str(vcf_path))

            variant_count = 0
            for variant in vcf:
                if variant_count >= max_variants:
                    break

                variant_count += 1
                validation_result = self._validate_single_variant(variant, filtered_bams)
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


class RealignmentBAMValidator:
    """
    Validate realignment filtered VCF against all relevant BAM files.
    
    This validator focuses on comprehensive validation of the final filtered VCF
    across all four samples: DNA_TUMOR, DNA_NORMAL, RNA_TUMOR (standard), and
    RNA_TUMOR_realign (realignment workflow).
    """
    
    def __init__(
        self,
        filtered_vcf_path: Path,
        bam_files: Dict[str, Path],
        reference_fasta: Optional[str] = None
    ):
        """
        Initialize validator with filtered VCF and BAM files.
        
        Args:
            filtered_vcf_path: Path to final realignment filtered VCF
            bam_files: Dictionary mapping sample names to BAM paths:
                {
                    "DNA_TUMOR": Path,
                    "DNA_NORMAL": Path,
                    "RNA_TUMOR": Path,
                    "RNA_TUMOR_realign": Path
                }
            reference_fasta: Optional path to reference FASTA file
        """
        if not BAM_SUPPORT:
            raise ImportError("pysam is required for BAM validation. Please install it.")
        
        self.filtered_vcf_path = Path(filtered_vcf_path)
        self.bam_files = {k: Path(v) for k, v in bam_files.items()}
        self.reference_fasta = reference_fasta
        
        # Validate that all required samples are provided
        required_samples = {"DNA_TUMOR", "DNA_NORMAL", "RNA_TUMOR", "RNA_TUMOR_realign"}
        provided_samples = set(self.bam_files.keys())
        
        if not required_samples.issubset(provided_samples):
            missing = required_samples - provided_samples
            print(f"Warning: Missing BAM files for samples: {missing}")
    
    def validate_all_samples(
        self,
        max_variants: int = 100
    ) -> pd.DataFrame:
        """
        Validate variants across all four samples.
        
        For each variant, check read support in:
        - DNA tumor
        - DNA normal
        - RNA tumor (standard)
        - RNA tumor (realignment)
        
        Args:
            max_variants: Maximum number of variants to validate (for performance)
        
        Returns:
            DataFrame with columns:
            - CHROM, POS, REF, ALT
            - DNA_TUMOR_support, DNA_TUMOR_VAF
            - DNA_NORMAL_support, DNA_NORMAL_VAF
            - RNA_TUMOR_support, RNA_TUMOR_VAF
            - RNA_TUMOR_realign_support, RNA_TUMOR_realign_VAF
            - RNA_VAF_improvement (realignment - standard)
            - Validation_status
        """
        if not BAM_SUPPORT:
            print("BAM validation not available (pysam not installed)")
            return pd.DataFrame()
        
        validation_results = []
        
        try:
            # Open VCF file
            vcf = pysam.VariantFile(str(self.filtered_vcf_path))
            
            variant_count = 0
            for variant in vcf:
                if variant_count >= max_variants:
                    break
                
                variant_count += 1
                result = self._validate_variant_all_samples(variant)
                validation_results.append(result)
            
            vcf.close()
            
        except Exception as e:
            print(f"Error validating variants in {self.filtered_vcf_path}: {e}")
            return pd.DataFrame()
        
        # Convert to DataFrame
        if not validation_results:
            return pd.DataFrame()
        
        df = pd.DataFrame(validation_results)
        
        # Calculate RNA VAF improvement
        if "RNA_TUMOR_VAF" in df.columns and "RNA_TUMOR_realign_VAF" in df.columns:
            df["RNA_VAF_improvement"] = df["RNA_TUMOR_realign_VAF"] - df["RNA_TUMOR_VAF"]
        
        # Determine validation status
        df["Validation_status"] = df.apply(self._determine_validation_status, axis=1)
        
        return df
    
    def _validate_variant_all_samples(self, variant: Any) -> Dict[str, Any]:
        """
        Validate a single variant against all four BAM files.
        
        Args:
            variant: pysam.VariantRecord object
        
        Returns:
            Dictionary with validation results for all samples
        """
        result = {
            "CHROM": variant.chrom,
            "POS": variant.pos,
            "REF": str(variant.ref) if variant.ref else "",
            "ALT": str(variant.alts[0]) if variant.alts else "",
            "QUAL": variant.qual,
            "FILTER": ";".join(variant.filter.keys()) if variant.filter else "PASS"
        }
        
        # Validate against each BAM file
        for sample_name, bam_path in self.bam_files.items():
            sample_result = self._validate_variant_in_bam(variant, bam_path)
            
            # Add sample-specific columns
            result[f"{sample_name}_ref_depth"] = sample_result["ref_depth"]
            result[f"{sample_name}_alt_depth"] = sample_result["alt_depth"]
            result[f"{sample_name}_total_depth"] = sample_result["total_depth"]
            result[f"{sample_name}_VAF"] = sample_result["vaf"]
            result[f"{sample_name}_support"] = sample_result["support"]
        
        return result
    
    def _validate_variant_in_bam(
        self,
        variant: Any,
        bam_path: Path
    ) -> Dict[str, Any]:
        """
        Validate a variant in a single BAM file.
        
        Args:
            variant: pysam.VariantRecord object
            bam_path: Path to BAM/CRAM file
        
        Returns:
            Dictionary with validation results for this sample
        """
        result = {
            "ref_depth": 0,
            "alt_depth": 0,
            "total_depth": 0,
            "vaf": 0.0,
            "support": "unknown"
        }
        
        try:
            # Open BAM file
            with pysam.AlignmentFile(str(bam_path), "rb") as bam:
                # Handle chromosome naming (chr prefix)
                chrom = self._normalize_chromosome_name(variant.chrom, bam.references)
                
                # Fetch reads at variant position
                reads = list(bam.fetch(
                    chrom, 
                    variant.pos - 1, 
                    variant.pos + len(variant.ref)
                ))
                
                # Process reads
                ref_depth = 0
                alt_depth = 0
                total_depth = 0
                
                ref_base = str(variant.ref)
                alt_base = str(variant.alts[0]) if variant.alts else ""
                
                for read in reads:
                    # Skip low quality or unmapped reads
                    if read.mapping_quality < 20:
                        continue
                    if read.is_unmapped:
                        continue
                    
                    total_depth += 1
                    
                    # Get base at variant position
                    try:
                        read_positions = read.get_reference_positions()
                        if variant.pos - 1 in read_positions:
                            idx = read_positions.index(variant.pos - 1)
                            if idx < len(read.query_sequence):
                                base = read.query_sequence[idx]
                                
                                # Count reference vs alternate
                                if base == ref_base:
                                    ref_depth += 1
                                elif base == alt_base:
                                    alt_depth += 1
                    except (IndexError, AttributeError):
                        continue
                
                # Calculate metrics
                result["ref_depth"] = ref_depth
                result["alt_depth"] = alt_depth
                result["total_depth"] = total_depth
                result["vaf"] = alt_depth / max(total_depth, 1)
                
                # Determine support level
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
    
    def _normalize_chromosome_name(self, chrom: str, bam_references: List[str]) -> str:
        """
        Normalize chromosome name to match BAM file references.
        
        Args:
            chrom: Chromosome name from VCF
            bam_references: List of reference names in BAM file
        
        Returns:
            Normalized chromosome name
        """
        # Check if BAM uses chr prefix
        has_chr_prefix = any(ref.startswith('chr') for ref in bam_references)
        
        if chrom.startswith('chr') and not has_chr_prefix:
            # Remove chr prefix
            return chrom[3:]
        elif not chrom.startswith('chr') and has_chr_prefix:
            # Add chr prefix
            return f"chr{chrom}"
        else:
            # Keep as is
            return chrom
    
    def _determine_validation_status(self, row: pd.Series) -> str:
        """
        Determine overall validation status for a variant.
        
        Args:
            row: DataFrame row with validation results
        
        Returns:
            Validation status string
        """
        # Check DNA tumor support (primary evidence)
        dna_tumor_support = row.get("DNA_TUMOR_support", "unknown")
        
        # Check RNA realignment improvement
        rna_improvement = row.get("RNA_VAF_improvement", 0.0)
        rna_realign_support = row.get("RNA_TUMOR_realign_support", "unknown")
        
        # Determine status
        if dna_tumor_support == "supported":
            if rna_improvement > 0.1 and rna_realign_support in ["supported", "weak_support"]:
                return "validated_with_realignment_improvement"
            elif rna_realign_support in ["supported", "weak_support"]:
                return "validated"
            else:
                return "validated_dna_only"
        elif dna_tumor_support in ["weak_support", "minimal_support"]:
            if rna_improvement > 0.1 and rna_realign_support == "supported":
                return "weak_dna_but_strong_realignment"
            else:
                return "weak_validation"
        else:
            return "unsupported"
    
    def get_realignment_improvements(self, validation_df: pd.DataFrame) -> pd.DataFrame:
        """
        Identify variants with improved support in realignment.
        
        Args:
            validation_df: DataFrame from validate_all_samples()
        
        Returns:
            DataFrame with only variants showing realignment improvement
        """
        if validation_df.empty:
            return pd.DataFrame()
        
        # Filter for variants with positive RNA VAF improvement
        improved = validation_df[
            (validation_df["RNA_VAF_improvement"] > 0.05) &
            (validation_df["RNA_TUMOR_realign_support"].isin(["supported", "weak_support"]))
        ].copy()
        
        # Sort by improvement
        improved = improved.sort_values("RNA_VAF_improvement", ascending=False)
        
        return improved
    
    def summarize_validation(self, validation_df: pd.DataFrame) -> Dict[str, Any]:
        """
        Create summary statistics for validation results.
        
        Args:
            validation_df: DataFrame from validate_all_samples()
        
        Returns:
            Dictionary with summary statistics
        """
        if validation_df.empty:
            return {}
        
        summary = {
            "total_variants": len(validation_df),
            "validation_status_counts": validation_df["Validation_status"].value_counts().to_dict(),
            "realignment_improvements": len(validation_df[validation_df["RNA_VAF_improvement"] > 0.05]),
            "sample_statistics": {}
        }
        
        # Statistics per sample
        for sample in ["DNA_TUMOR", "DNA_NORMAL", "RNA_TUMOR", "RNA_TUMOR_realign"]:
            if f"{sample}_VAF" in validation_df.columns:
                summary["sample_statistics"][sample] = {
                    "mean_VAF": validation_df[f"{sample}_VAF"].mean(),
                    "median_VAF": validation_df[f"{sample}_VAF"].median(),
                    "mean_depth": validation_df[f"{sample}_total_depth"].mean(),
                    "median_depth": validation_df[f"{sample}_total_depth"].median(),
                    "supported_count": len(validation_df[validation_df[f"{sample}_support"] == "supported"]),
                    "weak_support_count": len(validation_df[validation_df[f"{sample}_support"] == "weak_support"]),
                    "unsupported_count": len(validation_df[validation_df[f"{sample}_support"] == "unsupported"])
                }
        
        # RNA comparison
        if "RNA_TUMOR_VAF" in validation_df.columns and "RNA_TUMOR_realign_VAF" in validation_df.columns:
            summary["rna_comparison"] = {
                "mean_standard_VAF": validation_df["RNA_TUMOR_VAF"].mean(),
                "mean_realignment_VAF": validation_df["RNA_TUMOR_realign_VAF"].mean(),
                "mean_improvement": validation_df["RNA_VAF_improvement"].mean(),
                "variants_with_improvement": len(validation_df[validation_df["RNA_VAF_improvement"] > 0]),
                "variants_with_significant_improvement": len(validation_df[validation_df["RNA_VAF_improvement"] > 0.1])
            }
        
        return summary
    
    def export_validation_results(
        self,
        validation_df: pd.DataFrame,
        output_dir: Path,
        prefix: str = "realignment_validation"
    ):
        """
        Export validation results to files.
        
        Args:
            validation_df: DataFrame from validate_all_samples()
            output_dir: Directory to save results
            prefix: Prefix for output filenames
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        if validation_df.empty:
            print("No validation results to export")
            return
        
        # Export full results to Excel
        excel_path = output_dir / f"{prefix}_full_results.xlsx"
        with pd.ExcelWriter(excel_path, engine='openpyxl') as writer:
            validation_df.to_excel(writer, sheet_name="All_Variants", index=False)
            
            # Export realignment improvements
            improvements = self.get_realignment_improvements(validation_df)
            if not improvements.empty:
                improvements.to_excel(writer, sheet_name="Realignment_Improvements", index=False)
            
            # Export summary
            summary = self.summarize_validation(validation_df)
            if summary:
                summary_df = pd.DataFrame([
                    {"Metric": k, "Value": str(v)}
                    for k, v in summary.items()
                    if not isinstance(v, dict)
                ])
                summary_df.to_excel(writer, sheet_name="Summary", index=False)
        
        print(f"✓ Validation results exported to: {excel_path}")
        
        # Export CSV for easy analysis
        csv_path = output_dir / f"{prefix}_full_results.csv"
        validation_df.to_csv(csv_path, index=False)
        print(f"✓ Validation results exported to: {csv_path}")
        
        # Export improvements separately
        improvements = self.get_realignment_improvements(validation_df)
        if not improvements.empty:
            improvements_path = output_dir / f"{prefix}_improvements.csv"
            improvements.to_csv(improvements_path, index=False)
            print(f"✓ Realignment improvements exported to: {improvements_path}")