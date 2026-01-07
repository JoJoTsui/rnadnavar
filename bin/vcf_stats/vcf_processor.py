#!/usr/bin/env python3
"""
VCF Statistics Processor Module

Comprehensive VCF statistics extraction with stage-aware classification
and tiering-relevant INFO field parsing.

Supports:
- Stage-aware classification (caller-specific for normalized, FILTER-based for others)
- Tiering INFO field extraction (N_DNA_CALLERS_SUPPORT, N_RNA_CALLERS_SUPPORT, etc.)
- Database evidence extraction (gnomAD_AF, COSMIC_CNT, REDIportal, DARNED)
- FILTER_NORMALIZED_* field extraction for category-aware caller counting
"""

import re
import traceback
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict

# Import constants from main module
from . import TOOLS

# Import classification functions
from .classifiers import classify_by_filter, normalize_category

# Try to import optional dependencies
try:
    import matplotlib.pyplot as plt
    import plotly.express as px
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    VISUALIZATION_AVAILABLE = True
except ImportError:
    VISUALIZATION_AVAILABLE = False


# ============================================================================
# TIERING-RELEVANT INFO FIELDS
# ============================================================================

# Pattern for FILTER_NORMALIZED_* fields
FILTER_NORMALIZED_PATTERN = re.compile(r"FILTER_NORMALIZED_(\w+)_(\w+)")


class VCFStatisticsExtractor:
    """
    Extract comprehensive statistics from VCF files with biological classification,
    FILTER category tracking, and tiering-relevant INFO fields.
    """

    def __init__(self, vcf_path: Path, caller_name: str = None, stage: str = None):
        """
        Initialize VCF statistics extractor.

        Args:
            vcf_path: Path to VCF file
            caller_name: Name of the variant caller ('strelka', 'deepsomatic', 'mutect2')
            stage: VCF processing stage ('normalized', 'consensus', 'rescue', etc.)
        """
        self.vcf_path = Path(vcf_path)
        self.caller_name = caller_name
        self.stage = stage
        self.vcf = None
        self.sample_indices = None  # Will be detected when VCF is opened
        self.stats = {}

    def extract_basic_stats(self):
        """
        Extract basic variant statistics with stage-aware classification.

        Uses stage-aware classify_by_filter() for classification:
        - Normalized stage: caller-specific classification
        - Other stages: FILTER-based classification

        Returns statistics including:
        - Total variants, SNPs, INDELs, MNPs, complex variants
        - Chromosomes
        - Variant types
        - Classification by category (Somatic, Germline, Reference, Artifact, RNAedit, NoConsensus)
        - Category distribution percentages

        Note: Quality scores are not extracted as they are not always available.
        """
        try:
            from cyvcf2 import VCF

            self.vcf = VCF(str(self.vcf_path))

            # Detect sample indices once for normalized stage Strelka classification
            if self.stage == "normalized" and self.caller_name:
                from .classifiers import get_sample_indices

                self.sample_indices = get_sample_indices(self.vcf, self.caller_name)

            stats = {
                "total_variants": 0,
                "snps": 0,
                "indels": 0,
                "mnps": 0,
                "complex": 0,
                "chromosomes": set(),
                "variant_types": defaultdict(int),
                # Classification: count distribution by category
                "classification": {},
                # Category counts (for tracking distribution across categories)
                "category_distribution": {},
            }

            for variant in self.vcf:
                stats["total_variants"] += 1
                stats["chromosomes"].add(variant.CHROM)

                # Variant type
                if variant.is_snp:
                    stats["snps"] += 1
                    stats["variant_types"]["SNP"] += 1
                elif variant.is_indel:
                    stats["indels"] += 1
                    if variant.is_deletion:
                        stats["variant_types"]["DEL"] += 1
                    else:
                        stats["variant_types"]["INS"] += 1
                else:
                    stats["complex"] += 1
                    stats["variant_types"]["COMPLEX"] += 1

                # Stage-aware classification with sample indices for normalized stage
                try:
                    classification = classify_by_filter(
                        variant,
                        stage=self.stage,
                        caller_name=self.caller_name,
                        sample_indices=self.sample_indices,
                    )
                    # Normalize category name for consistency
                    classification = normalize_category(classification)

                    # Count this category
                    stats["classification"][classification] = (
                        stats["classification"].get(classification, 0) + 1
                    )

                except Exception:
                    # Fallback: default to Artifact if classification fails
                    stats["classification"]["Artifact"] = (
                        stats["classification"].get("Artifact", 0) + 1
                    )

            # Convert chromosomes to sorted list
            stats["chromosomes"] = sorted(list(stats["chromosomes"]))

            # Compute category distribution percentages
            total = stats["total_variants"]
            if total > 0:
                for cat, count in stats["classification"].items():
                    stats["category_distribution"][cat] = (count / total) * 100

            self.stats["basic"] = stats
            return stats

        except Exception as e:
            print(f"Error processing {self.vcf_path}: {e}")
            traceback.print_exc()
            return None

    def extract_tiering_fields(self) -> Dict[str, Any]:
        """
        Extract tiering-relevant INFO fields for tier assignment.

        Parses:
        - N_DNA_CALLERS_SUPPORT, N_RNA_CALLERS_SUPPORT
        - gnomAD_AF, COSMIC_CNT
        - REDIportal, DARNED
        - FILTER_NORMALIZED_* fields for category-aware counting

        Returns:
            Dictionary with tiering field statistics
        """
        try:
            from cyvcf2 import VCF

            self.vcf = VCF(str(self.vcf_path))

            tiering_stats = {
                "dna_caller_support": [],
                "rna_caller_support": [],
                "gnomad_af": [],
                "cosmic_cnt": [],
                "has_rediportal": 0,
                "has_darned": 0,
                "filter_normalized_fields": defaultdict(lambda: defaultdict(int)),
                "variants_with_db_support": 0,
            }

            variant_count = 0
            for variant in self.vcf:
                variant_count += 1

                try:
                    info = variant.INFO

                    # Caller support fields
                    dna_support = info.get("N_DNA_CALLERS_SUPPORT")
                    if dna_support is not None:
                        tiering_stats["dna_caller_support"].append(int(dna_support))

                    rna_support = info.get("N_RNA_CALLERS_SUPPORT")
                    if rna_support is not None:
                        tiering_stats["rna_caller_support"].append(int(rna_support))

                    # Database fields
                    gnomad_af = info.get("gnomAD_AF") or info.get("GNOMAD_AF")
                    if gnomad_af is not None:
                        try:
                            tiering_stats["gnomad_af"].append(float(gnomad_af))
                        except (ValueError, TypeError):
                            pass

                    cosmic_cnt = info.get("COSMIC_CNT") or info.get("COSMIC_COUNT")
                    if cosmic_cnt is not None:
                        try:
                            tiering_stats["cosmic_cnt"].append(int(cosmic_cnt))
                        except (ValueError, TypeError):
                            pass

                    # REDIportal / DARNED
                    if info.get("REDIportal") or info.get("REDIPORTAL"):
                        tiering_stats["has_rediportal"] += 1
                    if info.get("DARNED"):
                        tiering_stats["has_darned"] += 1

                    # Track database support
                    has_db = (
                        (
                            gnomad_af is not None and float(gnomad_af) > 0.001
                            if gnomad_af
                            else False
                        )
                        or (
                            cosmic_cnt is not None and int(cosmic_cnt) > 0
                            if cosmic_cnt
                            else False
                        )
                        or info.get("REDIportal")
                        or info.get("DARNED")
                    )
                    if has_db:
                        tiering_stats["variants_with_db_support"] += 1

                    # FILTER_NORMALIZED_* fields
                    info_dict = dict(info)
                    for key, value in info_dict.items():
                        match = FILTER_NORMALIZED_PATTERN.match(key)
                        if match:
                            caller = match.group(1)
                            modality = match.group(2)
                            if value:
                                normalized_val = normalize_category(str(value))
                                tiering_stats["filter_normalized_fields"][
                                    f"{caller}_{modality}"
                                ][normalized_val] += 1

                except Exception:
                    pass

                # Limit for efficiency
                if variant_count > 10000:
                    break

            # Compute summary statistics
            import numpy as np

            tiering_summary = {
                "total_variants_sampled": variant_count,
                "variants_with_db_support": tiering_stats["variants_with_db_support"],
                "variants_with_rediportal": tiering_stats["has_rediportal"],
                "variants_with_darned": tiering_stats["has_darned"],
            }

            if tiering_stats["dna_caller_support"]:
                tiering_summary["dna_support_mean"] = np.mean(
                    tiering_stats["dna_caller_support"]
                )
                tiering_summary["dna_support_median"] = np.median(
                    tiering_stats["dna_caller_support"]
                )
                tiering_summary["dna_support_distribution"] = dict(
                    zip(
                        *np.unique(
                            tiering_stats["dna_caller_support"], return_counts=True
                        )
                    )
                )

            if tiering_stats["rna_caller_support"]:
                tiering_summary["rna_support_mean"] = np.mean(
                    tiering_stats["rna_caller_support"]
                )
                tiering_summary["rna_support_median"] = np.median(
                    tiering_stats["rna_caller_support"]
                )
                tiering_summary["rna_support_distribution"] = dict(
                    zip(
                        *np.unique(
                            tiering_stats["rna_caller_support"], return_counts=True
                        )
                    )
                )

            if tiering_stats["gnomad_af"]:
                tiering_summary["gnomad_af_mean"] = np.mean(tiering_stats["gnomad_af"])
                tiering_summary["gnomad_af_max"] = np.max(tiering_stats["gnomad_af"])

            if tiering_stats["cosmic_cnt"]:
                tiering_summary["cosmic_cnt_mean"] = np.mean(
                    tiering_stats["cosmic_cnt"]
                )
                tiering_summary["cosmic_cnt_max"] = np.max(tiering_stats["cosmic_cnt"])

            # Convert filter_normalized_fields to regular dict
            tiering_summary["filter_normalized_fields"] = {
                k: dict(v) for k, v in tiering_stats["filter_normalized_fields"].items()
            }

            self.stats["tiering"] = tiering_summary
            return tiering_summary

        except Exception as e:
            print(f"Error extracting tiering fields from {self.vcf_path}: {e}")
            traceback.print_exc()
            return {}

    def extract_info_fields(self):
        """Extract INFO field statistics."""
        try:
            from cyvcf2 import VCF

            # Always reopen VCF to reset the iterator
            self.vcf = VCF(str(self.vcf_path))

            # Get available INFO fields from header
            info_fields = {}
            print("  [DEBUG] Starting header parsing...")
            try:
                for key in self.vcf.header_iter():
                    try:
                        # HREC objects: use dictionary-style access like notebook
                        record_type = key["HeaderType"]
                        if record_type == "INFO":
                            field_id = key["ID"]
                            try:
                                field_type = key["Type"]
                            except KeyError:
                                field_type = "String"

                            info_fields[field_id] = {
                                "type": field_type,
                                "values": [],
                            }
                    except (KeyError, AttributeError, TypeError):
                        # Skip this header entry if we can't parse it
                        continue
            except Exception as header_err:
                print(
                    f"  Error parsing header: {type(header_err).__name__}: {header_err}"
                )
                # Continue without header parsing, we'll collect from actual data

            print(f"  [DEBUG] Found {len(info_fields)} INFO fields in header")

            # Collect values
            variant_count = 0
            for variant in self.vcf:
                variant_count += 1
                for info_id in info_fields.keys():
                    try:
                        val = variant.INFO.get(info_id)
                        if val is not None:
                            info_fields[info_id]["values"].append(val)
                    except:
                        pass

                # Limit to first 10000 variants for efficiency
                if variant_count > 10000:
                    break

            print(
                f"  [DEBUG] Processed {variant_count} variants, calculating statistics..."
            )

            # Calculate statistics for numeric fields
            import numpy as np

            info_stats = {}
            for info_id, data in info_fields.items():
                if data["values"]:
                    try:
                        # Try to convert to numeric
                        numeric_vals = []
                        for v in data["values"]:
                            if isinstance(v, (list, tuple)):
                                numeric_vals.extend(
                                    [float(x) for x in v if x is not None]
                                )
                            else:
                                numeric_vals.append(float(v))

                        if numeric_vals:
                            info_stats[info_id] = {
                                "count": len(numeric_vals),
                                "mean": np.mean(numeric_vals),
                                "median": np.median(numeric_vals),
                                "std": np.std(numeric_vals),
                                "min": np.min(numeric_vals),
                                "max": np.max(numeric_vals),
                                "q25": np.percentile(numeric_vals, 25),
                                "q75": np.percentile(numeric_vals, 75),
                            }
                    except (ValueError, TypeError):
                        # Non-numeric field
                        info_stats[info_id] = {
                            "count": len(data["values"]),
                            "type": "categorical",
                        }

            print(f"  [DEBUG] Calculated statistics for {len(info_stats)} INFO fields")
            self.stats["info"] = info_stats
            return info_stats

        except Exception as e:
            print(f"Error extracting INFO fields from {self.vcf_path}:")
            print(f"  {type(e).__name__}: {str(e)}")
            traceback.print_exc()
            return {}

    # ...existing code...
    def extract_format_fields(self):
        """Extract FORMAT field statistics (sample-level)."""
        try:
            import numpy as np
            from cyvcf2 import VCF

            # Always reopen VCF to reset the iterator
            self.vcf = VCF(str(self.vcf_path))

            samples = self.vcf.samples
            format_stats = {sample: {} for sample in samples}

            # Common FORMAT fields to extract
            format_fields = ["DP", "AD", "AF", "GQ"]

            for sample in samples:
                for field in format_fields:
                    format_stats[sample][field] = []

            variant_count = 0
            for variant in self.vcf:
                variant_count += 1

                for i, sample in enumerate(samples):
                    # Depth
                    try:
                        dp = variant.format("DP")[i]
                        if dp is not None and dp[0] > 0:
                            format_stats[sample]["DP"].append(dp[0])
                    except:
                        pass

                    # Allelic depth
                    try:
                        ad = variant.format("AD")[i]
                        if ad is not None:
                            format_stats[sample]["AD"].append(ad)
                    except:
                        pass

                    # Allele frequency
                    try:
                        af = variant.format("AF")[i]
                        if af is not None and af[0] is not None:
                            format_stats[sample]["AF"].append(af[0])
                    except:
                        pass

                    # Genotype quality
                    try:
                        gq = variant.format("GQ")[i]
                        if gq is not None and gq[0] is not None:
                            format_stats[sample]["GQ"].append(gq[0])
                    except:
                        pass

                # Limit for efficiency
                if variant_count > 10000:
                    break

            # Calculate statistics
            format_summary = {}
            for sample, fields in format_stats.items():
                format_summary[sample] = {}
                for field, values in fields.items():
                    if values and field != "AD":
                        format_summary[sample][field] = {
                            "count": len(values),
                            "mean": np.mean(values),
                            "median": np.median(values),
                            "min": np.min(values),
                            "max": np.max(values),
                            "q25": np.percentile(values, 25),
                            "q75": np.percentile(values, 75),
                        }

            self.stats["format"] = format_summary
            return format_summary

        except Exception as e:
            print(f"Error extracting FORMAT fields from {self.vcf_path}: {e}")
            return {}

    # ...existing code...
    def extract_all_stats(self, verbose: bool = True, metadata: dict = None):
        """
        Extract all statistics from VCF file.

        Args:
            verbose: Whether to print progress information
            metadata: Optional metadata dict with stage, tool, sample, file_id
        """
        # Update stage and caller from metadata if provided
        if metadata:
            if metadata.get("stage") and not self.stage:
                self.stage = metadata.get("stage")
            if metadata.get("tool") and not self.caller_name:
                self.caller_name = metadata.get("tool")

        if verbose:
            print(f"\nProcessing: {self.vcf_path.name}")

            # Print metadata if available
            if metadata:
                print("  📋 Metadata:")
                if metadata.get("stage"):
                    print(f"     Stage: {metadata['stage']}")
                if metadata.get("tool"):
                    print(f"     Tool: {metadata['tool']}")
                if metadata.get("sample"):
                    print(f"     Sample: {metadata['sample']}")
                if metadata.get("file_id"):
                    print(f"     File ID: {metadata['file_id']}")

        # Extract basic statistics (with stage-aware classification)
        basic = self.extract_basic_stats()

        # Extract INFO field statistics
        info = self.extract_info_fields()

        # Extract FORMAT field statistics
        format_stats = self.extract_format_fields()

        # Extract tiering-relevant fields for rescue/annotation stages
        tiering = {}
        if self.stage in ("rescue", "cosmic_gnomad", "rna_editing", "filtered_rescue"):
            tiering = self.extract_tiering_fields()

        all_stats = {
            "basic": basic,
            "info": info,
            "format": format_stats,
            "tiering": tiering,
            "file_path": str(self.vcf_path),
            "caller_name": self.caller_name,
            "stage": self.stage,
        }

        if verbose and basic:
            print(f"  ✓ Total variants: {basic.get('total_variants', 0):,}")
            print(f"  ✓ SNPs: {basic.get('snps', 0):,}")
            print(f"  ✓ INDELs: {basic.get('indels', 0):,}")

            if basic.get("classification"):
                print(f"  ✓ Classification: {basic['classification']}")
            elif basic.get("filter_categories"):
                print(f"  ✓ Filter categories: {basic['filter_categories']}")

            if basic.get("chromosomes"):
                print(f"  ✓ Chromosomes: {len(basic['chromosomes'])}")

            if tiering and tiering.get("variants_with_db_support"):
                print(f"  ✓ DB Support: {tiering['variants_with_db_support']} variants")

        return all_stats


def process_all_vcfs(vcf_files_dict):
    """
    Process all VCF files and collect statistics with stage-aware classification.

    Args:
        vcf_files_dict: Dictionary mapping stage -> {name: {path, stage, tool, sample, file_id}}
                       Each entry contains VCF metadata including the file path

    Returns:
        Dictionary mapping stage -> {name: {path, stats, metadata}}
    """
    all_stats = {}

    for category, files in vcf_files_dict.items():
        if not files:
            continue

        print(f"\n{'=' * 80}")
        print(f"PROCESSING: {category.upper()}")
        print(f"{'=' * 80}")

        category_stats = {}

        for tool_modality, vcf_info in files.items():
            # Handle both old format (Path) and new format (dict with metadata)
            if isinstance(vcf_info, dict):
                vcf_path = vcf_info["path"]
                metadata = {
                    "stage": vcf_info.get("stage", category),
                    "tool": vcf_info.get("tool"),
                    "sample": vcf_info.get("sample"),
                    "file_id": vcf_info.get("file_id", tool_modality),
                }
            else:
                # Backward compatibility: treat as Path
                vcf_path = vcf_info
                metadata = {
                    "stage": category,
                    "tool": None,
                    "sample": None,
                    "file_id": tool_modality,
                }

            # Extract caller name from tool/modality string or metadata
            caller_name = metadata.get("tool")
            if not caller_name:
                for tool in TOOLS:
                    if tool in tool_modality.lower():
                        caller_name = tool
                        break

            # Get stage for stage-aware classification
            stage = metadata.get("stage", category)

            # Extract statistics with metadata and stage-aware classification
            extractor = VCFStatisticsExtractor(
                vcf_path, caller_name=caller_name, stage=stage
            )
            stats = extractor.extract_all_stats(verbose=True, metadata=metadata)

            category_stats[tool_modality] = {
                "path": vcf_path,
                "stats": stats,
                "metadata": metadata,
            }

        all_stats[category] = category_stats

    return all_stats


print("✓ VCF Statistics Extractor with stage-aware classification loaded")
