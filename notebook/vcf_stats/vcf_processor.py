#!/usr/bin/env python3
"""
VCF Statistics Processor Module

This file is a clean implementation extracted directly from the notebook
to serve as a working replacement for the refactored module with syntax errors.
"""

from pathlib import Path
from typing import Dict, Optional, Any, Tuple
import pandas as pd
import traceback
from collections import defaultdict

# Import classification functions
from .classifiers import (
    classify_variant,
    get_sample_indices
)

# Import constants from main module
from . import TOOLS, MODALITIES, CATEGORY_ORDER

# Try to import optional dependencies
try:
    import matplotlib.pyplot as plt
    import seaborn as sns
    import plotly.express as px
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    VISUALIZATION_AVAILABLE = True
except ImportError:
    VISUALIZATION_AVAILABLE = False


class VCFStatisticsExtractor:
    """
    Extract comprehensive statistics from VCF files with biological classification
    and FILTER category tracking.
    """

    def __init__(self, vcf_path: Path, caller_name: str = None):
        """
        Initialize VCF statistics extractor.

        Args:
            vcf_path: Path to VCF file
            caller_name: Name of the variant caller ('strelka', 'deepsomatic', 'mutect2')
        """
        self.vcf_path = Path(vcf_path)
        self.caller_name = caller_name
        self.vcf = None
        self.stats = {}

    # ...existing code...
    def _is_consensus_or_rescue(self) -> bool:
        """
        Determine if the current VCF is a consensus or rescue VCF.
        """
        p = str(self.vcf_path).lower()
        return (
            ".consensus.vcf.gz" in p
            or ".rescued.vcf.gz" in p
            or ".consensus.vcf" in p
            or ".rescued.vcf" in p
        )

    # ...existing code...
    def extract_basic_stats(self):
        """
        Extract basic variant statistics with biological classification or FILTER category.
        """
        try:
            from cyvcf2 import VCF
            self.vcf = VCF(str(self.vcf_path))

            # Get sample indices for classification
            sample_indices = get_sample_indices(self.vcf, self.caller_name)

            stats = {
                "total_variants": 0,
                "snps": 0,
                "indels": 0,
                "mnps": 0,
                "complex": 0,
                "passed": 0,
                "filtered": 0,
                "chromosomes": set(),
                "qualities": [],
                "variant_types": defaultdict(int),
                # Classification or category (based on FILTER for consensus/rescue)
                "classification": {},
            }

            use_filter_as_category = self._is_consensus_or_rescue()

            for variant in self.vcf:
                stats["total_variants"] += 1
                stats["chromosomes"].add(variant.CHROM)

                # Quality scores
                if variant.QUAL is not None and variant.QUAL > 0:
                    stats["qualities"].append(variant.QUAL)

                # Filter status
                if (
                    variant.FILTER is None
                    or variant.FILTER == "PASS"
                    or variant.FILTER == "."
                ):
                    stats["passed"] += 1
                else:
                    stats["filtered"] += 1

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

                # Classification or FILTER-based category
                try:
                    if use_filter_as_category:
                        # Normalize FILTER into unified categories
                        raw_filter = variant.FILTER if variant.FILTER else "PASS"
                        cat = "Artifact" if raw_filter == "NoConsensus" else raw_filter
                        stats["classification"][cat] = (
                            stats["classification"].get(cat, 0) + 1
                        )
                    else:
                        classification = classify_variant(
                            variant, self.caller_name, sample_indices
                        )
                        stats["classification"][classification] = (
                            stats["classification"].get(classification, 0) + 1
                        )
                except Exception:
                    # Fallback
                    fallback_filter = variant.FILTER if variant.FILTER else "Artifact"
                    fallback_cat = (
                        "Artifact"
                        if (use_filter_as_category and fallback_filter == "NoConsensus")
                        else fallback_filter
                    )
                    stats["classification"][fallback_cat] = (
                        stats["classification"].get(fallback_cat, 0) + 1
                    )

            stats["chromosomes"] = sorted(list(stats["chromosomes"]))

            self.stats["basic"] = stats
            return stats

        except Exception as e:
            print(f"Error processing {self.vcf_path}: {e}")
            traceback.print_exc()
            return None

    # ...existing code...
    def extract_info_fields(self):
        """Extract INFO field statistics."""
        try:
            from cyvcf2 import VCF
            # Always reopen VCF to reset the iterator
            self.vcf = VCF(str(self.vcf_path))

            # Get available INFO fields from header
            info_fields = {}
            print(f"  [DEBUG] Starting header parsing...")
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
                print(f"  Error parsing header: {type(header_err).__name__}: {header_err}")
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

            print(f"  [DEBUG] Processed {variant_count} variants, calculating statistics...")

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
            from cyvcf2 import VCF
            import numpy as np
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
    def extract_all_stats(self, verbose: bool = True):
        """
        Extract all statistics from VCF file.
        """
        if verbose:
            print(f"\nProcessing: {self.vcf_path.name}")

        # Extract basic statistics
        basic = self.extract_basic_stats()

        # Extract INFO field statistics
        info = self.extract_info_fields()

        # Extract FORMAT field statistics
        format_stats = self.extract_format_fields()

        all_stats = {
            "basic": basic,
            "info": info,
            "format": format_stats,
            "file_path": str(self.vcf_path),
            "caller_name": self.caller_name
        }

        if verbose and basic:
            print(f"  ✓ Total variants: {basic.get('total_variants', 0):,}")
            print(f"  ✓ SNPs: {basic.get('snps', 0):,}")
            print(f"  ✓ INDELs: {basic.get('indels', 0):,}")

            if basic.get('classification'):
                print(f"  ✓ Classification: {basic['classification']}")
            elif basic.get('filter_categories'):
                print(f"  ✓ Filter categories: {basic['filter_categories']}")

            if basic.get('chromosomes'):
                print(f"  ✓ Chromosomes: {len(basic['chromosomes'])}")

        return all_stats


def process_all_vcfs(vcf_files_dict):
    """
    Process all VCF files and collect statistics.
    """
    all_stats = {}

    for category, files in vcf_files_dict.items():
        if not files:
            continue

        print(f"\n{'=' * 80}")
        print(f"PROCESSING: {category.upper()}")
        print(f"{'=' * 80}")

        category_stats = {}

        for tool_modality, vcf_path in files.items():
            # Extract caller name from tool/modality string
            caller_name = None
            for tool in TOOLS:
                if tool in tool_modality.lower():
                    caller_name = tool
                    break

            # Extract statistics
            extractor = VCFStatisticsExtractor(vcf_path, caller_name)
            stats = extractor.extract_all_stats(verbose=True)

            category_stats[tool_modality] = {
                "path": vcf_path,
                "stats": stats
            }

        all_stats[category] = category_stats

    return all_stats


print("✓ VCF Statistics Extractor (Notebook Version) loaded successfully")