#!/usr/bin/env python3
"""
VCF File Discovery Module

Discovers and organizes all VCF files in the pipeline output directory.
"""

from pathlib import Path
from typing import Dict, List, Optional

# Import constants from main module
from . import TOOLS, MODALITIES


class VCFFileDiscovery:
    """Discover and organize all VCF files in the pipeline output"""

    def __init__(self, base_dir: Path):
        """
        Initialize VCF file discovery.

        Args:
            base_dir: Base directory of the pipeline output
        """
        self.base_dir = Path(base_dir)
        self.vcf_files = {
            "variant_calling": {},  # Raw tool outputs
            "normalized": {},  # Normalized VCFs
            "annotated": {},  # VEP annotated
            "consensus": {},  # Consensus VCFs
            "rescue": {},  # Rescued variants
        }
        self.bam_files = {}

    def discover_vcfs(self) -> Dict[str, Dict[str, Path]]:
        """Discover all VCF files"""
        # 1. Per-tool variant calling outputs
        for tool in TOOLS:
            for modality in MODALITIES:
                variant_dir = self.base_dir / "variant_calling" / tool / modality
                if variant_dir.exists():
                    vcf_files = list(variant_dir.glob("*.vcf.gz"))
                    # Filter out gVCF files
                    vcf_files = [f for f in vcf_files if ".g.vcf.gz" not in str(f)]
                    if vcf_files:
                        key = f"{tool}_{modality}"
                        self.vcf_files["variant_calling"][key] = vcf_files[0]

        # 2. Normalized VCFs
        normalized_dir = self.base_dir / "normalized"
        if normalized_dir.exists():
            for tool_modality_dir in normalized_dir.iterdir():
                if tool_modality_dir.is_dir():
                    for vcf_file in tool_modality_dir.glob("*.norm.vcf.gz"):
                        key = tool_modality_dir.name
                        self.vcf_files["normalized"][key] = vcf_file

        # 3. Annotated VCFs
        annotated_dir = self.base_dir / "annotated"
        if annotated_dir.exists():
            for tool_modality_dir in annotated_dir.iterdir():
                if tool_modality_dir.is_dir():
                    for vcf_file in tool_modality_dir.glob("*.vcf.gz"):
                        key = tool_modality_dir.name
                        self.vcf_files["annotated"][key] = vcf_file

        # 4. Consensus VCFs
        consensus_dir = self.base_dir / "consensus"
        if consensus_dir.exists():
            # Also check for consensus in subfolder
            if not consensus_dir.is_file():
                for consensus_subdir in consensus_dir.iterdir():
                    if consensus_subdir.is_dir():
                        for vcf_file in consensus_subdir.glob("*.consensus.vcf.gz"):
                            filename = vcf_file.name
                            # Try to extract tool and modality in different ways
                            # Handle different naming patterns
                            if "DNA_TUMOR_vs_DNA_NORMAL" in filename:
                                if "deepsomatic" in filename:
                                    key = "deepsomatic_DNA_TUMOR_vs_DNA_NORMAL"
                                elif "strelka" in filename:
                                    key = "strelka_DNA_TUMOR_vs_DNA_NORMAL"
                                elif "mutect2" in filename:
                                    key = "mutect2_DNA_TUMOR_vs_DNA_NORMAL"
                                else:
                                    continue
                            elif "RNA_TUMOR_vs_DNA_NORMAL" in filename:
                                if "deepsomatic" in filename:
                                    key = "deepsomatic_RNA_TUMOR_vs_DNA_NORMAL"
                                elif "strelka" in filename:
                                    key = "strelka_RNA_TUMOR_vs_DNA_NORMAL"
                                elif "mutect2" in filename:
                                    key = "mutect2_RNA_TUMOR_vs_DNA_NORMAL"
                                else:
                                    continue
                            else:
                                # Try generic pattern splitting
                                parts = filename.split('_')
                                if len(parts) >= 3:
                                    tool = parts[0]
                                    modality = '_'.join(parts[1:])
                                    key = f"{tool}_{modality}"
                                    self.vcf_files["consensus"][key] = vcf_file

        # 5. Rescue VCFs
        rescue_dir = self.base_dir / "rescue"
        if rescue_dir.exists():
            for vcf_file in rescue_dir.glob("*.rescued.vcf.gz"):
                # Try multiple filename patterns
                filename = vcf_file.name
                
                # Simple pattern recognition for rescue files
                if ".rescued.vcf.gz" in filename:
                    # Remove .rescued.vcf.gz to get base name
                    base_name = filename.replace(".rescued.vcf.gz", "")
                    
                    # Try to determine tool and modality
                    if "DNA" in base_name:
                        key = "DNA"
                    elif "RNA" in base_name:
                        key = "RNA"
                    else:
                        # Fall back to full name as key
                        key = base_name
                    
                    # Create rescue entry
                    self.vcf_files["rescue"][key] = vcf_file

        

        return self.vcf_files

    def discover_alignments(self) -> Dict[str, Path]:
        """Discover alignment files (CRAM/BAM)"""
        recal_dir = self.base_dir / "preprocessing" / "recalibrated"

        if recal_dir.exists():
            for sample_dir in recal_dir.iterdir():
                if sample_dir.is_dir():
                    sample_name = sample_dir.name

                    # Look for CRAM first, then BAM
                    cram_files = list(sample_dir.glob("*.cram"))
                    bam_files = list(sample_dir.glob("*.bam"))

                    if cram_files:
                        self.bam_files[sample_name] = cram_files[0]
                    elif bam_files:
                        self.bam_files[sample_name] = bam_files[0]

        return self.bam_files

    def print_summary(self):
        """Print discovery summary"""
        print("=" * 80)
        print("VCF FILE DISCOVERY SUMMARY")
        print("=" * 80)

        for category, files in self.vcf_files.items():
            if files:
                print(f"\n{category.upper()} VCFs ({len(files)} files):")
                for name, path in files.items():
                    print(f"  {name}: {path.name}")

        if self.bam_files:
            print(f"\nALIGNMENT FILES:")
            for sample, bam_path in self.bam_files.items():
                print(f"  {sample}: {bam_path.name}")

    def get_consensus_files(self) -> Dict[str, Path]:
        """Get only consensus VCF files"""
        return self.vcf_files.get("consensus", {})

    def get_variant_calling_files(self) -> Dict[str, Path]:
        """Get only raw variant calling output files"""
        return self.vcf_files.get("variant_calling", {})

    def get_rescue_files(self) -> Dict[str, Path]:
        """Get only rescue VCF files"""
        return self.vcf_files.get("rescue", {})

    def get_all_files_by_tool(self, tool: str) -> Dict[str, Dict[str, Path]]:
        """
        Get all files for a specific tool across all categories.

        Args:
            tool: Tool name (e.g., 'strelka', 'deepsomatic', 'mutect2')

        Returns:
            Dictionary of files organized by category
        """
        tool_files = {}

        for category, files in self.vcf_files.items():
            category_tool_files = {}
            for name, path in files.items():
                if tool.lower() in name.lower():
                    category_tool_files[name] = path

            if category_tool_files:
                tool_files[category] = category_tool_files

        return tool_files

    def get_all_files_by_modality(self, modality: str) -> Dict[str, Dict[str, Path]]:
        """
        Get all files for a specific modality across all categories.

        Args:
            modality: Modality (e.g., 'DNA_TUMOR', 'RNA_TUMOR')

        Returns:
            Dictionary of files organized by category
        """
        modality_files = {}

        for category, files in self.vcf_files.items():
            category_modality_files = {}
            for name, path in files.items():
                if modality in name:
                    category_modality_files[name] = path

            if category_modality_files:
                modality_files[category] = category_modality_files

        return modality_files

    def find_matching_files(self, tool: str, modality: str) -> Dict[str, Optional[Path]]:
        """
        Find files matching specific tool and modality across all categories.

        Args:
            tool: Tool name
            modality: Modality

        Returns:
            Dictionary with file paths or None if not found
        """
        matching_files = {}

        for category, files in self.vcf_files.items():
            key = f"{tool}_{modality}"
            matching_files[category] = files.get(key)

        return matching_files