#!/usr/bin/env python3
"""
VCF File Discovery Module

Discovers and organizes all VCF files in the pipeline output directory.
Refactored to support annotation stages: rescue → cosmic_gnomad → rna_editing → filtered_rescue.
"""

from pathlib import Path
from typing import Dict, List, Optional

# Import constants from main module
from . import TOOLS, MODALITIES, VCF_STAGE_ORDER


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
            "normalized": {},              # Normalized VCFs from standalone callers
            "consensus": {},               # Consensus VCFs (DNA & RNA)
            "rescue": {},                  # Rescue VCFs (DNA + RNA combined)
            "cosmic_gnomad": {},           # Cosmic/gnomAD annotated
            "rna_editing": {},             # RNA editing annotated
            "filtered_rescue": {},         # Final filtered rescue VCFs
        }
        self.bam_files = {}

    def discover_vcfs(self) -> Dict[str, Dict[str, Path]]:
        """
        Discover all VCF files organized by processing stage.
        
        Only discovers normalized, consensus, rescue, and annotated (cosmic_gnomad, 
        rna_editing, filtered_rescue) VCFs. Raw variant_calling VCFs are skipped.
        """
        # 1. Normalized VCFs - from standalone callers (primary input for analysis)
        self._discover_normalized_vcfs()
        
        # 2. Consensus VCFs - DNA and RNA consensus
        self._discover_consensus_vcfs()
        
        # 3. Rescue VCFs - DNA + RNA combined
        self._discover_rescue_vcfs()
        
        # 4. Annotation stage VCFs (in order):
        #    cosmic_gnomad → rna_editing → filtered_rescue
        self._discover_cosmic_gnomad_vcfs()
        self._discover_rna_editing_vcfs()
        self._discover_filtered_rescue_vcfs()

        return self.vcf_files

    def _discover_normalized_vcfs(self):
        """Discover normalized VCFs from standalone callers."""
        for tool in TOOLS:
            norm_tool_dir = self.base_dir / "normalized" / tool
            if not norm_tool_dir.exists():
                continue
                
            # Try modality-based structure
            for modality in MODALITIES:
                norm_dir = norm_tool_dir / modality
                if norm_dir.exists():
                    vcf_files = list(norm_dir.glob("*.norm.vcf.gz"))
                    if vcf_files:
                        key = f"{tool}_{modality}"
                        self.vcf_files["normalized"][key] = vcf_files[0]
            
            # Also try sample-pair structure
            if norm_tool_dir.exists():
                for subdir in norm_tool_dir.iterdir():
                    if subdir.is_dir() and "_vs_" in subdir.name:
                        vcf_files = list(subdir.glob("*.norm.vcf.gz"))
                        if vcf_files:
                            sample_pair = subdir.name
                            tumor_sample = sample_pair.split("_vs_")[0]
                            modality = self._infer_modality(tumor_sample)
                            
                            key = f"{tool}_{modality}"
                            self.vcf_files["normalized"][key] = vcf_files[0]

    def _discover_consensus_vcfs(self):
        """Discover consensus VCFs (DNA and RNA)."""
        consensus_dir = self.base_dir / "consensus" / "vcf"
        if not consensus_dir.exists():
            return
            
        # Try modality-based structure
        for modality in MODALITIES:
            vcf_dir = consensus_dir / modality
            if vcf_dir.exists():
                vcf_files = list(vcf_dir.glob("*.consensus.vcf.gz"))
                if vcf_files:
                    self.vcf_files["consensus"][modality] = vcf_files[0]
        
        # Also try sample-pair structure
        for subdir in consensus_dir.iterdir():
            if subdir.is_dir() and "_vs_" in subdir.name:
                vcf_files = list(subdir.glob("*.consensus.vcf.gz"))
                if vcf_files:
                    sample_pair = subdir.name
                    tumor_sample = sample_pair.split("_vs_")[0]
                    modality = self._infer_modality(tumor_sample)
                    
                    self.vcf_files["consensus"][modality] = vcf_files[0]

    def _discover_rescue_vcfs(self):
        """Discover rescue VCFs (DNA + RNA combined)."""
        rescue_dir = self.base_dir / "rescue"
        if not rescue_dir.exists():
            return
            
        for subdir in rescue_dir.iterdir():
            if subdir.is_dir():
                vcf_files = list(subdir.glob("*.rescued.vcf.gz"))
                if vcf_files:
                    self.vcf_files["rescue"][subdir.name] = vcf_files[0]

    def _discover_cosmic_gnomad_vcfs(self):
        """
        Discover Cosmic/gnomAD annotated VCFs.
        
        Expected path: annotation/cosmic_gnomad/{sample_pair}/{sample_pair}.rescue.cosmic_gnomad_annotated.final.vcf.gz
        """
        cosmic_dir = self.base_dir / "annotation" / "cosmic_gnomad"
        if not cosmic_dir.exists():
            return
            
        for sample_pair_dir in cosmic_dir.iterdir():
            if sample_pair_dir.is_dir():
                vcf_files = list(sample_pair_dir.glob("*.rescue.cosmic_gnomad_annotated.final.vcf.gz"))
                if vcf_files:
                    # Use sample pair as key
                    self.vcf_files["cosmic_gnomad"][sample_pair_dir.name] = vcf_files[0]

    def _discover_rna_editing_vcfs(self):
        """
        Discover RNA editing annotated VCFs.
        
        Expected path: annotation/rna_editing/{sample_pair}/{sample_pair}.rescue.rna_annotated.vcf.gz
        """
        rna_edit_dir = self.base_dir / "annotation" / "rna_editing"
        if not rna_edit_dir.exists():
            return
            
        for sample_pair_dir in rna_edit_dir.iterdir():
            if sample_pair_dir.is_dir():
                vcf_files = list(sample_pair_dir.glob("*.rescue.rna_annotated.vcf.gz"))
                if vcf_files:
                    self.vcf_files["rna_editing"][sample_pair_dir.name] = vcf_files[0]

    def _discover_filtered_rescue_vcfs(self):
        """
        Discover final filtered rescue VCFs.
        
        Expected path: filtered/{sample_pair}/{sample_pair}.filtered.vcf.gz
        """
        filtered_dir = self.base_dir / "filtered"
        if not filtered_dir.exists():
            return
            
        for sample_pair_dir in filtered_dir.iterdir():
            if sample_pair_dir.is_dir():
                vcf_files = list(sample_pair_dir.glob("*.filtered.vcf.gz"))
                if vcf_files:
                    self.vcf_files["filtered_rescue"][sample_pair_dir.name] = vcf_files[0]

    @staticmethod
    def _infer_modality(tumor_sample: str) -> str:
        """
        Infer modality from tumor sample name suffix.
        
        Args:
            tumor_sample: Sample name (e.g., 2374372DT, 2374372RT)
            
        Returns:
            Modality name (DNA_TUMOR_vs_DNA_NORMAL or RNA_TUMOR_vs_DNA_NORMAL)
        """
        if tumor_sample.endswith("DT"):
            return "DNA_TUMOR_vs_DNA_NORMAL"
        elif tumor_sample.endswith("RT"):
            return "RNA_TUMOR_vs_DNA_NORMAL"
        else:
            return tumor_sample

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

    def get_rescue_files(self) -> Dict[str, Path]:
        """Get only rescue VCF files"""
        return self.vcf_files.get("rescue", {})

    def get_annotation_files(self, stage: str) -> Dict[str, Path]:
        """
        Get VCF files from a specific annotation stage.
        
        Args:
            stage: Annotation stage ('cosmic_gnomad', 'rna_editing', 'filtered_rescue')
            
        Returns:
            Dictionary of VCF files for the stage
        """
        return self.vcf_files.get(stage, {})

    def get_all_annotation_stages(self) -> Dict[str, Dict[str, Path]]:
        """Get all annotation stage VCFs in order"""
        stages = {}
        for stage in VCF_STAGE_ORDER:
            if stage in self.vcf_files and self.vcf_files[stage]:
                stages[stage] = self.vcf_files[stage]
        return stages

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