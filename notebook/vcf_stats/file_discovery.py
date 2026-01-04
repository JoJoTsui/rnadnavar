#!/usr/bin/env python3
"""
VCF File Discovery Module

Discovers and organizes all VCF files in the pipeline output directory.
Refactored to support annotation stages: rescue → cosmic_gnomad → rna_editing → filtered_rescue.
Enhanced to support realignment workflow discovery.

Key Features:
    - Automatic VCF file discovery across all processing stages
    - Support for both standard and realignment workflows
    - Workflow-aware file organization
    - BAM/CRAM file discovery for validation
    - Flexible sample naming conventions (DT/DN/RT suffixes)

Supported Workflows:
    Standard Workflow:
        - Base path: pipeline_output/
        - Stages: normalized, consensus, rescue, cosmic_gnomad, rna_editing, filtered_rescue
        - Samples: DNA-tumor, DNA-normal, RNA-tumor
        
    Realignment Workflow:
        - Base path: pipeline_output/vcf_realignment/
        - Stages: normalized, consensus, rescue, cosmic_gnomad, rna_editing, filtered_rescue
        - Samples: RNA-tumor (realigned)

Usage Example:
    >>> from vcf_stats.file_discovery import VCFFileDiscovery
    >>> from vcf_stats.workflow import WorkflowManager
    >>> 
    >>> # Basic discovery (standard workflow only)
    >>> discovery = VCFFileDiscovery("/path/to/pipeline/output")
    >>> vcf_files = discovery.discover_vcfs()
    >>> print(f"Found {len(vcf_files['normalized'])} normalized VCFs")
    >>> 
    >>> # Workflow-aware discovery
    >>> manager = WorkflowManager("/path/to/pipeline/output")
    >>> discovery = VCFFileDiscovery("/path/to/pipeline/output", workflow_manager=manager)
    >>> 
    >>> # Discover files for all workflows
    >>> all_workflows = discovery.discover_all_workflows()
    >>> for workflow_name, stages in all_workflows.items():
    ...     print(f"{workflow_name}: {len(stages)} stages")
    >>> 
    >>> # Discover files for specific workflow
    >>> from vcf_stats.workflow import WorkflowType
    >>> realignment_config = manager.get_workflow_config(WorkflowType.REALIGNMENT)
    >>> realignment_vcfs = discovery.discover_vcfs_for_workflow(realignment_config)
    >>> 
    >>> # Discover BAM files
    >>> bam_files = discovery.discover_bam_files()
    >>> print(f"Found {len(bam_files)} BAM/CRAM files")

Processing Stages:
    1. normalized: VCFs from individual variant callers (DeepSomatic, Mutect2, Strelka)
    2. consensus: Consensus VCFs within modality (DNA consensus, RNA consensus)
    3. rescue: Cross-modality rescue VCFs (DNA rescued by RNA)
    4. cosmic_gnomad: COSMIC/gnomAD annotated VCFs
    5. rna_editing: RNA editing annotated VCFs
    6. filtered_rescue: Final filtered VCFs

Design Principles:
    - Workflow-agnostic: Same discovery logic for all workflow types
    - Configuration-driven: Uses WorkflowConfig for paths
    - Flexible naming: Handles various sample naming conventions
    - Comprehensive: Discovers VCFs and alignment files
"""

from pathlib import Path
from typing import Dict, List, Optional

# Import constants from main module
from . import TOOLS, MODALITIES, VCF_STAGE_ORDER
from .workflow import WorkflowManager, WorkflowConfig, WorkflowType


class VCFFileDiscovery:
    """
    Discover and organize all VCF files in the pipeline output.
    
    This class provides comprehensive VCF file discovery across all processing stages
    and workflow types. It supports both standard and realignment workflows, with
    automatic detection and organization of files.
    
    Key Features:
        - Multi-stage discovery: normalized, consensus, rescue, annotation stages
        - Multi-workflow support: standard and realignment workflows
        - Flexible naming: Handles various sample naming conventions
        - BAM/CRAM discovery: Finds alignment files for validation
        - Organized output: Files grouped by stage and workflow
        
    Discovery Process:
        1. Identify workflow type (standard or realignment)
        2. Scan each processing stage directory
        3. Match VCF files using naming patterns
        4. Organize files by stage and sample
        5. Return structured dictionary of file paths
        
    Usage Example:
        >>> # Basic discovery
        >>> discovery = VCFFileDiscovery("/path/to/pipeline/output")
        >>> vcf_files = discovery.discover_vcfs()
        >>> 
        >>> # Access discovered files
        >>> normalized_vcfs = vcf_files['normalized']
        >>> for name, path in normalized_vcfs.items():
        ...     print(f"{name}: {path}")
        >>> 
        >>> # Workflow-aware discovery
        >>> from vcf_stats.workflow import WorkflowManager, WorkflowType
        >>> manager = WorkflowManager("/path/to/pipeline/output")
        >>> discovery = VCFFileDiscovery("/path/to/pipeline/output", workflow_manager=manager)
        >>> 
        >>> # Discover all workflows
        >>> all_workflows = discovery.discover_all_workflows()
        >>> standard_vcfs = all_workflows['standard']
        >>> realignment_vcfs = all_workflows['realignment']
        >>> 
        >>> # Discover specific workflow
        >>> config = manager.get_workflow_config(WorkflowType.REALIGNMENT)
        >>> realignment_vcfs = discovery.discover_vcfs_for_workflow(config)
        
    Attributes:
        base_dir: Base directory of the pipeline output
        workflow_manager: Optional WorkflowManager for workflow-aware discovery
        vcf_files: Dictionary of discovered VCF files organized by stage
        bam_files: Dictionary of discovered BAM/CRAM files
    """

    def __init__(self, base_dir: Path, workflow_manager: Optional[WorkflowManager] = None):
        """
        Initialize VCF file discovery.

        Args:
            base_dir: Base directory of the pipeline output
            workflow_manager: Optional WorkflowManager for workflow-aware discovery
                If provided, enables multi-workflow discovery capabilities.
                If None, only standard workflow discovery is performed.
                
        Example:
            >>> # Standard discovery only
            >>> discovery = VCFFileDiscovery("/path/to/output")
            >>> 
            >>> # Workflow-aware discovery
            >>> manager = WorkflowManager("/path/to/output")
            >>> discovery = VCFFileDiscovery("/path/to/output", workflow_manager=manager)
        """
        self.base_dir = Path(base_dir)
        self.workflow_manager = workflow_manager or WorkflowManager(base_dir)
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

    @staticmethod
    def _suffix(sample: str) -> Optional[str]:
        """Return DT/DN/RT suffix if present."""
        if not sample:
            return None
        suf = sample[-2:]
        return suf if suf in {"DT", "DN", "RT"} else None

    @classmethod
    def _pair_to_suffix_key(cls, pair: str) -> Optional[str]:
        """
        Convert a sample pair name like '2374372RT_vs_2374372DN' to 'RT_vs_DN'.
        """
        if "_vs_" not in pair:
            return None
        t, n = pair.split("_vs_", 1)
        ts, ns = cls._suffix(t), cls._suffix(n)
        if ts and ns:
            return f"{ts}_vs_{ns}"
        return None

    @staticmethod
    def _modality_to_suffix_pair(modality: str) -> Optional[str]:
        """Map modality token to compact suffix pair DT/RT vs DN."""
        if modality.startswith("DNA_TUMOR"):
            return "DT_vs_DN"
        if modality.startswith("RNA_TUMOR"):
            return "RT_vs_DN"
        return None

    @classmethod
    def _rescued_pair_to_suffix_key(cls, rescued_pair: str) -> Optional[str]:
        """
        Convert '2374372DT_vs_2374372DN_rescued_2374372RT_vs_2374372DN' ->
        'DT_vs_DN_rescued_RT_vs_DN'.
        """
        if "_rescued_" not in rescued_pair:
            # Try simple pair
            return cls._pair_to_suffix_key(rescued_pair)
        left, right = rescued_pair.split("_rescued_", 1)
        left_key = cls._pair_to_suffix_key(left) or ""
        right_key = cls._pair_to_suffix_key(right) or ""
        if left_key and right_key:
            return f"{left_key}_rescued_{right_key}"
        return None

    def _discover_normalized_vcfs(self):
        """Discover normalized VCFs from standalone callers."""
        norm_dir = self.base_dir / "normalized"
        if norm_dir.exists():
            self.vcf_files["normalized"] = self._discover_stage_normalized(norm_dir, "standard")

    def _discover_consensus_vcfs(self):
        """Discover consensus VCFs (DNA and RNA)."""
        consensus_dir = self.base_dir / "consensus"
        if consensus_dir.exists():
            self.vcf_files["consensus"] = self._discover_stage_consensus(consensus_dir, "standard")

    def _discover_rescue_vcfs(self):
        """Discover rescue VCFs (DNA + RNA combined)."""
        rescue_dir = self.base_dir / "rescue"
        if rescue_dir.exists():
            self.vcf_files["rescue"] = self._discover_stage_rescue(rescue_dir, "standard")

    def _discover_cosmic_gnomad_vcfs(self):
        """
        Discover Cosmic/gnomAD annotated VCFs.
        
        Expected path: annotation/cosmic_gnomad/{sample_pair}/{sample_pair}.rescue.cosmic_gnomad_annotated.final.vcf.gz
        """
        cosmic_dir = self.base_dir / "annotation" / "cosmic_gnomad"
        if cosmic_dir.exists():
            self.vcf_files["cosmic_gnomad"] = self._discover_stage_annotation(
                cosmic_dir, "*.rescue.cosmic_gnomad_annotated.final.vcf.gz", "standard"
            )

    def _discover_rna_editing_vcfs(self):
        """
        Discover RNA editing annotated VCFs.
        
        Expected path: annotation/rna_editing/{sample_pair}/{sample_pair}.rescue.rna_annotated.vcf.gz
        """
        rna_edit_dir = self.base_dir / "annotation" / "rna_editing"
        if rna_edit_dir.exists():
            self.vcf_files["rna_editing"] = self._discover_stage_annotation(
                rna_edit_dir, "*.rescue.rna_annotated.vcf.gz", "standard"
            )

    def _discover_filtered_rescue_vcfs(self):
        """
        Discover final filtered rescue VCFs.
        
        Expected path: filtered/{sample_pair}/{sample_pair}.filtered.vcf.gz
        """
        filtered_dir = self.base_dir / "filtered"
        if filtered_dir.exists():
            self.vcf_files["filtered_rescue"] = self._discover_stage_filtered(filtered_dir, "standard")

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
        """Get only rescue VCF files (deduplicated to prefer suffix-based keys)"""
        rescue_files = self.vcf_files.get("rescue", {})
        # Deduplicate: prefer keys like 'DT_vs_DN_rescued_RT_vs_DN' over full sample names
        # If multiple keys point to the same file, keep only the suffix-based one
        seen_paths = {}
        result = {}
        for key, path in rescue_files.items():
            path_str = str(path)
            if path_str in seen_paths:
                # Already seen - prefer shorter/suffix key
                existing_key = seen_paths[path_str]
                if len(key) < len(existing_key):
                    # Replace with shorter key
                    del result[existing_key]
                    result[key] = path
                    seen_paths[path_str] = key
            else:
                result[key] = path
                seen_paths[path_str] = key
        return result

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

        return matching_files

    def discover_vcfs_for_workflow(
        self, 
        workflow_config: WorkflowConfig
    ) -> Dict[str, Dict[str, Path]]:
        """
        Discover VCF files for a specific workflow.
        
        Args:
            workflow_config: Configuration for the workflow to discover
            
        Returns:
            Dictionary organized by stage, with VCF paths:
            {
                "normalized": {name: path},
                "consensus": {name: path},
                "rescue": {name: path},
                "cosmic_gnomad": {name: path},
                "rna_editing": {name: path},
                "filtered_rescue": {name: path}
            }
        """
        workflow_vcfs = {}
        
        # Discover each stage using the workflow configuration
        for stage in workflow_config.stages:
            stage_path = workflow_config.get_stage_path(stage)
            if not stage_path or not stage_path.exists():
                workflow_vcfs[stage] = {}
                continue
                
            # Use stage-specific discovery logic
            if stage == "normalized":
                workflow_vcfs[stage] = self._discover_stage_normalized(stage_path, workflow_config.name)
            elif stage == "consensus":
                workflow_vcfs[stage] = self._discover_stage_consensus(stage_path, workflow_config.name)
            elif stage == "rescue":
                workflow_vcfs[stage] = self._discover_stage_rescue(stage_path, workflow_config.name)
            elif stage == "cosmic_gnomad":
                workflow_vcfs[stage] = self._discover_stage_annotation(
                    stage_path, "*.rescue.cosmic_gnomad_annotated.final.vcf.gz", workflow_config.name
                )
            elif stage == "rna_editing":
                workflow_vcfs[stage] = self._discover_stage_annotation(
                    stage_path, "*.rescue.rna_annotated.vcf.gz", workflow_config.name
                )
            elif stage == "filtered_rescue":
                workflow_vcfs[stage] = self._discover_stage_filtered(stage_path, workflow_config.name)
            else:
                workflow_vcfs[stage] = {}
                
        return workflow_vcfs

    def discover_all_workflows(self) -> Dict[str, Dict[str, Dict[str, Path]]]:
        """
        Discover VCF files for all detected workflows.
        
        Returns:
            Dictionary organized by workflow type and stage:
            {
                "standard": {
                    "normalized": {name: path},
                    "consensus": {name: path},
                    ...
                },
                "realignment": {
                    "normalized": {name: path},
                    "consensus": {name: path},
                    ...
                }
            }
        """
        all_workflows = {}
        
        # Get all detected workflow configurations
        workflow_configs = self.workflow_manager.get_all_configs()
        
        for workflow_type, config in workflow_configs.items():
            workflow_name = config.name
            all_workflows[workflow_name] = self.discover_vcfs_for_workflow(config)
            
        return all_workflows

    def discover_bam_files_for_workflow(
        self,
        workflow_config: WorkflowConfig
    ) -> Dict[str, Path]:
        """
        Discover BAM/CRAM files for a specific workflow.
        
        Args:
            workflow_config: Configuration for the workflow
            
        Returns:
            Dictionary mapping sample names to BAM/CRAM file paths
        """
        bam_files = {}
        
        # Look for preprocessing/recalibrated directory in the workflow base path
        recal_dir = workflow_config.base_path / "preprocessing" / "recalibrated"
        
        if not recal_dir.exists():
            return bam_files
            
        for sample_dir in recal_dir.iterdir():
            if not sample_dir.is_dir():
                continue
                
            sample_name = sample_dir.name
            
            # Add workflow prefix to sample name for realignment
            if workflow_config.name == "realignment":
                sample_name = f"realignment_{sample_name}"
            
            # Look for CRAM first, then BAM
            cram_files = list(sample_dir.glob("*.cram"))
            bam_files_list = list(sample_dir.glob("*.bam"))
            
            if cram_files:
                bam_files[sample_name] = cram_files[0]
            elif bam_files_list:
                bam_files[sample_name] = bam_files_list[0]
                
        return bam_files

    def _discover_stage_normalized(
        self, 
        stage_path: Path, 
        workflow_name: str
    ) -> Dict[str, Path]:
        """
        Discover normalized VCFs in a stage directory.
        
        Args:
            stage_path: Path to the normalized stage directory
            workflow_name: Name of the workflow ("standard" or "realignment")
            
        Returns:
            Dictionary of discovered VCF files
        """
        vcfs = {}
        
        for tool in TOOLS:
            tool_dir = stage_path / tool
            if not tool_dir.exists():
                continue
                
            # Try modality-based structure
            for modality in MODALITIES:
                modality_dir = tool_dir / modality
                if modality_dir.exists():
                    vcf_files = list(modality_dir.glob("*.norm.vcf.gz"))
                    if vcf_files:
                        suffix_pair = self._modality_to_suffix_pair(modality)
                        key = f"{tool}_{suffix_pair}" if suffix_pair else f"{tool}_{modality}"
                        # Add workflow prefix for realignment
                        if workflow_name == "realignment":
                            key = f"realignment_{key}"
                        vcfs[key] = vcf_files[0]
            
            # Also try sample-pair structure
            for subdir in tool_dir.iterdir():
                if subdir.is_dir() and "_vs_" in subdir.name:
                    vcf_files = list(subdir.glob("*.norm.vcf.gz"))
                    if vcf_files:
                        sample_pair = subdir.name
                        pair_key = self._pair_to_suffix_key(sample_pair)
                        modality_key = pair_key if pair_key else self._infer_modality(sample_pair.split("_vs_")[0])
                        key = f"{tool}_{modality_key}"
                        # Add workflow prefix for realignment
                        if workflow_name == "realignment":
                            key = f"realignment_{key}"
                        vcfs[key] = vcf_files[0]
                        
        return vcfs

    def _discover_stage_consensus(
        self, 
        stage_path: Path, 
        workflow_name: str
    ) -> Dict[str, Path]:
        """
        Discover consensus VCFs in a stage directory.
        
        Args:
            stage_path: Path to the consensus stage directory
            workflow_name: Name of the workflow ("standard" or "realignment")
            
        Returns:
            Dictionary of discovered VCF files
        """
        vcfs = {}
        vcf_dir = stage_path / "vcf"
        
        if not vcf_dir.exists():
            return vcfs
            
        # Try modality-based structure
        for modality in MODALITIES:
            modality_vcf_dir = vcf_dir / modality
            if modality_vcf_dir.exists():
                vcf_files = list(modality_vcf_dir.glob("*.consensus.vcf.gz"))
                if vcf_files:
                    suffix_pair = self._modality_to_suffix_pair(modality)
                    key = f"consensus_{suffix_pair}" if suffix_pair else modality
                    # Add workflow prefix for realignment
                    if workflow_name == "realignment":
                        key = f"realignment_{key}"
                    vcfs[key] = vcf_files[0]
        
        # Also try sample-pair structure
        for subdir in vcf_dir.iterdir():
            if subdir.is_dir() and "_vs_" in subdir.name:
                vcf_files = list(subdir.glob("*.consensus.vcf.gz"))
                if vcf_files:
                    sample_pair = subdir.name
                    pair_key = self._pair_to_suffix_key(sample_pair)
                    key = f"consensus_{pair_key}" if pair_key else sample_pair
                    # Add workflow prefix for realignment
                    if workflow_name == "realignment":
                        key = f"realignment_{key}"
                    vcfs[key] = vcf_files[0]
                    
        return vcfs

    def _discover_stage_rescue(
        self, 
        stage_path: Path, 
        workflow_name: str
    ) -> Dict[str, Path]:
        """
        Discover rescue VCFs in a stage directory.
        
        Args:
            stage_path: Path to the rescue stage directory
            workflow_name: Name of the workflow ("standard" or "realignment")
            
        Returns:
            Dictionary of discovered VCF files
        """
        vcfs = {}
        
        if not stage_path.exists():
            return vcfs
            
        for subdir in stage_path.iterdir():
            if subdir.is_dir():
                vcf_files = list(subdir.glob("*.rescued.vcf.gz"))
                if vcf_files:
                    key = self._rescued_pair_to_suffix_key(subdir.name) or subdir.name
                    # Add workflow prefix for realignment
                    if workflow_name == "realignment":
                        key = f"realignment_{key}"
                    vcfs[key] = vcf_files[0]
                    
        return vcfs

    def _discover_stage_annotation(
        self, 
        stage_path: Path, 
        pattern: str,
        workflow_name: str
    ) -> Dict[str, Path]:
        """
        Discover annotation VCFs in a stage directory.
        
        Args:
            stage_path: Path to the annotation stage directory
            pattern: Glob pattern for VCF files
            workflow_name: Name of the workflow ("standard" or "realignment")
            
        Returns:
            Dictionary of discovered VCF files
        """
        vcfs = {}
        
        if not stage_path.exists():
            return vcfs
            
        for sample_pair_dir in stage_path.iterdir():
            if sample_pair_dir.is_dir():
                vcf_files = list(sample_pair_dir.glob(pattern))
                if vcf_files:
                    key = self._rescued_pair_to_suffix_key(sample_pair_dir.name) or sample_pair_dir.name
                    # Add workflow prefix for realignment
                    if workflow_name == "realignment":
                        key = f"realignment_{key}"
                    vcfs[key] = vcf_files[0]
                    
        return vcfs

    def _discover_stage_filtered(
        self, 
        stage_path: Path, 
        workflow_name: str
    ) -> Dict[str, Path]:
        """
        Discover filtered VCFs in a stage directory.
        
        Args:
            stage_path: Path to the filtered stage directory
            workflow_name: Name of the workflow ("standard" or "realignment")
            
        Returns:
            Dictionary of discovered VCF files
        """
        vcfs = {}
        
        if not stage_path.exists():
            return vcfs
            
        for sample_pair_dir in stage_path.iterdir():
            if sample_pair_dir.is_dir():
                vcf_files = list(sample_pair_dir.glob("*.filtered.vcf.gz"))
                if not vcf_files:
                    continue
                # Only include combined rescued filtered result, omit single-modality filtered
                # i.e., keep keys that contain '_rescued_'
                if "_rescued_" not in sample_pair_dir.name:
                    continue
                key = self._rescued_pair_to_suffix_key(sample_pair_dir.name) or sample_pair_dir.name
                # Add workflow prefix for realignment
                if workflow_name == "realignment":
                    key = f"realignment_{key}"
                vcfs[key] = vcf_files[0]
                
        return vcfs