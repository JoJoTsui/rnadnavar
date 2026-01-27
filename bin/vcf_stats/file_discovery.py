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
from typing import Any, Dict, Optional

# Import constants from main module
from . import TOOLS, VCF_STAGE_ORDER
from .workflow import WorkflowConfig, WorkflowManager


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

    def __init__(
        self, base_dir: Path, workflow_manager: Optional[WorkflowManager] = None
    ):
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
            "normalized": {},  # Normalized VCFs from standalone callers
            "consensus": {},  # Consensus VCFs (DNA & RNA)
            "rescue": {},  # Rescue VCFs (DNA + RNA combined)
            "cosmic_gnomad": {},  # Cosmic/gnomAD annotated
            "rna_editing": {},  # RNA editing annotated
            "filtered_rescue": {},  # Final filtered rescue VCFs
        }
        self.bam_files = {}

    def discover_vcfs(self) -> Dict[str, Dict[str, Dict[str, Any]]]:
        """
        Discover all VCF files organized by processing stage.

        Only discovers normalized, consensus, rescue, and annotated (cosmic_gnomad,
        rna_editing, filtered_rescue) VCFs. Raw variant_calling VCFs are skipped.

        Returns:
            Dictionary mapping stage -> {name: {path, stage, tool, sample, file_id}}
            Each VCF entry contains metadata:
                - path: Path object to the VCF file
                - stage: Processing stage (e.g., "normalized", "consensus")
                - tool: Tool name (e.g., "deepsomatic", "consensus")
                - sample: Sample identifier (e.g., "DT_vs_DN", "RT_vs_DN")
                - file_id: Full file identifier key
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
        """
        Return DT/DN/RT suffix if present.
        Handles both sample IDs ending with suffixes (e.g., '2374372RT')
        and modality names (e.g., 'DNA_TUMOR', 'RNA_TUMOR', 'DNA_NORMAL').
        """
        if not sample:
            return None

        # Check modality names first
        if sample == "DNA_TUMOR":
            return "DT"
        elif sample == "RNA_TUMOR":
            return "RT"
        elif sample == "DNA_NORMAL":
            return "DN"

        # Otherwise check last 2 characters for suffix
        suf = sample[-2:]
        return suf if suf in {"DT", "DN", "RT"} else None

    @classmethod
    def _pair_to_suffix_key(cls, pair: str) -> Optional[str]:
        """
        Convert a sample pair name like '2374372RT_vs_2374372DN' or
        'RNA_TUMOR_realign_vs_DNA_NORMAL' to 'RT_vs_DN'.
        """
        if "_vs_" not in pair:
            return None
        t, n = pair.split("_vs_", 1)
        # Strip _realign suffix if present (for realignment workflow samples)
        t = t.replace("_realign", "")
        n = n.replace("_realign", "")
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

    @staticmethod
    def _extract_sample_pair(dir_name: str) -> Optional[str]:
        """
        Extract sample pair identifier using the *T*_vs_*N* pattern.

        This method handles various tumor vs normal naming conventions used in
        DNA-only mode and other workflows. It uses regex pattern matching to
        identify valid sample pairs.

        Supported naming formats:
            - Full sample names: "WES_LL_T_1_vs_WES_LL_N_1"
            - DT/DN suffix patterns: "2374372DT_vs_2374372DN"
            - RT/DN suffix patterns: "2374372RT_vs_2374372DN"
            - Mixed patterns: "Sample_T_vs_Sample_N"
            - Underscore variations: "SAMPLE_T_1_vs_SAMPLE_N_1"

        Pattern matching rules:
            1. Must contain "_vs_" separator
            2. Tumor part (before _vs_) must contain "T" indicator:
               - Ends with "DT" (DNA Tumor)
               - Ends with "RT" (RNA Tumor)
               - Contains "_T_" or "_T" followed by optional suffix
            3. Normal part (after _vs_) must contain "N" indicator:
               - Ends with "DN" (DNA Normal)
               - Contains "_N_" or "_N" followed by optional suffix

        Args:
            dir_name: Directory name representing a sample pair
                     (e.g., "WES_LL_T_1_vs_WES_LL_N_1")

        Returns:
            The sample pair identifier if pattern matches, None otherwise.
            Returns the original dir_name if it matches the pattern.

        Example:
            >>> VCFFileDiscovery._extract_sample_pair("WES_LL_T_1_vs_WES_LL_N_1")
            'WES_LL_T_1_vs_WES_LL_N_1'
            >>> VCFFileDiscovery._extract_sample_pair("2374372DT_vs_2374372DN")
            '2374372DT_vs_2374372DN'
            >>> VCFFileDiscovery._extract_sample_pair("invalid_name")
            None

        Validates: Requirements 2.3
        """
        import re

        if not dir_name or "_vs_" not in dir_name:
            return None

        # Split into tumor and normal parts
        parts = dir_name.split("_vs_", 1)
        if len(parts) != 2:
            return None

        tumor_part, normal_part = parts

        # Define regex patterns for tumor identification
        # Pattern 1: Ends with DT or RT (e.g., "2374372DT", "2374372RT")
        # Pattern 2: Contains _T_ or ends with _T followed by optional suffix
        #            (e.g., "WES_LL_T_1", "Sample_T", "Sample_T_1")
        tumor_patterns = [
            r".*[DR]T$",           # Ends with DT or RT
            r".*_T_.*",            # Contains _T_ in the middle
            r".*_T$",              # Ends with _T
            r".*_T\d+$",           # Ends with _T followed by digits (e.g., _T1)
        ]

        # Define regex patterns for normal identification
        # Pattern 1: Ends with DN (e.g., "2374372DN")
        # Pattern 2: Contains _N_ or ends with _N followed by optional suffix
        #            (e.g., "WES_LL_N_1", "Sample_N", "Sample_N_1")
        normal_patterns = [
            r".*DN$",              # Ends with DN
            r".*_N_.*",            # Contains _N_ in the middle
            r".*_N$",              # Ends with _N
            r".*_N\d+$",           # Ends with _N followed by digits (e.g., _N1)
        ]

        # Check if tumor part matches any tumor pattern
        tumor_match = any(
            re.match(pattern, tumor_part, re.IGNORECASE)
            for pattern in tumor_patterns
        )

        # Check if normal part matches any normal pattern
        normal_match = any(
            re.match(pattern, normal_part, re.IGNORECASE)
            for pattern in normal_patterns
        )

        # Return the sample pair if both parts match their respective patterns
        if tumor_match and normal_match:
            return dir_name

        return None

    def _discover_normalized_vcfs(self):
        """Discover normalized VCFs from standalone callers."""
        norm_dir = self.base_dir / "normalized"
        if norm_dir.exists():
            self.vcf_files["normalized"] = self._discover_stage_normalized(
                norm_dir, "standard"
            )

    def _discover_consensus_vcfs(self):
        """Discover consensus VCFs (DNA and RNA)."""
        consensus_dir = self.base_dir / "consensus"
        if consensus_dir.exists():
            self.vcf_files["consensus"] = self._discover_stage_consensus(
                consensus_dir, "standard"
            )

    def _discover_rescue_vcfs(self):
        """Discover rescue VCFs (DNA + RNA combined)."""
        rescue_dir = self.base_dir / "rescue"
        if rescue_dir.exists():
            self.vcf_files["rescue"] = self._discover_stage_rescue(
                rescue_dir, "standard"
            )

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
            self.vcf_files["filtered_rescue"] = self._discover_stage_filtered(
                filtered_dir, "standard"
            )

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
            print("\nALIGNMENT FILES:")
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

    def find_matching_files(
        self, tool: str, modality: str
    ) -> Dict[str, Optional[Path]]:
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
        self, workflow_config: WorkflowConfig
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

            # For annotation stages, don't check if stage_path exists since files are in rescue/
            is_annotation_stage = stage in ["cosmic_gnomad", "rna_editing"]

            if not stage_path or (not is_annotation_stage and not stage_path.exists()):
                workflow_vcfs[stage] = {}
                continue

            # Use stage-specific discovery logic
            if stage == "normalized":
                workflow_vcfs[stage] = self._discover_stage_normalized(
                    stage_path, workflow_config.name
                )
            elif stage == "consensus":
                workflow_vcfs[stage] = self._discover_stage_consensus(
                    stage_path, workflow_config.name
                )
            elif stage == "rescue":
                workflow_vcfs[stage] = self._discover_stage_rescue(
                    stage_path, workflow_config.name
                )
            elif stage == "cosmic_gnomad":
                workflow_vcfs[stage] = self._discover_stage_annotation(
                    stage_path,
                    "*.rescue.cosmic_gnomad_annotated.final.vcf.gz",
                    workflow_config.name,
                )
            elif stage == "rna_editing":
                workflow_vcfs[stage] = self._discover_stage_annotation(
                    stage_path, "*.rescue.rna_annotated.vcf.gz", workflow_config.name
                )
            elif stage == "filtered_rescue":
                workflow_vcfs[stage] = self._discover_stage_filtered(
                    stage_path, workflow_config.name
                )
            else:
                workflow_vcfs[stage] = {}

        return workflow_vcfs

    def discover_all_workflows(self) -> Dict[str, Dict[str, Dict[str, Path]]]:
        """
        Discover VCF files for all detected workflows.

        Supports three workflow types:
        - standard: Full RNA+DNA pipeline with normalized, consensus, rescue, annotation stages
        - realignment: RNA realignment workflow with same stages as standard
        - dna_only: DNA-only pipeline with variant_calling and consensus stages only

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
                },
                "dna_only": {
                    "variant_calling": {name: {path, stage, tool, sample, file_id}},
                    "consensus": {name: {path, stage, tool, sample, file_id}},
                }
            }

        Validates: Requirements 2.4
        """
        from .workflow import WorkflowType

        all_workflows = {}

        # Get all detected workflow configurations
        workflow_configs = self.workflow_manager.get_all_configs()

        for workflow_type, config in workflow_configs.items():
            workflow_name = config.name

            # Use DNA-only specific discovery for DNA_ONLY workflow
            if workflow_type == WorkflowType.DNA_ONLY:
                all_workflows[workflow_name] = self.discover_vcfs_for_dna_only(config)
            else:
                # Use standard discovery for STANDARD and REALIGNMENT workflows
                all_workflows[workflow_name] = self.discover_vcfs_for_workflow(config)

        return all_workflows

    def discover_bam_files_for_workflow(
        self, workflow_config: WorkflowConfig
    ) -> Dict[str, Path]:
        """
        Discover BAM/CRAM files for a specific workflow.

        Checks multiple locations:
        - preprocessing/recalibrated/ (standard workflow, *.recal.cram)
        - preprocessing/markduplicates/ (realignment workflow, *.md.cram)

        Args:
            workflow_config: Configuration for the workflow

        Returns:
            Dictionary mapping sample names to BAM/CRAM file paths
        """
        bam_files = {}

        # Try recalibrated directory first (standard workflow)
        recal_dir = workflow_config.base_path / "preprocessing" / "recalibrated"
        markdup_dir = workflow_config.base_path / "preprocessing" / "markduplicates"

        # Check which directory exists and use appropriate patterns
        search_dirs = []
        if recal_dir.exists():
            search_dirs.append((recal_dir, ["*.recal.cram", "*.cram", "*.bam"]))
        if markdup_dir.exists():
            search_dirs.append((markdup_dir, ["*.md.cram", "*.cram", "*.bam"]))

        if not search_dirs:
            return bam_files

        for base_dir, patterns in search_dirs:
            for sample_dir in base_dir.iterdir():
                if not sample_dir.is_dir():
                    continue

                sample_name = sample_dir.name

                # Skip if we already found this sample (prioritize recalibrated)
                if sample_name in bam_files:
                    continue

                # Add workflow prefix to sample name for realignment
                full_sample_name = sample_name
                if workflow_config.name == "realignment":
                    full_sample_name = f"realignment_{sample_name}"

                # Try each pattern in order
                for pattern in patterns:
                    found_files = list(sample_dir.glob(pattern))
                    if found_files:
                        bam_files[full_sample_name] = found_files[0]
                        break

        return bam_files

    def _discover_stage_normalized(
        self, stage_path: Path, workflow_name: str
    ) -> Dict[str, Path]:
        """
        Discover normalized VCFs in a stage directory.

        Args:
            stage_path: Path to the normalized stage directory
            workflow_name: Name of the workflow ("standard" or "realignment")

        Returns:
            Dictionary of discovered VCF files (deduplicated)
        """
        vcfs = {}
        seen_files = set()  # Track files we've already added

        for tool in TOOLS:
            tool_dir = stage_path / tool
            if not tool_dir.exists():
                continue

            # Scan subdirectories with _vs_ pattern (covers both modality and sample-pair structures)
            for subdir in tool_dir.iterdir():
                if not subdir.is_dir() or "_vs_" not in subdir.name:
                    continue

                vcf_files = list(subdir.glob("*.norm.vcf.gz"))
                if not vcf_files:
                    continue

                vcf_path = vcf_files[0]
                # Skip if we've already added this file
                if str(vcf_path) in seen_files:
                    continue

                seen_files.add(str(vcf_path))
                sample_pair = subdir.name

                # Extract suffix key
                pair_key = self._pair_to_suffix_key(sample_pair)
                if pair_key:
                    key = f"{tool}_{pair_key}"
                else:
                    # Fallback: try to infer from sample name
                    key = f"{tool}_{sample_pair}"

                # Add workflow prefix for realignment
                if workflow_name == "realignment":
                    key = f"realignment_{key}"

                # Store with metadata
                vcfs[key] = {
                    "path": vcf_path,
                    "stage": "normalized",
                    "tool": tool.lower(),
                    "sample": pair_key or sample_pair,
                    "file_id": key,
                }

        return vcfs

    def _discover_stage_consensus(
        self, stage_path: Path, workflow_name: str
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
        seen_files = set()  # Track discovered files to avoid duplicates

        if not stage_path.exists():
            return vcfs

        # Consensus files are in subdirectories with _vs_ pattern
        for subdir in stage_path.iterdir():
            if not subdir.is_dir() or "_vs_" not in subdir.name:
                continue

            vcf_files = list(subdir.glob("*.consensus.vcf.gz"))
            if not vcf_files:
                continue

            vcf_path = vcf_files[0]
            # Skip if we've already added this file
            if str(vcf_path) in seen_files:
                continue

            seen_files.add(str(vcf_path))
            sample_pair = subdir.name

            # Extract suffix key
            pair_key = self._pair_to_suffix_key(sample_pair)
            if pair_key:
                key = pair_key
            else:
                key = sample_pair

            # Add workflow prefix for realignment
            if workflow_name == "realignment":
                key = f"realignment_{key}"

            # Store with metadata
            vcfs[key] = {
                "path": vcf_path,
                "stage": "consensus",
                "tool": "consensus",
                "sample": pair_key or sample_pair,
                "file_id": key,
            }

        return vcfs

    def _discover_stage_rescue(
        self, stage_path: Path, workflow_name: str
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
                    rescued_key = self._rescued_pair_to_suffix_key(subdir.name)
                    key = rescued_key or subdir.name
                    # Add workflow prefix for realignment
                    if workflow_name == "realignment":
                        key = f"realignment_{key}"
                    # Store with metadata
                    vcfs[key] = {
                        "path": vcf_files[0],
                        "stage": "rescue",
                        "tool": "rescue",
                        "sample": rescued_key or subdir.name,
                        "file_id": key,
                    }

        return vcfs

    def _discover_stage_annotation(
        self, stage_path: Path, pattern: str, workflow_name: str
    ) -> Dict[str, Path]:
        """
        Discover annotation VCFs in a stage directory.

        For this pipeline, annotation files are actually in rescue subdirectories,
        not in a separate annotation directory.

        Args:
            stage_path: Path to the annotation stage directory (unused, we look in rescue)
            pattern: Glob pattern for VCF files
            workflow_name: Name of the workflow ("standard" or "realignment")

        Returns:
            Dictionary of discovered VCF files
        """
        vcfs = {}

        # Annotation files are in rescue subdirectories
        rescue_base = (
            stage_path.parent.parent
            if "annotation" in str(stage_path)
            else stage_path.parent
        )
        rescue_dir = rescue_base / "rescue"

        if not rescue_dir.exists():
            return vcfs

        for sample_pair_dir in rescue_dir.iterdir():
            if sample_pair_dir.is_dir() and "rescued" in sample_pair_dir.name:
                vcf_files = list(sample_pair_dir.glob(pattern))
                if vcf_files:
                    rescued_key = self._rescued_pair_to_suffix_key(sample_pair_dir.name)
                    key = rescued_key or sample_pair_dir.name
                    # Add workflow prefix for realignment
                    if workflow_name == "realignment":
                        key = f"realignment_{key}"

                    # Infer stage from pattern
                    stage_name = "annotation"
                    if "cosmic_gnomad" in pattern:
                        stage_name = "cosmic_gnomad"
                    elif "rna_annotated" in pattern:
                        stage_name = "rna_editing"

                    # Store with metadata
                    vcfs[key] = {
                        "path": vcf_files[0],
                        "stage": stage_name,
                        "tool": stage_name,
                        "sample": rescued_key or sample_pair_dir.name,
                        "file_id": key,
                    }

        return vcfs

    def _discover_stage_filtered(
        self, stage_path: Path, workflow_name: str
    ) -> Dict[str, Path]:
        """
        Discover filtered VCFs in a stage directory.

        Filtered files can be in both:
        - filtered/ directory subdirectories
        - rescue/ directory subdirectories

        Args:
            stage_path: Path to the filtered stage directory
            workflow_name: Name of the workflow ("standard" or "realignment")

        Returns:
            Dictionary of discovered VCF files
        """
        vcfs = {}

        # Check both filtered/ and rescue/ directories
        search_dirs = []
        if stage_path.exists():
            search_dirs.append(stage_path)

        # Also check rescue directory
        rescue_dir = stage_path.parent / "rescue"
        if rescue_dir.exists():
            search_dirs.append(rescue_dir)

        for search_dir in search_dirs:
            for sample_pair_dir in search_dir.iterdir():
                if not sample_pair_dir.is_dir():
                    continue

                vcf_files = list(sample_pair_dir.glob("*.filtered.vcf.gz"))
                if not vcf_files:
                    continue

                # Only include combined rescued filtered result, omit single-modality filtered
                # i.e., keep keys that contain '_rescued_'
                if "_rescued_" not in sample_pair_dir.name:
                    continue

                rescued_key = self._rescued_pair_to_suffix_key(sample_pair_dir.name)
                key = rescued_key or sample_pair_dir.name
                # Add workflow prefix for realignment
                if workflow_name == "realignment":
                    key = f"realignment_{key}"

                # Store with metadata
                vcfs[key] = {
                    "path": vcf_files[0],
                    "stage": "filtered_rescue",
                    "tool": "filtered_rescue",
                    "sample": rescued_key or sample_pair_dir.name,
                    "file_id": key,
                }

        return vcfs

    def _discover_variant_calling_vcfs(self, base_path: Path) -> Dict[str, Dict[str, Any]]:
        """
        Discover VCFs from variant_calling stage for DNA-only workflow.

        Scans variant_calling/deepsomatic/, variant_calling/mutect2/, and
        variant_calling/strelka/ subdirectories for VCF files.

        Expected directory structure:
            variant_calling/{caller}/{sample_pair}/*.vcf.gz

        Sample pair pattern: *T*_vs_*N* (tumor vs normal naming convention)

        File selection priority (per caller):
            - mutect2: Prefer *.filtered.vcf.gz over raw *.vcf.gz
            - strelka: Prefer *.variants.vcf.gz over other patterns
            - deepsomatic: Use *.vcf.gz (single output)

        Args:
            base_path: Base path of the pipeline output directory

        Returns:
            Dictionary mapping file_id to metadata dict containing:
                - path: Path object to the VCF file
                - stage: "variant_calling"
                - tool: Tool name (e.g., "deepsomatic", "mutect2", "strelka")
                - sample: Sample pair identifier (e.g., "WES_LL_T_1_vs_WES_LL_N_1")
                - file_id: Full file identifier key (e.g., "deepsomatic_WES_LL_T_1_vs_WES_LL_N_1")

        Example:
            >>> discovery = VCFFileDiscovery("/path/to/output")
            >>> vcfs = discovery._discover_variant_calling_vcfs(Path("/path/to/output"))
            >>> for file_id, metadata in vcfs.items():
            ...     print(f"{file_id}: {metadata['tool']} - {metadata['sample']}")
            deepsomatic_WES_LL_T_1_vs_WES_LL_N_1: deepsomatic - WES_LL_T_1_vs_WES_LL_N_1
            mutect2_WES_LL_T_1_vs_WES_LL_N_1: mutect2 - WES_LL_T_1_vs_WES_LL_N_1
            strelka_WES_LL_T_1_vs_WES_LL_N_1: strelka - WES_LL_T_1_vs_WES_LL_N_1

        Validates: Requirements 2.1
        """
        vcfs = {}
        seen_files = set()  # Track files we've already added to avoid duplicates

        variant_calling_dir = base_path / "variant_calling"
        if not variant_calling_dir.exists():
            return vcfs

        # Scan each caller subdirectory (deepsomatic, mutect2, strelka)
        for tool in TOOLS:
            tool_dir = variant_calling_dir / tool
            if not tool_dir.exists():
                continue

            # Scan sample pair subdirectories with _vs_ pattern
            for subdir in tool_dir.iterdir():
                if not subdir.is_dir() or "_vs_" not in subdir.name:
                    continue

                # Look for VCF files with tool-specific priority
                vcf_path = self._select_vcf_for_tool(subdir, tool)
                if vcf_path is None:
                    continue

                # Skip if we've already added this file
                if str(vcf_path) in seen_files:
                    continue

                seen_files.add(str(vcf_path))
                sample_pair = subdir.name

                # Extract suffix key using existing pattern (e.g., "DT_vs_DN")
                # For DNA-only, we may have different naming conventions like "WES_LL_T_1_vs_WES_LL_N_1"
                pair_key = self._pair_to_suffix_key(sample_pair)

                # Use the sample pair name directly if suffix extraction fails
                # This handles DNA-only naming conventions like "*T*_vs_*N*"
                sample_identifier = pair_key or sample_pair

                # Create file identifier key
                key = f"{tool}_{sample_identifier}"

                # Store with metadata
                vcfs[key] = {
                    "path": vcf_path,
                    "stage": "variant_calling",
                    "tool": tool.lower(),
                    "sample": sample_identifier,
                    "file_id": key,
                }

        return vcfs

    @staticmethod
    def _select_vcf_for_tool(subdir: Path, tool: str) -> Optional[Path]:
        """
        Select the appropriate VCF file for a given tool.

        Different variant callers produce different output files:
        - mutect2: Produces both raw (*.mutect2.vcf.gz) and filtered (*.mutect2.filtered.vcf.gz)
                   We prefer the filtered version as it contains quality-filtered variants.
        - strelka: Produces *.variants.vcf.gz for SNVs/indels
        - deepsomatic: Produces *.deepsomatic.vcf.gz

        Args:
            subdir: Directory containing VCF files for a sample pair
            tool: Tool name (e.g., "mutect2", "strelka", "deepsomatic")

        Returns:
            Path to the selected VCF file, or None if no suitable file found
        """
        # Tool-specific file selection patterns (in priority order)
        tool_patterns = {
            "mutect2": [
                "*.mutect2.filtered.vcf.gz",  # Prefer filtered output
                "*.filtered.vcf.gz",           # Alternative filtered pattern
            ],
            "strelka": [
                "*.strelka.variants.vcf.gz",  # SNV/indel variants
                "*.variants.vcf.gz",           # Alternative pattern
            ],
            "deepsomatic": [
                "*.deepsomatic.vcf.gz",       # DeepSomatic output
            ],
        }

        # Get patterns for this tool, fallback to generic if not specified
        patterns = tool_patterns.get(tool.lower(), [])

        # Try each pattern in priority order
        for pattern in patterns:
            vcf_files = list(subdir.glob(pattern))
            if vcf_files:
                return vcf_files[0]

        # Fallback: use any VCF file, but exclude index files (.tbi)
        # and prefer non-raw files when multiple exist
        all_vcfs = [f for f in subdir.glob("*.vcf.gz") if not f.name.endswith(".tbi")]
        if not all_vcfs:
            return None

        # If multiple VCFs exist, prefer filtered/processed over raw
        # Raw mutect2 files typically don't have "filtered" in the name
        for vcf in all_vcfs:
            if "filtered" in vcf.name or "variants" in vcf.name:
                return vcf

        # Last resort: return first VCF found
        return all_vcfs[0]

    def discover_vcfs_for_dna_only(
        self, workflow_config: WorkflowConfig
    ) -> Dict[str, Dict[str, Any]]:
        """
        Discover VCFs for DNA-only workflow.

        DNA-only mode runs variant calling with DNA-tumor and DNA-normal samples only,
        producing a simpler two-stage pipeline:
        - variant_calling/{caller}/{sample_pair}/*.vcf.gz
        - consensus/{sample_pair}/*.consensus.vcf.gz

        This method explicitly does NOT discover from normalized/, rescue/, filtered/,
        or annotation/ directories, as these are not present in DNA-only mode.

        Sample pair pattern: *T*_vs_*N* (tumor vs normal naming convention)

        Args:
            workflow_config: Configuration for the DNA_ONLY workflow type.
                            Must have stages=["variant_calling", "consensus"]

        Returns:
            Dictionary organized by stage, with VCF metadata:
            {
                "variant_calling": {
                    file_id: {
                        "path": Path,
                        "stage": "variant_calling",
                        "tool": str,
                        "sample": str,
                        "file_id": str
                    },
                    ...
                },
                "consensus": {
                    file_id: {
                        "path": Path,
                        "stage": "consensus",
                        "tool": "consensus",
                        "sample": str,
                        "file_id": str
                    },
                    ...
                }
            }

        Example:
            >>> from vcf_stats.workflow import WorkflowManager, WorkflowType
            >>> manager = WorkflowManager("/path/to/dna_only_output")
            >>> discovery = VCFFileDiscovery("/path/to/dna_only_output", workflow_manager=manager)
            >>> config = manager.get_workflow_config(WorkflowType.DNA_ONLY)
            >>> vcfs = discovery.discover_vcfs_for_dna_only(config)
            >>> print(f"Found {len(vcfs['variant_calling'])} variant_calling VCFs")
            >>> print(f"Found {len(vcfs['consensus'])} consensus VCFs")

        Validates: Requirements 2.2, 2.5
        """
        workflow_vcfs = {}

        # Only discover from the two DNA-only stages: variant_calling and consensus
        # Explicitly skip normalized/, rescue/, filtered/, annotation/ directories
        for stage in workflow_config.stages:
            stage_path = workflow_config.get_stage_path(stage)

            if not stage_path or not stage_path.exists():
                workflow_vcfs[stage] = {}
                continue

            if stage == "variant_calling":
                # Use the existing _discover_variant_calling_vcfs method
                workflow_vcfs[stage] = self._discover_variant_calling_vcfs(
                    workflow_config.base_path
                )
            elif stage == "consensus":
                # Discover consensus VCFs using the DNA-only pattern
                workflow_vcfs[stage] = self._discover_consensus_vcfs_dna_only(
                    stage_path
                )
            else:
                # Skip any other stages (should not happen for DNA_ONLY workflow)
                workflow_vcfs[stage] = {}

        return workflow_vcfs

    def _discover_consensus_vcfs_dna_only(
        self, stage_path: Path
    ) -> Dict[str, Dict[str, Any]]:
        """
        Discover consensus VCFs for DNA-only workflow.

        Expected directory structure:
            consensus/{sample_pair}/*.consensus.vcf.gz

        Sample pair pattern: *T*_vs_*N* (tumor vs normal naming convention)

        Args:
            stage_path: Path to the consensus stage directory

        Returns:
            Dictionary mapping file_id to metadata dict containing:
                - path: Path object to the VCF file
                - stage: "consensus"
                - tool: "consensus"
                - sample: Sample pair identifier (e.g., "WES_LL_T_1_vs_WES_LL_N_1")
                - file_id: Full file identifier key

        Example:
            >>> discovery = VCFFileDiscovery("/path/to/output")
            >>> vcfs = discovery._discover_consensus_vcfs_dna_only(Path("/path/to/output/consensus"))
            >>> for file_id, metadata in vcfs.items():
            ...     print(f"{file_id}: {metadata['sample']}")
            WES_LL_T_1_vs_WES_LL_N_1: WES_LL_T_1_vs_WES_LL_N_1

        Validates: Requirements 2.2
        """
        vcfs = {}
        seen_files = set()  # Track discovered files to avoid duplicates

        if not stage_path.exists():
            return vcfs

        # Consensus files are in subdirectories with _vs_ pattern
        for subdir in stage_path.iterdir():
            if not subdir.is_dir() or "_vs_" not in subdir.name:
                continue

            # Validate sample pair pattern using *T*_vs_*N* convention
            sample_pair = self._extract_sample_pair(subdir.name)
            if sample_pair is None:
                # If pattern doesn't match, use directory name as-is
                sample_pair = subdir.name

            # Look for consensus VCF files (*.consensus.vcf.gz)
            vcf_files = list(subdir.glob("*.consensus.vcf.gz"))
            if not vcf_files:
                continue

            vcf_path = vcf_files[0]

            # Skip if we've already added this file
            if str(vcf_path) in seen_files:
                continue

            seen_files.add(str(vcf_path))

            # Extract suffix key (e.g., "DT_vs_DN") if applicable
            pair_key = self._pair_to_suffix_key(sample_pair)

            # Use the sample pair name directly if suffix extraction fails
            # This handles DNA-only naming conventions like "*T*_vs_*N*"
            sample_identifier = pair_key or sample_pair

            # Create file identifier key
            key = sample_identifier

            # Store with metadata
            vcfs[key] = {
                "path": vcf_path,
                "stage": "consensus",
                "tool": "consensus",
                "sample": sample_identifier,
                "file_id": key,
            }

        return vcfs
