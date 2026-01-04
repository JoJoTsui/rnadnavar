#!/usr/bin/env python3
"""
Workflow Abstraction Module

Provides a unified interface for both standard and realignment workflows.
This module enables seamless handling of different workflow types without code duplication.
"""

from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Dict, List, Optional


class WorkflowType(Enum):
    """Workflow types supported by the system."""
    STANDARD = "standard"
    REALIGNMENT = "realignment"


@dataclass
class WorkflowConfig:
    """Configuration for a workflow type."""
    name: str                           # "standard" or "realignment"
    base_path: Path                     # Base directory for this workflow
    stages: List[str]                   # Processing stages
    stage_paths: Dict[str, str]         # Stage name -> directory path mapping
    
    def get_stage_path(self, stage: str) -> Optional[Path]:
        """
        Get the full path for a specific stage.
        
        Args:
            stage: Stage name (e.g., "normalized", "consensus")
            
        Returns:
            Full path to the stage directory, or None if stage not found
        """
        if stage not in self.stage_paths:
            return None
        return self.base_path / self.stage_paths[stage]


class WorkflowManager:
    """Manages workflow detection and configuration."""
    
    # Standard workflow stages and paths
    STANDARD_STAGES = [
        "normalized",
        "consensus",
        "rescue",
        "cosmic_gnomad",
        "rna_editing",
        "filtered_rescue"
    ]
    
    STANDARD_STAGE_PATHS = {
        "normalized": "normalized",
        "consensus": "consensus",
        "rescue": "rescue",
        "cosmic_gnomad": "annotation/cosmic_gnomad",
        "rna_editing": "annotation/rna_editing",
        "filtered_rescue": "filtered"
    }
    
    # Realignment workflow stages and paths
    REALIGNMENT_STAGES = [
        "normalized",
        "consensus",
        "rescue",
        "cosmic_gnomad",
        "rna_editing",
        "filtered_rescue"
    ]
    
    REALIGNMENT_STAGE_PATHS = {
        "normalized": "normalized",
        "consensus": "consensus",
        "rescue": "rescue",
        "cosmic_gnomad": "annotation/cosmic_gnomad",
        "rna_editing": "annotation/rna_editing",
        "filtered_rescue": "filtered"
    }
    
    def __init__(self, base_dir: Path):
        """
        Initialize workflow manager with pipeline output directory.
        
        Args:
            base_dir: Base directory of the pipeline output
        """
        self.base_dir = Path(base_dir)
        self._detected_workflows = None
        
    def detect_workflows(self) -> List[WorkflowType]:
        """
        Detect which workflows are present in the output directory.
        
        Checks for the existence of expected directory structures:
        - Standard workflow: base_dir/normalized, base_dir/consensus, etc.
        - Realignment workflow: base_dir/vcf_realignment/normalized, etc.
        
        Returns:
            List of detected workflow types
        """
        if self._detected_workflows is not None:
            return self._detected_workflows
            
        detected = []
        
        # Check for standard workflow
        # Look for at least one standard stage directory
        standard_indicators = [
            self.base_dir / "normalized",
            self.base_dir / "consensus",
            self.base_dir / "rescue"
        ]
        if any(path.exists() for path in standard_indicators):
            detected.append(WorkflowType.STANDARD)
            
        # Check for realignment workflow
        # Look for vcf_realignment directory with at least one stage
        realignment_base = self.base_dir / "vcf_realignment"
        if realignment_base.exists():
            realignment_indicators = [
                realignment_base / "normalized",
                realignment_base / "consensus",
                realignment_base / "rescue"
            ]
            if any(path.exists() for path in realignment_indicators):
                detected.append(WorkflowType.REALIGNMENT)
        
        self._detected_workflows = detected
        return detected
        
    def get_workflow_config(self, workflow_type: WorkflowType) -> WorkflowConfig:
        """
        Get configuration for a specific workflow.
        
        Args:
            workflow_type: Type of workflow to configure
            
        Returns:
            WorkflowConfig object with paths and stages
            
        Raises:
            ValueError: If workflow_type is not recognized
        """
        if workflow_type == WorkflowType.STANDARD:
            return WorkflowConfig(
                name="standard",
                base_path=self.base_dir,
                stages=self.STANDARD_STAGES.copy(),
                stage_paths=self.STANDARD_STAGE_PATHS.copy()
            )
        elif workflow_type == WorkflowType.REALIGNMENT:
            return WorkflowConfig(
                name="realignment",
                base_path=self.base_dir / "vcf_realignment",
                stages=self.REALIGNMENT_STAGES.copy(),
                stage_paths=self.REALIGNMENT_STAGE_PATHS.copy()
            )
        else:
            raise ValueError(f"Unknown workflow type: {workflow_type}")
            
    def get_all_configs(self) -> Dict[WorkflowType, WorkflowConfig]:
        """
        Get configurations for all detected workflows.
        
        Returns:
            Dictionary mapping workflow types to their configurations
        """
        detected = self.detect_workflows()
        return {
            workflow_type: self.get_workflow_config(workflow_type)
            for workflow_type in detected
        }
        
    def has_workflow(self, workflow_type: WorkflowType) -> bool:
        """
        Check if a specific workflow is present.
        
        Args:
            workflow_type: Type of workflow to check
            
        Returns:
            True if workflow is detected, False otherwise
        """
        return workflow_type in self.detect_workflows()
        
    def get_available_stages(self, workflow_type: WorkflowType) -> List[str]:
        """
        Get list of stages that actually exist for a workflow.
        
        Args:
            workflow_type: Type of workflow to check
            
        Returns:
            List of stage names that have corresponding directories
        """
        config = self.get_workflow_config(workflow_type)
        available = []
        
        for stage in config.stages:
            stage_path = config.get_stage_path(stage)
            if stage_path and stage_path.exists():
                available.append(stage)
                
        return available


# Workflow constants for easy access
WORKFLOW_TYPES = [WorkflowType.STANDARD, WorkflowType.REALIGNMENT]

# Stage names (shared across workflows)
STAGE_NORMALIZED = "normalized"
STAGE_CONSENSUS = "consensus"
STAGE_RESCUE = "rescue"
STAGE_COSMIC_GNOMAD = "cosmic_gnomad"
STAGE_RNA_EDITING = "rna_editing"
STAGE_FILTERED_RESCUE = "filtered_rescue"

ALL_STAGES = [
    STAGE_NORMALIZED,
    STAGE_CONSENSUS,
    STAGE_RESCUE,
    STAGE_COSMIC_GNOMAD,
    STAGE_RNA_EDITING,
    STAGE_FILTERED_RESCUE
]
