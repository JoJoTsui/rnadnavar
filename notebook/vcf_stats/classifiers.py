#!/usr/bin/env python3
"""
Variant Classifier Module

Unified classification system based on VCF FILTER field values.
Classifies variants into biological categories: Somatic, Germline, Reference, 
Artifact, RNA_Edit, and NoConsensus.

This module uses a simplified, unified approach:
- All classification is based solely on the FILTER field
- No complex rule-based logic or caller-specific classification
- Consistent behavior across all VCF types and processing stages
"""


def classify_by_filter(variant) -> str:
    """
    Classify variant based solely on FILTER field.
    
    This is the core classification function that maps FILTER values
    to biological categories. It does not use genotype, depth, or other
    fields - only the FILTER field as set by the pipeline.
    
    Args:
        variant: cyvcf2 Variant object
        
    Returns:
        Category name: 'Somatic', 'Germline', 'Reference', 'Artifact',
                      'RNA_Edit', or 'NoConsensus'
    """
    filter_val = variant.FILTER
    
    # Handle PASS or missing filter (None, ".")
    if filter_val in (None, ".", "PASS"):
        return "Somatic"
    
    # Convert to lowercase string for case-insensitive matching
    filter_str = str(filter_val).lower()
    
    # Check for specific filter values
    if "germline" in filter_str:
        return "Germline"
    if "reference" in filter_str:
        return "Reference"
    if "rna_edit" in filter_str or "rnaedit" in filter_str:
        return "RNA_Edit"
    if "noconsensus" in filter_str:
        return "NoConsensus"
    
    # Default to Artifact for any other filter value
    return "Artifact"


def classify_annotated_variant(variant) -> str:
    """
    Classify variants using unified FILTER-based classification.
    
    Uses the classify_by_filter() function for all classification.
    This provides consistent behavior across all VCF types and stages.
    
    Returns one of: Somatic, Germline, Reference, Artifact, RNA_Edit, NoConsensus
    
    Args:
        variant: cyvcf2 Variant object
        
    Returns:
        Category name based on FILTER value
    """
    return classify_by_filter(variant)


def get_sample_indices(vcf_obj, caller_name):
    """
    Determine tumor and normal sample indices from VCF samples.

    Args:
        vcf_obj: cyvcf2.VCF object
        caller_name (str): Name of the variant caller

    Returns:
        dict: {'tumor': int, 'normal': int} or None
    """
    samples = vcf_obj.samples

    if len(samples) < 2:
        return None

    indices = {}

    # Try to identify tumor and normal samples
    for i, sample in enumerate(samples):
        sample_lower = sample.lower()
        if "tumor" in sample_lower or "tumour" in sample_lower:
            indices["tumor"] = i
        elif "normal" in sample_lower:
            indices["normal"] = i

    # Fallback: assume first is tumor, second is normal
    if "tumor" not in indices and len(samples) >= 1:
        indices["tumor"] = 0
    if "normal" not in indices and len(samples) >= 2:
        indices["normal"] = 0
        indices["tumor"] = 1

    return indices


print("✓ Variant classification functions defined")