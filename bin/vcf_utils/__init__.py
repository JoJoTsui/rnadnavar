"""
VCF utilities package for consensus and rescue workflows.

This package provides reusable functions for VCF variant aggregation,
tagging, statistics generation, I/O operations, and filtering.

Modules:
    aggregation: Core variant aggregation logic for combining variants from
        multiple callers and across modalities
    tagging: Variant tagging and annotation functions for caller support,
        modality tracking, and rescue marking
    statistics: Statistics generation for consensus and rescue operations
    io_utils: VCF I/O utilities for reading and writing VCF files
    filters: Filter normalization and categorization functions
    classification: Variant classification functions for individual callers
    classification_config: Centralized classification configuration and thresholds
    variant_classifier_unified: Unified classifier for all classification scenarios
    chromosome_utils: Canonical chromosome filtering utilities

Example:
    Basic consensus workflow:

    >>> from vcf_utils.aggregation import read_variants_from_vcf, aggregate_variants
    >>> from vcf_utils.statistics import compute_consensus_statistics, print_statistics
    >>> from vcf_utils.io_utils import write_union_vcf
    >>>
    >>> # Read variants from multiple callers
    >>> mutect2_vars = read_variants_from_vcf('mutect2.vcf.gz', 'mutect2')
    >>> strelka_vars = read_variants_from_vcf('strelka.vcf.gz', 'strelka')
    >>>
    >>> # Aggregate variants
    >>> collections = [('mutect2', mutect2_vars, None), ('strelka', strelka_vars, None)]
    >>> aggregated = aggregate_variants(collections, snv_threshold=2, indel_threshold=2)
    >>>
    >>> # Generate statistics
    >>> stats = compute_consensus_statistics(aggregated, 2, 2)
    >>> print_statistics(stats, 'consensus')
    >>>
    >>> # Write output
    >>> write_union_vcf(aggregated, template_header, 'sample', 'output.vcf.gz',
    ...                 'vcf.gz', ['mutect2', 'strelka'])

    Cross-modality rescue workflow:

    >>> from vcf_utils.tagging import mark_rescued_variants, tag_variant_with_modality
    >>>
    >>> # Read DNA and RNA consensus VCFs
    >>> dna_vars = read_variants_from_vcf('dna_consensus.vcf.gz', 'consensus', modality='DNA')
    >>> rna_vars = read_variants_from_vcf('rna_consensus.vcf.gz', 'consensus', modality='RNA')
    >>>
    >>> # Aggregate across modalities
    >>> collections = [('consensus', dna_vars, 'DNA'), ('consensus', rna_vars, 'RNA')]
    >>> aggregated = aggregate_variants(collections, snv_threshold=2, indel_threshold=2)
    >>>
    >>> # Mark rescued variants
    >>> dna_keys = set(dna_vars.keys())
    >>> rna_keys = set(rna_vars.keys())
    >>> aggregated = mark_rescued_variants(aggregated, dna_keys, rna_keys)
    >>>
    >>> # Tag with modality
    >>> modality_map = {'consensus': 'DNA'}  # Build full map
    >>> for vkey, data in aggregated.items():
    ...     tag_variant_with_modality(data, modality_map)

    Unified classification workflow:

    >>> from vcf_utils.variant_classifier_unified import UnifiedVariantClassifier
    >>> from vcf_utils.classification_config import DEFAULT_THRESHOLDS
    >>> from vcf_utils.chromosome_utils import is_canonical_chromosome
    >>>
    >>> # Create classifier with custom thresholds
    >>> classifier = UnifiedVariantClassifier({
    ...     'consensus_snv_threshold': 2,
    ...     'consensus_indel_threshold': 2
    ... })
    >>>
    >>> # Classify at different stages
    >>> standalone_class = classifier.classify_standalone_variant(variant, 'strelka', sample_indices)
    >>> consensus_class = classifier.classify_consensus_variant(variant_data)
    >>> rescue_class = classifier.classify_rescue_variant(variant_data, modality_map)
    >>>
    >>> # Filter to canonical chromosomes
    >>> if is_canonical_chromosome(variant.CHROM):
    ...     # Process variant
    ...     pass
"""

__version__ = "1.1.0"

# Import key classes and functions for convenience
from .chromosome_utils import (
    get_canonical_chromosome_list,
    is_canonical_chromosome,
    should_skip_variant,
    standardize_chromosome_name,
)
from .classification_config import (
    BIOLOGICAL_CATEGORIES,
    CANONICAL_SET,
    CATEGORIES,
    DEFAULT_THRESHOLDS,
    FILTER_NON_CANONICAL_DEFAULT,
    get_config_with_overrides,
    validate_classification,
)
from .variant_classifier_unified import UnifiedVariantClassifier, default_classifier
