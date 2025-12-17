#!/usr/bin/env python3
"""
Multi-Modal Variant Classification Module

This module provides sophisticated variant classification logic using multi-modal
evidence from population databases, caller consensus, and cancer mutation catalogs.
It implements classification rules for germline, somatic, and artifact variants.

FEATURES:
- Multi-modal evidence evaluation (DNA + RNA callers)
- Population frequency-based germline classification
- Caller consensus-based somatic classification
- Cross-modality rescue classification
- COSMIC recurrence-based evidence weighting
- Comprehensive statistics generation

Requirements Satisfied: 3.1, 3.2, 3.3, 3.4, 3.5

Author: COSMIC/gnomAD Enhancement Pipeline
Date: 2025-12-17
"""

import logging
from typing import Dict, Any, Optional, List, Tuple
from dataclasses import dataclass

logger = logging.getLogger(__name__)


@dataclass
class ClassificationEvidence:
    """Container for variant classification evidence."""
    population_frequency: Optional[float] = None
    cosmic_recurrence: Optional[int] = None
    dna_caller_support: int = 0
    rna_caller_support: int = 0
    total_dna_callers: int = 0
    total_rna_callers: int = 0
    has_cross_modality: bool = False
    original_filter: str = "PASS"


@dataclass
class ClassificationResult:
    """Container for variant classification result."""
    classification: str
    confidence: str
    evidence_summary: str
    rescue_flag: bool = False
    statistics: Dict[str, Any] = None


class VariantClassifier:
    """
    Multi-modal variant classifier using population and cancer database evidence.
    
    This class implements sophisticated variant classification logic that considers:
    - Population frequency from gnomAD for germline classification
    - Caller consensus across DNA and RNA modalities
    - COSMIC recurrence for cancer-specific evidence
    - Cross-modality support for rescue classification
    """
    
    def __init__(self, 
                 germline_frequency_threshold: float = 0.001,
                 somatic_consensus_threshold: int = 3,
                 cosmic_recurrence_threshold: int = 5,
                 cross_modality_min_support: int = 1):
        """
        Initialize variant classifier with configurable thresholds.
        
        Args:
            germline_frequency_threshold: Population frequency threshold for germline classification
            somatic_consensus_threshold: Minimum caller support for high-confidence somatic
            cosmic_recurrence_threshold: Minimum COSMIC recurrence for rescue classification
            cross_modality_min_support: Minimum support required in each modality for cross-modality rescue
        """
        self.germline_freq_threshold = germline_frequency_threshold
        self.somatic_consensus_threshold = somatic_consensus_threshold
        self.cosmic_recurrence_threshold = cosmic_recurrence_threshold
        self.cross_modality_min_support = cross_modality_min_support
        
        # Statistics tracking
        self.stats = {
            'total_classified': 0,
            'germline_count': 0,
            'somatic_count': 0,
            'somatic_rescue_count': 0,
            'artifact_count': 0,
            'unchanged_count': 0,
            'classification_changes': {},
            'evidence_distribution': {
                'population_frequency': 0,
                'cosmic_recurrence': 0,
                'dna_consensus': 0,
                'cross_modality': 0
            }
        }
        
        logger.info("=== Multi-Modal Variant Classifier Initialized ===")
        logger.info(f"Germline frequency threshold: {self.germline_freq_threshold}")
        logger.info(f"Somatic consensus threshold: {self.somatic_consensus_threshold}")
        logger.info(f"COSMIC recurrence threshold: {self.cosmic_recurrence_threshold}")
        logger.info(f"Cross-modality minimum support: {self.cross_modality_min_support}")
    
    def extract_evidence(self, variant_info: Dict[str, Any]) -> ClassificationEvidence:
        """
        Extract classification evidence from variant INFO fields.
        
        Args:
            variant_info: Dictionary containing variant INFO field values
            
        Returns:
            ClassificationEvidence object with extracted evidence
        """
        evidence = ClassificationEvidence()
        
        # Extract population frequency from gnomAD
        # Try multiple gnomAD frequency fields in order of preference
        freq_fields = ['GNOMAD_FAF95', 'GNOMAD_AF', 'AF', 'faf95']
        for field in freq_fields:
            if field in variant_info and variant_info[field] is not None:
                try:
                    freq_value = variant_info[field]
                    if isinstance(freq_value, str):
                        # Handle comma-separated values (take maximum)
                        freq_values = [float(x) for x in freq_value.split(',') if x.strip()]
                        evidence.population_frequency = max(freq_values) if freq_values else None
                    else:
                        evidence.population_frequency = float(freq_value)
                    
                    if evidence.population_frequency is not None:
                        logger.debug(f"Found population frequency: {evidence.population_frequency} (field: {field})")
                        break
                except (ValueError, TypeError) as e:
                    logger.debug(f"Could not parse frequency field {field}: {e}")
                    continue
        
        # Extract COSMIC recurrence
        cosmic_fields = ['COSMIC_CNT', 'GENOME_SCREEN_SAMPLE_COUNT', 'COSMIC_ID']
        for field in cosmic_fields:
            if field in variant_info and variant_info[field] is not None:
                try:
                    cosmic_value = variant_info[field]
                    if isinstance(cosmic_value, str):
                        # Handle comma-separated values (take maximum)
                        cosmic_values = [int(x) for x in cosmic_value.split(',') if x.strip().isdigit()]
                        evidence.cosmic_recurrence = max(cosmic_values) if cosmic_values else None
                    else:
                        evidence.cosmic_recurrence = int(cosmic_value)
                    
                    if evidence.cosmic_recurrence is not None:
                        logger.debug(f"Found COSMIC recurrence: {evidence.cosmic_recurrence} (field: {field})")
                        break
                except (ValueError, TypeError) as e:
                    logger.debug(f"Could not parse COSMIC field {field}: {e}")
                    continue
        
        # Extract caller support information
        if 'N_DNA_CALLERS_SUPPORT' in variant_info:
            try:
                evidence.dna_caller_support = int(variant_info['N_DNA_CALLERS_SUPPORT'])
            except (ValueError, TypeError):
                evidence.dna_caller_support = 0
        
        if 'N_RNA_CALLERS_SUPPORT' in variant_info:
            try:
                evidence.rna_caller_support = int(variant_info['N_RNA_CALLERS_SUPPORT'])
            except (ValueError, TypeError):
                evidence.rna_caller_support = 0
        
        if 'N_DNA_CALLERS' in variant_info:
            try:
                evidence.total_dna_callers = int(variant_info['N_DNA_CALLERS'])
            except (ValueError, TypeError):
                evidence.total_dna_callers = 0
        
        if 'N_RNA_CALLERS' in variant_info:
            try:
                evidence.total_rna_callers = int(variant_info['N_RNA_CALLERS'])
            except (ValueError, TypeError):
                evidence.total_rna_callers = 0
        
        # Check for cross-modality support
        evidence.has_cross_modality = (
            evidence.dna_caller_support >= self.cross_modality_min_support and
            evidence.rna_caller_support >= self.cross_modality_min_support
        )
        
        # Extract original filter
        if 'FILTER' in variant_info:
            evidence.original_filter = variant_info['FILTER']
        elif 'original_filter' in variant_info:
            evidence.original_filter = variant_info['original_filter']
        
        return evidence
    
    def classify_variant(self, evidence: ClassificationEvidence) -> ClassificationResult:
        """
        Classify variant based on multi-modal evidence.
        
        Args:
            evidence: ClassificationEvidence object with variant evidence
            
        Returns:
            ClassificationResult with classification and supporting information
        """
        self.stats['total_classified'] += 1
        
        # Rule 1: Population frequency > threshold → Germline (highest priority)
        if (evidence.population_frequency is not None and 
            evidence.population_frequency > self.germline_freq_threshold):
            
            self.stats['germline_count'] += 1
            self.stats['evidence_distribution']['population_frequency'] += 1
            
            result = ClassificationResult(
                classification="Germline",
                confidence="High",
                evidence_summary=f"Population frequency {evidence.population_frequency:.4f} > {self.germline_freq_threshold}",
                rescue_flag=False
            )
            
            logger.debug(f"Classified as Germline: {result.evidence_summary}")
            return result
        
        # Rule 2: Unanimous DNA caller support → High-confidence Somatic
        if (evidence.total_dna_callers > 0 and 
            evidence.dna_caller_support >= self.somatic_consensus_threshold and
            evidence.dna_caller_support == evidence.total_dna_callers):
            
            self.stats['somatic_count'] += 1
            self.stats['evidence_distribution']['dna_consensus'] += 1
            
            result = ClassificationResult(
                classification="Somatic",
                confidence="High",
                evidence_summary=f"Unanimous DNA caller support ({evidence.dna_caller_support}/{evidence.total_dna_callers})",
                rescue_flag=False
            )
            
            logger.debug(f"Classified as Somatic: {result.evidence_summary}")
            return result
        
        # Rule 3: Cross-modality support + COSMIC recurrence → Somatic Rescue
        if (evidence.has_cross_modality and 
            evidence.cosmic_recurrence is not None and
            evidence.cosmic_recurrence >= self.cosmic_recurrence_threshold):
            
            self.stats['somatic_rescue_count'] += 1
            self.stats['evidence_distribution']['cross_modality'] += 1
            self.stats['evidence_distribution']['cosmic_recurrence'] += 1
            
            result = ClassificationResult(
                classification="Somatic",
                confidence="Medium",
                evidence_summary=f"Cross-modality support (DNA:{evidence.dna_caller_support}, RNA:{evidence.rna_caller_support}) + COSMIC recurrence {evidence.cosmic_recurrence}",
                rescue_flag=True
            )
            
            logger.debug(f"Classified as Somatic Rescue: {result.evidence_summary}")
            return result
        
        # Rule 4: No classification criteria met → Preserve original or assign based on quality
        self.stats['unchanged_count'] += 1
        
        # Determine classification based on original filter and evidence quality
        if evidence.original_filter in ["PASS", "Somatic"]:
            # If originally passed but doesn't meet new criteria, check evidence quality
            if (evidence.dna_caller_support > 0 or evidence.rna_caller_support > 0):
                classification = "Somatic"  # Some evidence, keep as somatic
                confidence = "Low"
            else:
                classification = "Artifact"  # No caller support
                confidence = "Medium"
        elif evidence.original_filter in ["Germline"]:
            classification = "Germline"
            confidence = "Low"
        elif evidence.original_filter in ["Artifact", "LowQuality", "LowDepth"]:
            classification = "Artifact"
            confidence = "Medium"
        else:
            classification = "Artifact"  # Default for unknown filters
            confidence = "Low"
        
        result = ClassificationResult(
            classification=classification,
            confidence=confidence,
            evidence_summary=f"No strong evidence, based on original filter: {evidence.original_filter}",
            rescue_flag=False
        )
        
        logger.debug(f"Preserved/assigned classification: {result.classification} ({result.evidence_summary})")
        return result
    
    def classify_variant_from_info(self, variant_info: Dict[str, Any]) -> ClassificationResult:
        """
        Classify variant directly from INFO field dictionary.
        
        Args:
            variant_info: Dictionary containing variant INFO field values
            
        Returns:
            ClassificationResult with classification and supporting information
        """
        evidence = self.extract_evidence(variant_info)
        result = self.classify_variant(evidence)
        
        # Track classification changes
        original_filter = evidence.original_filter
        new_classification = result.classification
        
        change_key = f"{original_filter} -> {new_classification}"
        if change_key not in self.stats['classification_changes']:
            self.stats['classification_changes'][change_key] = 0
        self.stats['classification_changes'][change_key] += 1
        
        return result
    
    def get_classification_statistics(self) -> Dict[str, Any]:
        """
        Get comprehensive classification statistics.
        
        Returns:
            Dictionary containing classification statistics and metrics
        """
        total = self.stats['total_classified']
        
        statistics = {
            'total_variants_classified': total,
            'classification_counts': {
                'Germline': self.stats['germline_count'],
                'Somatic': self.stats['somatic_count'],
                'Somatic_Rescue': self.stats['somatic_rescue_count'],
                'Artifact': self.stats['artifact_count'],
                'Unchanged': self.stats['unchanged_count']
            },
            'classification_rates': {},
            'evidence_usage': self.stats['evidence_distribution'].copy(),
            'classification_changes': self.stats['classification_changes'].copy(),
            'thresholds_used': {
                'germline_frequency_threshold': self.germline_freq_threshold,
                'somatic_consensus_threshold': self.somatic_consensus_threshold,
                'cosmic_recurrence_threshold': self.cosmic_recurrence_threshold,
                'cross_modality_min_support': self.cross_modality_min_support
            }
        }
        
        # Calculate classification rates
        if total > 0:
            for classification, count in statistics['classification_counts'].items():
                statistics['classification_rates'][classification] = (count / total) * 100
        
        return statistics
    
    def log_classification_summary(self) -> None:
        """Log comprehensive classification summary statistics."""
        stats = self.get_classification_statistics()
        total = stats['total_variants_classified']
        
        logger.info("=== Multi-Modal Variant Classification Summary ===")
        logger.info(f"Total variants classified: {total:,}")
        
        if total > 0:
            logger.info("Classification distribution:")
            for classification, count in stats['classification_counts'].items():
                rate = stats['classification_rates'][classification]
                logger.info(f"  {classification}: {count:,} ({rate:.1f}%)")
            
            logger.info("Evidence usage:")
            for evidence_type, count in stats['evidence_usage'].items():
                rate = (count / total) * 100 if total > 0 else 0
                logger.info(f"  {evidence_type}: {count:,} ({rate:.1f}%)")
            
            if stats['classification_changes']:
                logger.info("Classification changes:")
                for change, count in sorted(stats['classification_changes'].items(), 
                                          key=lambda x: x[1], reverse=True):
                    rate = (count / total) * 100
                    logger.info(f"  {change}: {count:,} ({rate:.1f}%)")
            
            logger.info("Thresholds used:")
            for threshold, value in stats['thresholds_used'].items():
                logger.info(f"  {threshold}: {value}")
        
        logger.info("=== Classification Summary Complete ===")


def create_variant_classifier(config: Optional[Dict[str, Any]] = None) -> VariantClassifier:
    """
    Create a variant classifier with optional configuration.
    
    Args:
        config: Optional configuration dictionary with threshold overrides
        
    Returns:
        Configured VariantClassifier instance
    """
    if config is None:
        config = {}
    
    # Default configuration with overrides
    classifier_config = {
        'germline_frequency_threshold': config.get('germline_frequency_threshold', 0.001),
        'somatic_consensus_threshold': config.get('somatic_consensus_threshold', 3),
        'cosmic_recurrence_threshold': config.get('cosmic_recurrence_threshold', 5),
        'cross_modality_min_support': config.get('cross_modality_min_support', 1)
    }
    
    return VariantClassifier(**classifier_config)


def classify_variants_batch(variants: List[Dict[str, Any]], 
                          classifier: Optional[VariantClassifier] = None) -> List[ClassificationResult]:
    """
    Classify a batch of variants using multi-modal evidence.
    
    Args:
        variants: List of variant dictionaries with INFO field data
        classifier: Optional pre-configured classifier (creates default if None)
        
    Returns:
        List of ClassificationResult objects
    """
    if classifier is None:
        classifier = create_variant_classifier()
    
    results = []
    for variant in variants:
        result = classifier.classify_variant_from_info(variant)
        results.append(result)
    
    return results


def main():
    """Main entry point for standalone usage and testing."""
    import argparse
    import sys
    import json
    
    parser = argparse.ArgumentParser(
        description="Multi-modal variant classifier",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Multi-Modal Variant Classification Engine:
  This module provides sophisticated variant classification logic using multi-modal
  evidence from population databases, caller consensus, and cancer mutation catalogs.

Key Features:
  - Multi-modal evidence evaluation (DNA + RNA callers)
  - Population frequency-based germline classification
  - Caller consensus-based somatic classification
  - Cross-modality rescue classification
  - COSMIC recurrence-based evidence weighting

Examples:
  # Test classification with sample data
  python variant_classifier.py --test
  
  # Classify variants from JSON file
  python variant_classifier.py -i variants.json -o results.json
        """
    )
    
    parser.add_argument(
        '-i', '--input',
        metavar='JSON_FILE',
        help='Input JSON file with variant data'
    )
    
    parser.add_argument(
        '-o', '--output',
        metavar='JSON_FILE',
        help='Output JSON file for classification results'
    )
    
    parser.add_argument(
        '--test',
        action='store_true',
        help='Run test classification with sample data'
    )
    
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Enable verbose logging'
    )
    
    args = parser.parse_args()
    
    # Set up logging
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    
    try:
        if args.test:
            # Run test with sample data
            logger.info("Running variant classification test...")
            
            # Sample test variants
            test_variants = [
                {
                    'GNOMAD_AF': 0.01,  # High frequency -> Germline
                    'N_DNA_CALLERS_SUPPORT': 2,
                    'N_DNA_CALLERS': 3,
                    'FILTER': 'PASS'
                },
                {
                    'GNOMAD_AF': 0.0001,  # Low frequency
                    'N_DNA_CALLERS_SUPPORT': 3,  # Unanimous -> Somatic
                    'N_DNA_CALLERS': 3,
                    'FILTER': 'PASS'
                },
                {
                    'GNOMAD_AF': 0.0001,  # Low frequency
                    'N_DNA_CALLERS_SUPPORT': 1,
                    'N_RNA_CALLERS_SUPPORT': 2,  # Cross-modality
                    'COSMIC_CNT': 10,  # High COSMIC -> Rescue
                    'FILTER': 'Artifact'
                }
            ]
            
            classifier = create_variant_classifier()
            results = classify_variants_batch(test_variants, classifier)
            
            logger.info("Test results:")
            for i, result in enumerate(results):
                logger.info(f"  Variant {i+1}: {result.classification} ({result.confidence}) - {result.evidence_summary}")
                if result.rescue_flag:
                    logger.info(f"    RESCUE FLAG SET")
            
            classifier.log_classification_summary()
            
        elif args.input and args.output:
            # Process input file
            logger.info(f"Processing variants from {args.input}...")
            
            with open(args.input, 'r') as f:
                variants = json.load(f)
            
            classifier = create_variant_classifier()
            results = classify_variants_batch(variants, classifier)
            
            # Convert results to serializable format
            output_data = []
            for result in results:
                output_data.append({
                    'classification': result.classification,
                    'confidence': result.confidence,
                    'evidence_summary': result.evidence_summary,
                    'rescue_flag': result.rescue_flag
                })
            
            with open(args.output, 'w') as f:
                json.dump(output_data, f, indent=2)
            
            logger.info(f"Results written to {args.output}")
            classifier.log_classification_summary()
            
        else:
            parser.print_help()
            sys.exit(1)
        
        logger.info("✓ Variant classification completed successfully!")
        sys.exit(0)
        
    except KeyboardInterrupt:
        logger.error("Variant classification interrupted by user")
        sys.exit(130)
    except Exception as e:
        logger.error(f"Variant classification failed: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()