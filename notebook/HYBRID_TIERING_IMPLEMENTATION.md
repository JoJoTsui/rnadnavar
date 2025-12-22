# Hybrid Tiering System Implementation Summary

## Overview

Successfully implemented a two-dimensional hybrid tiering system for variant quality assessment that combines:
1. **Caller Support (C1-C7)**: Category-aware DNA/RNA caller consensus
2. **Database Evidence (D0-D1)**: Presence in external variant databases

Final tiers use intuitive **CxDy notation** (e.g., C1D1, C2D0, C7D1) for clarity and maintainability.

## Implementation Date

December 21, 2025

## Files Created

### Core Modules (bin/common/)

1. **[tier_config.py](../bin/common/tier_config.py)** (418 lines)
   - Defines C1-C7 caller tier rules with lambda conditions
   - Defines D0-D1 database tier rules with thresholds
   - Generates all 14 CxDy tier combinations
   - Provides tier metadata (colors, quality scores, display names)
   - Includes backward compatibility mapping (T1-T8 → CxDy)

2. **[category_matcher.py](../bin/common/category_matcher.py)** (333 lines)
   - Implements category-aware caller counting
   - Parses `FILTER_NORMALIZED_<CALLER>_<MODALITY>` fields
   - Only counts callers voting for same category as final FILTER
   - Supports cyvcf2 VCF record extraction
   - Includes validation and debugging utilities

3. **[database_checker.py](../bin/common/database_checker.py)** (358 lines)
   - Detects variant presence in external databases
   - Checks gnomAD (AF > 0.001), COSMIC (CNT > 0), REDIportal, DARNED
   - Applies significance thresholds for D1 classification
   - Provides database details extraction
   - Includes validation utilities

### Tiering Engine (notebook/vcf_stats/)

4. **[tiering_engine.py](tiering_engine.py)** (418 lines)
   - Main `TieringEngine` class combining caller + database evidence
   - Methods: `compute_caller_tier()`, `compute_database_tier()`, `compute_final_tier()`
   - Batch processing utilities for VCF files
   - Supports both parsed data and cyvcf2 VCF records
   - Returns comprehensive tier information dictionaries

## Files Modified

### Integration Updates

5. **[tiering.py](tiering.py)** (2 functions updated)
   - Updated `tier_rule()` to use new CxDy system
   - Added `tier_rule_legacy()` for backward compatibility (T1-T8)
   - Integrated TieringEngine for consistent tier assignment

6. **[visualizer.py](visualizer.py)** (`plot_tier_distribution()`)
   - Updated to display all 14 CxDy tiers
   - Uses TIER_ORDER from tier_config for proper ordering
   - Enhanced title with CxDy format explanation
   - Rotated x-axis labels for readability

7. **[rescue_analyzer.py](rescue_analyzer.py)** (`analyze_rescue_vcf()`)
   - Added count labels to stacked bar plots
   - Labels display inside bars with white text
   - Improved hover templates with category counts
   - Better visibility of variant distribution

### Documentation

8. **[vcf_statistics_P2374372.ipynb](vcf_statistics_P2374372.ipynb)** (Cell 35 - markdown)
   - Complete rewrite with hybrid system documentation
   - Detailed explanation of C1-C7 caller tiers
   - Database tier D0-D1 definitions
   - Full 14-tier combination matrix table
   - Category-specific tier expectations
   - Technical implementation details

## Tier Structure

### Caller Support Tiers (C1-C7)

| Tier | DNA Callers | RNA Callers | Description |
|------|-------------|-------------|-------------|
| C1 | ≥2 | ≥2 | Both modalities strong (highest confidence) |
| C2 | ≥2 | 0-1 | DNA-strong, RNA weak/absent |
| C3 | 0-1 | ≥2 | RNA-strong, DNA weak/absent |
| C4 | 1 | 1 | Both modalities weak |
| C5 | 1 | 0 | DNA-only weak |
| C6 | 0 | 1 | RNA-only weak |
| C7 | 0 | 0 | No caller support (RNA_Edit, NoConsensus) |

### Database Evidence Tiers (D0-D1)

| Tier | Criteria |
|------|----------|
| D1 | Present in gnomAD (AF>0.001), COSMIC (CNT>0), REDIportal, or DARNED |
| D0 | Not present in any database (novel variant) |

### Final Tier Combinations (14 Total)

```
C1D1, C1D0  (Quality: 140, 130)
C2D1, C2D0  (Quality: 120, 110)
C3D1, C3D0  (Quality: 100, 90)
C4D1, C4D0  (Quality: 80, 70)
C5D1, C5D0  (Quality: 60, 50)
C6D1, C6D0  (Quality: 40, 30)
C7D1, C7D0  (Quality: 20, 10)
```

## Key Features

### Category-Aware Caller Counting

**Problem Solved**: Original tiering counted all callers regardless of their vote, which could assign high tiers to variants where callers disagreed on classification.

**Solution**: Only count callers whose `FILTER_NORMALIZED` field matches the final `FILTER` category.

**Example**:
```python
# Variant with final FILTER="Somatic"
FILTER_NORMALIZED_Strelka_DNA_TUMOR = "Somatic"   # ✓ Counts
FILTER_NORMALIZED_Mutect2_DNA_TUMOR = "Somatic"   # ✓ Counts
FILTER_NORMALIZED_Strelka_RNA_TUMOR = "Artifact"  # ✗ Excluded

Result: 2 DNA callers, 0 RNA callers → C2D0
```

### Database Evidence Integration

**Problem Solved**: Original system didn't consider external validation from population/mutation databases.

**Solution**: Add database dimension (D0/D1) providing orthogonal quality metric.

**Impact**:
- Germline variants with gnomAD support: Higher confidence
- Somatic variants in COSMIC: Recurrent mutations, known drivers
- RNA_Edit in REDIportal: C7D1 (no callers but database-validated)
- Novel variants: D0 tier requires additional validation

### Intuitive CxDy Notation

**Problem Solved**: T1-T12 numbering obscured the two-dimensional nature of tiering.

**Solution**: CxDy format explicitly shows caller tier (C1-C7) and database tier (D0-D1).

**Benefits**:
- Instant understanding: "C2D0" = DNA-strong, no database
- Easier filtering: Select all D1 tiers for database-validated variants
- Better maintenance: No need to memorize numeric mappings

## Testing Results

All core modules tested successfully:

```bash
✓ tier_config.py: 7 caller tiers, 2 database tiers, 14 final tiers
✓ category_matcher.py: Category-aware counting (3 DNA, 1 RNA) - PASS
✓ database_checker.py: gnomAD/COSMIC/RNA_Edit detection - PASS
✓ tiering_engine.py: Hybrid tier computation (C1D1, C2D0, C7D1) - PASS
```

## Usage Examples

### From Python Script

```python
from vcf_stats.tiering_engine import TieringEngine

engine = TieringEngine()

# Method 1: From parsed data
tier_info = engine.compute_tier(
    final_filter="Somatic",
    filter_normalized_fields={
        "FILTER_NORMALIZED_Strelka_DNA_TUMOR": "Somatic",
        "FILTER_NORMALIZED_Mutect2_DNA_TUMOR": "Somatic",
    },
    info_dict={"gnomAD_AF": 0.0, "COSMIC_CNT": 15}
)
print(tier_info["final_tier"])  # "C2D1"

# Method 2: From VCF record
import cyvcf2
vcf = cyvcf2.VCF("rescue.vcf.gz")
for variant in vcf:
    tier_info = engine.compute_tier_from_vcf_record(variant)
    print(f"{variant.CHROM}:{variant.POS} - {tier_info['final_tier']}")

# Method 3: Simple tier computation
tier = engine.compute_tier_simple(
    dna_caller_count=2,
    rna_caller_count=1,
    has_database_support=False
)
print(tier)  # "C2D0"
```

### From Jupyter Notebook

```python
# Already integrated - just run existing cells
visualizer.plot_tier_distribution()  # Now shows C1D1-C7D0 tiers
```

## Category-Specific Tier Expectations

| Category | Typical Tiers | Rationale |
|----------|---------------|-----------|
| Somatic | C1-C6 (any D) | Caller-supported mutations |
| Germline | C1D1, C2D1 | High caller + gnomAD validation |
| RNA_Edit | C7D1 | No callers (DNA/RNA differ), REDIportal DB |
| NoConsensus | C7D0, C7D1 | DNA/RNA disagree, no consensus |
| Artifact | C5D0, C6D0, C7D0 | Weak support, no database |
| Reference | C1D1, C2D1 | Strong consensus, gnomAD match |

## Integration Points

### For vcf_stats (Analysis)

- ✓ Tiering engine imported in tiering.py
- ✓ Visualizations updated for CxDy display
- ✓ Rescue analyzer shows category counts
- ✓ Notebook documentation complete

### For vcf_utils (Pipeline) - Future

Modules in `bin/common/` are ready for integration:
- Import `tier_config.py` for tier definitions
- Use `category_matcher.py` for caller counting during consensus
- Use `database_checker.py` for annotation checks
- Use `tiering_engine.py` for final tier assignment

## Backward Compatibility

### Legacy T1-T8 Mapping

`tier_rule_legacy()` function maintains old behavior:

```python
T1 → C1D0  # 2+ DNA + 2+ RNA
T2 → C2D0  # 2+ DNA + 1 RNA
T3 → C2D0  # 2+ DNA only
T4 → C4D0  # 1 DNA + 1+ RNA
T5 → C5D0  # 1 DNA only
T6 → C3D0  # 2+ RNA only
T7 → C6D0  # 1 RNA only
```

**Note**: Legacy mapping assumes D0 (no database) since original system didn't consider databases.

## Performance Considerations

- **Memory**: Tier config loaded once, reused across all variants
- **Speed**: Lambda conditions in C1-C7 rules enable fast tier assignment
- **Scalability**: Batch processing utilities handle large VCF files efficiently

## Future Enhancements

1. **Category-aware counting in pipeline**: Integrate category_matcher.py into vcf_utils consensus step
2. **Quality score filtering**: Use TIER_QUALITY_SCORES for automated variant prioritization
3. **Interactive tier selection**: Add UI controls to filter visualizations by tier
4. **Tier-based reports**: Generate separate reports for high-confidence (C1-C3) vs low-confidence (C5-C7) tiers
5. **Machine learning features**: Use tier as feature for variant effect prediction models

## Validation Checklist

- ✓ All 7 caller tiers correctly defined with mutually exclusive conditions
- ✓ Database thresholds match biological significance (gnomAD AF > 0.1%)
- ✓ Category-aware counting filters non-concordant callers
- ✓ CxDy notation generates all 14 expected combinations
- ✓ Tier quality scores monotonically decrease (C1D1=140 → C7D0=10)
- ✓ Colors provide visual quality gradient (blue → red)
- ✓ Visualizations display all tiers with proper ordering
- ✓ Rescue analyzer shows category counts in stacked bars
- ✓ Notebook documentation explains system comprehensively
- ✓ All test cases pass successfully

## References

- **Tier Config**: `/t9k/mnt/hdd/work/Vax/pipeline/rnadnavar/bin/common/tier_config.py`
- **Category Matcher**: `/t9k/mnt/hdd/work/Vax/pipeline/rnadnavar/bin/common/category_matcher.py`
- **Database Checker**: `/t9k/mnt/hdd/work/Vax/pipeline/rnadnavar/bin/common/database_checker.py`
- **Tiering Engine**: `/t9k/mnt/hdd/work/Vax/pipeline/rnadnavar/notebook/vcf_stats/tiering_engine.py`
- **Notebook Documentation**: Cell 35 in `vcf_statistics_P2374372.ipynb`

## Support

For questions or issues:
1. Check module docstrings for detailed API documentation
2. Run modules with `python3 <module>.py` to see configuration summaries
3. Use debugging utilities: `print_tier_summary()`, `print_caller_vote_summary()`, `print_database_summary()`

---

**Implementation Status**: ✅ Complete

All planned features implemented, tested, and documented.
