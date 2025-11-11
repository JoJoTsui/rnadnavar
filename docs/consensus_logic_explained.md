# Consensus Variant Calling Logic

## Overview

The `run_consensus.R` script performs multi-caller consensus variant calling by finding overlapping variants across different variant callers and applying threshold-based filtering rules. It can work with both VCF and MAF formats and supports hierarchical consensus (using previous consensus results as input).

## Input Requirements

- Multiple variant files (VCF or MAF format) from different callers
- Files must follow naming convention: `sample.caller.maf` or `sample.caller.vcf`
- Caller names are extracted from filenames using regex: `sub(".*(\\.|_)(.*?)\\.maf", "\\2", filename)`

## Core Algorithm

### Step 1: Header Parsing

For each input file, the script extracts:
```bash
zgrep -E '##|#version' <file>  # Metadata
zgrep 'Hugo_Symbol' <file>     # Column headers (MAF)
zgrep '#CHROM' <file>          # Column headers (VCF)
```

This preserves original file metadata for output reconstruction.

### Step 2: Data Loading and Genomic Coordinate Conversion

For each caller's variants:

1. Load variant data (skip headers):
   ```bash
   zgrep -v '#' <file>
   ```

2. Create standardized genomic coordinates:
   - `#CHROM` ← `Chromosome`
   - `POS` ← `Start_Position`
   - `REF` ← `Reference_Allele`
   - `ALT` ← `Tumor_Seq_Allele2`
   - `DNAchange` ← `chr:g.posREF>ALT` (unique variant identifier)

3. Calculate range endpoints:
   ```r
   start <- POS
   end <- ifelse(nchar(REF) > nchar(ALT),
                 POS + nchar(REF) - 1,      # Deletion
                 ifelse(nchar(REF) < nchar(ALT),
                        POS + nchar(ALT) - 1, # Insertion
                        POS))                 # SNV or MNV
   ```

4. Convert to GenomicRanges objects for efficient overlap detection

### Step 3: Pairwise Overlap Detection

Performs all-vs-all pairwise comparisons:

```r
for (caller1 in callers) {
    for (caller2 in callers) {
        if (caller1 != caller2) {
            hits <- findOverlaps(query = mutsGR[[caller1]], 
                                subject = mutsGR[[caller2]], 
                                maxgap = 0)
            # Record overlapping variants
        }
    }
}
```

For each overlapping variant, records:
- `DNAchange`: Unique variant identifier
- `caller`: Which caller(s) detected it
- `FILTER`: Filter status from each caller

### Step 4: Consensus Rules Application

Aggregates caller support per variant:
```r
overlapping.variants.count <- str_count(caller, "|") + 1
```

Applies different rules for SNVs vs Indels:

#### SNVs (Single Nucleotide Variants)
- Pattern: `[0-9](A|C|G|T)>(A|C|G|T)$`
- Requirement: **Exact position and allele match**
- Threshold: ≥ N callers (default: 2)
- Rationale: SNV positions should be consistent across callers

#### Indels (Insertions/Deletions)
- Pattern: Everything else
- Requirement: **Any overlap** (maxgap=0)
- Threshold: **ALL indels with overlap are rescued**
- Rationale: Different normalization strategies can shift indel positions slightly

### Step 5: Filter Consensus Determination

For each variant, calculates consensus filter status:

```r
simplified.filter <- sapply(filters, function(x) {
    filt.val <- strsplit(x, "|", fixed=TRUE)[[1]]
    filt.t <- table(filt.val == "PASS")
    ifelse(prop.table(filt.t)["FALSE"] > 0.5, "FAIL", "PASS")
})
```

Logic:
- If ≥50% of supporting callers mark variant as "PASS" → `PASS`
- Otherwise → `FAIL`

### Step 6: Parallel Annotation

Using multiple CPU cores, annotates each variant with consensus information:

```r
what.caller.called <- function(row, consensus, variants) {
    if (variant in consensus) {
        callers <- paste(supporting_callers, collapse="|")
        filters <- paste(filter_values, collapse="|")
    }
    return(list(callers=callers, filters=filters))
}
```

Adds columns to each variant:
- `callers`: Pipe-separated list (e.g., "deepsomatic|mutect2|strelka")
- `filters`: Corresponding FILTER values (e.g., "PASS|PASS|weak_evidence")
- `FILTER_consensus`: Aggregated filter status
- `isconsensus`: Boolean indicating if variant meets consensus criteria

### Step 7: Output Generation

Creates multiple output files:

#### Per-Caller Files
Format: `{prefix}_{caller}.maf` or `{prefix}_{caller}.vcf`

Contains:
- All variants from that caller
- Added consensus annotation columns
- Original metadata preserved

#### Final Consensus File
Format: `{prefix}.maf` or `{prefix}.vcf`

Contains:
- Only variants where `isconsensus == TRUE`
- Deduplicated by `DNAchange`
- Metadata includes consensus information:
  ```
  ##INFO=<ID=callers,Number=1,Type=String,Description="Variant callers that called this mutation, separated by |">
  ##INFO=<ID=filters,Number=1,Type=String,Description="Filters provided by each variant caller, separated by |">
  ##INFO=<ID=consensus_filter,Number=1,Type=String,Description="PASS if 50% or more of the callers give the mutation, otherwise FAIL.">
  ```

#### Visualization PDF
Format: `{prefix}.pdf`

Contains three plots:
- **Plot A**: Bar chart showing consensus vs non-consensus variants per caller, faceted by PASS/FAIL
- **Plot B**: UpSet plot showing all variant overlaps across callers
- **Plot C**: UpSet plot showing PASS-only variant overlaps

## Hierarchical Consensus

The script can use previous consensus results as input, enabling multi-level consensus:

### Special Handling of Consensus Input
```r
all.consensus.muts <- all.muts[stringi::stri_detect_regex(all.muts$Caller, "consensus", case_insensitive=TRUE),]
all.consensus.muts <- all.consensus.muts[!duplicated(all.consensus.muts$DNAchange),]
```

When a "consensus" caller is detected:
- Variants are deduplicated
- Treated as a single caller in overlap detection
- Represents aggregated evidence from previous consensus round

## Real-World Example Comparison

### Example 1: DNA Callers + RNA Consensus → 29,588 variants

**Inputs:**
```
27,703  DNA_TUMOR_vs_DNA_NORMAL.deepsomatic.maf    (individual caller)
   762  DNA_TUMOR_vs_DNA_NORMAL.mutect2.maf       (individual caller)
15,557  DNA_TUMOR_vs_DNA_NORMAL.strelka.maf       (individual caller)
 6,950  RNA_TUMOR_vs_DNA_NORMAL.consensus.maf     (pre-validated consensus)
```

**Logic:**
- 3 DNA individual callers with heavy overlap
- Deepsomatic is very sensitive (27K variants)
- When 2+ DNA callers agree → enters consensus
- RNA consensus adds orthogonal validation
- **Result: 29,588 variants** (DNA caller overlaps dominate)

### Example 2: DNA Consensus + RNA Callers → 20,057 variants

**Inputs:**
```
26,141  DNA_TUMOR_vs_DNA_NORMAL.consensus.maf      (pre-validated consensus)
13,724  RNA_TUMOR_vs_DNA_NORMAL.deepsomatic.maf    (individual caller)
   340  RNA_TUMOR_vs_DNA_NORMAL.mutect2.maf       (individual caller)
 8,740  RNA_TUMOR_vs_DNA_NORMAL.strelka.maf       (individual caller)
```

**Logic:**
- 1 DNA consensus (already filtered)
- 3 RNA individual callers with limited overlap
- RNA has lower coverage (expression-dependent)
- RNA_mutect2 very sparse (340 variants)
- **Result: 20,057 variants** (RNA caller overlap is limiting factor)

### Why the Difference?

**Example 1 (More Variants):**
- Multiple raw DNA callers create large overlap set
- DNA has uniform genome coverage
- RNA consensus adds validated variants
- Formula: `DNA_overlap + RNA_validated = 29,588`

**Example 2 (Fewer Variants):**
- RNA callers have less overlap (expression-biased coverage)
- DNA consensus is already filtered
- RNA detection limited by transcript expression
- Formula: `RNA_overlap + DNA_validated = 20,057`

## Best Practices

### Choosing Input Strategy

**Recommended: Multiple Individual Callers from Primary Platform + Consensus from Secondary Platform**

For DNA-primary analysis:
```
✓ DNA_caller1.maf + DNA_caller2.maf + DNA_caller3.maf + RNA_consensus.maf
✗ DNA_consensus.maf + RNA_caller1.maf + RNA_caller2.maf + RNA_caller3.maf
```

Rationale:
- DNA-seq has uniform coverage across genome
- Multiple DNA callers provide redundancy and confidence
- RNA consensus adds orthogonal validation without expression bias
- Captures comprehensive somatic variant landscape

### Parameter Tuning

**Threshold (`--thr`):**
- Default: 2 callers
- Increase for higher specificity (fewer false positives)
- Decrease for higher sensitivity (more variants, more false positives)

**CPU (`--cpu`):**
- Use multiple cores for large datasets
- Parallelizes annotation step (most time-consuming)

## Key Takeaways

1. **SNVs require exact matches** - position and alleles must be identical
2. **Indels are rescued liberally** - any overlap counts due to normalization issues
3. **Filter consensus uses majority vote** - ≥50% PASS → consensus PASS
4. **Hierarchical consensus is supported** - previous consensus can be input
5. **Input order matters** - affects final variant count significantly
6. **DNA-primary strategy recommended** - better coverage and reliability than RNA
