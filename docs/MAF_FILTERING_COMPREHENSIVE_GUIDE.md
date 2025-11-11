# MAF Filtering Comprehensive Guide

## Overview

The `filter_mutations.py` script applies multi-layered quality filters to somatic variants in MAF format, producing a filtered MAF with a `RaVeX_FILTER` column that documents which variants pass quality control or why they fail.

## Pipeline Integration

### Nextflow Workflow
```groovy
// STEP 7: FILTERING
MAF_FILTERING(maf_to_filter, fasta, input_sample, realignment)
filtered_maf = MAF_FILTERING.out.maf
```

### Module Execution
```bash
filter_mutations.py \
    -i ${input}.maf \
    --output ${prefix}.filtered.maf \
    --ref ${genome}.fa \
    ${additional_args}
```

## Complete Filtering Workflow

### Step 1: Input Loading and Preparation

**Read MAF file(s):**
```python
maf = pd.read_csv(input_file, sep="\t", comment="#", low_memory=False)
```

**Create unique variant identifier:**
```python
DNAchange = Chromosome:g.Start_Position + Reference_Allele + ">" + Tumor_Seq_Allele2
# Example: chr7:g.140453136A>T
```

**Expected input columns from consensus:**
- Basic variant info: `Hugo_Symbol`, `Chromosome`, `Start_Position`, `Reference_Allele`, `Tumor_Seq_Allele2`
- Read counts: `t_depth`, `t_ref_count`, `t_alt_count`
- Annotations: `Consequence`, `BIOTYPE`, `SYMBOL`, `MAX_AF`, `SOMATIC`
- Consensus info: `Caller`, `callers`, `filters`, `FILTER_consensus`, `isconsensus`

### Step 2: Load Filter Resources

**Whitelist (optional):**
```
BED format: CHROM  POS  END  REF  ALT
Purpose: Force PASS for specific variants (overrides all other filters)
```

**Blacklist (optional):**
```
BED format: CHROM  START  END  [REASON]
Purpose: Exclude variants in problematic regions
```

### Step 3: Apply Basic Filters

#### gnomAD Population Frequency Filter
```python
maf["ingnomAD"] = maf["MAX_AF"] >= gnomad_thr  # default: 0.0001 (0.01%)
```

**Purpose:** Remove common germline variants
**Logic:** Variants with population frequency ≥0.01% are likely germline polymorphisms

**Examples:**
```
chr7:g.140453136A>T  MAX_AF=0.00001 (0.001%)  → ingnomAD=False ✓ (rare, somatic)
chr7:g.12345678C>T   MAX_AF=0.15 (15%)        → ingnomAD=True  ✗ (common germline)
```

#### Whitelist Check
```python
maf["whitelist"] = maf["DNAchange"].isin(whitelist_variants)
```

**Purpose:** Force PASS for known pathogenic variants
**Priority:** Absolute - overrides all other filters

#### Blacklist Check
```python
maf = remove_muts_in_range(df=maf, blacklist=blacklist_regions)
```

**Purpose:** Exclude variants in problematic genomic regions
**Adds columns:** `blacklist` (True/False), `blk_reason` (explanation)

### Step 4: Annotate Variant Context

#### 4.1 Noncoding Classification

**Noncoding consequence types:**
```python
noncoding_list = [
    "intron_variant",
    "intergenic_variant", 
    "non_coding_transcript_variant",
    "non_coding_transcript_exon_variant",
    "mature_miRNA_variant",
    "regulatory_region_variant",
    "IGR", "INTRON", "RNA"
]
```

**Logic:**
```python
# Extract first consequence from VEP annotation
first_consequence = Consequence.split("&")[0].split(",")[0]
maf["noncoding"] = first_consequence in noncoding_list
```

**Examples:**
```
Consequence: "missense_variant&splice_region_variant"
  → First: "missense_variant" → noncoding=False ✓

Consequence: "intron_variant"
  → First: "intron_variant" → noncoding=True
```

#### 4.2 Immunoglobulin and Pseudogene Filter

**Flagged biotypes:**
```python
flagged_patterns = [
    "IG_C_gene", "IG_D_gene", "IG_J_gene", "IG_V_gene",
    "TR_C_gene", "TR_J_gene", "TR_V_gene", 
    "pseudogene"
]
maf["ig_pseudo"] = BIOTYPE.contains(flagged_patterns)
```

**Purpose:** Flag variants in immunoglobulin genes and pseudogenes (often artifacts)

**Examples:**
```
BRAF, BIOTYPE="protein_coding"  → ig_pseudo=False ✓
IGHV3-23, BIOTYPE="IG_V_gene"   → ig_pseudo=True
```

#### 4.3 Homopolymer Context Detection

**Algorithm:**
1. Extract genomic context (±10bp around variant position)
2. Substitute reference allele with alternate allele
3. Scan for 6+ consecutive identical bases

**Detailed process:**
```python
# Open reference genome
genome = pysam.FastaFile(reference_fasta)

# Extract context
context = genome.fetch(chrom, pos-11, pos+10)
# Example: "ATCGATCGATCAGTCGATCG"
#                     ^
#                   position

# For insertion: check if inserted base creates homopolymer
# For deletion: skip check (return False)
# For SNV: substitute and check

# Substitute ref with alt
context_with_alt = context.replace(ref, alt)

# Scan for homopolymer (6+ consecutive identical bases)
for window in sliding_windows(context_with_alt, size=6):
    if len(set(window)) == 1:  # All bases identical
        homopolymer = True
        break
```

**Examples:**
```
Normal context:
  "ATCGATCGATCAGTCGATCG"
  homopolymer=False ✓

Homopolymer context:
  "ATCGAAAAAAAAGTCG"
        ^^^^^^^^
  8 consecutive A's → homopolymer=True
```

### Step 5: Apply RaVeX Filter Logic

```python
maf = add_ravex_filters(
    maf=maf,
    filters=["PASS"],           # Allowed caller filters
    noncoding=False,            # Keep noncoding? (default: True)
    homopolymer=False,          # Keep homopolymers? (default: True)
    ig_pseudo=False,            # Keep IG/pseudo? (default: True)
    min_alt_reads=2,            # Minimum alt read count
    blacklist=blacklist,
    whitelist=whitelist
)
```

#### Filter Decision Tree

For each variant, the following checks are performed in order:

```
ravex_filter = []

# 1. Minimum alt reads
if t_alt_count <= min_alt_reads:
    ravex_filter.append("min_alt_reads")

# 2. gnomAD frequency
if ingnomAD:
    ravex_filter.append("gnomad")

# 3. Blacklist
if blacklist:
    ravex_filter.append("blacklist")

# 4. Noncoding (if filtering enabled)
if filter_noncoding and noncoding:
    ravex_filter.append("noncoding")

# 5. Homopolymer (if filtering enabled)
if filter_homopolymer and homopolymer:
    ravex_filter.append("homopolymer")

# 6. IG/Pseudogene (if filtering enabled)
if filter_ig_pseudo and ig_pseudo:
    ravex_filter.append("ig_pseudo")

# 7. Variant caller filter
if isconsensus:
    # Use consensus filter (aggregated from multiple callers)
    if FILTER_consensus not in allowed_filters:
        ravex_filter.append("vc_filter")
else:
    # Use individual caller filter
    if FILTER not in allowed_filters:
        ravex_filter.append("vc_filter")
    # Penalize single-caller variants
    ravex_filter.append("not_consensus")

# 8. Whitelist override (absolute priority)
if whitelist and variant in whitelist:
    ravex_filter = ["PASS"]

# 9. Final decision
if not ravex_filter:
    ravex_filter = ["PASS"]

# Join all filters
RaVeX_FILTER = ";".join(ravex_filter)
```

### Step 6: Variant Deduplication

**Purpose:** Remove duplicate variant calls while preserving the highest-priority caller

**Priority order (default):**
```python
vc_priority = ["mutect2", "sage", "strelka", "consensus"]
```

**Deduplication algorithm:**

1. **Separate multiallelic variants:**
   ```python
   multiallelic = maf[maf["FILTER"].contains("multiallelic")]
   multiallelic["Caller"] = Caller + "_multiallelic"
   other = maf[~maf["FILTER"].contains("multiallelic")]
   ```

2. **Sort by caller priority:**
   ```python
   sorted_variants = []
   for caller in vc_priority:
       sorted_variants.append(maf[maf["Caller"] == caller])
   ordered = pd.concat(sorted_variants)
   ```

3. **Remove duplicates (keep first by priority):**
   ```python
   deduplicated = ordered.drop_duplicates(
       subset=["DNAchange", "Tumor_Sample_Barcode"],
       keep="first"
   )
   ```

4. **Merge consensus info for same DNAchange:**
   ```python
   # If same variant appears in multiple consensus runs
   # (e.g., DNA consensus + RNA consensus)
   grouped = deduplicated.groupby("DNAchange")
   
   for group in grouped:
       if len(group) > 1:
           merged['callers'] = "|".join(group['callers'])
           merged['filters'] = "|".join(group['filters'])
           merged['Tumor_Sample_Barcode_consensus'] = group['Tumor_Sample_Barcode'].iloc[1]
   ```

**Example deduplication:**
```
Before:
  Row 1: chr7:g.140453136A>T, Sample=RNA_realign, Caller=consensus
  Row 2: chr7:g.140453136A>T, Sample=DNA, Caller=mutect2
  Row 3: chr7:g.140453136A>T, Sample=RNA, Caller=consensus

After (by priority):
  Row 2: chr7:g.140453136A>T, Sample=DNA, Caller=mutect2 (kept - highest priority)
  Row 1: chr7:g.140453136A>T, Sample=RNA_realign, Caller=consensus (kept - different sample)
  Row 3: Merged with Row 1 (same DNAchange, both consensus)
```

### Step 7: Write Output

**Output structure:**
```
#version 2.4
##fileformat=MAFv2.4
##source=consensus3Callers (deepsomatic,mutect2,strelka)
[Original MAF headers preserved]

Hugo_Symbol  Chromosome  Start_Position  ...  RaVeX_FILTER  [New columns]
BRAF         chr7        140453136       ...  PASS          ...
EGFR         chr7        55249071        ...  min_alt_reads;not_consensus  ...
```

**New columns added:**
- `RaVeX_FILTER`: Final filter status (PASS or semicolon-separated failure reasons)
- `whitelist`: Boolean indicating if variant is whitelisted
- `blacklist`: Boolean indicating if variant is blacklisted
- `blk_reason`: Reason for blacklisting (if applicable)
- `ingnomAD`: Boolean indicating if variant exceeds gnomAD threshold
- `noncoding`: Boolean indicating noncoding consequence
- `ig_pseudo`: Boolean indicating IG gene or pseudogene
- `homopolymer`: Boolean indicating homopolymer context
- `CONTEXT`: Genomic sequence context (±10bp)
- `Tumor_Sample_Barcode_consensus`: Sample barcode from merged consensus (if applicable)

## Filter Categories Reference

| Filter | Column | Criteria | Default | Effect |
|--------|--------|----------|---------|--------|
| **min_alt_reads** | `t_alt_count` | ≤ 2 reads | Always applied | FAIL |
| **gnomad** | `MAX_AF` | ≥ 0.0001 (0.01%) | Always applied | FAIL |
| **blacklist** | Position overlap | In blacklist BED | If BED provided | FAIL |
| **noncoding** | `Consequence` | First consequence in noncoding list | Optional (keep by default) | FAIL if enabled |
| **homopolymer** | Context | 6+ consecutive identical bases | Optional (keep by default) | FAIL if enabled |
| **ig_pseudo** | `BIOTYPE` | IG genes or pseudogenes | Optional (keep by default) | FAIL if enabled |
| **vc_filter** | `FILTER` or `FILTER_consensus` | Not in allowed filters | Always applied | FAIL |
| **not_consensus** | `isconsensus` | Single-caller only | Always applied | FAIL |
| **whitelist** | `DNAchange` | In whitelist BED | If BED provided | PASS (overrides all) |

## Real-World Example

### Command
```bash
filter_mutations.py \
    -i TCRBOA7-T-RNA_realign_vs_TCRBOA7-N_realign_with_TCRBOA7-T_vs_TCRBOA7-N_with_TCRBOA7-T-RNA_vs_TCRBOA7-N.consensus.maf \
    --output TCRBOA7-T-RNA_realign_vs_TCRBOA7-N_realign_with_TCRBOA7-T_vs_TCRBOA7-N_with_TCRBOA7-T-RNA_vs_TCRBOA7-N.filtered.maf \
    --ref GRCh38.d1.vd1.chr7.mini.fa
```

### Input File Name Breakdown

The complex filename reveals a **3-way rescue consensus** strategy:

```
TCRBOA7-T-RNA_realign_vs_TCRBOA7-N_realign  (Level 1: RNA realigned)
    with_TCRBOA7-T_vs_TCRBOA7-N             (Level 2: DNA original)
    with_TCRBOA7-T-RNA_vs_TCRBOA7-N         (Level 3: RNA original)
```

**Strategy:** Maximize variant detection by combining:
- Realignment improvements (better mapping quality)
- DNA evidence (comprehensive genome coverage)
- RNA evidence (expression validation)

### Example Variants

#### Variant 1: High-Quality Somatic Mutation (PASS)
```
DNAchange: chr7:g.140453136A>T (BRAF V600E)
Hugo_Symbol: BRAF
Variant_Classification: Missense_Mutation
t_depth: 150
t_alt_count: 75
MAX_AF: 0.00001
Consequence: missense_variant
BIOTYPE: protein_coding
Caller: consensus
callers: "deepsomatic|mutect2|strelka"
FILTER_consensus: PASS
isconsensus: True
homopolymer: False

Filter checks:
  ✓ t_alt_count (75) > 2
  ✓ MAX_AF (0.00001) < 0.0001
  ✓ FILTER_consensus = PASS
  ✓ isconsensus = True
  ✓ Not in homopolymer context
  ✓ Coding variant

→ RaVeX_FILTER: PASS
```

#### Variant 2: Low Coverage Variant (FAIL)
```
DNAchange: chr7:g.55249071G>A (EGFR L858R)
t_alt_count: 2
MAX_AF: 0.00001
Caller: strelka
callers: "strelka"
FILTER: PASS
isconsensus: False

Filter checks:
  ✗ t_alt_count (2) <= 2
  ✓ MAX_AF (0.00001) < 0.0001
  ✓ FILTER = PASS
  ✗ isconsensus = False

→ RaVeX_FILTER: min_alt_reads;not_consensus
```

#### Variant 3: Germline Variant (FAIL)
```
DNAchange: chr7:g.12345678C>T
t_alt_count: 50
MAX_AF: 0.15 (15% in gnomAD)
Caller: consensus
FILTER_consensus: PASS
isconsensus: True

Filter checks:
  ✓ t_alt_count (50) > 2
  ✗ MAX_AF (0.15) >= 0.0001
  ✓ FILTER_consensus = PASS
  ✓ isconsensus = True

→ RaVeX_FILTER: gnomad
```

#### Variant 4: Failed Caller Filter (FAIL)
```
DNAchange: chr7:g.100000000A>G
t_alt_count: 10
MAX_AF: 0.00001
Caller: consensus
callers: "deepsomatic|mutect2"
filters: "weak_evidence|PASS"
FILTER_consensus: FAIL (only 50% PASS)
isconsensus: True

Filter checks:
  ✓ t_alt_count (10) > 2
  ✓ MAX_AF (0.00001) < 0.0001
  ✗ FILTER_consensus (FAIL) not in ["PASS"]
  ✓ isconsensus = True

→ RaVeX_FILTER: vc_filter
```

#### Variant 5: Whitelisted Variant (PASS - Override)
```
DNAchange: chr3:g.11111T>C
t_alt_count: 1
MAX_AF: 0.5
FILTER: FAIL
whitelist: True

Filter checks:
  ✗ t_alt_count (1) <= 2
  ✗ MAX_AF (0.5) >= 0.0001
  ✗ FILTER = FAIL
  ✓ whitelist = True (OVERRIDES ALL)

→ RaVeX_FILTER: PASS
```

## Command-Line Arguments

### Required Arguments
```bash
-i, --input          Input MAF file(s) (can specify multiple)
--output             Output filtered MAF file
--ref                Reference genome FASTA (for homopolymer detection)
```

### Optional Arguments
```bash
--gnomad_thr         gnomAD frequency threshold (default: 0.0001)
--whitelist          BED file with variants to force PASS
--blacklist          BED file with regions to exclude
--filters            Allowed caller filter values (default: ["PASS"])
--vc_priority        Caller priority for deduplication (default: ["mutect2", "sage", "strelka", "consensus"])
```

### Example with All Options
```bash
filter_mutations.py \
    -i sample1.maf sample2.maf \
    --output merged.filtered.maf \
    --ref GRCh38.fa \
    --gnomad_thr 0.001 \
    --whitelist known_pathogenic.bed \
    --blacklist problematic_regions.bed \
    --filters PASS weak_evidence \
    --vc_priority mutect2 sage strelka consensus
```

## Integration with Consensus Module

The filtering step seamlessly integrates with consensus output:

**Consensus provides:**
```
callers: "deepsomatic|mutect2|strelka"
filters: "PASS|PASS|weak_evidence"
FILTER_consensus: "PASS"  (≥50% of callers say PASS)
isconsensus: True
```

**Filtering uses:**
```python
if isconsensus == True:
    # Use aggregated consensus filter
    check FILTER_consensus against allowed_filters
else:
    # Use individual caller filter
    check FILTER against allowed_filters
    # Add penalty for single-caller variants
    add "not_consensus" to RaVeX_FILTER
```

**Result:** Multi-caller supported variants are treated more leniently than single-caller variants.

## Key Design Principles

1. **Whitelist has absolute priority** - Overrides all other filters
2. **Consensus variants preferred** - Use aggregated FILTER_consensus
3. **Single-caller penalty** - Non-consensus variants marked with "not_consensus"
4. **Context-aware filtering** - Homopolymers checked against actual sequence
5. **Priority-based deduplication** - Higher-priority callers preferred
6. **Multiallelic preservation** - Tracked separately to prevent loss
7. **Consensus merging** - DNA and RNA consensus for same variant are merged
8. **Full traceability** - All filter decisions documented in RaVeX_FILTER

## Expected Output Statistics

```
Total variants in input: 15,234
├─ Deduplicated to: 12,456
├─ PASS: 8,234 (66%)
├─ gnomad: 2,100 (17%)
├─ min_alt_reads: 1,234 (10%)
├─ not_consensus: 567 (5%)
├─ vc_filter: 234 (2%)
└─ Multiple filters: 87 (1%)
```

## Downstream Usage

The filtered MAF can be used for:

1. **Strict filtering:** Select only `RaVeX_FILTER == "PASS"`
2. **Lenient filtering:** Include variants with specific acceptable filters
3. **Custom filtering:** Use individual filter columns for custom logic
4. **Quality assessment:** Analyze filter distribution to assess data quality

**Example downstream filtering:**
```python
# Strict: Only PASS variants
strict = maf[maf["RaVeX_FILTER"] == "PASS"]

# Lenient: PASS or only not_consensus
lenient = maf[maf["RaVeX_FILTER"].isin(["PASS", "not_consensus"])]

# Custom: PASS or single acceptable filter
custom = maf[
    (maf["RaVeX_FILTER"] == "PASS") |
    (maf["RaVeX_FILTER"] == "homopolymer")
]
```

The filtered MAF provides complete transparency for quality control decisions, enabling flexible downstream analysis based on project-specific requirements.
