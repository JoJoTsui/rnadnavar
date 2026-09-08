# Truth checks and SNP-only FN attribution

Supplement to [the FP/FN and INFO audit](SEQC2_FP_FN_INFO_AUDIT.md), 2026-09-08. Source baseline `152d304`. Read-only checks; original truth and workflow artifacts unchanged.

## SNP-only attribution

Using all 2205 truth SNPs and splitting comma-separated ALT alleles in original DNA caller VCFs removes the unsplit-multiallelic-key ambiguity for these SNP comparisons. MNV/haplotype representations can still differ.

| SNP FN category | DNA consensus | First rescue | Realignment rescue |
| --- | ---: | ---: | ---: |
| Absent output and all DNA callers | 1059 | 893 | 912 |
| Present output, absent all DNA callers | 0 | 151 | 132 |
| Rejected output, DNA callers present but none PASS | 103 | 101 | 101 |
| Rejected output despite at least one DNA caller PASS | 93 | 93 | 93 |
| Total FN | 1255 | 1238 | 1238 |

Totals match the native SNP FN metrics. This identifies 93 SNP truth alleles already passed by a DNA caller but lost from all three selected integration products. It does not establish that all 1059 DNA-caller-absent alleles have low BAM depth.

## BED boundary discrepancy

The additional truth deletion under default scoring is `chr4:6690569 CCGTGGTGATAGGGCGGCCTTGCCGAAACAAGGCCACAT>C`. It is included because of HC `-R` record overlap. `--regions-overlap 0` removes it; `--targets-overlap 0` changes nothing. Thus the exact diagnostic has 2299 records whereas native scoring has 2300. There is no truth duplicate or multiallelic record behind this discrepancy.

## Strong truth evidence conflicts with high population AF

Both rescue outputs give `chr14:106324560 C>A` `GNOMAD_AF=0.311307` and Somatic through `dna_consensus_only` with two DNA Somatic votes. The supplied SEQC2 truth marks it SOMATIC and records:

| Truth evidence | Value |
| --- | --- |
| Tumor VAF | .748 |
| Normal VAF | .001 |
| Passing / rejecting calls | 62 / 1 |
| BWA / Bowtie / Novo classification | Strong / Strong / Strong |
| PacBio tumor / normal VAF | .674 / 0 |
| Normal VAF interval in truth | .00003 to .0065 |

This is a population-annotation/truth conflict, not evidence that the truth is simply wrong. The supplied truth itself also reports cross-aligner inconsistent-VAF flags, so retain that context. Before turning population frequency into an unconditional veto, audit the database's source allele/reference representation and examine matched-normal evidence at the locus. The review does not relabel the truth or tune a locus-specific exception.

This locus explains why removing every rescue record with GNOMAD_AF at least .01 sacrifices one truth match while removing only 16 first-rescue or 7 realignment-rescue nonmatches. Population evidence alone cannot solve the measured RNA-only FP burden.
