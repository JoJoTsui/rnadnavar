# SEQC2 FP/FN and INFO audit

Date: 2026-09-08. Source baseline: `152d304e69b3618be4d6df0f982a1d229fea0749` on `seqc2-consolidated`.

This report supersedes the preliminary follow-up's requirement to regenerate outputs before review. The user requested review of this commit with existing full hybrid results. Historical output generation differs from current source; findings below distinguish actual artifacts, current-source defects and unproven performance hypotheses. No production code, input or pipeline output was changed.

## Accepted scope

- Workflow labels are selected by `FILTER=Somatic` only. Ignore obsolete `RaVeX_FILTER` completely.
- gnomAD comes from annotated rescue VCFs in this run, not per-caller VCFs. Its absence from caller files is expected. Caller FORMAT and caller diagnostics come from original caller VCFs; classification/support come from consensus/rescue; database evidence comes from the annotated stage.
- Preserve eight-column consensus/rescue VCFs and biological FILTER vocabulary. Extend per-caller evidence compatibly, including tumor and normal where available; audit duplicate aliases and conflicting meanings first.
- STAR MAPQ 60 should be parameter-backed and initially hybrid-scoped. Preserve original seq2neo FASTQ configuration defaults and inputs/outputs.
- Aim for higher F1 without lower precision against DNA DeepSomatic, separately evaluating all three integration products and reporting SNP/INDEL uncertainty and regressions. Earlier stricter per-type precision/recall superiority remains a separately reported target. No superiority is claimed by this review.
- Freeze numerical rules before held-out evaluation; chr1 remains excluded from tuning. Same-cell-line WES-IL/WGS-IL replication is robustness evidence, not independent biological validation.

## 1. Measured changes in errors

Sources: full results in `examples/seqc2/hybrid/output/seqc2.wes.ll.hybrid.realign.full`; historical normalized scores in `examples/seqc2/hybrid/comparison/comprehensive_realign/WES_LL_T_1_vs_WES_LL_N_1`.

| Normalized som.py records | TP | FP | FN | Precision | Recall | F1 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| DNA DeepSomatic | 1048 | 38 | 1252 | .9650 | .4557 | .6190 |
| DNA consensus | 985 | 107 | 1315 | .9020 | .4283 | .5808 |
| First rescue | 1002 | 539 | 1298 | .6502 | .4357 | .5217 |
| Realignment rescue | 1002 | 508 | 1298 | .6636 | .4357 | .5260 |

The following attribution is exact CHROM/POS/REF/ALT matching with merged HC and target BED point membership. **Nonmatch is not an independently verified biological FP.** It is diagnostic, not replacement normalized or haplotype-aware scoring.

| Transition | Added truth matches | Added nonmatches | Lost truth matches | Removed nonmatches |
| --- | ---: | ---: | ---: | ---: |
| DeepSomatic to DNA consensus | 8 | 80 | 71 | 9 |
| DNA consensus to first rescue | 17 | 432 | 0 | 0 |
| DNA consensus to realignment rescue | 17 | 401 | 0 | 0 |
| First to realignment rescue | 2 | 325 | 2 | 356 |

Equal TP totals between rescue rounds hide two gains and two losses. Realignment changes hundreds of nonmatching alleles; it does not simply remove 31 false calls.

Rescue admission attribution:

| Admission rule | First rescue added truth / nonmatch | Realignment rescue added truth / nonmatch |
| --- | ---: | ---: |
| RNA-only consensus pass-through | 17 / 430 | 17 / 397 |
| Explicit cross-modality promotion | 0 / 2 | 0 / 4 |

The principal measured rescue problem is RNA-only pass-through. Tightening explicit promotion alone cannot address it. Consensus adds 80 nonmatches while losing 71 DeepSomatic truth matches: the two-caller policy exchanges useful calls for a less precise set rather than merely trading recall for precision.

### FN attribution

For DNA consensus, the exact diagnostic has 1314 FN:

- 1106 keys absent from the output and absent as exact keys from all three original DNA caller VCFs.
- 104 present but rejected downstream despite at least one DNA caller PASS.
- 104 present but rejected with caller records present but no caller PASS.

The first category is **not a depth diagnosis**. Raw DeepSomatic and Mutect2 have 35 and 76 multiallelic records in this domain; representation must be resolved before converting exact-key absence into biological absence. Some uncalled alleles may have adequate reads but fail caller modeling or filters. BAM interrogation is required to separate these causes.

First rescue has 937 absent-output keys plus 360 present-but-rejected truth keys; realignment rescue has 958 plus 339. Both retain 104 rejected truth keys with at least one DNA caller PASS. These are directly actionable integration-review strata before attempting global depth or aligner changes.

### gnomAD attribution

Counts below are truth matches / nonmatches among selected rescue records in the exact diagnostic domain. These are strata, not proposed filtering cutoffs.

| Rescue GNOMAD_AF | First rescue | Realignment rescue |
| --- | ---: | ---: |
| Missing | 699 / 173 | 698 / 266 |
| Zero | 103 / 79 | 104 / 70 |
| Greater than zero, below .001 | 197 / 262 | 197 / 159 |
| .001 to below .01 | 2 / 8 | 2 / 5 |
| At least .01 | 1 / 16 | 1 / 7 |

A common-AF rule alone addresses only a small fraction of the rescue nonmatches and intersects one truth match. Missing annotation is unknown, not AF=0. Database evidence needs allele/version provenance and must not substitute for matched-normal evidence.

## 2. Truth and benchmark checks

The truth-count discrepancy is explained: default `bcftools view -R HC -T target` includes 2300 records/95 indels, whereas explicitly selecting point overlap for both BEDs includes 2299/94. The extra deletion is:

`chr4:6690569 CCGTGGTGATAGGGCGGCCTTGCCGAAACAAGGCCACAT>C`.

No duplicate or multiallelic truth records explain this difference. It is BED overlap semantics, not a missing truth call or normalization bug. The diagnostic's point-domain counts must not silently replace official domain counts.

Installed som.py normalizes then performs allele intersection. Retain `-N`, `-R`, `-T` and no `-P`; Somatic-to-PASS conversion belongs only in benchmark copies. `--scratch-prefix` retains intermediate TP/FN/FP VCFs; `--feature-table generic` provides allele/outcome tags, not arbitrary INFO/FORMAT. Join these to source evidence rather than expecting generic features to contain AD/DP/gnomAD. `--count-filtered-fn` requires `-P`; full workflow VCF joins provide rejected-versus-absent attribution without changing official selection.

Further metric limitations: native precision/recall confidence intervals are omitted by the aggregate CSV; records can include complex/other events beyond SNP+INDEL; reported FP region size uses chromosome lengths and is inappropriate for HC-target FP/Mb. These do not by themselves invalidate the reported precision/recall counts. A future benchmark needs retained normalized alleles and haplotype-aware adjudication for representation-sensitive events.

## 3. INFO correctness and evidence preservation

### Observed output contradictions, partly already repaired in source

Full-file counts, not restricted benchmark counts:

| Artifact | All records | Somatic and PASSES_CONSENSUS=NO |
| --- | ---: | ---: |
| DNA consensus | 188775 | 27 |
| RNA consensus | 835637 | 43 |
| Realigned RNA consensus | 273614 | 232 |
| First rescue | 998530 | 51 |
| Realignment rescue | 440008 | 234 |

Example `chr1:3858673 G>A`: Somatic in DNA consensus and both rescues, `CALLERS_SUPPORT=.` and tumor alternate count 2 for each caller. Applying the current classifier to those labels/counts yields NoConsensus. The eligible-support fix already exists in this commit; do not open a redundant implementation ticket claiming it is absent. Historical imported labels can persist into rescue; a later validation run must regenerate the whole consensus-to-rescue chain under the frozen code.

### Duplicate concepts and inconsistent names

`N_SUPPORT_CALLERS` / `CALLERS_SUPPORT` concern eligible support, whereas `DNA_SUPPORT`, `RNA_SUPPORT` and some modality counts concern observed records. In first rescue, observed DNA counts differ from eligible DNA support on 133518 records and RNA on 379809; second rescue differences are 133518 and 198074. This is incompatible terminology, not proof that either raw count is arithmetically corrupt. These fields must not be interchangeable model features or evidence gates.

No duplicate caller aliases were observed in `FILTERS_NORMALIZED` across the five scanned full VCFs. Alias collision and duplicate-record overwrite remain source risks, not measured duplicate corruption in this dataset. Prefer a canonical semantic dictionary with explicit compatibility aliases rather than deleting similarly named fields blindly.

Relevant writer locations: `bin/vcf_utils/io_utils.py:927`, `:977`, `:1189`, `:1259`.

### Sample and allele evidence gaps

Current source extracts tumor GT/DP/AD/VAF/GQ, but emits only tumor GT/DP/VAF/alternate counts by caller. Full reference/alternate AD, normal metrics, per-caller GQ/QUAL, sample identity and caller-native diagnostics are not retained in that interface. Standard rescue rereads original callers; consensus-only rescue does not decode INFO back into per-caller evidence.

Actual DeepSomatic has only the tumor sample: normal FORMAT is unavailable there, not merely lost in serialization. Mutect2 and Strelka provide paired evidence. Preserve their caller-reported AF separately from read-count-derived fractions: at `chr1:842015 G>A`, Mutect2 normal AD=20,0 and AF=.071; tumor AD=9,2 and AF=.200. Treating caller AF and AD-alt/DP as identical would fabricate contradictions.

Confirmed source risks needing focused regression cases:

- Independent classification/genotype sample resolvers disagree on ambiguous names; partially named samples can assign both roles to one index. Actual SEQC2 canonical names avoid the reproduced cases.
- Short FORMAT arrays silently fall back to row zero. Missing sentinels and numeric bounds are not consistently validated.
- Alternate count/VAF extraction assumes ALT1; require validated decomposed biallelic inputs or per-allele indexing.
- Strelka SGT is looked up in FORMAT although it is INFO. Avoid inventing a conventional GT from an incompatible categorical field.
- Duplicate allele dictionaries overwrite records; repeated caller identities can overwrite genotype evidence while lists still append.
- Generic reader dictionaries omit arbitrary INFO. Rescue-only database annotations must be joined at their own stage.

Source references: `aggregation.py:45`, `:183`, `:334`, `:348`, `:401`, `:748`, `:776`, `:905`; `classification.py:342`; `run_rescue_vcf.py:360`.

## 4. Classification of remedies

| Category | Evidence and appropriate action |
| --- | --- |
| Already-fixed defect | Historical eligible-support/FILTER mismatch: validate current source, regenerate later; do not reimplement the existing fix. |
| Current correctness defects | Sample resolution, missing-value/allele handling, evidence loss, duplicate overwrite risks: specify explicit invariants and regression cases. |
| Proven policy weakness in these artifacts | Losing DeepSomatic truth matches and RNA-only admission with low incremental yield: test DeepSomatic-preserving integration plus DNA-verified additions. |
| Database limitation | Common AF flags few nonmatches and a truth match; use allele-aware database evidence with normal evidence, not a universal truth override. |
| Coverage/caller limitation unresolved | Exact-key absence does not prove low depth; normalize candidates and measure original tumor/normal BAM evidence. |
| Alignment hypothesis | STAR already becomes MAPQ60 in inspected SplitNCigarReads output; upstream 60 is a scoped convention fix, not proven FP reduction. |

RNA realignment is a reassessment of the same underlying reads, not an independent vote. Preserve caller/sample/modality/round evidence; never sum per-caller DP into an apparent larger sequencing depth. Keep both caller-reported AF and count-derived allele fractions with their denominators.

## 5. Specification readiness and remaining experiments

The review supports specifications for canonical evidence semantics, validated paired-sample/allele extraction, lossless evidence propagation, explicit RNA nomination admission, retained benchmark outcomes and scoped STAR MAPQ configuration. It does **not** support invented numeric verification thresholds or a claim that redesigned integration already beats DeepSomatic.

Prioritize correctness and an auditable DeepSomatic-preserving control, then measure evidence-qualified additions/removals. Maintaining all DeepSomatic calls only improves precision if additions are at least as precise as the baseline; otherwise selected baseline errors must also be removed. Current RNA-only additions are far below that bar. Run individual ablations with frozen evaluation domains and separate SNP/INDEL results.

Remaining bounded experiments: normalized/haplotype-aware allele attribution; tumor/normal read evidence for unresolved FNs and RNA additions; editing/splice/repeat context strata; frozen-rule WES-IL/WGS-IL and held-out validation. No BAM-depth cause or optimized threshold is claimed from this VCF audit.

Reproduce the exact diagnostic with `.venv/bin/python docs/review/seqc2-evidence/audit_current_attributes.py`. Generated evidence is `seqc2-evidence/current_attributes.json`; duplicate keys are listed explicitly. Original VCFs remain read-only.
