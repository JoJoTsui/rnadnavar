# SEQC2 consensus and rescue logic evidence

Review date: 2026-09-07. Source: HEAD 931c476, original repository read-only. CodeGraph was consulted before source inspection. These are implementation findings and causal hypotheses; attributing a measured number of SEQC2 false positives/negatives requires the matched benchmark VCFs and truth-overlap tables. No workflow implementation was modified.

## Confirmed defects and ineffective controls

### 1. The alt-read support floor does not constrain the biological FILTER

`bin/vcf_utils/aggregation.py:404–430` defines eligible support as non-Artifact and at least `min_alt_support` tumor alt reads when AD exists. Lines 901–917 preserve every caller and normalized classification while collecting eligible `support_callers` separately; lines 930–938 compute `passes_consensus` from eligible support.

The final consensus classifier instead counts **all individual normalized labels** (`bin/vcf_utils/variant_classifier_unified.py:122–168`), without consulting `support_callers` or `passes_consensus`. `bin/vcf_utils/io_utils.py:1047–1059` writes that classification as FILTER, while line 1109 writes the separate support result into `PASSES_CONSENSUS`.

Consequently, two PASS callers with only two alt reads each yield FILTER=Somatic despite zero eligible supporters and PASSES_CONSENSUS=NO. This contradicts the stated support floor's purpose and matters because biological FILTER is the downstream model/benchmark selection contract. Adjusting `--min_alt_support` alone changes metadata without necessarily changing the selected somatic set.

The same gap exists in rescue promotion: `variant_classifier_unified.py:246–261,388–398` promotes from normalized label arrays, ignoring eligible support. One DNA caller and one RNA caller, each with two alt reads, can produce Somatic/rescue_promotion despite zero eligible support. Missing AD deliberately bypasses the floor (`aggregation.py:413–415,427–430`), an additional evidence policy requiring separate treatment rather than treating missing evidence as zero.

### 2. Standard process argument overrides are silently discarded

`modules/local/vcf_consensus/main.nf:21` and `modules/local/vcf_rescue/main.nf:25` assign `args = task.ext.args`, but the actual Python commands at lines 42–47 and 43–50 never interpolate `${args}`. Thus setting `ext.args` to `--min_alt_support`, `--disable_rescue_promotion`, `--rescue_min_dna_callers`, `--rescue_min_rna_callers`, `--rescue_veto`, or consensus exclusion flags does not affect these processes. Direct invocation of the Python scripts can still vary these controls.

Additionally, `conf/modules/consensus/vcf_consensus_workflow.config:10,20` sets `ext.thr`; the consensus module reads `ext.snv_thr` and `ext.indel_thr` instead (module lines 23–24). The setting happens to equal the defaults but is an ineffective knob. Rescue's `ext.snv_thr=2`, `ext.indel_thr=2` (config lines 30–31) are consumed correctly.

| Control | Defined default | Effective standard module behavior |
|---|---|---|
| SNV/indel threshold | 2/2 | Module ext.snv_thr/ext.indel_thr, otherwise 2/2 |
| Eligible tumor alt-read floor | 3 | Python default; ext.args override discarded; metadata/classification discrepancy above |
| Promotion | enabled | Python default; disable flag via ext.args discarded |
| Promotion minimum DNA/RNA | 1/1 | Python defaults; ext.args override discarded |
| Artifact veto | DNA | Python default; ext.args override discarded |
| Canonical chromosomes only | yes | Python default; ext.args override discarded |

Defaults are in `bin/vcf_utils/classification_config.py:27–53`. Executed `.command.sh` files remain authoritative for the historical benchmark; these are current source/module conclusions.

## Deliberate rules that can lose DeepSomatic true positives or add false positives

### 3. Equal categorical voting discards singleton and disputed DeepSomatic evidence

The classifier counts observed callers, requires at least two, then chooses the most frequent class and resolves all ties to Artifact (`variant_classifier_unified.py:129–171`). Therefore:

- A DeepSomatic-only Somatic call becomes NoConsensus even if its calibrated confidence is excellent.
- DeepSomatic Somatic plus one other caller's Artifact/Reference/Germline becomes Artifact on a 1:1 tie.
- DeepSomatic Somatic plus two negative labels can become a negative majority; no quality weighting protects it.
- Mutect2 and Strelka Somatic can outvote DeepSomatic Artifact/Reference/Germline.

This is a sensitivity/specificity policy, not proof of a programming bug. It cannot guarantee dominance over the strongest member: its success depends on truth-conditioned error overlap. The decision function never uses caller QUAL, read-count strength, normal evidence quality, variant context, or modality-specific calibration. Aggregated QUAL is a mean across caller scales (`io_utils.py:909`), not an ensemble posterior.

Raising the caller threshold to three also does **not** require three Somatic votes: it requires three observed records, followed by majority class. With three records, two Somatic plus one Artifact still passes. Threshold semantics need to be stated as observed-record threshold versus eligible-support threshold versus Somatic-vote threshold.

### 4. Normal-reference status and caller filter reasons are converted into biological labels

`classification.py:110–125` makes a failed Strelka call with NT=ref and normal DP>=2 into Reference. Reference normal genotype is compatible with a true tumor somatic variant; it does not by itself establish absence in tumor. This label then has the same voting weight as a genuine no-variant call.

`classification.py:215–225` maps Mutect2 haplotype to Germline and panel_of_normals/contamination/possible_numt to Reference; mixed technical+germline filters are classified Germline by priority. These are categorical heuristics, not equivalent measurements of the same latent biological state. `_counts_toward_support` excludes only Artifact, so these caller-rejected but relabeled Reference/Germline records can count as support if they meet AD (`aggregation.py:425–430`). The documentation phrase 'caller did not reject it' is therefore broader than the implemented non-Artifact predicate.

A reproduced example—DeepSomatic Somatic plus Strelka Reference, both AD=80,20—produces two eligible supporters/PASSES_CONSENSUS=YES but FILTER=Artifact because of a tie. This illustrates why support counts cannot be interpreted as Somatic votes.

`classification.py:150–151,190–205` also treats missing/unfiltered FILTER as Somatic. Standard pipeline inputs should be validated to ensure Mutect2 is post-FilterMutectCalls; if an unfiltered VCF enters a rerun, candidates become affirmative votes. This is a conditional input-integrity risk, not a claim that the measured benchmark used unfiltered VCFs.

### 5. RNA-only consensus bypasses a DNA evidence requirement

NoConsensus becomes absent (`variant_classifier_unified.py:235–245`). When only RNA has a substantive label, it is returned unchanged (`351–354`). RNA Somatic is therefore accepted with no DNA call or with DNA NoConsensus, even if the latter includes a DNA caller's rejected record. There is no DNA coverage/callability/normal-evidence gate on this branch.

The DNA Artifact veto is narrower than a DNA-support requirement: it only acts when an actual DNA consensus label is Artifact (`311–318`). It cannot protect an absent/no-consensus DNA site. RNA-only Somatic additions can expand recall for expressed DNA somatic variants but also add RNA editing or mapping-derived calls. Their contribution must be benchmarked as a distinct set, not pooled with DNA-supported rescue.

### 6. Veto and unanimity may block recovery of disputed true DNA variants

DNA Artifact vetoes RNA non-Artifact evidence (`311–318`), including Artifact caused by a categorical tie, not necessarily strong direct artifact evidence. Thus a genuine DeepSomatic DNA call tied against another caller can become permanently unavailable to RNA rescue.

When neither modality has consensus, promotion requires both within-modality classification sets to be internally uniform and equal (`368–398`). Any within-modality disagreement yields Artifact (`408–411`). Consequently the advertised 1-DNA + 1-RNA promotion is not generally 'at least one affirmative vote on each side': additional contradictory records can block promotion. Conversely a single caller/model reused in DNA and RNA can satisfy 1+1; there is no independent-caller-family requirement.

Both modalities holding different non-Artifact consensus labels yields Artifact if each has at least two observed individual callers (`290–305`). This condition uses caller presence, not eligible or Somatic support. These policies trade sensitivity for consistency and require truth-conditioned ablation before weakening them.

## Downstream interpretation traps

### 7. A file named filtered still preserves every biological FILTER and record

`bin/filter_rescue_vcf.py:314–337` copies the original FILTER, writes exclusion reasons into INFO/RaVeX_FILTER, and writes every record. `--min_alt_reads` defaults to 2 and `--gnomad_thr` to 0.0001 (lines 70–95). `unified_filters.py:182–207,310–319` tests maximum alt count across callers and gnomAD threshold. This stage cannot improve FILTER=Somatic-only precision by itself. A benchmark intended to evaluate these exclusions must explicitly select RaVeX_FILTER=PASS as well, with a separately documented selection rule. The stripped output chiefly changes representation/multiallelic inclusion; do not infer somatic QC from the filename.

The actual post-rescue order is COSMIC/gnomAD, then RNA-editing annotation, then filtering (`subworkflows/local/vcf_rescue_post_processing/main.nf:52–96`). Annotation is another classifier: DNA majority can reclassify to Somatic (`variant_classifier.py:377–411`) and cross-modality Somatic plus COSMIC recurrence can rescue (`413–453`). Current code checks common population frequency and prior artifact-veto rationale for the COSMIC rule. Results should isolate raw rescue, annotation changes, and INFO-filter selection rather than attributing every difference to rescue.

### 8. Exact allele matching is a prerequisite

Consensus joins chrom:pos:ref:alt, not haplotype equivalence (`io_utils.py:104–153`; `aggregation.py:833`). The workflow's preceding normalization needs to be verified on actual inputs, especially indels/MNVs. Residual equivalent-but-differently-represented calls split votes into singletons; a haplotype-aware benchmark can then count truth absent from consensus even when two callers represented the same event differently. This is a conditional representation issue, not evidence that normalization failed in the current run.

## Reproduced evidence

Executed with the repository `.venv/bin/python`, `PYTHONPATH=bin`, importing `aggregate_variants`, `compute_unified_classification_consensus`, and `compute_unified_classification_rescue`. In-memory single SNV dictionaries supplied explicit genotype AD/DP/VAF, normalized labels and classification. No fixture or source files were changed.

```text
consensus low alt support set() passes False classification
('Somatic', 'rule:majority|class:Somatic|votes:Somaticx2|callers:2|threshold:2')
rescue low alt support set() passes False classification
('Somatic', 'rule:rescue_promotion|class:Somatic|dna_consensus:none|rna_consensus:none|dna_votes:Somaticx1|rna_votes:Somaticx1')
consensus somatic+rejected support {'strelka', 'deepsomatic'} passes True classification
('Artifact', 'rule:majority_tie|class:Artifact|votes:Referencex1+Somaticx1|callers:2|threshold:2')
```

The first two cases used AD=98,2, DP=100, VAF=0.02 for both callers; the third used AD=80,20 and labels Somatic/Reference. These are actual decision-function reproductions, not counts of affected benchmark sites.

## Recommended ablations, in priority order

1. For the fixed historical artifacts and truth region, partition DeepSomatic PASS into retained, NoConsensus, majority-negative, tie, veto, annotation-demoted and representation-mismatch outcomes. Partition ensemble-only Somatic into DNA two-caller agreement, RNA-only consensus, cross-modality promotion, annotation rescue. Report TP/FP/FN and counts for each rule, separately SNV/indel and DNA-callable/expression strata.
2. Compare DeepSomatic alone; DNA 2-of-3; DeepSomatic-preserving backbone plus qualified additional sites; raw first rescue; RNA-only excluded; promotion off; promotion requiring two DNA or two RNA callers; veto variants. These are evaluation alternatives, not approved label-contract changes.
3. Fix the floor/FILTER inconsistency in an isolated implementation change before interpreting floor sweeps. Then vary floor 0/2/3/5 using actual command arguments; assess unknown-AD explicitly. Currently a floor sweep alone can be misleading.
4. Ensure process arguments are actually emitted before any Nextflow-based ablation. Capture command files, software versions/container digests, caller panel and thresholds.
5. Separate biological FILTER-only from FILTER plus INFO-QC selection. Measure each annotation transition; avoid making filtering appear effective by silently changing the benchmark selector.
6. Validate exact allele normalization; stratify false negatives by caller representation, VAF, normal depth, repeats, RNA splice proximity and callability. Only then consider learned/weighted integration; model confidence must be calibrated on held-out samples to avoid overfitting HCC1395.

The strongest conclusion from code alone is that the current integration has both information-losing policies and an unimplemented connection between eligible support and final classification. Neither proves that a revised ensemble will beat DeepSomatic on both precision and recall; that requires the rule-level truth attribution above.

Existing scoped regression verification completed: `.venv/bin/python -m pytest tests/vcf_utils/test_classifier_plumbing.py tests/vcf_utils/test_rescue_contract.py tests/vcf_utils/test_filter_info_contract.py -q` — **50 passed**, 15 NumPy scalar-conversion deprecation warnings, 32.79 seconds. Passing tests establish existing expected behaviors; they do not invalidate the synthetic support/FILTER counterexamples, which cover the missing connection between support eligibility and final class.
