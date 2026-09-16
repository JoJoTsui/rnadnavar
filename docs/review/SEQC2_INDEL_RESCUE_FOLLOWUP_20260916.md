# SEQC2 indel refinement and rescue-template audit

Completed development round after
[the fixed candidate experiment](SEQC2_CURRENT_NATIVE_EXPERIMENT_20260916.md).
No HG008 evidence was read. Production rules, original outputs and caller
caches remain unchanged. All experiments use existing VCFs; the rescue audit
only reads existing WES DNA BAMs. No alignment or calling was launched.

## Result

The refined consensus adds three WES UKB TP, one WES MedExome TP and one WGS
UKB TP to the previous candidate, with zero added FP and no removed TP in
any cell. WGS MedExome is unchanged. It becomes the next development comparator,
not a workflow default or a validated production label set.

| Dataset / target | DeepSomatic TP/FP/FN; F1 | Refined consensus TP/FP/FN; F1 | Refined + frozen realignment rescue TP/FP/FN; F1 |
|---|---|---|---|
| WES / UKB | 1048/38/1252; 0.619019 | 1052/35/1248; 0.621199 | 1063/37/1237; 0.625294 |
| WES / MedExome | 563/19/266; 0.798016 | 564/18/265; 0.799433 | 570/20/259; 0.803383 |
| WGS / UKB | 2168/19/132; 0.966347 | 2172/14/128; 0.968346 | 2172/14/128; 0.968346 |
| WGS / MedExome | 702/6/127; 0.913468 | 704/4/125; 0.916070 | 704/4/125; 0.916070 |

| Dataset / target | Refined precision / recall | Refined + rescue precision / recall |
|---|---|---|
| WES / UKB | 0.967801 / 0.457391 | 0.966364 / 0.462174 |
| WES / MedExome | 0.969072 / 0.680338 | 0.966102 / 0.687575 |
| WGS / UKB | 0.993596 / 0.944348 | 0.993596 / 0.944348 |
| WGS / MedExome | 0.994350 / 0.849216 | 0.994350 / 0.849216 |

All four aggregate F1 values beat DeepSomatic. Refined consensus also has
higher aggregate precision in all four cells. Rescue still adds two WES FP:
WES MedExome precision remains below DeepSomatic (0.966102 versus 0.967354).
The rescue gate was not tightened to eliminate those two loci.

## Variant-type results

SNP metrics are unchanged from the previous candidate. Indel rescue remains
disabled: new indels enter DNA consensus through caller evidence, not RNA.

| Dataset / target | Refined SNP TP/FP/FN; F1 | Refined indel TP/FP/FN; F1 | DeepSomatic indel F1 |
|---|---|---|---|
| WES / UKB | 1012/34/1193; 0.622578 | 40/1/55; 0.588235 | 0.585714 |
| WES / MedExome | 543/17/252; 0.801476 | 21/1/13; 0.750000 | 0.779661 |
| WGS / UKB | 2091/8/114; 0.971654 | 81/6/14; 0.890110 | 0.890052 |
| WGS / MedExome | 675/2/120; 0.917120 | 29/2/5; 0.892308 | 0.848485 |

WGS UKB indel F1 improvement is extremely small. WES MedExome indels still
trail DeepSomatic; do not claim superiority for every variant type. Detailed
per-type precision and recall are retained in each evaluation JSON.

## Failed hypotheses retained

Each row reports incremental TP/FP relative to the fixed previous candidate.
Selection uses caller evidence only; truth is used afterward to score it.

| Extension | WES UKB | WES MedExome | WGS UKB | WGS MedExome |
|---|---|---|---|---|
| Mutect2 PASS + DeepSomatic alternate reads | +2/+1 | +1/+0 | +0/+1 | 0/0 |
| Strelka PASS + Mutect2 corroboration | +1/+0 | 0/0 | +1/+1 | 0/0 |
| Low-QUAL DeepSomatic PASS + both DNA strands | +0/+1 | +0/+1 | +3/+0 | 0/0 |
| Union of the three | +3/+2 | +1/+1 | +4/+2 | 0/0 |

None qualifies as a shared TP-gaining, FP-nonincreasing replacement.
The reciprocal Mutect2 FP is the same deletion in WES and WGS:
chr2:197500318 CCGCT>C. Mutect2 phases it with a nearby deletion at 197500324;
ECNT=2. No nearby truth record was found in the normalized benchmark window.
This is a benchmark FP and a multi-event representation concern, not proof of
its biological cause. Merely raising TLOD would not solve it: TLOD is >20.

## Selected refinement

Keep the complete previous candidate, including its SNP rule and existing
indels. Add an indel if either branch below passes. Missing or malformed
required evidence fails the relevant branch. Multi-ALT/symbolic records are
not admitted by these extension assays.

```text
Existing fixed consensus ------------------------------------------ retain

New indel with clean Mutect2 normal evidence
  normal DP >=10, normal ALT count =0, GERMQ >=20, TLOD >0
  |
  +-- Reciprocal single-event branch
  |     Mutect2 PASS/., tumor ALT count >=3, ECNT =1
  |     DeepSomatic PASS/./RefCall, QUAL >0, tumor ALT count >=3
  |                                                               -> add
  |
  +-- Moderate-confidence, strand-corroborated branch
        DeepSomatic PASS/., QUAL >=10, tumor ALT count >=3
        Mutect2 tumor ALT count >=2
        Mutect2 alternate forward >=1 AND alternate reverse >=1
        Mutect2 filters only PASS/./contamination/weak_evidence
                                                                  -> add
```

The first branch is not restricted to DeepSomatic PASS. The second reduces
the earlier QUAL>=20 requirement only when both DNA strands corroborate it.
These are SEQC2-informed choices; their improvement is development evidence,
not held-out validation. Scores are not ensemble confidence probabilities.

The source Mutect2 header defines ECNT as the number of potential somatic
events in the assembly region; GERMQ is the Phred-scaled confidence against
germline origin; TLOD is the log10 likelihood ratio for the variant's presence.
The per-sample SB array supplies ref-forward/ref-reverse/alt-forward/alt-reverse
counts. Strand support is not the F1R2/F2R1 orientation requirement rejected
earlier for RNA. Callers observe shared reads; counts are not summed.

Single-event reciprocal additions recover two WES UKB TP (one MedExome TP)
and introduce no scored WGS additions. The moderate-confidence branch recovers
one WES UKB TP and one WGS UKB TP. No region-specific thresholds or locus
exceptions are used.

## Why the rescue gate remains unchanged

Read-only BAM inspection covered all 77 frozen WES rescue additions, both DNA
tumor and normal. Reads with unmapped/secondary/QC-fail/duplicate/supplementary
flags were excluded (0xF04), along with MAPQ<20 or BQ<20. Templates are grouped
by read group and query name; conflicting mate bases are not alternate support.

| Benchmark FP | Tumor alternate reads / templates | Normal alternate templates | Alternate read quality |
|---|---|---|---|
| chr17:76353649 C>T | 2 / 1 | 0 | MAPQ 60, BQ 41/32 |
| chr19:33206715 T>C | 2 / 1 | 0 | MAPQ 60, BQ 41/41 |

Both are proper-pair observations with one edit and no obvious low-quality
signature. Eight of eleven recovered TP also have only one alternate template.
Requiring two templates would remove both FP but also eight TP. RNA base
quality, mapping quality, position and normal evidence also overlap the TP.
This does not establish that either FP is biologically false or identify an
error mechanism; it rules out a simple quality/template-count separator.
Outside-HC additions remain unassessed, not FP.

## Remaining recoverable candidate scope

Using exact normalized alleles within the UKB benchmark domain:

| Dataset | Remaining indel FN | Present in at least one DNA caller | Absent from all three DNA caller VCFs |
|---|---:|---:|---:|
| WES | 55 | 7 | 48 |
| WGS | 14 | 10 | 4 |

Absent VCF records do not prove absent BAM coverage. Alternative haplotype
representations and RNA evidence may still exist. These counts bound simple
DNA-VCF selection opportunities, not biological recoverability.

## Reproduction and provenance

Root: `examples/seqc2/comparison/current_native_audit_20260916_v2/`.

- `{wes,wgs}_indel_followup/`: all failed/broader ablations, never overwritten.
- `{wes,wgs}_refined_indel/evaluation.json`: refinement ablations; select
  `cells["ukb/combined"]` and `cells["medexome/combined"]` for the new comparator.
- `{wes,wgs}_refined_indel/{ukb,medexome}/combined/query.vcf.gz`: consensus
  benchmark copies, not annotated training-label VCFs.
- `{wes,wgs}_refined_gate_replay/`: corresponding preserved realignment-gate
  replay; its historical candidate scope and annotations are inherited.
  It is not a fresh first-round or realignment-rescue driver execution.
- `wes_rescue_templates/templates.json`: all 77 sites, filtered read/template
  counts, alternate-read details, alignment identities and unchanged file stats.

The new assay scripts are `explore_indel_corroboration.py`,
`refine_indel_corroboration.py`, and `audit_retained_rescue_templates.py` under
`examples/seqc2/scripts/`. Gate replay uses `benchmark_frozen_gate_transfer.py`
with `--policy combined`. JSON artifacts include exact som.py commands,
source/script hashes, selection lists, and every variant-type metric.
Both indel experiments verified their VCF source hashes remained unchanged.

| Artifact under root | SHA-256 |
|---|---|
| wes_refined_indel/evaluation.json | b075622e62d7014c24ba2486bb32157ece58261b18ea6257b4689268f3982ba1 |
| wgs_refined_indel/evaluation.json | 272576efc98d6a7fda01145e5900befb64c8259cbf6dbeadcaea75da20901cc6 |
| wes_refined_gate_replay/evaluation.json | 7070414deb6a3ae56665b0bffb029ec42dc05ab754b4f2cee642e4ade02be426 |
| wgs_refined_gate_replay/evaluation.json | 64f5d8a18b669229c69883234dc5e1eb01f5a98b15bfb692de780c641895d6d7 |
| wes_rescue_templates/templates.json | 91856011b5e4bc1ee78fc8886d7c2a73f34bcbc267528c4d52c1178646bbe387 |

Thirty-two som.py comparisons completed in this round. Twelve focused tests,
Python compilation, and whitespace checks passed. No defaults changed.

Next prioritize reproducing the selected allele sets through the actual
consensus/rescue code with provenance intact before a frozen HG008 evaluation.
Keep WES MedExome indel recall and rescue precision as explicit unresolved
limitations; do not tune further on individual rescue FP coordinates.

## Standalone driver integration

The selected consensus rules are now available as the explicitly experimental
`bin/run_consensus_vcf.py --experimental-refined-native` path. Workflow configs,
default behavior and caller/alignment processes were not modified. Required
measurements come from native FORMAT rows using resolved sample roles; missing
normal DP is not replaced by AD-derived depth. The implementation does not
depend on SEQC2 sample names. Independent indel branches avoid an exploratory
helper quirk where missing ECNT could also suppress the strand branch.

The driver emits the existing biological FILTER vocabulary and per-caller
evidence. Admitted variants carry `rule:seqc2_refined_v2`, the branch, native
filters, QUAL, TLOD/GERMQ/ECNT and paired read measurements in their rationale.
It is not yet a released training-label policy.

`examples/seqc2/scripts/validate_refined_native_integration.py` runs the actual
consensus driver against existing normalized caller VCFs in each of the four
development cells, then compares Somatic allele sets with the frozen refined
queries. Outputs and the evolving report are under
`current_native_audit_20260916_v2/driver_parity_v1/`. Read `parity.json` for the
completion state; a partially populated report is not a four-cell validation.
WES/UKB completed with exactly 1,087 Somatic alleles and WES/MedExome with
582; both have zero extra and zero missing alleles. WGS/UKB and WGS/MedExome
also completed with exact matches (2,186 and 708 alleles). The report status is
`parity_passed`, with source and code integrity checks passed.
This is target-domain normalized-input parity, not a whole-genome raw-input,
annotation-stage, cohort or HG008 validation.

The rescue driver is not changed in this integration step. Inspection confirms
that its modality-disagreement and eligible-Somatic promotion rules are not
the historical replay's preserve-baseline/positive-nomination gate. Therefore
the combined rescue metrics above remain replay results, not current driver
results. Portable rescue integration must explicitly reconcile this contract
and verify both rescue rounds before promotion.

## Rescue eligibility reconstruction and benchmark validation

The opt-in experimental primitives in `bin/vcf_utils/refined_rescue_policy.py`
now reconstruct the gate from explicit caller identities and native evidence.
They are not wired into `run_rescue_vcf.py` or any Nextflow default. The audit
reads existing annotated **realignment rescue** candidates and the original
realigned RNA caller VCFs; it does not regenerate annotation or first rescue.

The historical gate has an important eligibility condition that was implicit
in `N_RNA_CALLERS_SOMATIC`: each RNA caller needs at least **three tumor
alternate reads**, not just native PASS. DNA nomination remains a separate
condition: a positive tumor alternate count from Mutect2/Strelka (including
filtered variants), or from native-PASS DeepSomatic. Neither RNA nor DNA read
counts are summed between callers. RNA callers are counted once within the
realignment round; first and realigned evidence are not independent votes.

Fresh native-PASS presence alone gave 78 WES and 142 WGS additions, versus the
historical 77 and 131. Restoring the three-read eligibility floor reproduces
all historical additions exactly, without reading truth labels. For example,
chr6:27402267 T>G has RNA DeepSomatic and Mutect2 PASS but each has AD=0,2;
the recorded zero eligible RNA votes is correct, not a missing-caller bug.

The audit additionally evaluates the complete Somatic rescue candidate set
before excluding the historical DNA baseline. This guards against silently
re-admitting SNPs removed by the new consensus. Combining these freshly
evaluated candidates with the actual refined consensus-driver outputs and
rerunning som.py reproduces **all SNP, indel and aggregate metrics** above:

| Dataset / target | TP | FP | FN | F1 |
|---|---:|---:|---:|---:|
| WES / UKB | 1063 | 37 | 1237 | 0.625294 |
| WES / MedExome | 570 | 20 | 259 | 0.803383 |
| WGS / UKB | 2172 | 14 | 128 | 0.968346 |
| WGS / MedExome | 704 | 4 | 125 | 0.916070 |

UKB query allele sets exactly match the prior replay. MedExome queries include
147 extra WES and 246 extra WGS SNPs because reconstruction uses the full rescue
candidate set; all are outside the MedExome target BED, and none lie in the
scored HC-plus-target intersection. Thus this is exact benchmark-domain
reproduction, not whole-query allele equality in MedExome.

Artifacts under `current_native_audit_20260916_v2/`:

- `rescue_portability_v1/` and `rescue_portability_v2/`: retained raw-PASS
  ablations showing why eligibility must not be dropped.
- `rescue_portability_v3/audit.json`: recorded versus raw versus eligible
  comparisons, full accepted candidate list, source/code hashes and integrity.
- `portable_rescue_benchmark_v1/evaluation.json`: all four fresh som.py
  commands, metrics, comparison metrics and query-set differences.
- `portable_rescue_benchmark_v1/{wes,wgs}/{ukb,medexome}/query.vcf.gz`:
  benchmark-only PASS queries, not training-label VCFs.

Audit SHA-256:
`27db22acaebc9e1374eba6fdf60ba8a9c6cc5d98c10e5a0a98bb20b7ab1b1214`.
Benchmark report SHA-256:
`0427bc293b3fce9eb3d56f6c4dfe0cce8a807393689f49277528aa72a307f630`.
All input checksum checks passed. No HG008 evidence was used, no original
outputs were edited, and no alignment or variant calling was launched.

Remaining integration work is an evidence-preserving rescue output adapter
(including negative labels and final annotation consistency), first-round
validation, and subsequent frozen HG008 evaluation. These gate primitives and
benchmark queries are not a completed production rescue-driver integration.
