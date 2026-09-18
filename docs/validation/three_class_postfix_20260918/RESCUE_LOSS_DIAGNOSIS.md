# Completed Somatic-loss diagnosis

Reviewed 2026-09-19. This supersedes the running-status snapshot in FOLLOWUP.md,
but does not overwrite the original validation outputs.

Both rescue rounds and the causal replay completed for all three datasets.
The replay matched the actual v1 transitions and benchmark metrics, with
original-source integrity checks. Comparisons use the same truth, high-confidence
mask, target and reference; HG008 uses the recommended **tumorvariants** truth.

## Cause, not a depth-induced loss

The v1 three-class rescue applied additional population-frequency and inherited
biological-label vetoes to an existing DNA Somatic baseline. The old labels were
not independent normal genotypes. Consequently a legacy Artifact or Reference
label could suppress a valid established Somatic call. The common-population-AF
veto also removed true positives. No rejected/inconclusive DNA-verification value
contributed to these baseline losses. These are deterministic rule transitions,
not evidence that the same alleles disappeared through sequencing depth changes.

UKB aggregate loss attribution, relative to each v1 DNA consensus:

| Dataset | Round | TP removed | Inherited-label-only TP | AF-veto TP, including overlaps | FP removed |
| --- | --- | ---: | ---: | ---: | ---: |
| SEQC2 WES_LL | First | 50 | 46 | 4 | 11 |
| SEQC2 WES_LL | Realignment | 40 | 36 | 4 | 11 |
| SEQC2 WGS_IL | First | 86 | 32 | 54 | 7 |
| SEQC2 WGS_IL | Realignment | 83 | 29 | 54 | 7 |
| HG008 WGS | First | 39 | 7 | 32 | 6 |
| HG008 WGS | Realignment | 38 | 6 | 32 | 6 |

The two TP-cause columns are disjoint. Records without a matching benchmark
partition allele are representation-unresolved, not automatically FP. The v1
rescue also changed new-addition admission by requiring confirmed verification;
it made no new Somatic additions in these runs. Preserving only DNA membership
would therefore be insufficient to recover the established rescue policy:
each output needs its own stage-specific Somatic baseline.

There was also a consensus-level interaction: native negative conflict handling
changed three established Somatic alleles to NoConsensus in each WGS dataset.
The separated policy restores their baseline membership without asserting that
the conflicting model class is correct. Full native decision traces remain in
`THREE_CLASS_NATIVE_RATIONALE`, including conflicts that resolved to NoConsensus
rather than to an admitted Germline/Reference candidate.

| Dataset | Allele (GRCh38) | Native-v1 abstention reason |
| --- | --- | --- |
| SEQC2 WGS_IL | chrX:71663958 C>T | Somatic/Germline nomination conflict |
| SEQC2 WGS_IL | chrX:135426830 AAAGCCCT>A | Somatic/Germline nomination conflict |
| SEQC2 WGS_IL | chrX:155767914 GTTTTGTT>G | Somatic/Germline nomination conflict |
| HG008 WGS | chr7:76911089 TC>T | Mutect2 normal-reference evidence conflicts with native Germline |
| HG008 WGS | chr10:100494825 AT>A | Native RefCall with substantial tumor ALT evidence |
| HG008 WGS | chr13:34910529 AT>A | Native RefCall with substantial tumor ALT evidence |

The HG008 three alleles are restored TPs in the UKB comparison (500 to 503 TP,
17 FP unchanged). A RefCall tag with substantial ALT reads is not a valid
Reference training label. Somatic preservation and negative abstention can
therefore coexist without accepting the contradictory negative nomination.

## Corrective policy and limits

[Separated three-class v2](../../SEPARATED_THREE_CLASS_V2.md) retains the exact
established Somatic allele set separately for DNA consensus, first rescue and
realignment rescue. It independently adds native DNA negative nominations only
outside that set. Legacy rescue classes cannot nominate negatives. Annotation
or native-class conflicts are retained as review flags, not silently resolved
using truth labels. This preserves performance; it does not certify every
retained Somatic call as biologically correct.

The negative evidence collector additionally excluded supplementary alignments
and recorded its read-filter settings. Reference remains subject to the approved
1%/95% bound (299 zero-ALT observations per DNA sample under the idealized model).
Germline support remains provisional; negative indels require haplotype-aware
validation. No candidate is automatically approved for model training.

## Evidence locations

All paths below are relative to
`examples/seqc2/comparison/three_class_postfix_20260918/`:

- `seqc2_wes_ll/rescue_loss_audit_v3/audit.json`
- `seqc2_wgs_il/rescue_loss_audit/audit.json`
- `hg008_wgs/rescue_loss_audit/audit.json`
- Each dataset's `rescue_validation/validation.json` contains both v1 rounds.

The corrected revalidation is in the separate
`examples/seqc2/comparison/separated_three_class_v2_20260919/` namespace. Earlier
`separated_three_class_v2_20260918` and `separated_three_class_v2_retry_20260918`
attempts are failed/interrupted development evidence, not accepted validation.
