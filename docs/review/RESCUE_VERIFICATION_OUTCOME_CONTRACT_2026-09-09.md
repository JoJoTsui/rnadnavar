# Rescue verification and outcome contract

This records the approved Q6/Q8 changes to `classification.py`,
`variant_statistics.py`, and `run_rescue_vcf.py`. The rules apply equally to
first rescue and realignment rescue. They change downstream interpretation,
not alignment, variant calling, caller inputs, executor settings, or cache.

## Verification precedence

Rescue first computes the ordinary biological classification with its normal
consensus, promotion, and veto rules. An optional DNA verification manifest
can then withhold an otherwise-Somatic, RNA-dependent result. It cannot erase
Artifact, Germline, Reference, RNAedit, or NoConsensus outcomes. In particular,
an inconclusive manifest entry cannot turn a DNA Artifact veto into NoConsensus.

A Somatic result is independently supported by DNA if the DNA consensus passed
as Somatic, or the eligible DNA caller votes alone classify as Somatic under
the configured SNV/indel threshold and consensus majority rules. A single DNA
vote in a cross-modality promotion is not sufficient at a threshold of two.
Observed but ineligible DNA records cannot bypass verification.

Without independent DNA support, a supplied verification status other than
`confirmed` produces NoConsensus. An absent variant in a supplied manifest is
inconclusive; omitting the manifest preserves ordinary rescue classification.
The final rationale retains the original decision trace under a prior-rationale
prefix, and records the verification status and final class. `RESCUE_PROMOTED`
and `RESCUED` are cleared when the final class is not Somatic. Passed-input
consensus and cross-modality flags remain descriptions of the inputs.

## Descriptive outcomes and denominators

Statistics are computed after output classification and promotion flags are
finalized. The former `rescue_rate` divided promoted-plus-overlap records by
consensus overlap, so it could exceed 100%. It is removed in favor of
`rescued_union_fraction = final Somatic records flagged rescued / all union
records`. The fraction lies in [0, 1]; an empty union is reported as zero with
explicit 0/0 counts. Printed percentages are merely this fraction times 100.

For exact allele keys, let D be DNA input records with FILTER=Somatic and O be
final output records with FILTER=Somatic. Report:

| Field | Definition |
| --- | --- |
| `dna_somatic_baseline` | size of D |
| `output_somatic` | size of O |
| `somatic_retained_from_dna` | size of intersection of O and D |
| `somatic_new_vs_dna` | size of O minus D |
| `somatic_lost_from_dna` | size of D minus O |

Retained + lost equals the DNA baseline; retained + new equals output Somatic.
Counts are summed across chromosomes and the fraction is recomputed from the
summed numerator and denominator. These counts are not precision, recall, F1,
or evidence of improved truth accuracy. Statistics require labeled input record
dictionaries and finalized output labels, rather than treating unlabeled key
sets as if they were Somatic evidence.

## Regression coverage

`tests/vcf_utils/test_rescue_verification_outcomes.py` exercises negative-label
preservation, DNA Artifact veto, RNA-only withholding, confirmed and absent
verification, cross-modality promotion, DNA-independent support and eligibility,
final flag consistency, bounded fractions, baseline reconciliation, empty union,
and missing-label rejection. Existing rescue CLI tests check promotion and veto
against written VCF records.
