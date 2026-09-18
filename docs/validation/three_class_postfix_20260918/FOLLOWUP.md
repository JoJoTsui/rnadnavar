# Rescue validation and negative-label yield follow-up

Date: 2026-09-18. Status: validation in progress; cohort rerun not authorized
by these results. Production classification and calling/cache definitions are
unchanged. Original inputs remain read-only.

## Committed checkpoint

- `f28ede3b`: provisional negative-evidence gate, tests and the approved 1% limit.
- `f4fab73f`: retain evidence in candidate exports without training approval.
- `790a499a`: corrected-consensus validation tooling and lightweight checkpoint.

## What usable negative-label yield means

Here negative means non-Somatic: Germline and Reference are separate classes,
not interchangeable negative genotypes. Keep four counts per class:

1. Candidates available in the VCF.
2. Candidates actually assessed with usable paired-DNA measurements.
3. Candidates passing the provisional evidence gate.
4. Candidates approved for training after biological and export/QC review.

Evidence yield is count 3 divided by count 2; report assessment coverage (count
2 divided by count 1) separately. Final training yield is count 4 divided by
count 1. Approval is currently absent, so count 4 is zero; that is not a claim
of zero biological accuracy. Read support is not independent truth.

The bounded 32-SNP-per-class pilots supported 16/32 Germline candidates in
WES_LL, 15/32 in WGS_IL and 17/32 in HG008. Reference support was 0/32 in each.
These are pilot evidence-retention rates, not precision estimates or a forecast
for all 66 samples. Unsupported includes insufficient coverage, observed ALT,
and unresolved context; it does not mean the variant is necessarily Somatic.
The strict Reference requirement remains 299 zero-ALT observations in **each**
DNA sample (1% detection limit, 95% confidence per sample, idealized model).

## WES: both rescue rounds completed

Same truth, HC mask, FASTA and target BED within each comparison. Aggregate
record metrics below are som.py measurements, not raw VCF record counts.

| Region | Stage | TP | FP | FN | F1 |
| --- | --- | ---: | ---: | ---: | ---: |
| UKB | Corrected three-class DNA consensus | 1052 | 35 | 1248 | 0.621199 |
| UKB | Three-class first rescue | 1002 | 24 | 1298 | 0.602526 |
| UKB | Three-class realignment rescue | 1012 | 24 | 1288 | 0.606715 |
| MedExome | Corrected three-class DNA consensus | 564 | 18 | 265 | 0.799433 |
| MedExome | Three-class first rescue | 541 | 15 | 288 | 0.781227 |
| MedExome | Three-class realignment rescue | 545 | 15 | 284 | 0.784737 |

The original source/code integrity checks passed. A full-query benchmark replay
matched the baseline metrics before attributing losses to its TP/FP partitions.
Read-only replay of every baseline Somatic transition matched the actual output.

UKB first rescue loses 50 TPs: 46 through inherited rescue labels alone (41
Artifact, 5 Reference), and 4 through population-AF vetoes (one also has an
Artifact label). Realignment loses 40 TPs: 36 through inherited labels alone
(31 Artifact, 5 Reference), and the same 4 with AF vetoes. Both remove 11 FPs:
10 with AF vetoes, one with an inherited Artifact label alone. These categories
partition the records; overlapping reasons must not be summed twice.

Across the complete VCF, first rescue withholds 1278 baseline Somatic records:
1225 have population-AF vetoes, 53 additional records have inherited-label-only
vetoes. Realignment withholds 1265 (1225 plus 40). Most are not scored in these
regional comparisons. Unmatched partition records are reported as unresolved,
not falsely counted as FP. No rejected/inconclusive DNA-verification value
contributed to these baseline losses. No newly verified Somatic additions were
made by either experimental round.

The behavior follows the provisional conflict-abstention policy, but it fails
the intended Somatic-retention validation. A legacy biological label derived
from failed Somatic calling is not independent native negative evidence. The
policy needs to distinguish inherited review flags from supported class
conflicts; do not silently tune thresholds or restore calls using truth status.
Population-AF conflicts need separate, allele-specific biological review.

## Collector review finding

Synthetic BAM tests confirm exclusion of duplicate, secondary, QC-failed and
low-MAPQ/BQ reads, and count an overlapping read pair once. They also reproduce
a gap: the current pileup defaults retain supplementary alignments, so two
alignments of one molecule can contribute two observations. This weakens the
independence assumption behind the Reference bound. The characterization test
passes by documenting that behavior, **not** by certifying it safe.

Recommended correction before expanding negative-evidence validation: explicitly
exclude supplementary alignments and record filter settings/provenance. Even
after correction, 299 reads do not prove full biological independence or remove
mapping/systematic-error uncertainty. No collector or policy change has been
made as part of this read-only review.

## Artifacts and reproduction

Root: `examples/seqc2/comparison/three_class_postfix_20260918/` (heavy ignored outputs).

- `seqc2_wes_ll/rescue_validation/validation.json`: both rounds completed.
  SHA256: `65895ca17ce431ce7bfb2b75bf5108ac0ce42257b59786abf94c2788e0bfc9e1`.
- `seqc2_wes_ll/rescue_loss_audit_v2/audit.json`: completed metric-parity and
  causal replay; SHA256 `70fc836bc099f8d85bb2a05f4c3e3c725647f35fc754d6f5cf01feb2962da857`.
  `seqc2_wes_ll/rescue_loss_audit_v3/audit.json` also completed, additionally
  checking hashes of all original source inputs; source/code integrity passed.
- `seqc2_wgs_il/rescue_validation/validation.json`: validation running at this
  checkpoint; first and realignment rounds are scheduled sequentially.
- `hg008_wgs/rescue_validation/validation.json`: likewise running, using the
  recommended tumorvariants truth. Historical truth is not substituted.

The initial loss-audit attempt stopped on an audit-script schema comparison
error (wrapped metrics versus their values). The numerical benchmark scores
were identical. The parser was corrected; failed artifacts remain separate.

```bash
.venv/bin/python examples/seqc2/scripts/audit_three_class_rescue_losses.py --validation examples/seqc2/comparison/three_class_postfix_20260918/seqc2_wes_ll/rescue_validation/validation.json --outdir /path/to/new/loss_audit
.venv/bin/python -m pytest tests/seqc2/test_three_class_rescue_losses.py tests/seqc2/test_negative_bam_collection_review.py -q
```

## Sequence and stop gates

Finish the two remaining dataset validations; retain all failure evidence.
Resolve the rescue-policy and evidence-collection findings before expanding
negative-label validation. Independent normal evidence, indel/haplotype context
and candidate-to-export checks remain required. Then consider a small cohort
pilot from existing outputs. Do not start a 66-sample regeneration or produce
an approved training manifest merely because the execution/test suite passes.
