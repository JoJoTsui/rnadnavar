# Corrected three-class validation checkpoint — 2026-09-18

This lightweight snapshot records completed corrected-consensus screens and
paired-DNA pilot evidence for SEQC2 WES_LL, SEQC2 WGS_IL and HG008 WGS.
See [the evidence policy](../../NEGATIVE_LABEL_EVIDENCE_GATE.md) for definitions.

The approved Reference requirement is a 1% allele-fraction detection limit at
95% confidence **per DNA sample**: zero ALT observations and at least 299 usable
observations in each sample, under the idealized independent-observation model.
Other-allele evidence is also withheld. This is not biological training approval.

| Dataset | Germline SNPs with paired-read support | Reference SNPs passing the evidence gate |
| --- | ---: | ---: |
| SEQC2 WES_LL | 16 / 32 | 0 / 32 |
| SEQC2 WGS_IL | 15 / 32 | 0 / 32 |
| HG008 WGS | 17 / 32 | 0 / 32 |

These bounded candidate pilots are not precision estimates. Unmeasured records
remain withheld; all exported candidates remain training-ineligible. Somatic
truth overlap screens cannot establish Germline or Reference accuracy.

`summary.json` records source paths, metrics, integrity checks and copied-file
SHA256 values. `evidence/` preserves the detailed lightweight reports. HG008's
recommended-truth comparison is separate from its historical-truth screen.
The corrected SEQC2 consensus screens restore the frozen Somatic metrics in
both UKB and MedExome; this does not establish negative-label validity.

Rescue reports are **snapshots**, not completed validation claims: WES rescue
validation was running when captured; WGS and HG008 rescue plans were prepared
but not executed. Consult the original paths in the reports for later progress.
No full workflow or 66-sample cohort rerun was launched by this validation.
Original inputs and calling caches remain untouched.

Reproduce a new snapshot (use a new destination):

```bash
.venv/bin/python examples/seqc2/scripts/record_three_class_postfix.py --root examples/seqc2/comparison/three_class_postfix_20260918 --outdir /path/to/new/checkpoint
```

Heavy VCF/BAM-derived outputs remain under the ignored comparison directory;
only lightweight evidence and reproducible code are intended for version control.
