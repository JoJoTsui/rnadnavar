# Rescue VCF Documentation Consistency Notes

## Goal
Track documentation consistency for rescue behavior without changing algorithms or workflow code.

## Status
As of this update, rescue documentation has been synchronized to current implementation in:
- ../bin/vcf_utils/variant_classifier_unified.py
- ../bin/vcf_utils/io_utils.py
- ../subworkflows/local/vcf_consensus_workflow/main.nf
- ../subworkflows/local/second_rescue/main.nf

## Corrected documentation mismatches
1. Cross-modality disagreement handling
- Correct behavior: disagreement between non-Artifact DNA/RNA consensus labels is support-aware and can map to Artifact, a modality label, or NoConsensus depending on support.
- Documentation no longer states a fixed DNA-priority override in all disagreement cases.

2. Rescue output FILTER semantics
- Correct behavior: FILTER is unified biological class from rescue classification.
- Caller-level origin details remain in INFO fields.

3. First rescue versus second rescue naming
- Correct behavior: first rescue and realignment-based second rescue use distinct naming patterns.
- Documentation now includes concrete examples from output artifacts.

## Evidence anchor for naming and outputs
Verified rescue output examples in COO8801.shared:
- First rescue:
  - ../sequencing/aim_exp/rdv_test/output/COO8801.shared/rescue/COO8801DT_vs_COO8801DN_rescued_COO8801RT_vs_COO8801DN/COO8801DT_vs_COO8801DN_rescued_COO8801RT_vs_COO8801DN.rescued.vcf.gz
- Second rescue after realignment:
  - ../sequencing/aim_exp/rdv_test/output/COO8801.shared/vcf_realignment/rescue/COO8801DT_vs_COO8801DN_rescued_COO8801RT_realign_vs_COO8801DN/COO8801DT_vs_COO8801DN_rescued_COO8801RT_realign_vs_COO8801DN.rescued.vcf.gz

## Run configuration evidence
From ../sequencing/aim_exp/rdv_test/output/COO8801.shared/pipeline_info/params_2026-03-18_20-24-15.json:
- tools include rescue and realignment
- rescue_snv_thr = 2
- rescue_indel_thr = 2
- realignment_mode = vcf

## Remaining technical debt
No algorithmic or code changes were made in this work.
If future logic changes occur, this page should be updated together with the canonical guide.
