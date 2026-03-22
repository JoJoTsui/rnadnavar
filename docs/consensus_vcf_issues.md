# Consensus VCF Documentation Consistency Notes

## Goal
Track documentation consistency for consensus behavior without changing algorithms or code.

## Status
As of this update, consensus documentation has been synchronized to current implementation in:
- ../bin/vcf_utils/variant_classifier_unified.py
- ../bin/vcf_utils/classification.py
- ../bin/vcf_utils/io_utils.py

## Corrected documentation mismatches
1. Tie handling in majority vote
- Correct behavior: tie at top class -> Artifact.
- Documentation now reflects this and no longer uses priority tie-break language.

2. FILTER meaning in consensus outputs
- Correct behavior: output FILTER is unified biological class.
- Caller-origin filter strings are retained in INFO fields.

3. PASSES_CONSENSUS interpretation
- Correct behavior: informational field in INFO.
- FILTER assignment remains controlled by unified classification logic.

## Evidence anchor for naming and outputs
Verified consensus output examples in COO8801.shared:
- ../sequencing/aim_exp/rdv_test/output/COO8801.shared/consensus/COO8801DT_vs_COO8801DN/COO8801DT_vs_COO8801DN.consensus.vcf.gz
- ../sequencing/aim_exp/rdv_test/output/COO8801.shared/consensus/COO8801RT_vs_COO8801DN/COO8801RT_vs_COO8801DN.consensus.vcf.gz

## Remaining technical debt
No algorithmic changes were made in this work.
Potential performance or refactoring tasks, if needed, should be tracked separately from documentation maintenance.
