# 12: Rerun smoke + cohort QC handoff

**What to build:** The end-to-end verification that the fixed pipeline plus the QC gate deliver clean labels. Run the rerun driver on 1–2 real samples from existing caller VCFs; assert new consensus + rescue VCFs are produced in the new location and every input checksum is unchanged; run `label_qc.py` over the smoke outputs and confirm the verdict artifacts are produced. Deliver a cohort runbook: the exact commands to rerun all 66 samples and gate them for training, with the PASS/WARN/FAIL table as the training-inclusion contract.

**Blocked by:** 08: Consensus+rescue-only rerun path; 10: label_qc.py Tier A.

**Status:** ready-for-agent

- [ ] Smoke rerun on 1–2 samples produces new consensus + rescue VCFs
- [ ] Input checksums before/after are identical (originals untouched)
- [ ] No alignment or per-caller calling processes appear in the rerun trace
- [ ] `label_qc.py` verdicts are produced over the smoke outputs
- [ ] Cohort runbook documents the full 66-sample rerun + gating commands
