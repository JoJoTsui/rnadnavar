# 04: gnomAD guard fix (C3)

**What to build:** The population-frequency germline guard in rescue filtering actually fires. The annotation step writes a canonical gnomAD AF field and the filter probes it case-insensitively (or both standardize on one name), so an annotated common polymorphism (AF ≥ threshold) is flagged instead of silently passing with an effective AF of 0.0.

**Blocked by:** None (can start immediately).

**Status:** ready-for-agent

- [ ] Field-name lookup matches the annotated field (case-insensitive or standardized), verified on a real annotated VCF fixture
- [ ] Regression test: an annotated common-AF record is flagged by the gnomAD rule; before the fix it was not
- [ ] Threshold drift between the VCF path and MAF path is documented in the audit report (reconciliation itself may be deferred to the guidelines)
