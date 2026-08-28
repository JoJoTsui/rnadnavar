# 01: Audit report + fix specs

**What to build:** A committed audit report under `dev_docs/audit/` capturing the full adversarial review: every finding (critical C1–C3, major M1–M10, minors) with file:line evidence, the empirical reproduction transcripts (synthetic-VCF runs proving the normal-sample extraction bug and the 100%-bogus filter flags), severity, and a concrete fix spec per finding. Also records the three abnormal samples' signatures from the external TruthQC run as calibration data.

**Blocked by:** None (can start immediately).

**Status:** ready-for-agent

- [ ] Report exists under `dev_docs/audit/` and covers C1–C3, M1–M10, and the minor findings
- [ ] Each finding has file:line evidence, severity, and a fix spec
- [ ] The three abnormal samples' TruthQC signatures (self-contradiction rate, RNA-only fraction, normal alt-VAF, count inflation) are recorded as calibration data
- [ ] A "what looks sound" section documents verified-correct behavior worth preserving
