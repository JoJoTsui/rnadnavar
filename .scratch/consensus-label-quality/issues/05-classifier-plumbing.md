# 05: Classifier plumbing (M4, M9, tie-fallback, min-alt floor)

**What to build:** The classification layer's plumbing tells the truth. CLI consensus thresholds (`--snv_thr/--indel_thr`) actually reach the rescue classifier instead of being silently dropped to defaults. Caller-support counts include only records the caller itself did not reject (non-Artifact), so two callers' rejected calls at the same site no longer count as support. The documented tie-to-Artifact rule is applied consistently in the unified-filter fallback (no Somatic-favoring tiebreak). A configurable minimum tumor alt-read floor exists so a 1–2-alt-read caller PASS cannot count as a full Somatic vote.

**Blocked by:** 02: Tumor-sample-aware genotype extraction (C1) (support/floor logic depends on tumor-derived counts).

**Status:** ready-for-agent

- [ ] Rescue classification behavior changes when CLI thresholds change — wiring test
- [ ] Caller-rejected records do not increment caller support — regression test
- [ ] Ties resolve to Artifact in every fallback path, matching the documented consensus rules
- [ ] Min-alt-read floor is configurable and applied to support votes
- [ ] pytest regression tests pass
