# 03: Truthful RaVeX/filter flags (C2)

**What to build:** The filtering stages stop emitting garbage flags. Consensus output carries per-caller tumor alt-count INFO so the filtering stage can evaluate `min_alt_reads` from real evidence instead of a nonexistent FORMAT column; `vc_filter` no longer fires merely because consensus FILTER holds a biological class instead of PASS. The same fix applies to the rescue filtering path. On real inputs, a 2-caller Somatic consensus record no longer accumulates `min_alt_reads,vc_filter` flags.

**Blocked by:** 02: Tumor-sample-aware genotype extraction (C1) (needs tumor-derived alt counts written at consensus time).

**Status:** ready-for-agent

- [ ] Consensus VCFs carry per-caller tumor alt-count INFO
- [ ] `min_alt_reads` fires only when real tumor alt support is below threshold — regression test reproducing the original 100%-firing bug
- [ ] Biological-class FILTER values are exempt from the caller-rejection check
- [ ] Both the consensus filtering path and the rescue filtering path are fixed
- [ ] pytest regression tests pass
