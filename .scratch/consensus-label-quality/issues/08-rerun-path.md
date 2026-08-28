# 08: Consensus+rescue-only rerun path (M5 gating + driver)

**What to build:** A rerun driver that regenerates consensus and rescue VCFs for the whole cohort from the existing, read-only per-caller VCFs, writing to a new output location. The entry path guarantees normalized input (normalization can no longer be skipped when entering at the consensus step). FASTQ→BAM alignment and per-caller variant calling are structurally unreachable from this path — the driver builds its input manifest from existing caller VCFs and never references raw reads. Original outputs are verified untouched (checksums).

**Blocked by:** 03: Truthful RaVeX/filter flags (C2); 05: Classifier plumbing; 06: Rescue contract redesign; 07: FILTER/INFO contract (the rerun must run the fixed code end-to-end).

**Status:** ready-for-agent

- [ ] Entering at the consensus step always runs normalization first — nf-test or workflow-level check
- [ ] Driver consumes a manifest of existing per-caller VCFs and produces new consensus + rescue VCFs per sample
- [ ] No alignment or variant-calling process can be triggered by the rerun path
- [ ] Input files are opened read-only; a checksum guard proves originals are untouched
- [ ] Rerun outputs land in a new, separate output location
