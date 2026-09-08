# 01: Reproduce one complete SEQC2 error comparison

**What to build:** A reproducible benchmark boundary that selects `FILTER=Somatic`, ignores obsolete RaVeX fields, retains TP/FP/FN source attributes, separates SNP, INDEL, and records metrics, and records region/normalization semantics and uncertainty.

**Blocked by:** None (can start immediately).

**Status:** ready-for-agent

- [ ] Benchmark selection and region semantics are explicit and reproducible.
- [ ] TP/FP/FN records retain source-linked attributes and uncertainty.
- [ ] SNP, INDEL, and records precision, recall, and F1 are emitted.
- [ ] Existing benchmark fixtures and the HC-overlap edge case are covered.
