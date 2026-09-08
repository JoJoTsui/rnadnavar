# 04: Run an opt-in DeepSomatic-preserving DNA consensus

**What to build:** Add an opt-in DNA consensus strategy that retains the matched DeepSomatic baseline unless an explicit evidence-qualified rule justifies a change, while leaving the existing default strategy unchanged.

**Blocked by:** 01: Reproduce one complete SEQC2 error comparison; 02: Preserve trustworthy paired SNP evidence through consensus and rescue.

**Status:** ready-for-agent

- [ ] Opt-in configuration is validated and provenance is emitted.
- [ ] DeepSomatic baseline calls are retained by default in the experiment.
- [ ] Any promotion or rejection has a named evidence rationale.
- [ ] Legacy default configuration and output labels remain unchanged.
- [ ] Comparison reports per-type precision, recall, F1, and paired transitions.
