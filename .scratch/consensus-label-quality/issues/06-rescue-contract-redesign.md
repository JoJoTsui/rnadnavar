# 06: Rescue contract redesign (M1/M2/M3 + small rules)

**What to build:** The rescue stage gets a real contract. (a) Promotion: when at least one DNA caller and at least one RNA caller (configurable) agree on Somatic at a site that failed within-modality consensus, the variant is rescued as Somatic and tagged as such in INFO — cross-modality agreement no longer collapses to NoConsensus. (b) Veto: a DNA Artifact label outranks RNA non-Artifact evidence; RNA can no longer flip DNA-flagged artifacts to Somatic. (c) Truthful flags: `RESCUED`, `CROSS_MODALITY`, and `PASSES_CONSENSUS_*` are computed from records that actually passed as Somatic, never from mere presence in the union file. Behavior changes are param-gated and documented in the rescue rules doc.

**Blocked by:** 05: Classifier plumbing (M4, M9, tie-fallback, min-alt floor).

**Status:** ready-for-agent

- [ ] A fixture variant with agreeing DNA+RNA Somatic callers but no within-modality consensus is rescued and INFO-tagged — regression test reproducing the original NoConsensus outcome
- [ ] A DNA-Artifact / RNA-Somatic conflict resolves against promotion — regression test
- [ ] `RESCUED=YES` never appears on a record whose per-modality outcomes were NoConsensus/Artifact — regression test
- [ ] The rescue rules documentation matches implemented behavior
- [ ] pytest regression tests pass; nf-test updated if the rescue module interface changed
