# 07: FILTER/INFO contract + header fix (M8a, INFO enrichment)

**What to build:** The output contract is cleaned up without breaking the downstream model label contract. The FILTER column keeps the biological-class vocabulary (both model tryouts read it as labels). The bogus output sample column (a normal sample's name with empty FORMAT) is fixed — either a meaningful tumor sample column or none. INFO is enriched so every record carries the evidence needed to audit its label: per-caller tumor AD/DP/VAF and a classification-rationale field, making each FILTER value derivable from the record's own INFO.

**Blocked by:** 02: Tumor-sample-aware genotype extraction (C1) (enriched INFO carries tumor-derived values).

**Status:** ready-for-agent

- [ ] Output headers no longer carry a misleading sample column
- [ ] Every emitted record's FILTER is derivable from its own INFO fields — checked by a fixture test
- [ ] INFO includes per-caller tumor AD/DP/VAF and classification rationale
- [ ] FILTER vocabulary is unchanged (Somatic/Germline/Reference/Artifact/NoConsensus/RNAedit), verified against the model tryouts' expectations
- [ ] pytest tests pass
