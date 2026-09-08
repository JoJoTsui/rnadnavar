# 09: Select and freeze an evidence-guided candidate

**What to build:** Run isolated current-code controls and finite non-chr1 ablations, attribute BAM-level causes, and freeze one evidence-guided candidate or an explicit no-candidate result with code, parameters, and source versions.

**Blocked by:** 01: Reproduce one complete SEQC2 error comparison; 03: Extend the evidence path to indels and multiallelic records; 04: Run an opt-in DeepSomatic-preserving DNA consensus; 06: Admit DNA-verified indel nominations into first rescue; 07: Apply verification through final realignment rescue.

**Status:** ready-for-agent

- [ ] Tuning excludes chr1 and keeps HC/target region definitions fixed.
- [ ] Controls and ablations are checksum- and parameter-recorded.
- [ ] DNA consensus, first rescue, and realignment rescue are compared with DeepSomatic.
- [ ] A frozen candidate or no-candidate outcome is written with rationale.
