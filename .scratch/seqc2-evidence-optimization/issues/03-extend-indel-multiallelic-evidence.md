# 03: Extend the evidence path to indels and multiallelic records

**What to build:** Extend the validated evidence contract to insertions, deletions, and multiallelic records without allele swapping, silent truncation, duplicate overwrite, or representation loss.

**Blocked by:** 02: Preserve trustworthy paired SNP evidence through consensus and rescue.

**Status:** ready-for-agent

- [ ] Indel and multiallelic caller evidence is allele-specific.
- [ ] Duplicate and conflicting records are diagnosed deterministically.
- [ ] Missing FORMAT values produce explicit unavailable states.
- [ ] Consensus and rescue round trips preserve indel evidence and biological FILTERs.
