# 02: Preserve trustworthy paired SNP evidence through consensus and rescue

**What to build:** Carry validated tumor/normal, caller, modality, allele, and round evidence from caller VCFs through DNA consensus and both rescue interfaces, with explicit missingness and compatible canonical INFO fields.

**Blocked by:** None (can start immediately).

**Status:** ready-for-agent

- [ ] Sample roles resolve consistently or fail clearly.
- [ ] Full available AD, DP, GT, AF, quality, and diagnostics remain linked to the caller and allele.
- [ ] Zero, absent, unavailable, and invalid evidence remain distinguishable.
- [ ] Consensus-only and full rescue round trips preserve evidence and labels.
- [ ] Legacy consumers remain readable through compatible aliases.
