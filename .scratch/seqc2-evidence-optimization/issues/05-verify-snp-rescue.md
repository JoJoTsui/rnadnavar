# 05: Admit DNA-verified SNP nominations into first rescue

**What to build:** Verify RNA nominations against original DNA tumor and normal evidence, including nominations absent from DNA caller PASS, and emit confirmed, rejected, or inconclusive outcomes before first rescue labeling.

**Blocked by:** 04: Run an opt-in DeepSomatic-preserving DNA consensus.

**Status:** ready-for-agent

- [ ] Original BAM/BAI inputs are checked and never modified.
- [ ] SNP reference/alternate evidence is counted with read-filter provenance.
- [ ] Normal evidence can veto or make a nomination inconclusive.
- [ ] Only confirmed nominations can receive the Somatic rescue label.
- [ ] Benchmark output attributes each change and uncertainty.
