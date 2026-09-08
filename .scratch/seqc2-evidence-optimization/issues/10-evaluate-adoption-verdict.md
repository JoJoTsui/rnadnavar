# 10: Evaluate the frozen candidate and issue an adoption verdict

**What to build:** Evaluate the frozen candidate on available held-out or replicate SEQC2 products, report uncertainty and per-type metrics, and issue an explicit pass, fail, or inconclusive adoption verdict.

**Blocked by:** 09: Select and freeze an evidence-guided candidate. If evaluating the STAR change, also blocked by 08: Enable upstream MAPQ60 for hybrid STAR runs.

**Status:** ready-for-agent

- [ ] DNA consensus, first rescue, and realignment rescue are evaluated separately.
- [ ] SNP, INDEL, and records precision, recall, and F1 are reported.
- [ ] Precision regressions and uncertainty are visible.
- [ ] Original FASTQ compatibility is checked before adoption.
- [ ] Verdict and artifacts are reproducible from frozen inputs and parameters.
