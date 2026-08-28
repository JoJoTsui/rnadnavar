# 11: label_qc.py Tier B (BAM-based verification via samtools)

**What to build:** The optional deep-check tier of `label_qc.py`. Given the DN/DT/RT BAMs, it verifies normal contamination (alt-VAF in the DNA normal at truth sites, via samtools mpileup subprocess — no pysam) and strand/orientation bias, folding the evidence into the sample verdict and per-site flags. A Rust/PyO3 backend beside the existing stats_core precedent is the designated extension point if profiling shows samtools is the bottleneck.

**Blocked by:** 10: label_qc.py Tier A (VCF-only rules).

**Status:** ready-for-agent

- [ ] Normal-contamination check reproduces the known 4081/4255 signature (normal alt-VAF ≥5% at most truth sites → FAIL)
- [ ] Tier B is strictly optional: Tier A output is unchanged when BAMs are absent
- [ ] No pysam dependency; samtools invoked as a subprocess
- [ ] pytest suite with a tiny synthetic BAM fixture
