# 02: Tumor-sample-aware genotype extraction (C1)

**What to build:** Consensus output VCFs carry genotype, AD, DP, and VAF extracted from the **tumor** sample of each caller's VCF — never blindly from sample index 0 (which is the normal for Strelka and Mutect2 in this pipeline). Tumor resolution follows each caller's convention: Strelka by the TUMOR sample name, Mutect2 via the normal-sample metadata, DeepSomatic by its sample ordering. All downstream aggregates (mean VAF/DP, per-modality VAF means, RNA-editing low-DNA-VAF inputs) consequently reflect tumor biology.

**Blocked by:** None (can start immediately).

**Status:** ready-for-agent

- [ ] A synthetic paired caller-VCF fixture (normal first, tumor VAF 0.40) yields consensus records whose VAF/genotype fields match the tumor, not the normal — regression test reproducing the original bug
- [ ] Strelka's max-read-count row selection is replaced by explicit tumor resolution, so DP and VAF never come from different samples
- [ ] Per-modality VAF aggregates and RNA-editing tier inputs use tumor values
- [ ] pytest suite covering the extraction seam passes
