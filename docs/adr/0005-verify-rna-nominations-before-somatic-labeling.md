# Verify RNA nominations before assigning Somatic training labels

Status: accepted policy; implementation pending.

RNA-only consensus nominates a site for verification using DNA tumor and matched-normal evidence; it does not by itself establish a Somatic training label. The user accepted this boundary after the SEQC2 review found 17 truth matches and 606 nonmatches among RNA-only rescue additions. Unresolved nominations remain outside the Somatic label set, accepting possible recall loss to avoid propagating RNA/DNA discordance as validated DNA somatic truth.

This supersedes RNA-only pass-through as the desired label policy, while preserving the existing biological FILTER vocabulary. Direct read evidence from original caller-ready DNA tumor/normal alignments may establish verification even when no DNA caller passed the site. Insufficient coverage and ambiguous evidence remain inconclusive and outside Somatic. Realignment reassesses RNA evidence and does not count as independent confirmation. Apply this admission policy to the final second-round rescue output and its selected re-consensus/QC derivatives, which preserve the realigned-RNA training lineage; first-round outputs are diagnostic controls. Numerical verifier thresholds remain to be established through development analysis and frozen before held-out evaluation; this decision does not claim that the existing workflow implements verification or that it improves both benchmark metrics.

Evidence and interview decisions: [SEQC2 comprehensive review](../review/SEQC2_COMPREHENSIVE_REVIEW.md), especially sections 2 and 8.
