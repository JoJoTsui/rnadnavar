# 10: label_qc.py Tier A (VCF-only rules)

**What to build:** A standalone Python CLI (`label_qc.py`, following the repo's bin/ conventions) that gates cohort label quality without needing BAMs. It ports the external TruthQC's proven rules R1–R6 and expands them: sample-level gates (somatic-count outlier vs cohort, FILTER/INFO self-contradiction rate, RNA-only Somatic fraction, modality completeness, spectrum sanity — Ti/Tv, indel fraction, VAF distribution, per-caller contribution, coverage floors) and variant-level checks (gnomAD AF tiers, COSM_ hotspot-vs-frequency cross-checks, low-complexity/homopolymer context, clustered variants, RNA-editing overlap), with variant-level counterparts of sample rules where meaningful. Thresholds are cohort-adaptive (median/MAD) with absolute biological floors, all in a bundled overridable config. Output is the 4-part contract: human report, machine-readable summary + per-sample PASS/WARN/FAIL table, per-site flagged detail, cleaned VCFs only under an explicit apply flag. Inputs are never modified.

**Blocked by:** None (can start immediately; runs on existing cohort outputs).

**Status:** ready-for-agent

- [ ] Reproduces the external TruthQC baseline on the 66-sample cohort: the three known abnormal samples (4081, 4255, PRJNA298330_4032) are flagged by the expected rules
- [ ] Every rule is independently toggleable and threshold-configurable
- [ ] 4-part output contract emitted; cleaned VCFs only under explicit apply; inputs untouched
- [ ] pytest suite over the CLI seam with synthetic VCF fixtures per rule
- [ ] Design doc under `dev_docs/implementation/`
