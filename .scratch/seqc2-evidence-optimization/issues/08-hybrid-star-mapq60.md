# 08: Enable upstream MAPQ60 for hybrid STAR runs

**What to build:** Add a parameter-backed unique-MAPQ60 setting for the hybrid STAR path, validate its bounded preprocessing effect, and preserve the original seq2neo command when the setting is unset.

**Blocked by:** None (can start immediately).

**Status:** ready-for-agent

- [ ] Hybrid configuration enables unique MAPQ 60 explicitly.
- [ ] Original FASTQ workflow leaves the STAR option unset and unchanged.
- [ ] SplitNCigarReads and downstream BAM checks still pass.
- [ ] Mapping provenance records the effective MAPQ convention.
