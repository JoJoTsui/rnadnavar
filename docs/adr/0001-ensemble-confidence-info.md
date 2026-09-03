# Ensemble confidence as Wilson CI on caller support fraction, INFO-only

Status: superseded by ADR-0003

The decision below is retained as historical context.

Consensus VCFs now annotate every record with an ensemble-confidence interval:
`ENS_SUPPORT` (`k/n`), `ENS_CONF_LO`, `ENS_CONF_HI` — a 95% Wilson score
interval on the caller support fraction, where k is the number of supporting
callers (existing support semantics: non-Artifact classification plus the
min-alt-read floor) and n is **all** callers configured for the consensus
invocation. A caller with no record at a site counts as a non-support vote,
because absence means the caller evaluated the locus and did not call it —
treating absence as missing data would flatter low-support variants. The
annotation is deliberately INFO-only: the FILTER vocabulary is the frozen
training-label contract, and no categorical tier is derived yet — cutoffs
will be set only after SEQC2 benchmarking shows where precision actually
separates. Per-modality rescue-side variants (`*_DNA`/`*_RNA`) are deferred
until an RNA+DNA benchmark exists; consensus is per-modality de facto because
DNA and RNA run as separate `run_consensus_vcf.py` invocations.

## Considered Options

- Categorical tier from support count (HIGH/MEDIUM/LOW): simpler, but hard
  cutoffs with n=3 hide how little evidence three callers provide; the CI
  keeps that uncertainty explicit. The offline C1–C7 tier scheme in
  `bin/common/tier_config.py` was not reused because it is stats-only and
  reads INFO fields the writer never emitted.
- VAF binomial CI as the primary confidence: rejected — it measures sampling
  noise of the evidence, not ensemble agreement.

## Consequences

- With n=3 callers, `ENS_CONF_LO` takes only a handful of discrete values;
  downstream consumers should treat it as an ordinal signal, not a calibrated
  probability.
- The field names avoid bare "confidence" because `bin/label_qc.py` already
  uses that word for QC action levels.
