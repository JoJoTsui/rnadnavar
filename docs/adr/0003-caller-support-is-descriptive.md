# Keep caller support descriptive

Status: accepted

Supersedes: ADR-0001

## Context

The Wilson interval previously written as `ENS_CONF_LO`/`ENS_CONF_HI` treats
the configured callers as independent binomial trials.  Mutect2, Strelka2, and
DeepSomatic are a fixed, correlated panel, so those bounds are not calibrated
probabilities of label correctness.  With three callers they are also only a
monotonic transformation of the support count.

The implementation derived the denominator from discovered VCF files.  A
missing caller could therefore turn intended support of 2/3 into 2/2.

## Decision

Consensus VCFs retain `ENS_SUPPORT=k/n` as descriptive caller-detection
support.  They no longer declare or write `ENS_CONF_LO` or `ENS_CONF_HI`.
The denominator is the explicit expected caller panel for that consensus
sample/modality, never the set of files that happened to be discovered.
Missing, duplicate, and unexpected caller inputs are fatal.

`ENS_SUPPORT` is distinct from agreement with the final biological class.  A
future probability-like confidence field requires calibration and validation
against held-out truth data.

## Consequences

Consumers keep a compact ordinal support signal and can distinguish 2/3 from
an incomplete panel.  The two Wilson fields are an intentional output-schema
removal; the frozen FILTER vocabulary is unchanged.
