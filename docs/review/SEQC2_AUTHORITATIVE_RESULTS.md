# Authoritative SEQC2 manual-method comparison

For consensus VCFs, exact caller BAM/CRAM inputs, source workflow outputs,
and archived runs, use the [verified output bundle](SEQC2_VERIFIED_OUTPUTS_AND_ARCHIVE.md).

Verified 2026-09-14. These are historical manual native-consensus and gated
second-rescue results, not proof that deployed workflow defaults are equivalent.

All paths below are relative to the repository root. The result root is:

```text
examples/seqc2/comparison/historical_manual_replay_20260914/
```

Use only these consolidated files for current comparisons:

- `AUTHORITATIVE_RESULTS.json`: explicit allowlist of 12 metric files, their
  SHA-256 identities, full metrics, and retirement inventory.
- `AUTHORITATIVE_COMPARISON.csv`: 36 rows, covering SNP, indel, and records
  with TP/FP/FN, precision, recall, and F1.

Metric paths below are relative to that result root; `{target}` is `ukb` or
`medexome`.

| Dataset / method | Correct metrics path |
| --- | --- |
| WES-LL DeepSomatic | `wes_ll/{target}/deepsomatic.metrics.json` |
| WES-LL native | `wes_ll/{target}/historical_native.metrics.json` |
| WES-LL native + gated rescue | `wes_ll/{target}/historical_native_gated.metrics.json` |
| WGS-IL DeepSomatic | `wgs_il/{target}/deepsomatic.metrics.json` |
| WGS-IL native | `wgs_historical_scope/{target}/historical_native.metrics.json` |
| WGS-IL native + gated rescue | `wgs_historical_scope/{target}/historical_native_gated.metrics.json` |

Use the corresponding `.vcf.gz` query under `wes_ll/` or
`wgs_historical_scope/`. Commands and input identities are in `provenance.json`
and `wgs_historical_scope/provenance.json`. These benchmark query copies have
PASS labels and are not final VEP-annotated training VCFs.

## Cleanup

Thirty superseded paths were moved out of active comparisons into
`.artifacts/retired_seqc2_benchmarks_20260914/`, retaining recovery copies:

- `common_policy_20260914/policy_experiments/` (about 508 MB);
- its obsolete `policy_validation.json` and ten exact-key transition summaries;
- unrestricted WGS native/gated queries, metrics, and logs that omitted the
  recovered historical candidate restriction.

`retirement_manifest.json` in the archive maps every original path to its
recovery location. The original combined provenance is also archived. The
active provenance excludes the unrestricted WGS metrics. Caller controls,
truth/reference inputs, verified replay queries, original workflow outputs,
and Nextflow/Conda caches remain intact. Recovery copies do not free disk
space; they remove obsolete results from the active comparison namespace.

## Next action

Investigate the three WGS UKB rescue-only FPs (one also in MedExome), alongside
the true-positive WES rescue additions, before modifying biological thresholds.
Compare DNA tumor/normal AD, DP and VAF; caller FILTER, QUAL, TLOD and GERMQ;
gnomAD annotations; RNA support and missing coverage; editing and alignment
context. Realigned and first-pass RNA evidence must not count as independent
samples. Missing RNA evidence must not veto an accepted DNA baseline.

Document the candidate-universe restriction and manual-versus-deployed
differences now. Then test only evidence-supported changes to the gate across
both datasets and both regions, preserving the native baseline and reporting
per-type metrics. Do not remove named FP sites or tune separately by dataset.
Indel policy remains unresolved: historical parity depends on DeepSomatic-derived
indels and does not itself justify an indel consensus rule. Do not enable new
production defaults or rerun mapping/calling based on this comparison alone.

See the [verified replay report](SEQC2_HISTORICAL_NATIVE_GATE_REPLAY.md) for
counts and limitations. HG008 remains independent validation after policy freeze.
