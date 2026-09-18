# Cohort preparation verification

Date: 2026-09-19. Policy: `separated_three_class_v2`.
The validation archive was committed as `acd1975f` **before** cohort integration.

## Result

- New-policy read-only preflight passed for all 66 samples.
- Completed validation, archived evidence and current policy code hashes matched.
- All DNA DeepSomatic inputs expose native GQ/AD/DP schema.
- Source VCFs had readable headers, BGZF EOF markers and indexes recovering the
  tested first record; reference dictionaries matched. Full-record validation
  remains an execution check, not a claim of this preflight.
- 267 input checks found at least one index timestamp older than its VCF. These
  are warnings, not proof of corruption; tested indexed records were recovered.
  No original index was rewritten.
- Two 16 GiB worker budgets plus 8 GiB reserve fit the 83,751,862,272-byte cgroup
  limit. Estimated disk requirement: 662,851,262,912 bytes; filesystem free-space
  checks passed. This does not guarantee a separate storage quota.
- The new shared output and work directories did not exist after preflight.
  **No real cohort sample or real cohort pilot was generated.**
- Focused regression suite: **392 passed**, with 66 existing NumPy deprecation
  warnings. Includes synthetic execution of the old and new cohort modes,
  integrity/resume tests, stale-gate rejection, manifest and Parquet provenance,
  baseline-fallback rejection, and wrapper argument checks. Bash syntax passed.

## Provenance

Source manifest SHA256:
`1672f3ebb727bf0f1f743304ffd6af7fc61881cfab7b67c80fedf6a3430afd7f`.

Accepted validation summary SHA256:
`e3acf1bf3611b9b6d435ab963578f29a2a11a68c8c62d1c7dee872227061ee12`.

Detailed preflight (ignored generated artifact):
`examples/seqc2/comparison/separated_three_class_v2_20260919/cohort_candidate_preflight.json`.
SHA256: `9433e35a92f0f0172a793f2c3107d8d156bb464c8553d24abd7ded713094b02f`.

Executed command (preparation only):

```bash
.venv/bin/python examples/seq2neo/scripts/run_refined_cohort.py --config examples/seq2neo/config/separated_three_class_v2_cohort.json --plan examples/seqc2/comparison/separated_three_class_v2_20260919/cohort_candidate_preflight.json
```

See [the runbook](../../THREE_CLASS_COHORT_RERUN.md) for shared destinations,
pilot/full candidate-generation commands, and candidate export commands.
Biological training approval remains open; the new manifest and Parquet cannot
automatically designate these candidates as approved truth.
