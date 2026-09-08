# SEQC2 artifact provenance output

The bounded label and benchmark tools may emit a JSON provenance manifest using the `seqc2-artifact-provenance.v1` schema. The manifest identifies the artifact path, size and SHA-256 digest; optional reference and region inputs receive the same identity fields. It records stage, sample, modality, library, selectors, effective arguments, commands, model, databases and tool metadata. Values that cannot be observed from execution must be recorded as `unknown` or a limitation; caller-supplied metadata is not independent proof of executable identity.

The manifest is an audit input for training-artifact selection and benchmark reports. It does not replace successful second-pass execution evidence, caller/BAM/input identities, or final-label acceptance. Existing labels and historical outputs remain unchanged.
