# Hybrid realignment design

Status: all eight design decisions and the consolidated design confirmed on 2026-09-06. The user requested specification and tickets only; implementation and a new run are explicitly deferred.

## Established requirements

The SEQC2 hybrid analysis uses caller-ready WES_LL normal and tumor DNA alignments plus pooled RNA-tumor FASTQ inputs. Enable RNA realignment while preserving the original seq2neo DN/DT/RT FASTQ workflow, existing BAM validation, the SEQC2 reference and BED configuration, original input files, and existing output directories. The intended final artifact is the annotated second-round rescue VCF. Performance interpretation is deferred until this branch is validated.

Terminology is defined in [CONTEXT.md](../CONTEXT.md). Hybrid ingress remains governed by [ADR 0004](adr/0004-converge-hybrid-inputs-at-caller-ready-boundary.md).

## Observed implementation

- `examples/seqc2/seqc2.hybrid.shared.config` omits `realignment` from `tools` and sets `realignment = false`.
- `workflows/rnadnavar.nf` gates the second pass on membership of `realignment` in `tools`. Its VCF route extracts RNA reads, realigns with HISAT2, reuses first-pass normal DNA, calls realigned RNA, and performs second-round rescue.
- `examples/seq2neo/seq2neo.shared.config` already selects VCF realignment and includes `realignment` in `tools`.
- The second-pass caller intervals come from candidate-region BED output, rather than directly reusing the first-pass global intervals.
- Effective candidates are all FILTER classes in RT consensus: `modules/local/vcf2bed/main.nf` does not filter records. Although the main workflow mixes rescue VCFs into the candidate channel, its RNA-status filter excludes current rescue metadata, which lacks `status`. Adding rescue candidates would change current behavior.
- The previous `.scratch/hybrid-ingress/spec.md` explicitly deferred realignment acceptance; its tickets 01–07 do not contain a dedicated enablement ticket.
- Hybrid inherits `enable_rna_annotation = false` and `enable_cosmic_gnomad_annotation = true` from the DNA-only SEQC2 config. Annotation policy therefore needs explicit consideration for the hybrid second pass.
- `examples/seqc2/hybrid/config_pooling_fix.yaml` checks first-pass caller/consensus/rescue artifacts; it cannot prove that second-round rescue completed.

## Accepted decisions — round 1

The user accepted all five recommendations on 2026-09-06. Entries below follow the interview question order.

1. Scope: reuse the existing VCF realignment route with all three callers, first-pass DNA consensus, and second-round rescue. Keep enablement in a hybrid-specific configuration; make shared changes only where needed for correctness and verify both ingress modes.
2. Candidate policy: initially preserve all RT-consensus FILTER classes as candidates. Adding DNA-only or rescue-only candidates is a separate behavior change requiring defined merging and sample matching, since preparation currently joins by patient alone.
3. Regions: keep the SEQC2 reference and first-pass BED, preserve candidate-derived second-pass intervals, and compare outputs using the same high-confidence `-R` and UKB `-T` benchmark regions. Record candidate BED provenance and overlap with the configured BED.
4. Annotation: enable RNA-editing annotation alongside COSMIC/gnomAD and VEP for the final second-round rescue artifact. Preserve the biological FILTER vocabulary; benchmark PASS conversions remain separate derivatives. This changes annotation policy relative to the completed hybrid baseline and must be recorded separately from realignment effects.
5. Completion: fail on missing expected inputs or failed second-pass tasks; report a verified zero-candidate case explicitly instead of substituting first-round rescue as final. Require trace and exact second-round annotated VCF/index evidence in the runner.

## Follow-up evidence needed

Verify candidate selection and sample joins, HISAT2 index/reference compatibility, annotation resource availability, empty-candidate behavior, final publication names, and regression tests exercising actual realignment for both hybrid ingress and the original seq2neo FASTQ configuration. A fresh output/state/log location is required for validation.

## Resource and annotation findings

The seq2neo-configured HISAT2 directory contains all eight `Homo_sapiens_assembly38.*.ht2` files and `genes.splice_sites.txt`. The configured REDIportal VCF and index exist, as do the VEP and gnomAD directories. Presence alone does not establish reference compatibility or usable cache versions.

`VCF_ANNOTATE` executes when either `vep` is selected or its realignment input is true. Second-round rescue passes true, so its VEP annotation does not require enabling VEP across every first-round caller. Its input is the filtered, stripped second-round rescue VCF. Acceptance must verify preserved biological FILTER and evidence INFO fields through annotation and identify the exact final VCF/index from the run.

## Accepted decisions — round 2

The user accepted Q6–Q8 on 2026-09-06.

- Q6, annotation comparison: apply the accepted RNA-editing and COSMIC/gnomAD policy consistently to both rescue rounds in the new hybrid run. Preserve both outputs so realignment can be compared under matched annotation settings. The completed historical run remains a separate baseline with RNA-editing annotation disabled. Use the existing second-round VEP execution without unnecessarily adding first-round per-caller VEP work.
- Q7, compatibility acceptance: require routing/contract tests, bounded real-execution tests for both hybrid inputs and a representative DN/DT/RT FASTQ triplet using seq2neo configuration, followed by the full SEQC2 hybrid realignment validation. Cover realignment disabled/enabled, pooled RT identity, missing/duplicate pairing, empty candidates, final-artifact validation, biological label/evidence preservation, and input/output isolation. Full seq2neo cohort rerunning is outside this validation scope; successful smoke tests are bounded evidence, not proof for every cohort sample.
- Q8, final-artifact filtering: the current rescue filter is invoked with `--filter_multiallelic`, and the VEP input strips FORMAT/sample columns. Retain this existing filtering behavior, designate `vcf_realignment/rescue/<id>/<id>.rescue.filtered.stripped.vep.vcf.gz` plus index as final, retain pre-filter rescue and filtered pre-VEP outputs for auditing, and account for every record removal. Verify retained records' biological FILTER and evidence INFO through VEP rather than asserting that pre-filter and final record sets are identical.

Existing seq2neo acceptance scaffolding includes `tests/config/seq2neo_regression.config`, `tests/run_seq2neo_regression.sh`, `tests/seq2neo/validate_regression.py`, `tests/csv/3.0/fastq_triplet_local.csv`, and `tests/data/seq2neo_chr7_smoke.bed`. Reuse and verify these fixtures during implementation; their existence is not a passing test result.

## Implementation sequence and acceptance

1. Verify resources and baseline contracts: resolve reference/index/cache compatibility, candidate identity and interval handling, and the existing final annotation contract. Record baseline input identities and output provenance.
2. Add isolated hybrid enablement and execution configuration: explicitly select VCF realignment, reuse validated first-pass DNA, configure annotations, and use fresh output/state/log paths. Preserve original seq2neo configuration behavior and existing BAM audit/normalization.
3. Fix demonstrated shared routing or validation defects: validate patient/sample pairing and pooled RNA identity, distinguish empty candidate input from missing inputs, and prevent silent branch loss. Keep all RT-consensus classes as candidates; do not introduce DNA/rescue-only candidates.
4. Require final-artifact completion evidence: three second-pass caller outputs, RT consensus, second rescue, requested annotation stages, indexed final VCF, and successful execution trace. A verified zero-candidate outcome must be distinct from successful production of the designated final artifact; first rescue cannot satisfy second-rescue completion.
5. Run targeted tests and bounded hybrid/seq2neo executions. Verify realignment-enabled and disabled behavior, caller/normal pairing, biological labels and evidence, intentional record removal, and input/output isolation. Compare unaffected first-pass DNA calls by variant content rather than compressed-file bytes or volatile headers.
6. Run the full SEQC2 hybrid validation after the bounded checks pass, with the same original inputs and new publication/state/log locations. Reuse valid cache entries where possible without changing files underlying prior published outputs. Record actual cache reuse rather than assuming every first-pass task is reusable after annotation/configuration changes.
7. Produce a validation handoff listing final artifacts, provenance, execution evidence, regression results, and remaining limitations. Subsequent benchmark comparison must use matched regions and annotation policy; performance interpretation remains deferred until branch validation completes.

No unresolved design choices remain. Resource verification and execution results are implementation evidence to collect, not questions for the user to answer. Final shared-understanding confirmation was received; current authorization covers specification and ticket preparation only.
