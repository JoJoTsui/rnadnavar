# Seq2neo regression

`tests/run_seq2neo_regression.sh` is the bounded DN/DT/RT acceptance test for
the pipeline's primary seq2neo use case. Its overlay includes the real
`examples/seq2neo/seq2neo.shared.config`, replacing only cohort-sized paths and
resource limits with the repository's provenance-preserving TCRBOA7 test
triplet and mini-reference assets.

The validator requires the seq2neo topology (trimming, DNA/RNA alignment, all
three active callers, normalization, consensus, rescue, realignment, RNA
editing, COSMIC/gnomAD, and VEP), rejects any failed task, and proves that the
external-alignment normalization processes remain unreachable from FASTQ
input. Consensus headers must contain `ENS_SUPPORT` and must not contain the
removed Wilson fields.

Run from the repository root:

```bash
bash tests/run_seq2neo_regression.sh
```

This bounded run precedes the production acceptance run:

```bash
PROJECT=PRJNA298376 PATIENT=4060 bash examples/seq2neo/run_single.sh
```
