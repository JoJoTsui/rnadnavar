# Hybrid STAR failure investigation — 2026-09-07

## Findings

Two questions must be distinguished: why initial RNA alignment executed again,
and why that execution reported a malformed FASTQ record. The first is explained;
the underlying mechanism of the second is not yet established.

The latest successful baseline without realignment is `hybrid.pooling_fix`,
not the older `hybrid.modality` run. Its STAR command is byte-identical to the
failed full-run command, including 16 threads. Enabling downstream HISAT2
realignment did not change this initial STAR command.

| Execution | Work directory (relative to repository) | Threads | Result |
| --- | --- | --- | --- |
| Older modality baseline | `work/b5/5a22f570e93a6ea062b0805ed3255d` | 16 | Successful |
| Latest pooling-fix baseline | `work/84/98d0ad6937526cca52bf205f576e4d` | 16 | Successful, Sep 6 01:20:02 |
| Full-run failed attempt | `examples/seqc2/hybrid/work/27/e3bff6aec362b491e05e4316e60dfb` | 16 | Exit 104, Sep 6 22:42:03 |
| Full-run retry | `examples/seqc2/hybrid/work/82/1a77b1671dfa4a497edb75f3f29078` | 32 | Successful, Sep 7 00:31:22 |

Times are Asia/Shanghai. Evidence is retained in each task's `.command.sh`,
`.command.err`, `.exitcode`, staging symlinks, and STAR `Log.final.out`.
All three successful executions report 30,200,567 input read pairs and
25,706,589 uniquely mapped pairs. The older modality command differs in its
library read-group tag; that difference is already absent in the pooling-fix baseline.

## Why STAR was not cached

The successful baseline used repository-root `.nextflow` and `work`, with session
`09dfe166-b11f-432e-a648-28ef18c4ba9b`. The full realignment launch used
`examples/seqc2/hybrid/.nextflow` and its `work`, with session
`b9bca534-2e03-4159-83f1-e8664147ea62`. These are recorded in the respective
`.nextflow/history` files. Bare `-resume` in the latter launch context did not
select the former session's cache.

[Nextflow's cache documentation](https://docs.seqera.io/nextflow/cache-and-resume)
identifies session ID as a task-hash input and locates the default cache relative
to the launch directory. This is concrete evidence explaining the rerun;
the previous attribution to arbitrary global parameter or publish-directory
changes was not established.

## What failed inside STAR

Failed task 27 read the original RNA FASTQs through symlinks and `zcat`.
For R2 record `SRR9134727.3515511`, STAR logged a 73-base sequence but a
70-character quality string, ending `...A///<2`. A fresh targeted read of the
same source record returns 73 bases and 73 quality characters, ending
`...A///</</A`. Earlier full gzip integrity checks also passed.

Failed task 27 and successful retry 82 used the same source FASTQ targets,
GTF target, STAR index target, and STAR 2.7.11b executable. Their commands
differ in thread count (16 versus 32), so the retry is not a controlled
identical-command reproduction.

The evidence supports an intermittent discrepancy in the bytes STAR consumed
or interpreted. It does not establish which layer caused it. Storage reads,
decompression, STAR parsing/concurrency, and runtime memory remain hypotheses.
Ceph being the source filesystem is not proof of a Ceph fault. Present-day
integrity checks cannot reconstruct the failed execution's input stream.
There is no demonstrated deterministic FASTQ defect or direct causal link
from enabling realignment to this parser failure.

## Controlled reproduction results

The approved comparison completed September 7, 09:56:30 Asia/Shanghai.
Harness and raw evidence remain in `/tmp/hybrid-star-diagnostic.I7BmAg/`:
`run.py`, `status.json`, `driver.log`, and each arm's command, STAR logs,
exit code, and `zcat_status.log`. Temporary artifacts are not durable storage.

Both source FASTQs and independent copies were fully decompressed and checked
record by record for FASTQ structure and equal sequence/quality lengths.
Each mate contained 30,200,567 valid records. Compressed and decompressed SHA-256
checksums matched between source and copy. Source compressed checksums were
checked again after both STAR runs and remained unchanged:

- R1: `940849ea3e500fe8c83cff589c5d4add9776a93c81a9479727816f3dd49a5a37`
- R2: `6594fd82440d50f7e50f75e002fd09ea7e9ad8f7e6e72392fc431a85e66026de`

| Arm | Threads | Exit | Input pairs | Uniquely mapped pairs | Mismatch reproduced |
| --- | --- | --- | --- | --- | --- |
| Original source FASTQ paths | 16 | 0 | 30,200,567 | 25,706,589 | No |
| Verified copies on `/tmp` filesystem | 16 | 0 | 30,200,567 | 25,706,589 | No |

Both arms used the failed task's unchanged `.command.sh`, STAR 2.7.11b
environment, genome index target, and GTF target. A common PATH wrapper around
the real `zcat` recorded stderr and exit status while forwarding decompressed
output to STAR. All eight decompressor invocations (two mates, two passes,
two arms) exited 0. Both runs used separate output directories on `/tmp`.
Original pipeline inputs and outputs were not edited.

This comparison did not reproduce the failure. It demonstrates that the current
input succeeds at 16 threads from either location; it does not establish that
storage or thread count caused the historical failure. The wrapper, local
output filesystem, warmed caches, and absence of the original concurrent
workload mean these tests do not recreate every historical runtime condition.
Do not describe 32 threads, a STAR upgrade, or copying FASTQs as a proven fix.
No production code change or regression test claiming a resolved cause was made.

## Further diagnostic options

If the error recurs, preserve the decompressed stream and runtime diagnostics
from that failing execution before assigning a cause. Further tests could vary
output filesystem or concurrent load independently; the present comparison
does not justify a broad production change or blind thread-count tuning.

The later realignment tuple/join failure is a separate downstream incident and
does not explain the earlier STAR FASTQ error.
