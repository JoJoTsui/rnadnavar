#!/usr/bin/env python3
"""Check QC-DROP sites for 4032 against rerun FILTER assignments."""
import gzip
from collections import Counter

flagged = "/t9k/mnt/WorkSpace/data/ngs/zhanlingmin/ClairS_Train/truth_qc_out/flagged_sites.tsv.gz"
rerun = "/t9k/mnt/hdd/work/Vax/pipeline/rnadnavar/examples/seq2neo/output_reconsensus/PRJNA298330_4032/rescue/PRJNA298330_4032DT_vs_PRJNA298330_4032DN_rescued_PRJNA298330_4032RT_realign_vs_PRJNA298330_4032DN/PRJNA298330_4032DT_vs_PRJNA298330_4032DN_rescued_PRJNA298330_4032RT_realign_vs_PRJNA298330_4032DN.filtered.vcf.stripped.vcf.gz"

drop = {}
with gzip.open(flagged, "rt") as fh:
    header = fh.readline()
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if p[0] == "PRJNA298330_4032" and p[8] == "DROP":
            drop["%s\t%s\t%s\t%s" % (p[1], p[2], p[3], p[4])] = p[6]  # old filter
print("DROP sites:", len(drop))

found = Counter()
missing = 0
with gzip.open(rerun, "rt") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        p = line.split("\t", 7)
        key = "\t".join([p[0], p[1], p[3], p[4]])
        if key in drop:
            found[(drop[key], p[6])] += 1
            del drop[key]
missing = len(drop)
print("(old_filter, rerun_filter) counts:", found.most_common())
print("DROP sites absent from rerun:", missing)
