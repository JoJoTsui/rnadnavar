#!/usr/bin/env python3
"""For 4255/4081: check clean-side COSMIC_RESCUE on clean-Somatic records that became NoConsensus in rerun."""
import gzip, sys
from collections import Counter

clean_path, rerun_path, sample = sys.argv[1], sys.argv[2], sys.argv[3]

# clean: Somatic records -> COSMIC_RESCUE / GNOMAD_RESCUE / FILTERS_ORIGINAL
clean = {}
with gzip.open(clean_path, "rt") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        p = line.split("\t", 8)
        if p[6] != "Somatic":
            continue
        info = {}
        for f in p[7].split(";"):
            if "=" in f:
                k, v = f.split("=", 1)
                if k in ("COSMIC_RESCUE", "GNOMAD_RESCUE", "COSMIC_CNT", "GNOMAD_AF", "RESCUED"):
                    info[k] = v
        clean["\t".join(p[:5])] = info

counts = Counter()
with gzip.open(rerun_path, "rt") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        p = line.split("\t", 8)
        if p[6] not in ("NoConsensus", "Artifact"):
            continue
        key = "\t".join(p[:5])
        c = clean.get(key)
        if c is None:
            continue
        counts[(p[6], c.get("COSMIC_RESCUE", "<absent>"), "cosmic_cnt>0" if c.get("COSMIC_CNT", "0").split(",")[0] not in ("0", "<absent>") else "no_cosmic")] += 1

print(sample, counts.most_common())
