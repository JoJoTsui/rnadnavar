#!/usr/bin/env python3
"""Compare clean-side vs rerun-side caller-level classifications for Germline->NoConsensus sites."""
import gzip, sys, json
from collections import Counter

clean_path, rerun_path, sample, out_path = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]

def get_info(s, keys):
    out = {}
    for f in s.split(";"):
        if "=" in f:
            k, v = f.split("=", 1)
            if k in keys:
                out[k] = v
    return out

CK = {"VC_CALLERS", "FILTERS_ORIGINAL", "FILTERS_NORMALIZED", "CALLERS", "N_SUPPORT_CALLERS"}
# clean side: keep only Germline records
clean = {}
with gzip.open(clean_path, "rt") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        p = line.split("\t", 8)
        if p[6] == "Germline":
            clean["\t".join(p[:5])] = get_info(p[7], CK)

print(f"clean Germline records: {len(clean)}", file=sys.stderr)

pair_vc = Counter()   # (clean VC_CALLERS, rerun VC_CALLERS)
pair_fo = Counter()
examples = []
with gzip.open(rerun_path, "rt") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        if p[6] != "NoConsensus":
            continue
        key = "\t".join(p[:5])
        c = clean.get(key)
        if c is None:
            continue
        r = get_info(p[7], CK | {"CLASSIFICATION_RATIONALE"})
        rat = r.get("CLASSIFICATION_RATIONALE", "")
        if "insufficient_modality_support" not in rat and "below_rescue_promotion_threshold" not in rat:
            continue
        cv = c.get("VC_CALLERS", "<abs>")
        rv = r.get("VC_CALLERS", "<abs>")
        pair_vc[(cv, rv)] += 1
        pair_fo[(c.get("FILTERS_ORIGINAL", "<abs>"), r.get("FILTERS_ORIGINAL", "<abs>"))] += 1
        if len(examples) < 12:
            examples.append({"site": key.replace("\t", ":"), "clean_VC": cv, "rerun_VC": rv,
                             "clean_FO": c.get("FILTERS_ORIGINAL", "<abs>"),
                             "rerun_FO": r.get("FILTERS_ORIGINAL", "<abs>"),
                             "rationale": rat})

out = {"sample": sample,
       "vc_pairs_top": [(a, b, n) for (a, b), n in pair_vc.most_common(15)],
       "fo_pairs_top": [(a, b, n) for (a, b), n in pair_fo.most_common(15)],
       "examples": examples}
json.dump(out, open(out_path, "w"), indent=1)
print(json.dumps(out, indent=1)[:4000])
