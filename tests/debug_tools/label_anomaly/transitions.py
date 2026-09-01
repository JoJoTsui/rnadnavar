#!/usr/bin/env python3
"""Per-record clean-vs-rerun FILTER transition analysis (streaming, stdlib only).

Reads clean VCF.gz into a dict key->FILTER, then streams rerun VCF.gz,
joining on (CHROM,POS,REF,ALT) and tallying transitions plus rerun-side
INFO details per transition class.
"""
import gzip, json, sys
from collections import Counter

INFO_KEYS = [
    "CLASSIFICATION_RATIONALE", "RESCUE_PROMOTED", "RESCUED",
    "PASSES_CONSENSUS", "PASSES_CONSENSUS_DNA", "PASSES_CONSENSUS_RNA",
    "N_SUPPORT_CALLERS", "CALLERS_SUPPORT", "N_DNA_CALLERS_SUPPORT",
    "N_RNA_CALLERS_SUPPORT", "GNOMAD_AF", "gnomAD_AF", "COSMIC_CNT",
    "ALT_COUNT_BY_CALLER", "ALT_COUNT_MAX", "VAF_BY_CALLER", "VC_CONSENSUS",
    "MODALITIES", "CROSS_MODALITY", "UNIFIED_FILTER",
    "GNOMAD_RESCUE", "COSMIC_RESCUE",
]
WANT = set(INFO_KEYS)


def parse_info(s):
    d = {}
    for f in s.split(";"):
        if "=" in f:
            k, v = f.split("=", 1)
            if k in WANT:
                d[k] = v
        elif f in WANT:
            d[f] = "FLAG"
    return d


def load_clean(path):
    m = {}
    n = 0
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.split("\t", 7)
            m["\t".join(p[:5])] = p[6]
            n += 1
    return m, n


def main(clean_path, rerun_path, out_path):
    clean, n_clean = load_clean(clean_path)
    trans = Counter()
    only_rerun = Counter()
    # per-transition rerun-side detail aggregates
    rationale = {}        # (cf,rf) -> Counter
    rescue_prom = {}      # (cf,rf) -> Counter (value YES/NO/absent)
    rescued = {}
    nsupport = {}         # (cf,rf) -> Counter of N_SUPPORT_CALLERS
    callers = {}          # (cf,rf) -> Counter of CALLERS_SUPPORT
    vc_cons = {}
    passes_dna = {}
    passes_rna = {}
    gnomad = {}           # (cf,rf) -> [n_present, n_ge_1e-3, n_ge_1e-2, n]
    cosmic = {}           # (cf,rf) -> [n_cosmic, n]
    altmax = {}           # (cf,rf) -> Counter bucketed ALT_COUNT_MAX
    modalities = {}
    cross = {}
    n_rerun = 0

    def bucket(pair, store, val):
        c = store.get(pair)
        if c is None:
            c = store[pair] = Counter()
        c[val] += 1

    with gzip.open(rerun_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            n_rerun += 1
            p = line.rstrip("\n").split("\t")
            key = "\t".join(p[:5])
            rf = p[6]
            cf = clean.pop(key, None)
            if cf is None:
                only_rerun[rf] += 1
                cf = "ONLY_IN_RERUN"
            pair = (cf, rf)
            trans[pair] += 1
            info = parse_info(p[7]) if len(p) > 7 else {}
            bucket(pair, rationale, info.get("CLASSIFICATION_RATIONALE", "<absent>"))
            bucket(pair, rescue_prom, info.get("RESCUE_PROMOTED", "<absent>"))
            bucket(pair, rescued, info.get("RESCUED", "<absent>"))
            bucket(pair, nsupport, info.get("N_SUPPORT_CALLERS", "<absent>"))
            bucket(pair, callers, info.get("CALLERS_SUPPORT", "<absent>"))
            bucket(pair, vc_cons, info.get("VC_CONSENSUS", "<absent>"))
            bucket(pair, passes_dna, info.get("PASSES_CONSENSUS_DNA", "<absent>"))
            bucket(pair, passes_rna, info.get("PASSES_CONSENSUS_RNA", "<absent>"))
            bucket(pair, modalities, info.get("MODALITIES", "<absent>"))
            bucket(pair, cross, info.get("CROSS_MODALITY", "<absent>"))
            af = info.get("GNOMAD_AF") or info.get("gnomAD_AF")
            g = gnomad.setdefault(pair, [0, 0, 0, 0])
            g[3] += 1
            if af is not None:
                g[0] += 1
                try:
                    v = float(af.split(",")[0])
                    if v >= 1e-3:
                        g[1] += 1
                    if v >= 1e-2:
                        g[2] += 1
                except ValueError:
                    pass
            c = cosmic.setdefault(pair, [0, 0])
            c[1] += 1
            cc = info.get("COSMIC_CNT")
            if cc is not None:
                try:
                    if int(cc.split(",")[0]) > 0:
                        c[0] += 1
                except ValueError:
                    pass
            am = info.get("ALT_COUNT_MAX")
            if am is not None:
                try:
                    v = int(am)
                    b = "0" if v == 0 else "1-2" if v < 3 else "3-4" if v < 5 else "5-9" if v < 10 else "10-19" if v < 20 else ">=20"
                except ValueError:
                    b = "?"
            else:
                b = "<absent>"
            bucket(pair, altmax, b)

    only_clean = Counter(clean.values())
    result = {
        "clean_path": clean_path,
        "rerun_path": rerun_path,
        "n_clean_records": n_clean,
        "n_rerun_records": n_rerun,
        "n_only_in_clean": sum(only_clean.values()),
        "transitions": [
            {
                "clean": k[0], "rerun": k[1], "n": v,
                "rationale_top": rationale[k].most_common(8),
                "rescue_promoted": dict(rescue_prom[k]),
                "rescued": dict(rescued[k]),
                "n_support_callers": dict(nsupport[k]),
                "callers_support_top": callers[k].most_common(6),
                "vc_consensus": dict(vc_cons[k]),
                "passes_consensus_dna": dict(passes_dna[k]),
                "passes_consensus_rna": dict(passes_rna[k]),
                "modalities": dict(modalities[k]),
                "cross_modality": dict(cross[k]),
                "gnomad[af_present,af>=1e-3,af>=1e-2,n]": gnomad[k],
                "cosmic[cosmic>0,n]": cosmic[k],
                "alt_count_max_buckets": dict(altmax[k]),
            }
            for k, v in trans.most_common()
        ],
        "only_in_clean_by_filter": dict(only_clean.most_common()),
        "only_in_rerun_by_filter": dict(only_rerun.most_common()),
    }
    with open(out_path, "w") as fh:
        json.dump(result, fh, indent=1)
    print(f"done {out_path}: clean={n_clean} rerun={n_rerun} only_clean={sum(only_clean.values())} only_rerun={sum(only_rerun.values())}", file=sys.stderr)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
