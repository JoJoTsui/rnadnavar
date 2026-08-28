#!/usr/bin/env python3
"""label_qc.py — Tier A (VCF-only) label QC gate for EnsembleVar consensus/rescue VCFs.

Ports and expands the external TruthQC R1-R6 rules (reference implementation,
not imported) and adds cohort-adaptive sample-level gates. Pure stdlib: no
pysam, no third-party dependencies.

Variant-level rules (evaluated on records whose FILTER contains the target
label, default "Somatic"):
    R1  common polymorphism   gnomAD AF >= af_common (tiered: af_strong -> MID)
    R2  zero DNA endorsement  no DNA_* caller judged Somatic in FILTERS_CATEGORY
    R3  hotspot contradiction COSMIC hit AND common gnomAD AF
    R4  category conflict     co-located Germline/Reference record at same site
    R6  self-contradiction    FILTER=Somatic but UNIFIED_FILTER in {Reference, Germline}
    R7  RNA-editing overlap   REDIportal evidence on a Somatic-labelled record
    R8  clustered variants    >= min_records Somatic calls within window_bp
    R9  low caller agreement  N_SUPPORT_CALLERS <= max_support_callers

Disposition tiers (TruthQC-compatible):
    HIGH -> DROP candidate, MID -> RELABEL candidate, LOW -> report only.

Sample-level gates (cohort-adaptive median/MAD outliers with absolute floors;
all thresholds in the bundled label_qc_config.json, overridable via --config):
    S0  flag burden           (HIGH+MID) share of Somatic records
    S1  somatic-count outlier vs cohort
    S2  self-contradiction rate (R6-class share of Somatic records)
    S3  RNA-only Somatic fraction at common gnomAD AF (germline leakage)
    S4  normal contamination  (Tier B: normal alt-VAF at truth sites via
                               batched `samtools mpileup`; needs --verify-bam)
    S5  modality completeness (DNA callers participated at all)
    S6  spectrum sanity       (Ti/Tv, indel fraction, het-like VAF, caller share)
    S7  coverage floors       (fraction of Somatic records below DP floor)
    S8  strand bias           (Tier B, optional WARN gate; off by default)
Verdict per sample: PASS / WARN / FAIL.

Tier B (strictly optional, --verify-bam): BAM-based verification using samtools
as a subprocess (no pysam). Given the DNA-normal BAM it measures alt support in
the normal at the sample's Somatic truth sites (S4); given the DNA-tumor BAM it
computes per-site strand counts of alt reads and flags strand-biased sites
(B1, variant-level REPORT_ONLY). Without BAMs the Tier A output contract is
byte-identical (samples_qc.tsv / flagged_sites.tsv.gz unchanged).

Output contract (4 parts): report.md; summary.json + samples_qc.tsv;
flagged_sites.tsv.gz; cleaned VCFs only under an explicit --apply flag.
Inputs are never modified.

Usage:
    label_qc.py --samples manifest.tsv --out label_qc_out
    label_qc.py --samples manifest.tsv --out label_qc_out --apply
    label_qc.py --samples manifest.tsv --out label_qc_out --verify-bam
    label_qc.py --truth-vcf sample.vcf.gz --sample-id S1 --out label_qc_out \
        --verify-bam --normal-bam dn.bam --tumor-bam dt.bam
"""

import argparse
import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from bisect import bisect_left, bisect_right
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor

VERSION = "0.2.0"

DEFAULT_CONFIG_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                   "label_qc_config.json")

COLUMN_ALIASES = {
    "sample_id": ("sample_id", "sample", "id", "sampleid"),
    "truth_vcf": ("truth_vcf", "truth", "vcf", "truth_vcf_fn", "truthvcf"),
    "normal_bam": ("normal_bam", "dn_bam", "normal", "nbam"),
    "tumor_bam": ("tumor_bam", "dt_bam", "tumor", "tbam"),
}
REQUIRED_COLUMNS = ("sample_id", "truth_vcf")

STRONG_RULES = ("R1", "R2", "R3", "R7")
WEAK_RULES = ("R8", "R9")

TRANSITIONS = frozenset(("AG", "GA", "CT", "TC"))

FLAGGED_HEADER = ("sample", "chrom", "pos", "ref", "alt", "filter", "rules",
                  "confidence", "action", "gnomad_af", "cosmic_id",
                  "unified_filter", "filters_category", "vaf_mean", "dp", "note")

SAMPLES_QC_HEADER = ("sample", "verdict", "n_records", "n_somatic", "n_flagged",
                     "n_high", "n_mid", "n_low", "actioned_rate",
                     "contradiction_rate", "rna_only_fraction",
                     "rna_only_common_fraction",
                     "dna_participation", "ti_tv", "indel_fraction",
                     "median_vaf", "het_vaf_fraction", "median_dp",
                     "low_dp_fraction", "gates_failed", "gates_warned", "notes")

# Appended to samples_qc.tsv only when Tier B actually measured at least one
# BAM; otherwise the Tier A output contract is unchanged.
TIER_B_COLUMNS = ("normal_contam_fraction", "n_normal_sites_eval",
                  "n_strand_bias_sites")


# ------------------------------------------------------------------ helpers

def open_text(path):
    """Transparently open plain-text or gzip-compressed text."""
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def truthy(value):
    """Lenient truth test for INFO boolean/string values."""
    if value is True:
        return True
    if value is None:
        return False
    return str(value).strip().upper() not in ("", "NO", "0", ".", "FALSE", "NONE")


def parse_info(text):
    """Parse a VCF INFO column into a dict; bare flags map to True."""
    out = {}
    for tok in text.split(";"):
        if not tok:
            continue
        if "=" in tok:
            key, val = tok.split("=", 1)
            out[key] = val
        else:
            out[tok] = True
    return out


def parse_floats(value):
    """Parse possibly multi-valued numeric INFO fields into floats."""
    out = []
    for item in re.split(r"[,|]", str(value)):
        item = item.strip()
        if not item or item == ".":
            continue
        try:
            out.append(float(item))
        except ValueError:
            pass
    return out


def sanitize_sid(sid):
    return re.sub(r"[^A-Za-z0-9._-]", "_", str(sid))


def natural_chrom_key(chrom):
    m = re.match(r"^(?:chr)?(\d+)$", chrom, re.I)
    if m:
        return (0, int(m.group(1)), "")
    m = re.match(r"^(?:chr)?([XYM])$", chrom, re.I)
    if m:
        return (1, {"X": 0, "Y": 1, "M": 2}[m.group(1).upper()], "")
    return (2, 0, chrom)


def median(values):
    xs = sorted(values)
    n = len(xs)
    if n == 0:
        return None
    return xs[n // 2] if n % 2 else (xs[n // 2 - 1] + xs[n // 2]) / 2.0


def mad(values, med):
    return median([abs(x - med) for x in values])


def fmt_num(x, digits=4):
    if x is None:
        return ""
    return "%.*g" % (digits, x)


def deep_merge(base, override):
    """Recursively merge override into base (override wins on scalars)."""
    out = dict(base)
    for key, val in override.items():
        if key in out and isinstance(out[key], dict) and isinstance(val, dict):
            out[key] = deep_merge(out[key], val)
        else:
            out[key] = val
    return out


def load_config(config_path):
    with open(DEFAULT_CONFIG_PATH, "rt") as fh:
        cfg = json.load(fh)
    if config_path:
        with open(config_path, "rt") as fh:
            cfg = deep_merge(cfg, json.load(fh))
    return cfg


# ------------------------------------------------------------------ inputs

class Sample(object):
    __slots__ = ("sid", "truth_vcf", "normal_bam", "tumor_bam")

    def __init__(self, sid, truth_vcf, normal_bam=None, tumor_bam=None):
        self.sid = sid
        self.truth_vcf = truth_vcf
        self.normal_bam = normal_bam or None
        self.tumor_bam = tumor_bam or None


def _detect_columns(cells):
    mapping = {}
    for i, header in enumerate(cells):
        header = header.strip().lstrip("#").strip().lower()
        for canon, names in COLUMN_ALIASES.items():
            if header in names and canon not in mapping:
                mapping[canon] = i
    if all(k in mapping for k in REQUIRED_COLUMNS):
        return mapping
    return None


def _sid_from_path(path):
    base = os.path.basename(path)
    for suffix in (".vcf.gz", ".vcf", ".bcf"):
        if base.endswith(suffix):
            base = base[: -len(suffix)]
            break
    return sanitize_sid(base)


def load_samples(args):
    """Load the sample manifest / list / single-sample input."""
    samples = []
    if args.samples:
        with open(args.samples, "rt") as fh:
            lines = [ln.rstrip("\n") for ln in fh if ln.strip()]
        if not lines:
            sys.exit("[ERROR] empty sample manifest: %s" % args.samples)
        colmap, data_lines = None, lines
        if lines[0].lstrip().startswith("#"):
            colmap = _detect_columns(lines[0].split("\t"))
            if colmap is None:
                sys.exit("[ERROR] unrecognized manifest header: %r" % lines[0])
            data_lines = lines[1:]
        for ln in data_lines:
            cells = ln.split("\t")
            if colmap is not None:
                def cell(key):
                    idx = colmap.get(key)
                    if idx is None or idx >= len(cells):
                        return ""
                    return cells[idx].strip()
                sid, vcf = cell("sample_id"), cell("truth_vcf")
                nbam, tbam = cell("normal_bam"), cell("tumor_bam")
            elif len(cells) >= 4:  # ClairS 4-col convention: id/normal/tumor/truth
                sid, vcf = cells[0].strip(), cells[3].strip()
                nbam, tbam = cells[1].strip(), cells[2].strip()
            elif len(cells) == 2:
                sid, vcf = cells[0].strip(), cells[1].strip()
                nbam = tbam = ""
            else:
                sys.exit("[ERROR] unparseable manifest row: %r" % ln)
            if not sid or not vcf:
                sys.exit("[ERROR] row missing sample_id or truth_vcf: %r" % ln)
            samples.append(Sample(sid, vcf, nbam, tbam))
    elif args.vcf_list:
        with open(args.vcf_list, "rt") as fh:
            for ln in fh:
                path = ln.strip()
                if path:
                    samples.append(Sample(_sid_from_path(path), path))
    elif args.truth_vcf:
        sid = args.sample_id or _sid_from_path(args.truth_vcf)
        samples.append(Sample(sid, args.truth_vcf,
                              args.normal_bam, args.tumor_bam))
    else:
        sys.exit("[ERROR] one of --samples / --vcf-list / --truth-vcf is required")

    seen = set()
    for s in samples:
        if s.sid in seen:
            sys.exit("[ERROR] duplicate sample id: %s" % s.sid)
        seen.add(s.sid)
        if not os.path.isfile(s.truth_vcf):
            sys.exit("[ERROR] truth_vcf not found: %s (%s)" % (s.truth_vcf, s.sid))
    return samples


# ------------------------------------------------------------------ variant rules

def get_gnomad_af(info, fields):
    """Max AF across the probed population-frequency fields (or None)."""
    vals = []
    for field in fields:
        if field in info:
            vals.extend(parse_floats(info[field]))
    return max(vals) if vals else None


def has_cosmic_hit(info, fields):
    return any(truthy(info.get(f)) for f in fields)


def dna_categories(filters_category, dna_prefix):
    """'DNA_strelka:Somatic|RNA_mutect2:Artifact' -> ['SOMATIC', ...] (DNA only)."""
    cats = []
    for item in str(filters_category).split("|"):
        if ":" in item:
            caller, cat = item.split(":", 1)
            if caller.strip().upper().startswith(dna_prefix.upper()):
                cats.append(cat.strip().upper())
    return cats


def has_dna_entry(info, dna_prefix):
    """True when any DNA_* caller actually called this record.

    Only per-record detection fields count: CALLERS lists every aggregated
    caller sample-wide and is deliberately ignored.
    """
    prefix = dna_prefix.upper()
    for field in ("FILTERS_CATEGORY", "FILTERS_NORMALIZED", "FILTERS_ORIGINAL",
                  "CALLERS_SUPPORT"):
        value = info.get(field)
        if not value:
            continue
        for item in str(value).split("|"):
            caller = item.split(":", 1)[0].strip().upper()
            if caller.startswith(prefix):
                return True
    return False


def evaluate_record(info, conflict, cluster_count, cfg):
    """Evaluate variant-level rules for one target record.

    Returns None when no rule fires, else dict with rules/confidence/action/af.
    """
    vcfg = cfg["variant_rules"]
    target = vcfg["target_filter"].upper()

    af = get_gnomad_af(info, vcfg["gnomad_af_fields"])
    r1_cfg = vcfg["R1"]
    r1 = (r1_cfg["enabled"] and af is not None and af >= r1_cfg["af_common"])

    r2 = False
    if vcfg["R2"]["enabled"] and "FILTERS_CATEGORY" in info:
        r2 = target not in dna_categories(info["FILTERS_CATEGORY"],
                                          vcfg["dna_prefix"])

    r3 = (vcfg["R3"]["enabled"] and r1
          and has_cosmic_hit(info, vcfg["cosmic_id_fields"]))

    r4 = vcfg["R4"]["enabled"] and conflict

    r6 = False
    r6_cfg = vcfg["R6"]
    if r6_cfg["enabled"]:
        unified = str(info.get("UNIFIED_FILTER", "")).strip()
        if unified and unified != ".":
            contra = {c.upper() for c in r6_cfg["contradict_unified"]}
            r6 = unified.upper() in contra

    r7 = False
    r7_cfg = vcfg["R7"]
    if r7_cfg["enabled"]:
        evidence = str(info.get("REDI_EVIDENCE", "")).strip().upper()
        if evidence in {e.upper() for e in r7_cfg["evidence_levels"]}:
            r7 = True
        elif truthy(info.get("REDI_ACCESSION")):
            canonical = str(info.get("REDI_CANONICAL", "")).strip().upper()
            r7 = (not r7_cfg["require_canonical"]) or canonical == "YES"

    r8_cfg = vcfg["R8"]
    r8 = (r8_cfg["enabled"] and cluster_count >= r8_cfg["min_records"])

    r9 = False
    r9_cfg = vcfg["R9"]
    if r9_cfg["enabled"] and "N_SUPPORT_CALLERS" in info:
        vals = parse_floats(info["N_SUPPORT_CALLERS"])
        if vals:
            r9 = vals[0] <= r9_cfg["max_support_callers"]

    fired = {"R1": r1, "R2": r2, "R3": r3, "R4": r4, "R6": r6,
             "R7": r7, "R8": r8, "R9": r9}
    rules = [name for name in ("R1", "R2", "R3", "R4", "R6", "R7", "R8", "R9")
             if fired[name]]
    if not rules:
        return None

    n_strong = sum(1 for name in STRONG_RULES if fired[name])
    if fired["R4"] or fired["R6"] or n_strong >= 2:
        confidence, action = "high", "DROP"
    elif ((fired["R1"] and af >= r1_cfg["af_strong"])
          or (fired["R7"] and n_strong == 1)):
        confidence = "mid"
        action = "RELABEL_" + vcfg["mid_relabel_to"].upper()
    else:
        confidence, action = "low", "REPORT_ONLY"

    return {"rules": rules, "confidence": confidence, "action": action, "af": af}


# ------------------------------------------------------------------ sample scan

def _variant_class(ref, alt):
    if len(ref) == 1 and len(alt) == 1:
        return "ti" if (ref + alt).upper() in TRANSITIONS else "tv"
    return "indel"


def scan_sample(task):
    """Stream one truth VCF; evaluate rules; collect metrics and flagged entries."""
    sample, cfg = task
    vcfg = cfg["variant_rules"]
    target = vcfg["target_filter"]
    conflict_filters = {c.upper() for c in vcfg["conflict_filters"]}

    res = {
        "sid": sample.sid, "error": None, "notes": [], "tier_b": None,
        "n_records": 0, "n_somatic": 0,
        "dna_records": 0, "rna_only_somatic": 0, "rna_only_common": 0,
        "n_ti": 0, "n_tv": 0, "n_indel": 0, "n_snv": 0,
        "vafs": [], "dps": [], "n_het_vaf": 0, "n_low_dp": 0,
        "caller_support": defaultdict(int),
        "rule_hits": defaultdict(int),
        "conf_counts": defaultdict(int),
        "flagged": [], "entries": [],
    }

    try:
        somatic_records = []   # (parts, info) for target records, in file order
        conflict_keys = set()  # (chrom, pos, ref, alt) with Germline/Reference FILTER

        with open_text(sample.truth_vcf) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 8:
                    continue
                chrom, pos_s, ref, alt = parts[0], parts[1], parts[3], parts[4]
                try:
                    int(pos_s)
                except ValueError:
                    continue
                filt = parts[6]
                filt_set = set(filt.split(";"))
                filt_set_upper = {f.upper() for f in filt_set}
                info = parse_info(parts[7]) if parts[7] != "." else {}
                res["n_records"] += 1

                if has_dna_entry(info, vcfg["dna_prefix"]):
                    res["dna_records"] += 1
                if filt_set_upper & conflict_filters:
                    conflict_keys.add((chrom, pos_s, ref, alt))

                if target not in filt_set:
                    continue
                res["n_somatic"] += 1
                somatic_records.append((parts, info))

                # spectrum metrics
                cls = _variant_class(ref, alt)
                if cls == "ti":
                    res["n_ti"] += 1
                    res["n_snv"] += 1
                elif cls == "tv":
                    res["n_tv"] += 1
                    res["n_snv"] += 1
                else:
                    res["n_indel"] += 1

                vaf_vals = parse_floats(info.get("VAF_MEAN", ""))
                if vaf_vals:
                    res["vafs"].append(vaf_vals[0])
                dp = None
                for field in ("DP_DNA_MEAN", "DP_MEAN", "DP_MIN"):
                    vals = parse_floats(info.get(field, ""))
                    if vals:
                        dp = vals[0]
                        break
                if dp is not None:
                    res["dps"].append(dp)

                if not has_dna_entry(info, vcfg["dna_prefix"]):
                    res["rna_only_somatic"] += 1
                    af_vals = get_gnomad_af(info, vcfg["gnomad_af_fields"])
                    if (af_vals is not None
                            and af_vals >= vcfg["R1"]["af_common"]):
                        res["rna_only_common"] += 1

                support = info.get("CALLERS_SUPPORT") or info.get("CALLERS") or ""
                for item in str(support).split("|"):
                    caller = item.strip()
                    if caller:
                        res["caller_support"][caller] += 1

        # ---- R8: clustered Somatic records (sliding window per chromosome)
        positions = defaultdict(list)  # chrom -> [(pos, idx)]
        for idx, (parts, _info) in enumerate(somatic_records):
            positions[parts[0]].append((int(parts[1]), idx))
        cluster_of = [1] * len(somatic_records)
        window = vcfg["R8"]["window_bp"]
        for plist in positions.values():
            plist.sort()
            coords = [p for p, _ in plist]
            for i, (pos, idx) in enumerate(plist):
                lo = bisect_left(coords, pos - window)
                hi = bisect_right(coords, pos + window)
                cluster_of[idx] = hi - lo

        # ---- evaluate rules
        for idx, (parts, info) in enumerate(somatic_records):
            key = (parts[0], parts[1], parts[3], parts[4])
            conflict = key in conflict_keys
            ev = evaluate_record(info, conflict, cluster_of[idx], cfg)
            if ev is None:
                continue
            for name in ev["rules"]:
                res["rule_hits"][name] += 1
            res["conf_counts"][ev["confidence"]] += 1
            dp = None
            for field in ("DP_DNA_MEAN", "DP_MEAN", "DP_MIN"):
                vals = parse_floats(info.get(field, ""))
                if vals:
                    dp = vals[0]
                    break
            res["flagged"].append([
                sample.sid, parts[0], parts[1], parts[3], parts[4], parts[6],
                "+".join(ev["rules"]), ev["confidence"], ev["action"],
                fmt_num(ev["af"]),
                str(info.get("COSMIC_ID", "")),
                str(info.get("UNIFIED_FILTER", "")),
                str(info.get("FILTERS_CATEGORY", "")),
                str(info.get("VAF_MEAN", "")),
                fmt_num(dp), "",
            ])
            res["entries"].append({"parts": parts, "ev": ev})

        # het-like VAF / low-DP fractions need config thresholds
        s6 = cfg["sample_gates"]["S6_spectrum"]
        res["n_het_vaf"] = sum(1 for v in res["vafs"]
                               if s6["het_vaf_lo"] <= v <= s6["het_vaf_hi"])
        s7 = cfg["sample_gates"]["S7_coverage"]
        res["n_low_dp"] = sum(1 for d in res["dps"] if d < s7["dp_floor"])

        # ---- Tier B: optional BAM-based verification (samtools subprocess)
        if cfg.get("tier_b", {}).get("enabled"):
            res["tier_b"] = run_tier_b(sample, somatic_records, cfg, res)

    except Exception as exc:  # noqa: BLE001
        res["error"] = "%s: %s" % (type(exc).__name__, exc)

    res["caller_support"] = dict(res["caller_support"])
    res["rule_hits"] = dict(res["rule_hits"])
    res["conf_counts"] = dict(res["conf_counts"])
    return res


# ------------------------------------------------------------------ tier B (BAM)

def scan_mpileup_bases(bstr):
    """Parse an mpileup base string (run without -f: bases are actual letters,
    uppercase = forward strand, lowercase = reverse).

    Returns (base_chars, indel_counts{('+'|'-', len): count}).
    Ported from TruthQC's scan_mpileup_bases.
    """
    bases, indels = [], defaultdict(int)
    i, n = 0, len(bstr)
    while i < n:
        c = bstr[i]
        if c == "^":                      # ^MQ read start
            i += 2
            continue
        if c == "$":                      # read end
            i += 1
            continue
        if c in "+-":                     # indel: [+|-]<len><seq>
            i += 1
            num = ""
            while i < n and bstr[i].isdigit():
                num += bstr[i]
                i += 1
            ln = int(num) if num else 0
            i += ln
            if ln:
                indels[(c, ln)] += 1
            continue
        bases.append(c)
        i += 1
    return bases, indels


def alt_vaf_from_pileup(bases, indels, ref, alt, depth):
    """SNV exact alt-base match; indel approximated by +/-length. Denominator
    is the mpileup depth when available. Ported from TruthQC."""
    covered = sum(1 for b in bases if b.upper() in "ACGT.,*#")
    den = depth if depth > 0 else covered
    if den <= 0:
        return None
    if len(ref) == 1 and len(alt) == 1:
        au = alt.upper()
        cnt = sum(1 for b in bases if b.upper() == au)
    elif len(alt) > len(ref):
        cnt = indels.get(("+", len(alt) - len(ref)), 0)
    elif len(ref) > len(alt):
        cnt = indels.get(("-", len(ref) - len(alt)), 0)
    else:
        return None
    return cnt / float(den)


def run_mpileup(samtools, bam, sites_path, tb_cfg):
    """Single batched `samtools mpileup -l sites` pass over one BAM.

    Returns {(chrom, pos): (depth, base_string)}. Raises RuntimeError when
    samtools exits non-zero or times out.
    """
    cmd = [samtools, "mpileup", "-l", sites_path,
           "-q", str(tb_cfg["min_mapq"]), "-Q", str(tb_cfg["min_baseq"]), bam]
    try:
        proc = subprocess.run(cmd, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE, universal_newlines=True,
                              timeout=tb_cfg["timeout_sec"])
    except subprocess.TimeoutExpired:
        raise RuntimeError("mpileup timed out after %ds"
                           % tb_cfg["timeout_sec"])
    except OSError as exc:
        raise RuntimeError("mpileup failed to start: %s" % exc)
    if proc.returncode != 0:
        raise RuntimeError("mpileup exit %d: %s"
                           % (proc.returncode, proc.stderr.strip()[-200:]))
    out = {}
    for ln in proc.stdout.splitlines():
        f = ln.rstrip("\n").split("\t")
        if len(f) < 5:
            continue
        try:
            out[(f[0], f[1])] = (int(f[3]), f[4])
        except ValueError:
            continue
    return out


def _unique_sites(somatic_records):
    """Ordered unique (chrom, pos, ref, alt, filter) of the Somatic records."""
    seen, sites = set(), []
    for parts, _info in somatic_records:
        key = (parts[0], parts[1], parts[3], parts[4])
        if key not in seen:
            seen.add(key)
            sites.append((parts[0], parts[1], parts[3], parts[4], parts[6]))
    return sites


def _write_sites_file(sites):
    fh = tempfile.NamedTemporaryFile(mode="wt", prefix="label_qc_sites_",
                                     suffix=".tsv", delete=False)
    with fh:
        for chrom, pos, _ref, _alt, _filt in sites:
            fh.write("%s\t%s\n" % (chrom, pos))
    return fh.name


def measure_normal_contamination(pileup, sites, nc_cfg):
    """Fraction of truth sites with normal alt-VAF >= alt_vaf_min."""
    n_evaluable, n_contam = 0, 0
    for chrom, pos, ref, alt, _filt in sites:
        hit = pileup.get((chrom, pos))
        if hit is None:
            continue
        depth, bstr = hit
        if depth < nc_cfg["min_normal_dp"]:
            continue
        bases, indels = scan_mpileup_bases(bstr)
        vaf = alt_vaf_from_pileup(bases, indels, ref, alt, depth)
        if vaf is None:
            continue
        n_evaluable += 1
        if vaf >= nc_cfg["alt_vaf_min"]:
            n_contam += 1
    return {"n_sites": len(sites), "n_evaluable": n_evaluable,
            "n_contaminated": n_contam,
            "fraction": (n_contam / float(n_evaluable)
                         if n_evaluable else None)}


def measure_strand_bias(pileup, sites, sb_cfg):
    """Per-SNV-site strand counts of tumor alt reads (letter case = strand).

    A site is strand-biased when it has >= min_alt_reads alt reads and one
    strand carries >= max_strand_share of them.
    """
    n_evaluable, biased = 0, []
    for chrom, pos, ref, alt, filt in sites:
        if len(ref) != 1 or len(alt) != 1:
            continue  # indel strand bookkeeping is unreliable; SNVs only
        hit = pileup.get((chrom, pos))
        if hit is None:
            continue
        depth, bstr = hit
        if depth <= 0:
            continue
        bases, _indels = scan_mpileup_bases(bstr)
        au = alt.upper()
        fwd = sum(1 for b in bases if b == au)
        rev = sum(1 for b in bases if b == au.lower())
        tot = fwd + rev
        if tot < sb_cfg["min_alt_reads"]:
            continue
        n_evaluable += 1
        if max(fwd, rev) / float(tot) >= sb_cfg["max_strand_share"]:
            biased.append((chrom, pos, ref, alt, filt, fwd, rev))
    return {"n_sites": len(sites), "n_evaluable": n_evaluable,
            "n_strand_bias": len(biased),
            "fraction": (len(biased) / float(n_evaluable)
                         if n_evaluable else None),
            "biased_sites": biased}


def _flag_strand_bias(res, biased_sites):
    """Fold strand-bias sites into the flagged list as rule B1 (report-only)."""
    if not biased_sites:
        return
    idx = {}
    for i, row in enumerate(res["flagged"]):
        idx.setdefault((row[1], row[2], row[3], row[4]), i)
    for chrom, pos, ref, alt, filt, fwd, rev in biased_sites:
        note = "strand_bias: tumor alt reads %d fwd / %d rev" % (fwd, rev)
        key = (chrom, pos, ref, alt)
        if key in idx:
            row = res["flagged"][idx[key]]
            if "B1" not in row[6].split("+"):
                row[6] += "+B1"
            row[15] = (row[15] + "; " + note) if row[15] else note
        else:
            res["flagged"].append([
                res["sid"], chrom, pos, ref, alt, filt, "B1", "low",
                "REPORT_ONLY", "", "", "", "", "", "", note,
            ])
            res["conf_counts"]["low"] += 1
        res["rule_hits"]["B1"] += 1


def run_tier_b(sample, somatic_records, cfg, res):
    """Run the Tier B mpileup measurements for one sample.

    Returns {"normal": {...}|None, "tumor": {...}|None, "error": str|None};
    a measurement is None when its BAM was not provided, or {"error": ...}
    when the BAM was provided but could not be measured (never raises).
    """
    tb_cfg = cfg["tier_b"]
    out = {"normal": None, "tumor": None, "error": None}
    samtools = tb_cfg.get("samtools", "samtools")
    if shutil.which(samtools) is None and not os.path.isfile(samtools):
        out["error"] = "samtools not found: %s" % samtools
        return out

    sites = _unique_sites(somatic_records)
    if not sites:
        out["error"] = "no Somatic sites to verify"
        return out

    plans = []
    if tb_cfg.get("normal_contamination", {}).get("enabled", True) \
            and sample.normal_bam:
        plans.append(("normal", sample.normal_bam))
    if tb_cfg.get("strand_bias", {}).get("enabled", True) \
            and sample.tumor_bam:
        plans.append(("tumor", sample.tumor_bam))
    if not plans:
        return out

    sites_path = _write_sites_file(sites)
    try:
        for kind, bam in plans:
            if not os.path.isfile(bam):
                out[kind] = {"error": "BAM not found: %s" % bam}
                continue
            try:
                pileup = run_mpileup(samtools, bam, sites_path, tb_cfg)
            except RuntimeError as exc:
                out[kind] = {"error": str(exc)}
                continue
            if kind == "normal":
                out["normal"] = measure_normal_contamination(
                    pileup, sites, tb_cfg["normal_contamination"])
            else:
                tumor = measure_strand_bias(pileup, sites,
                                            tb_cfg["strand_bias"])
                _flag_strand_bias(res, tumor.pop("biased_sites"))
                out["tumor"] = tumor
    finally:
        try:
            os.remove(sites_path)
        except OSError:
            pass
    return out


# ------------------------------------------------------------------ sample gates

def _robust_z(x, values, min_mad):
    med = median(values)
    m = mad(values, med)
    denom = 1.4826 * max(m, min_mad)
    return (x - med) / denom if denom > 0 else 0.0


def _gate(gate, status, metric, detail):
    return {"gate": gate, "status": status, "metric": metric, "detail": detail}


def evaluate_gates(res, all_results, cfg):
    """Evaluate sample-level gates; returns a list of gate outcome dicts."""
    gates_cfg = cfg["sample_gates"]
    ok_results = [r for r in all_results if r["error"] is None]
    cohort_mode = len(ok_results) >= gates_cfg["min_cohort_for_outliers"]
    min_mad = gates_cfg["min_mad"]

    n_som = max(res["n_somatic"], 1)
    high = res["conf_counts"].get("high", 0)
    mid = res["conf_counts"].get("mid", 0)
    actioned_rate = (high + mid) / float(n_som)
    contradiction_rate = res["rule_hits"].get("R6", 0) / float(n_som)
    rna_only_fraction = res["rna_only_somatic"] / float(n_som)
    rna_only_common_fraction = res["rna_only_common"] / float(n_som)
    dna_participation = (res["dna_records"] / float(res["n_records"])
                         if res["n_records"] else 0.0)
    indel_fraction = (res["n_indel"] / float(n_som))
    ti_tv = (res["n_ti"] / float(res["n_tv"])) if res["n_tv"] else None
    median_vaf = median(res["vafs"])
    het_fraction = (res["n_het_vaf"] / float(len(res["vafs"]))
                    if res["vafs"] else None)
    median_dp = median(res["dps"])
    low_dp_fraction = (res["n_low_dp"] / float(len(res["dps"]))
                       if res["dps"] else None)
    total_support = sum(res["caller_support"].values())
    top_caller_share = (max(res["caller_support"].values()) / float(total_support)
                        if total_support else None)

    res["metrics"] = {
        "actioned_rate": actioned_rate,
        "contradiction_rate": contradiction_rate,
        "rna_only_fraction": rna_only_fraction,
        "rna_only_common_fraction": rna_only_common_fraction,
        "dna_participation": dna_participation,
        "indel_fraction": indel_fraction,
        "ti_tv": ti_tv,
        "median_vaf": median_vaf,
        "het_vaf_fraction": het_fraction,
        "median_dp": median_dp,
        "low_dp_fraction": low_dp_fraction,
        "top_caller_share": top_caller_share,
    }

    gates = []

    # S0: flag burden (port of TruthQC's abnormal-sample criterion)
    g = gates_cfg["S0_flag_burden"]
    if g["enabled"]:
        if res["n_somatic"] < g["min_somatic"]:
            gates.append(_gate("S0", "SKIP", actioned_rate,
                               "too few Somatic records (<%d)" % g["min_somatic"]))
        elif actioned_rate >= g["fail"]:
            gates.append(_gate("S0", "FAIL", actioned_rate,
                               "flag burden %.1f%% >= %.1f%%"
                               % (actioned_rate * 100, g["fail"] * 100)))
        elif actioned_rate >= g["warn"]:
            gates.append(_gate("S0", "WARN", actioned_rate,
                               "flag burden %.1f%% >= %.1f%%"
                               % (actioned_rate * 100, g["warn"] * 100)))
        else:
            gates.append(_gate("S0", "PASS", actioned_rate, ""))

    # S1: somatic-count outlier vs cohort, with absolute floors
    g = gates_cfg["S1_count_outlier"]
    if g["enabled"]:
        n = res["n_somatic"]
        z = (_robust_z(n, [r["n_somatic"] for r in ok_results], min_mad)
             if cohort_mode else 0.0)
        if n >= g["abs_fail"] or (cohort_mode and z >= g["z_fail"]):
            gates.append(_gate("S1", "FAIL", n,
                               "somatic count %d (z=%.1f, abs_fail=%d)"
                               % (n, z, g["abs_fail"])))
        elif n >= g["abs_warn"] or (cohort_mode and z >= g["z_warn"]):
            gates.append(_gate("S1", "WARN", n,
                               "somatic count %d (z=%.1f)" % (n, z)))
        else:
            gates.append(_gate("S1", "PASS", n, ""))

    # S2: FILTER/INFO self-contradiction rate
    g = gates_cfg["S2_contradiction"]
    if g["enabled"]:
        r = contradiction_rate
        if r >= g["fail"]:
            gates.append(_gate("S2", "FAIL", r,
                               "contradiction rate %.1f%% >= %.1f%%"
                               % (r * 100, g["fail"] * 100)))
        elif r >= g["warn"]:
            gates.append(_gate("S2", "WARN", r,
                               "contradiction rate %.1f%% >= %.1f%%"
                               % (r * 100, g["warn"] * 100)))
        else:
            gates.append(_gate("S2", "PASS", r, ""))

    # S3: RNA-only Somatic fraction at common population frequency
    # (germline leakage through the RNA branch). Plain RNA-only calls are the
    # norm for rescued records and are reported but not gated.
    g = gates_cfg["S3_rna_only"]
    if g["enabled"]:
        r = rna_only_common_fraction
        if r >= g["fail"]:
            gates.append(_gate("S3", "FAIL", r,
                               "RNA-only common-AF fraction %.1f%% >= %.1f%%"
                               % (r * 100, g["fail"] * 100)))
        elif r >= g["warn"]:
            gates.append(_gate("S3", "WARN", r,
                               "RNA-only common-AF fraction %.1f%% >= %.1f%%"
                               % (r * 100, g["warn"] * 100)))
        else:
            gates.append(_gate("S3", "PASS", r, ""))

    # S4: normal contamination — Tier B (BAM-based) verification
    g = gates_cfg["S4_normal_contamination"]
    tb = cfg.get("tier_b", {})
    if g.get("enabled") and tb.get("enabled"):
        tb_res = res.get("tier_b") or {}
        nc_cfg = tb.get("normal_contamination", {})
        nc = tb_res.get("normal")
        if tb_res.get("error"):
            gates.append(_gate("S4", "SKIP", None,
                               "Tier B error: %s" % tb_res["error"]))
        elif nc is None:
            gates.append(_gate("S4", "SKIP", None, "no normal BAM provided"))
        elif nc.get("error"):
            gates.append(_gate("S4", "SKIP", None,
                               "normal BAM not measurable: %s" % nc["error"]))
        elif nc["n_evaluable"] < nc_cfg["min_sites"]:
            gates.append(_gate("S4", "SKIP", nc["n_evaluable"],
                               "too few evaluable sites (%d < %d; zero coverage "
                               "may mean VCF/BAM contig naming mismatch)"
                               % (nc["n_evaluable"], nc_cfg["min_sites"])))
        else:
            frac = nc["fraction"] or 0.0
            detail = ("normal alt-VAF >= %.2f at %d/%d sites (%.1f%%)"
                      % (nc_cfg["alt_vaf_min"], nc["n_contaminated"],
                         nc["n_evaluable"], frac * 100))
            if frac >= nc_cfg["fail_fraction"]:
                gates.append(_gate("S4", "FAIL", frac, detail))
            elif frac >= nc_cfg["warn_fraction"]:
                gates.append(_gate("S4", "WARN", frac, detail))
            else:
                gates.append(_gate("S4", "PASS", frac, ""))
    elif g.get("enabled"):
        gates.append(_gate("S4", "SKIP", None,
                           "Tier B not enabled (use --verify-bam with BAMs)"))

    # S8: strand bias burden — optional WARN gate (Tier B; off by default)
    g = gates_cfg.get("S8_strand_bias", {})
    if g.get("enabled") and tb.get("enabled"):
        tm = (res.get("tier_b") or {}).get("tumor")
        if tm is None or tm.get("error"):
            gates.append(_gate("S8", "SKIP", None,
                               "no tumor BAM strand measurement"))
        elif tm["n_evaluable"] < g["min_sites"]:
            gates.append(_gate("S8", "SKIP", tm["n_evaluable"],
                               "too few evaluable SNV sites"))
        else:
            frac = tm["fraction"] or 0.0
            if frac >= g["warn"]:
                gates.append(_gate("S8", "WARN", frac,
                                   "strand-biased sites %d/%d (%.1f%%)"
                                   % (tm["n_strand_bias"], tm["n_evaluable"],
                                      frac * 100)))
            else:
                gates.append(_gate("S8", "PASS", frac, ""))

    # S5: modality completeness (DNA callers participated at all)
    g = gates_cfg["S5_modality"]
    if g["enabled"]:
        r = dna_participation
        z = (-_robust_z(r, [x["dna_records"] / float(max(x["n_records"], 1))
                            for x in ok_results], min_mad)
             if cohort_mode else 0.0)
        if r <= g["fail_floor"]:
            gates.append(_gate("S5", "FAIL", r,
                               "no DNA-caller participation in any record"))
        elif cohort_mode and z >= g["z_warn"]:
            gates.append(_gate("S5", "WARN", r,
                               "DNA participation %.3f is a cohort low outlier (z=%.1f)"
                               % (r, z)))
        else:
            gates.append(_gate("S5", "PASS", r, ""))

    # S6: spectrum sanity (WARN-level)
    g = gates_cfg["S6_spectrum"]
    if g["enabled"]:
        problems = []
        if res["n_snv"] >= g["min_snvs"] and ti_tv is not None:
            if not (g["ti_tv_min"] <= ti_tv <= g["ti_tv_max"]):
                problems.append("Ti/Tv=%.2f outside [%.2f, %.2f]"
                                % (ti_tv, g["ti_tv_min"], g["ti_tv_max"]))
        if indel_fraction > g["indel_frac_max"]:
            problems.append("indel fraction %.2f > %.2f"
                            % (indel_fraction, g["indel_frac_max"]))
        if het_fraction is not None and het_fraction > g["het_vaf_frac_max"]:
            problems.append("het-like VAF fraction %.2f > %.2f (germline-like)"
                            % (het_fraction, g["het_vaf_frac_max"]))
        if (top_caller_share is not None and res["n_somatic"] >= g["min_snvs"]
                and top_caller_share > g["caller_share_max"]):
            problems.append("single caller contributes %.1f%% of support"
                            % (top_caller_share * 100))
        status = "WARN" if problems else "PASS"
        gates.append(_gate("S6", status, None, "; ".join(problems)))

    # S7: coverage floors (WARN-level)
    g = gates_cfg["S7_coverage"]
    if g["enabled"]:
        if res["dps"] and len(res["dps"]) >= g["min_dp_records"]:
            if low_dp_fraction > g["low_dp_frac_max"]:
                gates.append(_gate("S7", "WARN", low_dp_fraction,
                                   "%.1f%% of Somatic records below DP floor %d"
                                   % (low_dp_fraction * 100, g["dp_floor"])))
            else:
                gates.append(_gate("S7", "PASS", low_dp_fraction, ""))
        else:
            gates.append(_gate("S7", "SKIP", None, "no usable DP fields in INFO"))

    failed = [x["gate"] for x in gates if x["status"] == "FAIL"]
    warned = [x["gate"] for x in gates if x["status"] == "WARN"]
    verdict = "FAIL" if failed else ("WARN" if warned else "PASS")
    res["gates"] = gates
    res["verdict"] = verdict
    res["gates_failed"] = failed
    res["gates_warned"] = warned
    return res


# ------------------------------------------------------------------ writers

def write_flagged_tsv(path, all_rows):
    rows = sorted(all_rows, key=lambda r: (natural_chrom_key(str(r[1])),
                                           int(r[2]), str(r[0]), r[3], r[4]))
    with gzip.open(path, "wt") as fh:
        fh.write("\t".join(FLAGGED_HEADER) + "\n")
        for row in rows:
            fh.write("\t".join(str(x) for x in row) + "\n")
    return len(rows)


def write_samples_qc(path, results, tier_b_active=False):
    header = SAMPLES_QC_HEADER + (TIER_B_COLUMNS if tier_b_active else ())
    with open(path, "wt") as fh:
        fh.write("\t".join(header) + "\n")
        for r in sorted(results, key=lambda x: x["sid"]):
            m = r.get("metrics", {})
            row = [
                r["sid"], r.get("verdict", "ERROR"),
                str(r["n_records"]), str(r["n_somatic"]),
                str(len(r["flagged"])),
                str(r["conf_counts"].get("high", 0)),
                str(r["conf_counts"].get("mid", 0)),
                str(r["conf_counts"].get("low", 0)),
                fmt_num(m.get("actioned_rate")),
                fmt_num(m.get("contradiction_rate")),
                fmt_num(m.get("rna_only_fraction")),
                fmt_num(m.get("rna_only_common_fraction")),
                fmt_num(m.get("dna_participation")),
                fmt_num(m.get("ti_tv")),
                fmt_num(m.get("indel_fraction")),
                fmt_num(m.get("median_vaf")),
                fmt_num(m.get("het_vaf_fraction")),
                fmt_num(m.get("median_dp")),
                fmt_num(m.get("low_dp_fraction")),
                ";".join(r.get("gates_failed", [])),
                ";".join(r.get("gates_warned", [])),
                ";".join(r["notes"]) or (r["error"] or ""),
            ]
            if tier_b_active:
                tb = r.get("tier_b") or {}
                nc = tb.get("normal") or {}
                tm = tb.get("tumor") or {}
                row += [fmt_num(nc.get("fraction")),
                        str(nc.get("n_evaluable", "")),
                        str(tm.get("n_strand_bias", ""))]
            fh.write("\t".join(row) + "\n")


def build_summary(cfg, results, run_id, params_fp):
    return {
        "run_id": run_id,
        "version": VERSION,
        "mode": "apply" if cfg["apply"] else "dry-run",
        "config_fingerprint": params_fp,
        "created": time.strftime("%Y-%m-%d %H:%M:%S"),
        "n_samples": len(results),
        "n_failed_samples": sum(1 for r in results if r["error"] is not None),
        "verdicts": {
            v: sum(1 for r in results if r.get("verdict") == v)
            for v in ("PASS", "WARN", "FAIL")
        },
        "rule_hits": _sum_counters(r["rule_hits"] for r in results),
        "samples": [_summary_sample(r)
                    for r in sorted(results, key=lambda x: x["sid"])],
    }


def _summary_sample(r):
    out = {
        "sample": r["sid"],
        "verdict": r.get("verdict", "ERROR"),
        "error": r["error"],
        "n_records": r["n_records"],
        "n_somatic": r["n_somatic"],
        "n_flagged": len(r["flagged"]),
        "conf_counts": r["conf_counts"],
        "rule_hits": r["rule_hits"],
        "metrics": r.get("metrics", {}),
        "gates": r.get("gates", []),
        "notes": r["notes"],
    }
    tb = r.get("tier_b")
    if tb and (tb.get("normal") or tb.get("tumor") or tb.get("error")):
        out["tier_b"] = tb
    return out


def _sum_counters(counters):
    out = defaultdict(int)
    for c in counters:
        for k, v in c.items():
            out[k] += v
    return dict(out)


def write_report(path, cfg, results, run_id, n_flagged, elapsed):
    ok = [r for r in results if r["error"] is None]
    bad = [r for r in results if r["error"] is not None]
    conf_all = _sum_counters(r["conf_counts"] for r in ok)
    rule_all = _sum_counters(r["rule_hits"] for r in ok)

    L = []
    tier_b_on = any(r.get("tier_b") for r in results)
    L.append("# LabelQC report (Tier A%s)\n"
             % (" + Tier B BAM verification" if tier_b_on
                else ", VCF-only"))
    L.append("- run id: `%s` | version: v%s" % (run_id, VERSION))
    L.append("- mode: **%s**%s"
             % ("APPLY" if cfg["apply"] else "dry-run",
                " — cleaned VCFs written (inputs untouched)" if cfg["apply"]
                else " — no files modified; use `--apply` to write cleaned VCFs"))
    L.append("- samples: %d ok / %d failed | elapsed: %.1fs"
             % (len(ok), len(bad), elapsed))
    L.append("")
    L.append("## Verdicts\n")
    L.append("| verdict | samples |")
    L.append("|---|---|")
    for v in ("PASS", "WARN", "FAIL"):
        n = sum(1 for r in ok if r.get("verdict") == v)
        L.append("| %s | %d |" % (v, n))
    L.append("")
    L.append("## Global flag statistics\n")
    L.append("| metric | value |")
    L.append("|---|---|")
    L.append("| scanned records | %s |"
             % format(sum(r["n_records"] for r in ok), ","))
    L.append("| target records (FILTER=%s) | %s |"
             % (cfg["variant_rules"]["target_filter"],
                format(sum(r["n_somatic"] for r in ok), ",")))
    L.append("| flagged total | %s |" % format(n_flagged, ","))
    L.append("| high (-> DROP) | %s |" % format(conf_all.get("high", 0), ","))
    L.append("| mid (-> RELABEL) | %s |" % format(conf_all.get("mid", 0), ","))
    L.append("| low (report only) | %s |" % format(conf_all.get("low", 0), ","))
    L.append("")
    L.append("### Rule hits\n")
    L.append("| rule | hits | meaning |")
    L.append("|---|---|---|")
    rule_desc = {
        "R1": "common gnomAD polymorphism (AF >= threshold)",
        "R2": "no DNA caller judged Somatic",
        "R3": "COSMIC hit with common gnomAD AF (hotspot contradiction)",
        "R4": "co-located Germline/Reference record at same site",
        "R6": "FILTER=Somatic but UNIFIED_FILTER in {Reference, Germline}",
        "R7": "RNA-editing (REDIportal) overlap",
        "R8": "clustered Somatic calls",
        "R9": "low caller agreement",
        "B1": "Tier B: strand-biased alt reads in the tumor BAM",
    }
    for name in sorted(rule_all):
        L.append("| %s | %s | %s |"
                 % (name, format(rule_all[name], ","), rule_desc.get(name, "")))
    L.append("")
    L.append("## Per-sample QC (FAIL first, then WARN, then PASS)\n")
    L.append("| sample | verdict | somatic | high | mid | low | actioned | "
             "contradiction | RNA-only common-AF | gates failed | gates warned |")
    L.append("|---|---|---|---|---|---|---|---|---|---|---|")
    order = {"FAIL": 0, "WARN": 1, "PASS": 2, "ERROR": 3}
    for r in sorted(ok, key=lambda x: (order.get(x.get("verdict"), 9),
                                       -x["metrics"]["actioned_rate"])):
        m = r["metrics"]
        L.append("| %s | %s | %d | %d | %d | %d | %.1f%% | %.1f%% | %.1f%% | %s | %s |"
                 % (r["sid"], r["verdict"], r["n_somatic"],
                    r["conf_counts"].get("high", 0),
                    r["conf_counts"].get("mid", 0),
                    r["conf_counts"].get("low", 0),
                    m["actioned_rate"] * 100,
                    m["contradiction_rate"] * 100,
                    m["rna_only_common_fraction"] * 100,
                    ";".join(r["gates_failed"]) or "-",
                    ";".join(r["gates_warned"]) or "-"))
    L.append("")
    tb_rows = [r for r in ok
               if (r.get("tier_b") or {}).get("normal")
               or (r.get("tier_b") or {}).get("tumor")]
    if tb_rows:
        nc_cfg = cfg.get("tier_b", {}).get("normal_contamination", {})
        L.append("## Tier B (BAM verification)\n")
        L.append("| sample | normal alt-VAF>=%.2f sites | evaluable | "
                 "contam fraction | strand-bias sites |"
                 % nc_cfg.get("alt_vaf_min", 0.05))
        L.append("|---|---|---|---|---|")
        for r in sorted(tb_rows, key=lambda x: x["sid"]):
            tb = r["tier_b"]
            nc, tm = tb.get("normal") or {}, tb.get("tumor") or {}
            L.append("| %s | %s | %s | %s | %s |"
                     % (r["sid"],
                        nc.get("n_contaminated", "-"),
                        nc.get("n_evaluable", "-"),
                        ("%.1f%%" % (nc["fraction"] * 100))
                        if nc.get("fraction") is not None else "-",
                        tm.get("n_strand_bias", "-")))
        L.append("")
    L.append("## Notes / warnings\n")
    any_note = False
    for r in ok:
        for gate in r.get("gates", []):
            if gate["status"] in ("WARN", "FAIL", "SKIP") and gate["detail"]:
                L.append("- [%s] %s %s: %s"
                         % (r["sid"], gate["gate"], gate["status"], gate["detail"]))
                any_note = True
    for r in bad:
        L.append("- **FAILED** [%s]: %s" % (r["sid"], r["error"]))
        any_note = True
    if not any_note:
        L.append("- none")
    L.append("")
    L.append("## Output files\n")
    L.append("```")
    L.append("report.md              this file")
    L.append("summary.json           machine-readable run summary")
    L.append("samples_qc.tsv         per-sample verdicts and metrics")
    L.append("flagged_sites.tsv.gz   per-site flagged detail (%d rows)" % n_flagged)
    if cfg["apply"]:
        L.append("cleaned_vcf/           cleaned VCFs (HIGH dropped, MID relabelled)")
        L.append("samples_cleaned.tsv    manifest pointing at cleaned VCFs")
    L.append("```")
    with open(path, "wt") as fh:
        fh.write("\n".join(L) + "\n")


# ------------------------------------------------------------------ apply mode

def write_cleaned_vcf(sample, res, cfg):
    """Second streaming pass: DROP high, RELABEL mid, pass everything else."""
    out_dir = cfg["out_cleaned"]
    os.makedirs(out_dir, exist_ok=True)
    final_path = os.path.join(out_dir, sanitize_sid(sample.sid) + ".vcf.gz")
    tmp_path = "%s.tmp%d" % (final_path, os.getpid())

    fmap = defaultdict(list)
    for e in res["entries"]:
        fmap[tuple(e["parts"])].append(e)

    relabel_to = cfg["variant_rules"]["mid_relabel_to"]
    with open_text(sample.truth_vcf) as src, gzip.open(tmp_path, "wt") as dst:
        for line in src:
            if line.startswith("#"):
                dst.write(line)
                continue
            parts = line.rstrip("\n").split("\t")
            queue = fmap.get(tuple(parts)) if len(parts) >= 8 else None
            if queue:
                e = queue.pop(0)
                if not queue:
                    fmap.pop(tuple(parts), None)
                ev = e["ev"]
                if ev["action"] == "DROP":
                    continue
                if ev["action"].startswith("RELABEL_"):
                    new_parts = list(parts)
                    new_parts[6] = relabel_to
                    tag = "LABEL_QC=%s[%s,%s]" % (ev["action"],
                                                  "+".join(ev["rules"]),
                                                  ev["confidence"])
                    new_parts[7] = (new_parts[7] + ";" + tag
                                    if new_parts[7] != "." else tag)
                    dst.write("\t".join(new_parts) + "\n")
                    continue
            dst.write(line)
    os.replace(tmp_path, final_path)
    return final_path


def write_samples_cleaned(path, samples, cleaned_paths):
    with open(path, "wt") as fh:
        fh.write("#sample_id\ttruth_vcf\n")
        for s in samples:
            fh.write("%s\t%s\n" % (s.sid, cleaned_paths.get(s.sid, s.truth_vcf)))


# ------------------------------------------------------------------ main

def parse_args(argv=None):
    ap = argparse.ArgumentParser(
        prog="label_qc.py",
        description="Tier A (VCF-only) label QC gate for consensus/rescue truth VCFs",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""examples:
  label_qc.py --samples samples.tsv --out label_qc_out
  label_qc.py --samples samples.tsv --out label_qc_out --apply
  label_qc.py --truth-vcf sample.vcf.gz --sample-id S1 --out label_qc_out
""")
    g_in = ap.add_argument_group("input (choose one)")
    g_in.add_argument("--samples",
                      help="sample manifest TSV (header aliases auto-detected; "
                           "or headerless 2-col id/vcf or 4-col ClairS convention)")
    g_in.add_argument("--vcf-list", help="file with one truth VCF path per line")
    g_in.add_argument("--truth-vcf", help="single-sample mode: truth VCF")
    g_in.add_argument("--sample-id", help="single-sample mode: sample id override")

    g_out = ap.add_argument_group("output and actions")
    g_out.add_argument("--out", default="label_qc_out",
                       help="output directory (default: %(default)s)")
    g_out.add_argument("--apply", action="store_true",
                       help="write cleaned VCFs (HIGH dropped, MID relabelled); "
                            "default is dry-run — inputs are never modified")

    g_cfg = ap.add_argument_group("configuration")
    g_cfg.add_argument("--config",
                       help="JSON config overriding the bundled label_qc_config.json")
    g_cfg.add_argument("--af-common", type=float,
                       help="override R1 common-polymorphism AF threshold")
    g_cfg.add_argument("--af-strong", type=float,
                       help="override R1 strong-AF (MID) threshold")
    g_cfg.add_argument("--relabel-to", choices=["Germline", "Reference"],
                       help="override MID relabel target class")

    g_tb = ap.add_argument_group("tier B (optional BAM verification)")
    g_tb.add_argument("--verify-bam", action="store_true",
                      help="enable Tier B: batched samtools mpileup checks "
                           "(S4 normal contamination, B1 strand bias)")
    g_tb.add_argument("--normal-bam",
                      help="single-sample mode: DNA normal BAM (for S4)")
    g_tb.add_argument("--tumor-bam",
                      help="single-sample mode: DNA tumor BAM (for B1)")
    g_tb.add_argument("--samtools",
                      help="samtools executable (default: tier_b.samtools "
                           "config value, else PATH)")

    g_run = ap.add_argument_group("run control")
    g_run.add_argument("--threads", type=int, default=8,
                       help="parallel samples (default: %(default)s)")
    g_run.add_argument("--version", action="version",
                       version="label_qc v" + VERSION)
    return ap.parse_args(argv)


def main(argv=None):
    t0 = time.time()
    args = parse_args(argv)
    if sum([bool(args.samples), bool(args.vcf_list), bool(args.truth_vcf)]) > 1:
        sys.exit("[ERROR] --samples / --vcf-list / --truth-vcf are mutually exclusive")

    cfg = load_config(args.config)
    if args.af_common is not None:
        cfg["variant_rules"]["R1"]["af_common"] = args.af_common
    if args.af_strong is not None:
        cfg["variant_rules"]["R1"]["af_strong"] = args.af_strong
    if args.relabel_to:
        cfg["variant_rules"]["mid_relabel_to"] = args.relabel_to
    if args.verify_bam:
        cfg.setdefault("tier_b", {})["enabled"] = True
        # --verify-bam implies the S4 gate (off by default in Tier A config)
        cfg["sample_gates"]["S4_normal_contamination"]["enabled"] = True
    if args.samtools:
        cfg.setdefault("tier_b", {})["samtools"] = args.samtools

    out_abs = os.path.abspath(args.out)
    cfg["apply"] = bool(args.apply)
    cfg["out"] = out_abs
    cfg["out_cleaned"] = os.path.join(out_abs, "cleaned_vcf")

    samples = load_samples(args)

    params_fp = hashlib.sha1(
        json.dumps(cfg["variant_rules"], sort_keys=True).encode()
        + json.dumps(cfg["sample_gates"], sort_keys=True).encode()
        + json.dumps(cfg.get("tier_b", {}), sort_keys=True).encode()
    ).hexdigest()[:10]
    run_id = "%s_%s" % (time.strftime("%Y%m%d_%H%M%S"), params_fp)

    tasks = [(s, cfg) for s in samples]
    workers = min(max(1, args.threads), len(tasks))
    if workers > 1:
        with ProcessPoolExecutor(max_workers=workers) as ex:
            results = list(ex.map(scan_sample, tasks))
    else:
        results = [scan_sample(t) for t in tasks]

    ok = [r for r in results if r["error"] is None]
    if not ok:
        for r in results:
            sys.stderr.write("[FAIL] %s: %s\n" % (r["sid"], r["error"]))
        sys.exit(2)

    for r in ok:
        evaluate_gates(r, ok, cfg)

    os.makedirs(out_abs, exist_ok=True)
    all_flagged = [row for r in ok for row in r["flagged"]]
    n_flagged = write_flagged_tsv(os.path.join(out_abs, "flagged_sites.tsv.gz"),
                                  all_flagged)
    tier_b_active = any((r.get("tier_b") or {}).get("normal")
                        or (r.get("tier_b") or {}).get("tumor")
                        for r in results)
    write_samples_qc(os.path.join(out_abs, "samples_qc.tsv"), results,
                     tier_b_active)

    summary = build_summary(cfg, results, run_id, params_fp)
    with open(os.path.join(out_abs, "summary.json"), "wt") as fh:
        json.dump(summary, fh, indent=2)

    if cfg["apply"]:
        cleaned_paths = {}
        for sample, r in zip(samples, results):
            if r["error"] is None:
                cleaned_paths[sample.sid] = write_cleaned_vcf(sample, r, cfg)
        write_samples_cleaned(os.path.join(out_abs, "samples_cleaned.tsv"),
                              samples, cleaned_paths)

    write_report(os.path.join(out_abs, "report.md"), cfg, results, run_id,
                 n_flagged, time.time() - t0)

    verdicts = summary["verdicts"]
    print("=" * 62)
    print("label_qc v%s done | run_id=%s" % (VERSION, run_id))
    print("  samples: %d ok / %d failed | verdicts: PASS=%d WARN=%d FAIL=%d"
          % (len(ok), len(results) - len(ok),
             verdicts["PASS"], verdicts["WARN"], verdicts["FAIL"]))
    print("  somatic records: %s | flagged: %s"
          % (format(sum(r["n_somatic"] for r in ok), ","),
             format(n_flagged, ",")))
    if cfg["apply"]:
        print("  cleaned VCFs: %s" % cfg["out_cleaned"])
    print("  report: %s" % os.path.join(out_abs, "report.md"))
    print("=" * 62)
    return 0


if __name__ == "__main__":
    sys.exit(main())
