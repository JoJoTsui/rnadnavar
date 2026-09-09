"""
Variant aggregation functions for consensus and rescue workflows.

This module provides reusable functions for aggregating variants
from multiple callers and across modalities. It handles reading VCF files,
extracting genotype information, and combining variant data with support
for both within-modality consensus and cross-modality rescue operations.

Functions:
    extract_genotype_info: Extract genotype, depth, and VAF from a variant
    aggregate_genotypes: Aggregate genotype information across callers
    read_variants_from_vcf: Read variants from a single VCF file
    aggregate_variants: Aggregate variants from multiple collections

Example:
    Reading and aggregating variants from multiple callers:

    >>> from vcf_utils.aggregation import read_variants_from_vcf, aggregate_variants
    >>>
    >>> # Read variants from each caller
    >>> mutect2_vars = read_variants_from_vcf('mutect2.vcf.gz', 'mutect2')
    >>> strelka_vars = read_variants_from_vcf('strelka.vcf.gz', 'strelka')
    >>> deepsomatic_vars = read_variants_from_vcf('deepsomatic.vcf.gz', 'deepsomatic',
    ...                                           exclude_refcall=True, exclude_germline=True)
    >>>
    >>> # Prepare collections for aggregation
    >>> collections = [
    ...     ('mutect2', mutect2_vars, None),
    ...     ('strelka', strelka_vars, None),
    ...     ('deepsomatic', deepsomatic_vars, None)
    ... ]
    >>>
    >>> # Aggregate with consensus thresholds
    >>> aggregated = aggregate_variants(collections, snv_threshold=2, indel_threshold=2)
    >>>
    >>> # Access aggregated data
    >>> for vkey, data in aggregated.items():
    ...     print(f"Variant: {vkey}")
    ...     print(f"  Supported by: {data['callers']}")
    ...     print(f"  Passes consensus: {data['passes_consensus']}")
    ...     print(f"  Mean VAF: {data['gt_aggregated']['vaf_mean']}")
"""


def sanitize_genotype(info):
    """Normalize unavailable/invalid measurements without turning them into zero.

    Keep invalid-field provenance; valid zero depth and zero AF remain valid.
    This is also used at the public aggregation boundary for in-memory callers.
    """
    import math

    result = dict(info or {})
    invalid = set(result.get("INVALID_FIELDS", []))
    for key in ("DP", "GQ", "VAF"):
        value = result.get(key)
        if value is None or (isinstance(value, str) and value in ("", ".")):
            result[key] = None
            continue
        try:
            number = float(value)
            valid = math.isfinite(number) and number >= 0
            valid = valid and (number <= 1 if key == "VAF" else number.is_integer())
        except (ValueError, TypeError, OverflowError):
            valid = False
        if valid:
            result[key] = number if key == "VAF" else int(number)
        else:
            result[key] = None
            invalid.add(key)
    ad = result.get("AD")
    if ad is None or (isinstance(ad, str) and ad in ("", ".")):
        result["AD"] = None
    else:
        try:
            values = str(ad).split(",") if isinstance(ad, str) else list(ad)
            counts = [int(v) for v in values]
            if len(counts) < 2 or any(c < 0 or float(v) != c for c, v in zip(counts, values)):
                raise ValueError("invalid allele depths")
            result["AD"] = ",".join(map(str, counts))
        except (ValueError, TypeError, OverflowError):
            result["AD"] = None
            invalid.add("AD")
    if result.get("VAF") is None:
        result["VAF_SOURCE"] = None
    result["INVALID_FIELDS"] = sorted(invalid)
    return result


def resolve_tumor_sample_index(samples, caller, normal_sample=None):
    """
    Resolve the tumor sample index in a caller VCF.

    Genotype metrics (GT/AD/DP/VAF) must be extracted from the tumor sample,
    never blindly from sample index 0, which is the NORMAL for Mutect2 and
    Strelka in this pipeline.

    Resolution order:
        1. Single-sample (tumor-only) VCF: the only sample is used
        2. Mutect2 ``##normal_sample=<name>`` header line (ground truth
           written by the caller): the tumor is the other sample
        3. Sample named *tumor*/*tumour* (Strelka 'TUMOR',
           DeepSomatic '<id>_tumor'): matched by name
        4. Sample named *normal* (Mutect2 with a named normal): the tumor is
           the other sample
        5. Caller ordering conventions: Strelka uses a fixed NORMAL,TUMOR
           order and Mutect2 is invoked with the normal CRAM first
           (subworkflows/local/bam_variant_calling/main.nf), so the tumor is
           the last sample; DeepSomatic lists the tumor first
        6. Fallback: index 0

    Args:
        samples (list): Sample names from the VCF header (vcf.samples)
        caller (str): Caller name (e.g., 'mutect2', 'strelka', 'deepsomatic')
        normal_sample (str, optional): Value of the VCF header's
            ``##normal_sample`` line, if present (Mutect2 writes this for
            paired runs). Default: None

    Returns:
        int: Index of the tumor sample in `samples`
    """
    if not samples or len(samples) == 1:
        return 0

    # Caller-declared normal (Mutect2 ##normal_sample header) is ground truth
    if normal_sample and normal_sample in samples:
        return (samples.index(normal_sample) + 1) % len(samples)

    # Name-based resolution: tumor name wins outright
    for i, sample in enumerate(samples):
        sample_lower = sample.lower()
        if "tumor" in sample_lower or "tumour" in sample_lower:
            return i

    # A named normal identifies the tumor as the other sample
    for i, sample in enumerate(samples):
        if "normal" in sample.lower():
            return (i + 1) % len(samples)

    # Caller ordering conventions when names are uninformative
    if caller and caller.lower() in ("strelka", "mutect2"):
        # Normal sample is listed first for both callers in this pipeline
        return len(samples) - 1

    # DeepSomatic lists the tumor first; unknown callers keep the legacy default
    return 0


def resolve_normal_sample_index(samples, caller, normal_sample=None):
    """Resolve the normal sample index without inventing tumor-only evidence."""
    if not samples or len(samples) < 2:
        return None
    if normal_sample and normal_sample in samples:
        return samples.index(normal_sample)
    for i, sample in enumerate(samples):
        if "normal" in sample.lower():
            return i
    if len(samples) == 2:
        tumor_idx = resolve_tumor_sample_index(samples, caller, normal_sample)
        return 1 - tumor_idx
    return None


def _normal_sample_from_header(raw_header):
    """
    Extract the sample name declared by a ``##normal_sample=<name>`` header
    line (written by Mutect2 for paired tumor/normal runs), if present.

    Args:
        raw_header (str): Full VCF header text (cyvcf2 ``VCF.raw_header``)

    Returns:
        str or None: The declared normal sample name, or None
    """
    if not raw_header:
        return None
    for line in raw_header.splitlines():
        if line.startswith("##normal_sample="):
            return line.split("=", 1)[1].strip() or None
    return None


def extract_genotype_info(variant, caller, sample_idx=0):
    """
    Extract genotype, depth, and VAF information from variant.

    This function extracts genotype information from a cyvcf2 Variant object,
    with special handling for caller-specific formats (particularly Strelka).
    It attempts to extract GT, DP, AD, VAF, and GQ fields using multiple
    strategies to handle different VCF format variations.

    Args:
        variant (cyvcf2.Variant): Variant object from cyvcf2 VCF reader
        caller (str): Caller name for caller-specific parsing. Special handling
            is provided for 'strelka' which uses non-standard format fields
            (TAR, TIR, AU, CU, GU, TU, SGT)
        sample_idx (int): Index of the sample to extract from. This should be
            the tumor sample index as resolved by resolve_tumor_sample_index().
            Default: 0 (legacy behavior for direct callers)

    Returns:
        dict: Genotype information dictionary with keys:
            - GT (str or None): Genotype string (e.g., '0/1', '1|1', '0|1')
                where '/' indicates unphased and '|' indicates phased
            - DP (int or None): Total read depth at the variant position
            - AD (str or None): Allele depths as comma-separated string (e.g., '50,25')
            - VAF (float or None): Variant allele frequency (0.0 to 1.0)
            - GQ (int or None): Genotype quality score

    Notes:
        - For Strelka variants, the function uses TAR/TIR (indels) or AU/CU/GU/TU
          (SNVs) to compute AD and VAF when standard fields are missing
        - VAF is calculated from AD if not directly available
        - DP is calculated from AD if not directly available
        - Returns None for fields that cannot be extracted
        - Prints warnings to stdout if extraction errors occur

    Example:
        >>> from cyvcf2 import VCF
        >>> vcf = VCF('variants.vcf.gz')
        >>> for variant in vcf:
        ...     gt_info = extract_genotype_info(variant, 'mutect2')
        ...     print(f"GT: {gt_info['GT']}, DP: {gt_info['DP']}, VAF: {gt_info['VAF']}")
        ...     break
        GT: 0/1, DP: 150, VAF: 0.25
    """
    info = {
        "GT": None,
        "DP": None,
        "AD": None,
        "VAF": None,
        "VAF_SOURCE": None,
        "GQ": None,
        "ALT_INDICES": [],
    }

    try:

        def fmt(key):
            try:
                return variant.format(key)
            except Exception:
                return None

        def pick(arr):
            # Per-sample FORMAT rows: select the tumor sample row, falling
            # back to row 0 if the field has fewer rows than samples
            return arr[sample_idx] if len(arr) > sample_idx else arr[0]

        # Extract genotype
        try:
            if variant.genotypes and len(variant.genotypes) > 0:
                gt = pick(variant.genotypes)
                if len(gt) >= 2:
                    a1 = "." if gt[0] == -1 else str(gt[0])
                    a2 = "." if gt[1] == -1 else str(gt[1])
                    phased = len(gt) > 2 and bool(gt[2])
                    sep = "|" if phased else "/"
                    info["GT"] = f"{a1}{sep}{a2}"
                    info["ALT_INDICES"] = sorted({a for a in (gt[0], gt[1]) if isinstance(a, int) and a > 0})
        except Exception:
            pass

        # Depth
        dp = fmt("DP")
        if dp is not None and len(dp) > 0:
            try:
                row = pick(dp)
                info["DP"] = int(row) if row is not None else None
            except Exception:
                pass

        # Allele depth
        ad = fmt("AD")
        if ad is not None and len(ad) > 0:
            row = pick(ad)
            if row is not None:
                try:
                    info["AD"] = ",".join(map(str, row))
                except Exception:
                    pass

        # VAF/AF
        if info["VAF"] is None:
            for vaf_field in ["VAF", "AF", "FREQ", "FA"]:
                vaf = fmt(vaf_field)
                if vaf is None or len(vaf) == 0:
                    continue
                row = pick(vaf)
                if row is None:
                    continue
                val = row
                if isinstance(val, str) and "%" in val:
                    val = float(val.replace("%", "")) / 100.0
                try:
                    info["VAF"] = float(val)
                    info["VAF_SOURCE"] = "reported"
                except Exception:
                    pass
                break

        info = sanitize_genotype(info)

        # Strelka-specific parsing when AD/AF are missing
        if caller.lower() == "strelka":
            # Extract from the tumor sample row explicitly (Strelka writes
            # samples in fixed NORMAL,TUMOR order); never mix rows across
            # samples, so DP and VAF always come from the same sample.
            row_idx = sample_idx
            if info["AD"] is None:
                tir = fmt("TIR")
                tar = fmt("TAR")
                if (
                    tir is not None
                    and tar is not None
                    and len(tir) > 0
                    and len(tar) > 0
                ):
                    try:
                        t_alt = tir[row_idx]
                        t_ref = tar[row_idx]
                        a1 = (
                            t_alt[0]
                            if hasattr(t_alt, "__len__") and len(t_alt) > 0
                            else t_alt
                        )
                        r1 = (
                            t_ref[0]
                            if hasattr(t_ref, "__len__") and len(t_ref) > 0
                            else t_ref
                        )
                        alt_c = int(a1) if a1 is not None else None
                        ref_c = int(r1) if r1 is not None else None
                        if ref_c is not None and alt_c is not None:
                            info["AD"] = f"{ref_c},{alt_c}"
                            if info["DP"] is None:
                                info["DP"] = ref_c + alt_c
                            if info["VAF"] is None and (ref_c + alt_c) > 0:
                                info["VAF"] = alt_c / (ref_c + alt_c)
                                info["VAF_SOURCE"] = "derived"
                    except Exception:
                        pass
                else:
                    au = fmt("AU")
                    cu = fmt("CU")
                    gu = fmt("GU")
                    tu = fmt("TU")
                    if (
                        au is not None
                        and cu is not None
                        and gu is not None
                        and tu is not None
                    ):
                        base_map = {"A": au, "C": cu, "G": gu, "T": tu}
                        ref_base = variant.REF if hasattr(variant, "REF") else None
                        alt_base = None
                        try:
                            if (
                                hasattr(variant, "ALT")
                                and variant.ALT
                                and len(variant.ALT) > 0
                            ):
                                alt_base = variant.ALT[0]
                        except Exception:
                            pass
                        try:
                            ref_pair = base_map.get(str(ref_base))
                            alt_pair = base_map.get(str(alt_base))
                            if (
                                ref_pair is not None
                                and alt_pair is not None
                                and len(ref_pair) > 0
                                and len(alt_pair) > 0
                            ):
                                ref_counts = ref_pair[row_idx]
                                alt_counts = alt_pair[row_idx]
                                r1 = (
                                    ref_counts[0]
                                    if hasattr(ref_counts, "__len__")
                                    and len(ref_counts) > 0
                                    else ref_counts
                                )
                                a1 = (
                                    alt_counts[0]
                                    if hasattr(alt_counts, "__len__")
                                    and len(alt_counts) > 0
                                    else alt_counts
                                )
                                ref_c = int(r1) if r1 is not None else None
                                alt_c = int(a1) if a1 is not None else None
                                info["AD"] = f"{ref_c},{alt_c}"
                                if info["DP"] is None:
                                    info["DP"] = ref_c + alt_c
                                if info["VAF"] is None and (ref_c + alt_c) > 0:
                                    info["VAF"] = alt_c / (ref_c + alt_c)
                                    info["VAF_SOURCE"] = "derived"
                        except Exception:
                            pass

            # Use SGT when GT is absent
            if info["GT"] is None:
                sgt = fmt("SGT")
                if sgt is not None and len(sgt) > 0:
                    row = pick(sgt)
                    if row is not None:
                        try:
                            info["GT"] = str(row)
                        except Exception:
                            pass

        info = sanitize_genotype(info)

        # Calculate VAF from AD if still not available
        if info["VAF"] is None and info["AD"] is not None:
            try:
                ad_values = [int(x) for x in info["AD"].split(",")]
                if len(ad_values) >= 2 and sum(ad_values) > 0:
                    indices = info.get("ALT_INDICES") or [1]
                    selected = [ad_values[i] for i in indices if i < len(ad_values)]
                    info["VAF"] = sum(selected) / sum(ad_values) if selected else None
                    info["VAF_SOURCE"] = "derived"
            except Exception:
                pass

        if info["DP"] is None and info["AD"] is not None:
            try:
                ad_values = [int(x) for x in info["AD"].split(",")]
                info["DP"] = sum(ad_values)
            except Exception:
                pass

        # Genotype quality
        gq = fmt("GQ")
        if gq is not None and len(gq) > 0:
            try:
                row = pick(gq)
                info["GQ"] = int(row) if row is not None else None
            except Exception:
                pass

    except (KeyError, IndexError, ValueError, TypeError) as e:
        print(f"Warning: Data format error extracting genotype info from {caller}: {e}")
    except Exception as e:
        print(f"Error: Unexpected error extracting genotype info from {caller}: {e}")
        # Don't re-raise for genotype extraction as it's not critical

    return sanitize_genotype(info)


def tumor_alt_count_from_genotype(genotype_info):
    """
    Tumor alt-read count from an extract_genotype_info() dict.

    The AD field is stored as a comma-separated "ref,alt[,alt2...]" string
    extracted from the tumor sample (see resolve_tumor_sample_index); the
    alt count is the second value.

    Args:
        genotype_info (dict or None): Genotype info dict from
            extract_genotype_info()

    Returns:
        int or None: Tumor alt-read count, or None if no AD evidence exists
    """
    if not genotype_info:
        return None
    ad = genotype_info.get("AD")
    if not ad:
        return None
    try:
        values = [int(x) for x in str(ad).split(",")]
    except (ValueError, TypeError):
        return None
    if len(values) < 2 or any(value < 0 for value in values):
        return None
    indices = genotype_info.get("ALT_INDICES") or [1]
    if any(not isinstance(i, int) or i <= 0 or i >= len(values) for i in indices):
        return None
    selected = [values[i] for i in indices]
    return max(selected) if selected else None


def _counts_toward_support(variant_data, min_alt_support):
    """
    Decide whether a caller's record votes toward consensus support.

    A record counts only if:
    - the caller itself did not reject it (per-record classification is not
      Artifact — audit finding M9), and
    - the tumor alt-read evidence meets the floor (audit minor note: a caller
      PASS with 1-2 alt reads must not count as a full Somatic vote). The
      positive floor requires available allele-specific AD evidence.
      Consensus labels are handled separately by the rescue classifier and
      never substitute for an independent caller vote.

    Args:
        variant_data (dict): Per-caller record from read_variants_from_vcf()
        min_alt_support (int): Minimum tumor alt reads for a support vote
            (0 disables the floor)

    Returns:
        bool: True if the record counts toward caller support
    """
    if variant_data.get("classification") == "Artifact":
        return False
    if variant_data.get("is_multiallelic"):
        return False
    genotype = variant_data.get("genotype") or {}
    if len(genotype.get("ALT_INDICES") or []) > 1:
        # A joint multi-ALT record is retained for provenance but cannot cast a
        # scalar support vote without decomposition into allele-specific rows.
        return False
    alt_count = tumor_alt_count_from_genotype(variant_data.get("genotype"))
    if min_alt_support and (alt_count is None or alt_count < min_alt_support):
        return False
    return True


def aggregate_genotypes(genotypes_by_caller, callers_order):
    """
    Aggregate genotype information across callers.

    This function combines genotype information from multiple variant callers,
    computing consensus genotypes and summary statistics for depth and VAF.
    The consensus genotype is determined by majority vote across callers.

    Args:
        genotypes_by_caller (dict): Dictionary mapping caller names to genotype
            info dictionaries (as returned by extract_genotype_info). Each
            genotype info dict should contain GT, DP, VAF, AD, and GQ fields.
        callers_order (list): Ordered list of caller names. This order is used
            to maintain consistent ordering in the output lists (gt_by_caller,
            dp_by_caller, vaf_by_caller).

    Returns:
        dict: Aggregated genotype statistics with keys:
            - consensus_gt (str or None): Most common genotype across callers
            - gt_list (list): List of all non-None genotypes
            - dp_values (list): List of all non-None depth values
            - dp_mean (float or None): Mean depth across callers
            - dp_min (int or None): Minimum depth across callers
            - dp_max (int or None): Maximum depth across callers
            - vaf_values (list): List of all non-None VAF values
            - vaf_mean (float or None): Mean VAF across callers
            - vaf_min (float or None): Minimum VAF across callers
            - vaf_max (float or None): Maximum VAF across callers
            - ad_values (list): List of all non-None allele depth strings
            - gq_values (list): List of all non-None genotype quality values
            - gt_by_caller (list): Ordered list of genotypes (or '.' if None) by caller
            - dp_by_caller (list): Ordered list of depths (or None) by caller
            - vaf_by_caller (list): Ordered list of VAFs (or None) by caller
            - alt_count_by_caller (list): Ordered list of tumor alt-read counts
                (or None) by caller, derived from the tumor-sample AD
            - alt_count_max (int or None): Maximum tumor alt-read count across callers

    Example:
        >>> genotypes = {
        ...     'mutect2': {'GT': '0/1', 'DP': 100, 'VAF': 0.25, 'AD': '75,25', 'GQ': 99},
        ...     'strelka': {'GT': '0/1', 'DP': 120, 'VAF': 0.27, 'AD': '88,32', 'GQ': 95},
        ...     'deepsomatic': {'GT': '0/1', 'DP': 110, 'VAF': 0.26, 'AD': '81,29', 'GQ': 98}
        ... }
        >>> callers = ['mutect2', 'strelka', 'deepsomatic']
        >>> agg = aggregate_genotypes(genotypes, callers)
        >>> print(f"Consensus GT: {agg['consensus_gt']}")
        Consensus GT: 0/1
        >>> print(f"Mean DP: {agg['dp_mean']:.1f}, Mean VAF: {agg['vaf_mean']:.3f}")
        Mean DP: 110.0, Mean VAF: 0.260
    """
    from collections import defaultdict
    from statistics import mean

    agg = {
        "consensus_gt": None,
        "gt_list": [],
        "dp_values": [],
        "dp_mean": None,
        "dp_min": None,
        "dp_max": None,
        "vaf_values": [],
        "vaf_mean": None,
        "vaf_min": None,
        "vaf_max": None,
        "ad_values": [],
        "gq_values": [],
        "gt_by_caller": [],
        "dp_by_caller": [],
        "vaf_by_caller": [],
        "alt_count_by_caller": [],
        "alt_count_max": None,
    }

    gt_counts = defaultdict(int)

    for caller in callers_order:
        info = sanitize_genotype(genotypes_by_caller.get(caller, {}))
        # Collect GTs
        if info and "GT" in info and info["GT"] is not None:
            agg["gt_list"].append(info["GT"])
            gt_counts[info["GT"]] += 1
            agg["gt_by_caller"].append(info["GT"])
        else:
            agg["gt_by_caller"].append(".")

        # Collect DPs
        if info and "DP" in info and info["DP"] is not None:
            agg["dp_values"].append(info["DP"])
            agg["dp_by_caller"].append(info["DP"])
        else:
            agg["dp_by_caller"].append(None)

        # Collect VAFs
        if info and "VAF" in info and info["VAF"] is not None:
            agg["vaf_values"].append(info["VAF"])
            agg["vaf_by_caller"].append(info["VAF"])
        else:
            agg["vaf_by_caller"].append(None)

        # Collect ADs
        if info and "AD" in info and info["AD"] is not None:
            agg["ad_values"].append(info["AD"])

        # Collect GQs
        if info and "GQ" in info and info["GQ"] is not None:
            agg["gq_values"].append(info["GQ"])

        # Collect tumor alt counts (from AD, tumor sample — see
        # resolve_tumor_sample_index / extract_genotype_info)
        alt_count = tumor_alt_count_from_genotype(info)
        agg["alt_count_by_caller"].append(alt_count)
        if alt_count is not None and (
            agg["alt_count_max"] is None or alt_count > agg["alt_count_max"]
        ):
            agg["alt_count_max"] = alt_count

    # Determine consensus genotype (most common)
    if gt_counts:
        agg["consensus_gt"] = max(gt_counts, key=gt_counts.get)

    # Calculate DP statistics
    if agg["dp_values"]:
        agg["dp_mean"] = mean(agg["dp_values"])
        agg["dp_min"] = min(agg["dp_values"])
        agg["dp_max"] = max(agg["dp_values"])

    # Calculate VAF statistics
    if agg["vaf_values"]:
        agg["vaf_mean"] = mean(agg["vaf_values"])
        agg["vaf_min"] = min(agg["vaf_values"])
        agg["vaf_max"] = max(agg["vaf_values"])

    return agg


def read_variants_from_vcf(
    vcf_path,
    caller_name,
    modality=None,
    exclude_refcall=False,
    exclude_germline=False,
    classify_variants=True,
    include_non_canonical=False,
    chrom=None,
    alignment_round="unknown",
):
    """
    Read variants from a single VCF file with biological classification.

    This function reads all variants from a VCF file and extracts relevant
    information including position, alleles, filters, quality, and genotype
    data. It supports optional filtering to exclude reference calls and
    germline variants, and can tag variants with modality information.
    Additionally, it classifies variants into biological categories.

    Args:
        vcf_path (str): Path to VCF file (can be .vcf, .vcf.gz, or .bcf)
        caller_name (str): Name of the variant caller (e.g., 'mutect2', 'strelka')
        modality (str, optional): Modality tag to add to variants ('DNA' or 'RNA').
            If None, no modality tag is added. Default: None
        exclude_refcall (bool, optional): If True, exclude variants with 'RefCall'
            in the FILTER field. Useful for DeepSomatic output. Default: False
        exclude_germline (bool, optional): If True, exclude variants with 'GERMLINE'
            in the FILTER field. Default: False
        classify_variants (bool, optional): If True, classify variants into
            biological categories (Somatic, Germline, Reference, Artifact).
            Default: True
        include_non_canonical (bool, optional): If True, include variants on
            non-canonical chromosomes. Default: False (only chr1-22, X, Y, M)
        chrom (str, optional): If set, only read variants on this chromosome
            (any naming style — matched after normalize_chromosome, with M/MT
            equivalence). Uses a tabix/csi region query when an index exists;
            falls back to a full scan with early skipping otherwise. Used by
            per-chromosome streaming drivers to bound memory. Default: None
            (read the whole file)

    Returns:
        dict: Dictionary mapping variant_key to variant_data. Each variant_data
            dict contains:
            - CHROM (str): Chromosome name
            - POS (int): 1-based position
            - REF (str): Reference allele
            - ALT (str): Comma-separated alternate allele(s)
            - is_snv (bool): True if variant is a single nucleotide variant
            - caller (str): Caller name
            - modality (str, optional): Modality tag (only if modality arg provided)
            - filter_original (str): Original FILTER field value
            - filter_normalized (str): Normalized filter string (unified)
            - filter_category (str): Filter category
            - classification (str): Biological classification (if classify_variants=True)
            - quality (float or None): QUAL score
            - genotype (dict): Genotype information from extract_genotype_info()
            - id (str or None): Variant ID from ID field

    Notes:
        - Variant keys are generated using variant_key() function which creates
          unique identifiers in format: "chrom:pos:ref:alt"
        - Filter normalization and categorization use functions from vcf_utils.filters
        - Genotype extraction uses extract_genotype_info() with caller-specific handling
        - Classification uses vcf_utils.classification module for unified categorization
        - Filter values are unified: PASS, GERMLINE, RefCall, LowQuality

    Example:
        >>> # Read DNA variants from Mutect2 with classification
        >>> dna_vars = read_variants_from_vcf('mutect2.vcf.gz', 'mutect2',
        ...                                    modality='DNA', exclude_germline=True)
        >>> print(f"Read {len(dna_vars)} DNA variants from Mutect2")
        >>>
        >>> # Read RNA variants from Strelka
        >>> rna_vars = read_variants_from_vcf('strelka.vcf.gz', 'strelka', modality='RNA')
        >>>
        >>> # Access variant data with classification
        >>> for vkey, data in list(dna_vars.items())[:3]:
        ...     print(f"{vkey}: {data['classification']} -> {data['filter_normalized']}")
    """
    from cyvcf2 import VCF

    from vcf_utils.filters import normalize_filter
    from vcf_utils.io_utils import is_snv, variant_key

    # Import classification functions if needed
    if classify_variants:
        from vcf_utils.classification import (
            classify_variant_from_record,
            get_sample_indices,
            normalize_filter_value,
        )

    variants = {}

    vcf = VCF(vcf_path)

    # Get sample indices for classification (needed for Strelka)
    sample_indices = None
    if classify_variants:
        sample_indices = get_sample_indices(vcf, caller_name)

    # Resolve the tumor sample once per VCF so genotype metrics (GT/AD/DP/VAF)
    # are extracted from the tumor sample, not blindly from sample index 0
    # (which is the NORMAL for Mutect2 and Strelka in this pipeline). Mutect2's
    # ##normal_sample header is ground truth and takes priority over heuristics.
    normal_sample = _normal_sample_from_header(vcf.raw_header)
    tumor_sample_idx = resolve_tumor_sample_index(vcf.samples, caller_name, normal_sample)
    normal_sample_idx = resolve_normal_sample_index(vcf.samples, caller_name, normal_sample)

    # Import chromosome filtering
    from vcf_utils.chromosome_utils import is_canonical_chromosome

    # Per-chromosome streaming: restrict iteration to the requested contig.
    # Region queries need a tabix/csi index; without one, fall back to a full
    # scan with an early chromosome skip (slower, same result).
    record_iter = vcf
    if chrom is not None:
        import os

        from vcf_utils.io_utils import normalize_chromosome

        target = normalize_chromosome(chrom)
        match = {"M", "MT"} if target in ("M", "MT") else {target}
        seqname = next(
            (n for n in vcf.seqnames if normalize_chromosome(n) in match), None
        )
        if seqname is None:
            vcf.close()
            return variants
        if os.path.exists(vcf_path + ".tbi") or os.path.exists(vcf_path + ".csi"):
            record_iter = vcf(seqname)
        else:
            record_iter = (
                v for v in vcf if normalize_chromosome(v.CHROM) in match
            )

    for variant in record_iter:
        # Skip non-canonical chromosomes if filtering is enabled
        if not include_non_canonical:
            if not is_canonical_chromosome(variant.CHROM):
                continue

        # Get filter and check exclusions
        filter_str = variant.FILTER if variant.FILTER else "PASS"

        if exclude_refcall and "RefCall" in filter_str:
            continue
        if exclude_germline and "GERMLINE" in filter_str:
            continue

        vkey = variant_key(variant, use_cyvcf2=True)

        # Classify variant before creating data dict
        classification = None
        if classify_variants:
            # Optimization: For consensus callers, FILTER field already contains biological category
            if "consensus" in caller_name.lower():
                # Use FILTER field directly for consensus VCFs
                if filter_str in ["Somatic", "Germline", "Reference", "Artifact"]:
                    classification = filter_str
                elif filter_str == "NoConsensus":
                    # NoConsensus variants shouldn't be in consensus VCFs being rescued
                    # but if they are, classify as Artifact
                    classification = "NoConsensus"
                else:
                    # Fallback to classification function
                    try:
                        classification = classify_variant_from_record(
                            variant, caller_name, sample_indices
                        )
                    except Exception:
                        classification = "Artifact"
            else:
                # Regular caller - classify using caller-specific logic
                try:
                    classification = classify_variant_from_record(
                        variant, caller_name, sample_indices
                    )
                except Exception:
                    classification = "Artifact"

        # Create variant data
        source_evidence = {}
        for evidence_key in (
            "GT_BY_CALLER", "DP_BY_CALLER", "AD_BY_CALLER", "VAF_BY_CALLER",
            "VAF_SOURCE_BY_CALLER", "NORMAL_GT_BY_CALLER",
            "NORMAL_DP_BY_CALLER", "NORMAL_AD_BY_CALLER", "NORMAL_VAF_BY_CALLER",
            "NORMAL_VAF_SOURCE_BY_CALLER",
        ):
            try:
                value = variant.INFO.get(evidence_key)
            except Exception:
                value = None
            if value not in (None, "", "."):
                source_evidence[evidence_key] = str(value)

        data = {
            "CHROM": variant.CHROM,
            "POS": variant.POS,
            "REF": variant.REF,
            "ALT": ",".join(variant.ALT) if variant.ALT else ".",
            "is_snv": is_snv(variant.REF, variant.ALT),
            "is_multiallelic": bool(variant.ALT and len(variant.ALT) > 1),
            "caller": caller_name,
            "filter_original": filter_str,
            "filter_normalized": normalize_filter_value(classification)
            if classification
            else normalize_filter(filter_str),
            "quality": float(variant.QUAL) if variant.QUAL is not None else None,
            "genotype": extract_genotype_info(variant, caller_name, tumor_sample_idx),
            "normal_genotype": (
                extract_genotype_info(variant, caller_name, normal_sample_idx)
                if normal_sample_idx is not None
                else None
            ),
            "source_evidence": source_evidence,
            "id": variant.ID if variant.ID else None,
        }

        # Add classification
        if classification:
            data["classification"] = classification

        # Filter category is the same as normalized filter (both are biological categories)
        data["filter_category"] = data["filter_normalized"]

        # Add modality if provided
        if modality:
            data["modality"] = modality

        from .caller_evidence import evidence_from_record, read_info
        data["caller_evidence"] = evidence_from_record(
            read_info(variant), caller_name, vkey,
            genotype=data["genotype"], normal_genotype=data["normal_genotype"],
            modality=modality or "unknown", alignment_round=alignment_round,
            sample_id=vcf.samples[tumor_sample_idx] if vcf.samples else "unknown",
            normal_sample_id=vcf.samples[normal_sample_idx] if normal_sample_idx is not None else "unknown",
        )

        variants[vkey] = data

    return variants


def aggregate_variants(
    variant_collections, snv_threshold=2, indel_threshold=2, min_alt_support=None,
    preserve_baseline_callers=None,
):
    """
    Aggregate variants from multiple collections.

    This function combines variants from multiple callers and/or modalities,
    merging information for variants at the same genomic position. It applies
    consensus thresholds to determine which variants have sufficient support,
    and aggregates genotype information across all supporting callers.

    Args:
        variant_collections (list): List of tuples, each containing:
            - caller_name (str): Name of the caller
            - variant_dict (dict): Dictionary of variants from read_variants_from_vcf()
            - modality (str or None): Modality tag ('DNA', 'RNA', or None)
        snv_threshold (int, optional): Minimum number of callers required for
            a SNV to pass consensus. Default: 2
        indel_threshold (int, optional): Minimum number of callers required for
            an indel to pass consensus. Default: 2
        min_alt_support (int, optional): Minimum tumor alt reads for a caller's
            record to count toward consensus support (audit minor note: a 1-2
            alt-read caller PASS must not count as a full Somatic vote).
            Records without alt-read evidence are not floored (legacy
            behavior); 0 disables the floor. Default: None, which uses
            DEFAULT_THRESHOLDS["consensus_min_alt_support"]

    Returns:
        dict: Dictionary mapping variant_key to aggregated variant_data. Each
            aggregated variant_data dict contains:
            - CHROM (str): Chromosome name
            - POS (int): 1-based position
            - REF (str): Reference allele
            - ALT (str): Comma-separated alternate allele(s)
            - is_snv (bool): True if variant is SNV
            - callers (list): List of all caller names that detected this variant
            - modalities (list): List of modalities (if provided in collections)
            - caller_modality_map (dict): Maps caller name to modality
            - filters_original (list): Original filter strings from each caller
            - filters_normalized (list): Normalized filter strings
            - filters_category (list): Filter categories
            - qualities (list): Quality scores from each caller
            - genotypes (dict): Maps caller name to genotype info dict
            - ids (list): Variant IDs from each caller
            - support_callers (set): Set of unique callers whose record counts
                toward consensus support (caller did not itself reject the
                record; min_alt_support floor applied when AD is available)
            - passes_consensus (bool): True if variant meets consensus threshold
            - gt_aggregated (dict): Aggregated genotype statistics from
                aggregate_genotypes()

    Notes:
        - Variants are matched by genomic position (chrom:pos:ref:alt)
        - The consensus threshold is applied based on variant type (SNV vs indel)
        - Genotype aggregation is performed automatically for all variants
        - Modality information is preserved when provided in variant_collections

    Example:
        >>> # Aggregate within-modality consensus (DNA only)
        >>> mutect2_vars = read_variants_from_vcf('mutect2.vcf.gz', 'mutect2')
        >>> strelka_vars = read_variants_from_vcf('strelka.vcf.gz', 'strelka')
        >>> collections = [
        ...     ('mutect2', mutect2_vars, None),
        ...     ('strelka', strelka_vars, None)
        ... ]
        >>> aggregated = aggregate_variants(collections, snv_threshold=2, indel_threshold=2)
        >>>
        >>> # Count consensus variants
        >>> consensus_count = sum(1 for v in aggregated.values() if v['passes_consensus'])
        >>> print(f"Variants passing consensus: {consensus_count}")
        >>>
        >>> # Aggregate cross-modality rescue
        >>> dna_vars = read_variants_from_vcf('dna_consensus.vcf.gz', 'consensus', 'DNA')
        >>> rna_vars = read_variants_from_vcf('rna_consensus.vcf.gz', 'consensus', 'RNA')
        >>> collections = [
        ...     ('consensus', dna_vars, 'DNA'),
        ...     ('consensus', rna_vars, 'RNA')
        ... ]
        >>> rescue_aggregated = aggregate_variants(collections, snv_threshold=1, indel_threshold=1)
    """
    from collections import defaultdict

    preserve_baseline_callers = {str(c).lower() for c in (preserve_baseline_callers or set())}

    if min_alt_support is None:
        from vcf_utils.classification_config import DEFAULT_THRESHOLDS

        min_alt_support = DEFAULT_THRESHOLDS["consensus_min_alt_support"]

    aggregated = defaultdict(
        lambda: {
            "CHROM": None,
            "POS": None,
            "REF": None,
            "ALT": None,
            "is_snv": None,
            "is_multiallelic": False,
            "callers": [],
            "modalities": [],
            "caller_modality_map": {},
            "filters_original": [],
            "filters_normalized": [],
            "filters_category": [],
            "qualities": [],
            "genotypes": {},
            "normal_genotypes": {},
            "source_evidence": {},
            "caller_evidence": [],
            "ids": [],
            "support_callers": set(),
        }
    )

    # Process each variant collection
    for caller_name, variant_dict, modality in variant_collections:
        for vkey, variant_data in variant_dict.items():
            data = aggregated[vkey]

            # Store variant location info (first time)
            if data["CHROM"] is None:
                data["CHROM"] = variant_data["CHROM"]
                data["POS"] = variant_data["POS"]
                data["REF"] = variant_data["REF"]
                data["ALT"] = variant_data["ALT"]
                data["is_snv"] = variant_data["is_snv"]
                data["is_multiallelic"] = variant_data.get("is_multiallelic", False)

            # Store caller-specific information. A caller's record only counts
            # toward consensus support if the caller itself did not reject it
            # (non-Artifact classification, audit M9) and its tumor alt-read
            # evidence meets the min_alt_support floor (when AD is available).
            data["callers"].append(caller_name)
            if _counts_toward_support(variant_data, min_alt_support):
                data["support_callers"].add(caller_name)

            # Store modality information
            if modality:
                data["modalities"].append(modality)
                data["caller_modality_map"][caller_name] = modality

            # Store filters
            data["filters_original"].append(variant_data["filter_original"])
            data["filters_normalized"].append(variant_data["filter_normalized"])
            data["filters_category"].append(variant_data["filter_category"])

            # Store quality
            if variant_data["quality"] is not None:
                data["qualities"].append(variant_data["quality"])

            # Store ID
            if variant_data["id"]:
                data["ids"].append(variant_data["id"])

            # Store genotype information
            data["genotypes"][caller_name] = sanitize_genotype(variant_data["genotype"])
            normal = variant_data.get("normal_genotype")
            data["normal_genotypes"][caller_name] = sanitize_genotype(normal) if normal is not None else None
            from .caller_evidence import evidence_from_record, bind_evidence
            observations = variant_data.get("caller_evidence")
            if observations is None:
                observations = evidence_from_record(
                    variant_data.get("source_evidence", {}), caller_name, vkey,
                    genotype=data["genotypes"][caller_name],
                    normal_genotype=data["normal_genotypes"][caller_name], modality=modality or "unknown",
                )
            data["caller_evidence"].extend(bind_evidence(observations, modality=modality or "unknown"))
            for key, value in variant_data.get("source_evidence", {}).items():
                if key not in data["source_evidence"]:
                    data["source_evidence"][key] = value
                elif data["source_evidence"][key] != value:
                    data["source_evidence"][key] += "||" + value

    # Apply consensus thresholds and aggregate genotypes
    for vkey, data in aggregated.items():
        n_support = len(data["support_callers"])

        # Determine if variant passes consensus threshold
        if data["is_snv"]:
            data["passes_consensus"] = n_support >= snv_threshold
        else:
            data["passes_consensus"] = n_support >= indel_threshold

        # Aggregate genotype information
        data["gt_aggregated"] = aggregate_genotypes(data["genotypes"], data["callers"])
        data["preserve_baseline"] = any(
            str(c).lower() in preserve_baseline_callers
            and c in data["support_callers"]
            and i < len(data["filters_normalized"])
            and data["filters_normalized"][i] == "Somatic"
            for i, c in enumerate(data["callers"])
        )

    return dict(aggregated)
