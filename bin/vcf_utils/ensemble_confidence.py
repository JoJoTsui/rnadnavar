"""
Ensemble confidence statistics for consensus variant calling.

Per-variant ensemble confidence is the Wilson score confidence interval on
the caller support fraction: k supporting callers out of n callers configured
for the invocation. A caller with no record at the site counts as a
non-support vote. Written to consensus VCF INFO as ENS_SUPPORT / ENS_CONF_LO /
ENS_CONF_HI (see docs/consensus_vcf_rules.md).
"""

import math

# z-score for a 95% confidence interval
DEFAULT_Z = 1.96


def wilson_interval(k, n, z=DEFAULT_Z):
    """
    Wilson score confidence interval for a binomial proportion.

    Args:
        k (int): Number of supporting callers (successes).
        n (int): Total number of callers configured for the invocation
            (trials); callers absent at the site count as non-support.
        z (float): z-score for the confidence level. Default: 1.96 (95%).

    Returns:
        tuple: (lower, upper) bounds in [0, 1]. The lower bound is 0 when
            k == 0; the upper bound is 1 when k == n.
        None: when n <= 0 (no callers configured — nothing to annotate).
    """
    if n <= 0:
        return None
    k = max(0, min(k, n))
    z2 = z * z
    p_hat = k / n
    denom = 1 + z2 / n
    center = (p_hat + z2 / (2 * n)) / denom
    half_width = (z / denom) * math.sqrt(p_hat * (1 - p_hat) / n + z2 / (4 * n * n))
    return (max(0.0, center - half_width), min(1.0, center + half_width))
