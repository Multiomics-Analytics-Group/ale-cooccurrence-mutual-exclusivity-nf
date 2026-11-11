"""
pairwise_tests.py

The three analysis methods, each a separate function operating on the SAME
shared binary matrix (rows = units, columns = samples):

    run_discover(matrix, tested_units, ...)  -> DISCOVER co-occurrence + mutual exclusivity
    run_fisher(matrix, tested_units, ...)    -> Fisher exact (two-sided + greater + less)
    run_pearson(matrix, tested_units)        -> Pearson/phi coefficient (no p-values)

Plus pair_contingency() for the method-agnostic 2x2 fields and correlation_matrix()
for the Pearson/phi matrix of the tested units. None of these functions know what a
"unit" represents. DISCOVER's p-values are used and its q-values are recomputed over
the testable pairs (no pi0 * qvalue rescaling).
"""
import os
from contextlib import contextmanager
from itertools import combinations

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from scipy.stats import fisher_exact, hypergeom
from statsmodels.stats.multitest import multipletests

import discover
from _discover import fdr as _discover_fdr

# A direction label needs a difference larger than numerical error. The expected overlap is
# computed from DISCOVER's background, which is fitted iteratively and so is only approximate.
# A pair whose observed and expected overlaps differ by less than DISCOVER_DIRECTION_TOL (in
# endpoints) has the direction "none". A Pearson r closer to zero than PEARSON_DIRECTION_TOL
# is zero up to rounding and has the direction "none" too.
DISCOVER_DIRECTION_TOL = 1e-4
PEARSON_DIRECTION_TOL = 1e-12


@contextmanager
def _silence_c_stdout():
    """Silence the DISCOVER maxent solver's C-level stdout (fd 1) only."""
    saved = os.dup(1)
    devnull = os.open(os.devnull, os.O_WRONLY)
    try:
        os.dup2(devnull, 1)
        yield
    finally:
        os.dup2(saved, 1)
        os.close(devnull)
        os.close(saved)


def _key(a, b):
    return (a, b) if str(a) <= str(b) else (b, a)


def _bh(pvals):
    pvals = list(pvals)
    if not pvals:
        return []
    return list(multipletests(pvals, method="fdr_bh")[1])


def _bh_subset(pvals, keep):
    """BH over ONLY the pairs flagged True in `keep`.

    A pair is kept when the direction being tested could, in principle, have produced a
    significant result: the p-value at the most extreme overlap the margins permit is
    below the threshold. Pairs that fail can never reject whatever the data show, so
    charging the surviving tests for them costs power and buys nothing. The criterion uses
    only (N, c1, c2) and, for DISCOVER, the fitted background with its structurally forced
    entries, never the observed overlap.
    Excluded pairs get q = NaN and are reported as not significant.
    """
    q = np.full(len(pvals), np.nan)
    idx = [i for i, k in enumerate(keep) if k]
    if not idx:
        return q
    sub = [pvals[i] for i in idx]
    vals = _bh(sub)
    for i, v in zip(idx, vals):
        q[i] = v
    return q


def _discover_dbh_subset(pvals, keep, rows, bg, lowertail):
    """DISCOVER's discrete Benjamini-Hochberg over ONLY the pairs flagged True in `keep`.

    Source: NKI-CCB/DISCOVER, python/src/fdr.f90, function `mutex` and subroutine `updateQ`.
    For every test, the steps of its null CDF (the Poisson-binomial with probabilities
    bg[gene1] * bg[gene2]) are added to an accumulator at the sorted observed p-value levels,
    so after a cumulative sum each level holds the sum of all tests' null CDFs there. That sum
    is divided by the number of observed p-values at or below the level, capped at 1, and made
    monotone by a cumulative minimum from the top. The per-test step is DISCOVER's own compiled
    `fdr.updatemultiq`, called here once per kept pair, so the family is the testable pairs
    instead of every pair. With every pair kept this reproduces DISCOVER's own DBH q-values.

    rows[i] is (row of gene1, row of gene2) in `bg` for pair i. lowertail is 1 for mutual
    exclusivity and 0 for co-occurrence.
    """
    q = np.full(len(pvals), np.nan)
    idx = [i for i, k in enumerate(keep) if k]
    if not idx:
        return q
    levels, inv = np.unique([pvals[i] for i in idx], return_inverse=True)
    acc = np.zeros(levels.size)
    expected_p0 = np.zeros(())   # DISCOVER's pi0 bookkeeping, not used here
    for i in idx:
        r1, r2 = rows[i]
        _discover_fdr.updatemultiq(bg[[r1]], bg[[r2]], lowertail, levels, acc, expected_p0)
    counts = np.cumsum(np.bincount(inv, minlength=levels.size))
    acc = np.minimum(np.cumsum(acc) / counts, 1.0)
    acc = np.minimum.accumulate(acc[::-1])[::-1]
    for i, v in zip(idx, acc[inv]):
        q[i] = v
    return q


def _structurally_fixed(X):
    """Entries of the binary matrix X (genes x endpoints) that its totals force.

    DISCOVER's background keeps every gene's and every endpoint's total. An entry that is 0
    (or 1) in every binary matrix with those totals is exactly 0 (or 1) in the background, but
    the maximum-entropy fit only approaches that value, so the fitted numbers cannot be relied
    on to show it (a forced 0 can come out as a small positive number). An entry is forced
    iff no alternating cycle passes through it: in the bipartite digraph with an edge
    gene -> endpoint for each 0 and endpoint -> gene for each 1, the gene and the endpoint lie
    in different strongly connected components. This covers genes mutated in every endpoint,
    endpoints with every gene mutated, and every entry the totals force through them, such as
    the other genes at an endpoint whose only mutations are in genes mutated everywhere.
    Returns a boolean mask of X's shape. The forced value is X itself.
    """
    X = np.asarray(X, dtype=bool)
    m, n = X.shape
    r0, c0 = np.nonzero(~X)
    r1, c1 = np.nonzero(X)
    src = np.concatenate([r0, m + c1])
    dst = np.concatenate([m + c0, r1])
    graph = csr_matrix((np.ones(src.size), (src, dst)), shape=(m + n, m + n))
    _, label = connected_components(graph, directed=True, connection="strong")
    return label[:m, None] != label[None, m:]


def _tested(matrix, tested_units):
    units = [u for u in tested_units if u in matrix.index]
    return matrix.loc[units], units


def pair_contingency(matrix, tested_units):
    """Method-agnostic 2x2 contingency fields for every tested pair."""
    sub, units = _tested(matrix, tested_units)
    n = sub.shape[1]
    arr = {u: sub.loc[u].values.astype(bool) for u in units}
    rows = []
    for g1, g2 in combinations(units, 2):
        x1, x2 = arr[g1], arr[g2]
        n11 = int(np.sum(x1 & x2)); n10 = int(np.sum(x1 & ~x2))
        n01 = int(np.sum(~x1 & x2)); n00 = int(np.sum(~x1 & ~x2))
        c1, c2 = n11 + n10, n11 + n01
        union = n11 + n10 + n01
        u1, u2 = _key(g1, g2)
        rows.append({
            "unit1": u1, "unit2": u2,
            "count_unit1": c1 if u1 == g1 else c2,
            "count_unit2": c2 if u1 == g1 else c1,
            "only_unit1": (n10 if u1 == g1 else n01),
            "only_unit2": (n01 if u1 == g1 else n10), "neither": n00,
            "observed_overlap": n11,
            "freq_unit1": (c1 if u1 == g1 else c2) / n,
            "freq_unit2": (c2 if u1 == g1 else c1) / n,
            "jaccard_index": (n11 / union) if union else 0.0,
        })
    return pd.DataFrame(rows)


def run_discover(matrix, tested_units, significance_threshold, fdr_method):
    """DISCOVER co-occurrence ('greater') and mutual exclusivity ('less').

    fdr_method: 'BH' (standard Benjamini-Hochberg) or 'DBH' (DISCOVER's
    discrete Benjamini-Hochberg, less conservative for discrete tests). Either one is
    applied within each direction over only the testable pairs.
    """
    if fdr_method not in ("BH", "DBH"):
        raise ValueError(f"fdr_method must be 'BH' or 'DBH', got {fdr_method!r}")
    with _silence_c_stdout():
        events = discover.DiscoverMatrix(matrix)
        mask = np.array([u in set(tested_units) for u in matrix.index])
        sub = events[mask]
        # Only the p-values are used. DISCOVER's own q-values cover every pair rather than
        # the testable ones, so the cheap BH setting is requested and q is recomputed below.
        co = discover.pairwise_discover_test(sub, alternative="greater", fdr_method="BH")
        me = discover.pairwise_discover_test(sub, alternative="less", fdr_method="BH")
    # DISCOVER's pi0 is not reported: it is 1.0 under BH, and neither correction below uses it
    # (DISCOVER's pi0 * qvalue rule is not applied).

    def pq(result, names):
        """Every pair's raw p-value, read from the triangular matrix.

        NOT via significant_pairs(q_threshold=1.0): that filter is a strict '<', so pairs
        whose q equals exactly 1.0 are silently dropped. Harmless when only significant
        pairs are wanted, but this function needs the complete set to re-correct over a
        filtered pool.
        """
        P = np.asarray(result.pvalues.values if hasattr(result.pvalues, "values")
                       else result.pvalues, dtype=float)
        d = {}
        for i in range(len(names)):
            for j in range(i + 1, len(names)):
                v = P[i, j]
                if np.isnan(v):
                    v = P[j, i]
                d[_key(names[i], names[j])] = (float(v), np.nan)
        return d
    _names = list(sub.rownames)
    co_pq, me_pq = pq(co, _names), pq(me, _names)

    bg = pd.DataFrame(sub.bg, index=sub.rownames, columns=sub.colnames)
    X = pd.DataFrame(sub.events, index=sub.rownames, columns=sub.colnames).astype(bool)

    # Background entries forced to exactly 0 or 1 by the totals, read from the BINARY MATRIX
    # rather than from fitted values: the maximum-entropy fit only approaches 0 or 1, so a
    # test on the fitted numbers misses them. The check uses the full matrix, since the
    # background was fitted on all of its genes.
    _full = np.asarray(matrix.values) > 0
    _fixed = _structurally_fixed(_full)
    _row = {u: i for i, u in enumerate(matrix.index)}
    _log_thresh = np.log(significance_threshold)
    _pos = {u: i for i, u in enumerate(_names)}

    rows, p_co_l, p_me_l, keep_co, keep_me, pair_rows = [], [], [], [], [], []
    for g1, g2 in combinations(list(sub.rownames), 2):
        k = _key(g1, g2)
        pair_rows.append((_pos[g1], _pos[g2]))
        qj = bg.loc[g1].values * bg.loc[g2].values
        # q_j is fixed when both entries are forced, or when either is forced to 0. Substitute
        # the exact value the totals imply for the optimizer's approximation to it. This is not
        # rounding: q_j is 0 or 1 by construction.
        f1, f2 = _fixed[_row[g1]], _fixed[_row[g2]]
        x1, x2 = _full[_row[g1]], _full[_row[g2]]
        exact = (f1 & f2) | (f1 & ~x1) | (f2 & ~x2)
        qj = np.where(exact, (x1 & x2).astype(float), np.clip(qj, 1e-15, 1 - 1e-15))
        exp = float(qj.sum())
        var = float((qj * (1.0 - qj)).sum())
        obs = int(np.sum(X.loc[g1].values & X.loc[g2].values))
        p_co, _ = co_pq.get(k, (np.nan, np.nan))
        p_me, _ = me_pq.get(k, (np.nan, np.nan))
        # Testability, per direction, at the Poisson-binomial's TRUE support boundaries:
        #   x_max = n - #{j : q_j = 0}   ->  P(X >= x_max) = prod over q_j > 0 of q_j
        #   x_min =     #{j : q_j = 1}   ->  P(X <= x_min) = prod over q_j < 1 of (1 - q_j)
        # Endpoints with q_j = 1 drop out of the exclusivity product because they contribute
        # their success with probability 1 and cannot affect P(no success elsewhere), and
        # endpoints with q_j = 0 drop out of the co-occurrence product the same way. An empty
        # product is 1: the pair cannot move in that direction and is not testable.
        # Both are exact products, evaluated in LOG SPACE so a very small probability is
        # never mistaken for an impossible outcome through underflow.
        _rest = qj[qj < 1.0]
        log_floor_me = float(np.log1p(-_rest).sum())
        log_floor_co = float(np.log(qj[qj > 0.0]).sum())
        floor_co, floor_me = float(np.exp(log_floor_co)), float(np.exp(log_floor_me))
        rows.append({
            "unit1": k[0], "unit2": k[1],
            "discover_expected_overlap": exp,
            "discover_observed_minus_expected": obs - exp,
            "discover_enrichment_ratio": obs / max(exp, 1e-12),
            "discover_z": (obs - exp) / np.sqrt(max(var, 1e-12)),
            "discover_p_co": p_co,
            "discover_p_mutex": p_me,
            "discover_floor_co": floor_co, "discover_floor_mutex": floor_me,
            "discover_testable_co": bool(log_floor_co < _log_thresh),
            "discover_testable_mutex": bool(log_floor_me < _log_thresh),
            "discover_direction": ("co" if obs - exp >= DISCOVER_DIRECTION_TOL else
                                   "mutex" if exp - obs >= DISCOVER_DIRECTION_TOL else "none"),
        })
        p_co_l.append(p_co); p_me_l.append(p_me)
        keep_co.append(bool(log_floor_co < _log_thresh) and np.isfinite(p_co))
        keep_me.append(bool(log_floor_me < _log_thresh) and np.isfinite(p_me))
    df = pd.DataFrame(rows)
    if len(df):
        # Correct within each direction over only the pairs testable in that direction.
        if fdr_method == "DBH":
            df["discover_q_co_local"] = _discover_dbh_subset(p_co_l, keep_co, pair_rows, bg.values, 0)
            df["discover_q_mutex_local"] = _discover_dbh_subset(p_me_l, keep_me, pair_rows, bg.values, 1)
        else:
            df["discover_q_co_local"] = _bh_subset(p_co_l, keep_co)
            df["discover_q_mutex_local"] = _bh_subset(p_me_l, keep_me)
        df["discover_sig_co"] = df["discover_q_co_local"] < significance_threshold
        df["discover_sig_mutex"] = df["discover_q_mutex_local"] < significance_threshold
        df["discover_sig_any"] = df["discover_sig_co"] | df["discover_sig_mutex"]
    return df


def run_fisher(matrix, tested_units, significance_threshold):
    """Fisher exact: two-sided + one-sided greater (co) and less (mutex), local per direction.

    q-values are standard Benjamini-Hochberg within each direction over only the testable pairs.
    """
    sub, units = _tested(matrix, tested_units)
    n = sub.shape[1]
    arr = {u: sub.loc[u].values.astype(int) for u in units}
    recs, p_co, p_me = [], [], []
    keep_co, keep_me = [], []
    for g1, g2 in combinations(units, 2):
        x1, x2 = arr[g1], arr[g2]
        n11 = int(np.sum((x1 == 1) & (x2 == 1))); n10 = int(np.sum((x1 == 1) & (x2 == 0)))
        n01 = int(np.sum((x1 == 0) & (x2 == 1))); n00 = int(np.sum((x1 == 0) & (x2 == 0)))
        table = [[n11, n10], [n01, n00]]
        odds, p_two = fisher_exact(table, alternative="two-sided")
        pc = fisher_exact(table, alternative="greater")[1]
        pm = fisher_exact(table, alternative="less")[1]
        c1, c2 = n11 + n10, n11 + n01
        expected = c1 * c2 / n if n else np.nan
        log2or = float(np.log2(odds)) if (odds and np.isfinite(odds) and odds > 0) else np.nan
        # Testability, per direction: the p-value at the most extreme overlap the margins
        # permit. ceiling = min(c1,c2) for co-occurrence, floor = max(0,c1+c2-N) for
        # exclusivity. For Fisher these are the hypergeometric's own support boundaries.
        hi, lo = min(c1, c2), max(0, c1 + c2 - n)
        floor_co = float(hypergeom.sf(hi - 1, n, c1, c2)) if c1 and c2 else 1.0
        floor_me = float(hypergeom.cdf(lo, n, c1, c2)) if c1 and c2 else 1.0
        k = _key(g1, g2)
        recs.append({
            "unit1": k[0], "unit2": k[1],
            "fisher_expected_overlap": expected,
            "fisher_observed_minus_expected": n11 - expected,
            "fisher_odds_ratio": odds, "fisher_log2_odds_ratio": log2or,
            "fisher_p_two_sided": p_two, "fisher_p_co": pc, "fisher_p_mutex": pm,
            "fisher_floor_co": floor_co, "fisher_floor_mutex": floor_me,
            "fisher_testable_co": floor_co < significance_threshold,
            "fisher_testable_mutex": floor_me < significance_threshold,
            "fisher_direction": "co" if odds > 1 else ("mutex" if odds < 1 else "none"),
        })
        p_co.append(pc); p_me.append(pm)
        keep_co.append(floor_co < significance_threshold); keep_me.append(floor_me < significance_threshold)
    df = pd.DataFrame(recs)
    if len(df):
        # Correct within each direction over only the pairs testable in that direction.
        df["fisher_q_co_bh_local"] = _bh_subset(p_co, keep_co)
        df["fisher_q_mutex_bh_local"] = _bh_subset(p_me, keep_me)
        df["fisher_sig_co"] = df["fisher_q_co_bh_local"] < significance_threshold
        df["fisher_sig_mutex"] = df["fisher_q_mutex_bh_local"] < significance_threshold
        df["fisher_sig_any"] = df["fisher_sig_co"] | df["fisher_sig_mutex"]
    return df


def run_pearson(matrix, tested_units):
    """Pearson/phi correlation coefficient as a DESCRIPTIVE effect size only.

    Only the coefficient (= phi for binary data) and its sign are reported, without
    p-values. The sign is "none" when r is closer to zero than PEARSON_DIRECTION_TOL.
    The parametric Pearson p-value is unreliable for the small, sparse
    binary tables here: r reaches +/-1 by chance among rare mutations, which the
    normal approximation reads as extreme significance. Significance is decided
    by Fisher's exact test and DISCOVER, and pearson_r is kept purely as an
    effect-size metric.
    """
    sub, units = _tested(matrix, tested_units)
    arr = {u: sub.loc[u].values.astype(float) for u in units}
    recs = []
    for g1, g2 in combinations(units, 2):
        x1, x2 = arr[g1], arr[g2]
        k = _key(g1, g2)
        if x1.std() == 0 or x2.std() == 0:
            recs.append({"unit1": k[0], "unit2": k[1], "pearson_r": np.nan, "pearson_direction": "none"})
            continue
        r = float(np.corrcoef(x1, x2)[0, 1])
        recs.append({"unit1": k[0], "unit2": k[1], "pearson_r": r,
                     "pearson_direction": ("positive" if r >= PEARSON_DIRECTION_TOL else
                                           "negative" if r <= -PEARSON_DIRECTION_TOL else "none")})
    return pd.DataFrame(recs)


def correlation_matrix(matrix, tested_units):
    """Full Pearson/phi correlation matrix over the tested units (descriptive).

    Rows and columns are the tested units. Each cell is the Pearson correlation
    (= phi for binary presence/absence) between two units across samples, with a
    1.0 diagonal. Intended for a per-experiment heatmap / clustering view. Units
    with no variation yield NaN entries.
    """
    sub, units = _tested(matrix, tested_units)
    units = list(units)
    if len(units) < 1:
        return pd.DataFrame()
    X = sub.loc[units].values.astype(float)
    with np.errstate(invalid="ignore", divide="ignore"):
        C = np.atleast_2d(np.asarray(np.corrcoef(X), dtype=float))
    return pd.DataFrame(C, index=units, columns=units)
