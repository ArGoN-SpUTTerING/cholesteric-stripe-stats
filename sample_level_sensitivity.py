#!/usr/bin/env python3
"""
How much do the reported p-values depend on the choice of statistical unit?

The tests in stats_stripe_spacing.py treat each peak-to-peak distance as one
observation. Those distances are not independent: several are read off each
line, several lines are drawn per image, and all of the measurements for one
condition come from roughly 5-6 samples. This script asks what the same tests
give when the statistical unit is the sample rather than the measurement.

TWO PARTS
---------
1. Evidence of clustering, taken from the data alone.
   The lag-1 autocorrelation of each column. If the measurements were in
   random order this would be ~0 (one standard error is 1/sqrt(n), which
   differs by column because n does).

2. Sensitivity analysis.
   Samples were measured one at a time, so all measurements from one sample are
   contiguous in the file, and the samples contributed unequal numbers of
   measurements (both confirmed by the author; the first is also consistent
   with the autocorrelation in part 1). The positions of the boundaries are not
   recorded.

   Rather than guess them, the boundaries are drawn at random: each column is
   split into K contiguous blocks of random unequal size (each at least MIN
   measurements), the median of each block stands in for one sample, and the
   whole analysis is repeated N times. This reports the range of p-values
   consistent with what is actually known, instead of depending on one
   arbitrary choice of boundaries. Fitting the boundaries to the data would
   instead maximise the apparent differences between samples.

   K is swept because the number of samples per condition is known only
   approximately (5-6). The random seed is fixed, so the output is
   reproducible.

Run:  python3 sample_level_sensitivity.py
"""
import os
import random
import statistics as st

from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data", "pitch_series.txt")

COMP = ["DLPC", "1:0.43", "1:1", "1:2.3", "DOPC"]
CONC = ["2.8", "5", "10"]

MIN_BLOCK = 8      # fewest measurements a sample may contribute
N_DRAWS = 2000     # random partitions per K
SEED = 20260721    # fixed, so the output is reproducible


def load():
    rows = [ln.rstrip("\n").split("\t") for ln in open(DATA) if ln.strip()]
    cols = []
    for i in range(15):
        v = []
        for r in rows:
            if i < len(r):
                t = r[i].strip()
                if t and t != "--":
                    try:
                        v.append(float(t))
                    except ValueError:
                        pass
        cols.append(v)
    return {(c, k): cols[i * 3 + j]
            for i, c in enumerate(COMP) for j, k in enumerate(CONC)}


def acf1(v):
    m = st.mean(v)
    num = sum((v[i] - m) * (v[i + 1] - m) for i in range(len(v) - 1))
    return num / sum((x - m) ** 2 for x in v)


def random_blocks(v, K, min_block=None):
    """Median of each of K contiguous blocks of random, unequal size."""
    MIN_BLOCK_ = MIN_BLOCK if min_block is None else min_block
    n = len(v)
    if n < K * MIN_BLOCK_:
        return None
    extra = n - K * MIN_BLOCK_
    cuts = sorted(random.randint(0, extra) for _ in range(K - 1))
    edges = [0] + cuts + [extra]
    out, p = [], 0
    for i in range(K):
        size = MIN_BLOCK_ + edges[i + 1] - edges[i]
        out.append(st.median(v[p:p + size]))
        p += size
    return out


def main():
    if not os.path.exists(DATA):
        raise SystemExit(f"Measurements not found: {os.path.relpath(DATA, HERE)}")
    G = load()

    print("=" * 78)
    print("1.  ARE THE MEASUREMENTS INDEPENDENT?  lag-1 autocorrelation per column")
    print("=" * 78)
    print("  Each column against its own n. One standard error is 1/sqrt(n),")
    print("  so the multiple below is what matters, not the raw correlation.\n")
    print(f"  {'composition':<10}{'CB15':>7}{'n':>6}{'acf':>8}{'1/sqrt(n)':>11}{'x SE':>8}")
    mult = []
    for c in COMP:
        for k in CONC:
            v = G[(c, k)]
            a, se = acf1(v), 1 / len(v) ** 0.5
            mult.append(a / se)
            print(f"  {c:<10}{k:>7}{len(v):>6}{a:>8.3f}{se:>11.3f}{a/se:>8.2f}")
    print(f"\n  Autocorrelations span {min(acf1(G[(c, k)]) for c in COMP for k in CONC):.3f}"
          f"-{max(acf1(G[(c, k)]) for c in COMP for k in CONC):.3f};")
    print(f"  every column is at least {min(mult):.1f} standard errors from zero "
          f"(largest {max(mult):.1f}).")
    print("  The measurements are therefore strongly clustered in file order.")

    print("\n" + "=" * 78)
    print("2.  SENSITIVITY:  p-values when the unit is the sample, not the measurement")
    print("=" * 78)
    print(f"  Samples were measured one at a time and contributed unequal numbers")
    print(f"  of measurements, so each column is split into K contiguous blocks of")
    print(f"  random unequal size (>= {MIN_BLOCK} each). {N_DRAWS} random partitions per K.\n")
    random.seed(SEED)
    print(f"  {'K':>3}  {'test':<12}{'median p':>11}{'5th pct':>10}{'95th pct':>10}"
          f"{'p < 0.05':>10}{'at min':>9}")
    print("  " + "-" * 60)
    for K in range(4, 9):
        acc = {n: [] for n in ["KW worst", "MW 2.8", "MW 5", "MW 10"]}
        for _ in range(N_DRAWS):
            b = {key: random_blocks(v, K) for key, v in G.items()}
            if any(x is None for x in b.values()):
                break
            acc["KW worst"].append(
                max(stats.kruskal(*[b[(c, k)] for k in CONC])[1] for c in COMP))
            for k, nm in zip(CONC, ["MW 2.8", "MW 5", "MW 10"]):
                acc[nm].append(stats.mannwhitneyu(
                    b[("DLPC", k)], b[("DOPC", k)], alternative="two-sided")[1])
        if not acc["KW worst"]:
            print(f"  {K:>3}  (too few measurements for {K} blocks of {MIN_BLOCK})")
            continue
        for nm, vals in acc.items():
            vals.sort()
            q = lambda fr: vals[int(fr * len(vals))]
            frac = sum(1 for x in vals if x < 0.05) / len(vals)
            atmin = sum(1 for x in vals if x <= vals[0] * (1 + 1e-12)) / len(vals)
            print(f"  {K if nm=='KW worst' else '':>3}  {nm:<12}"
                  f"{st.median(vals):>11.4f}{q(.05):>10.4f}{q(.95):>10.4f}{frac:>9.0%}"
                  f"{atmin:>9.1%}")
        print()

    kw = max(stats.kruskal(*[G[(c, k)] for k in CONC])[1] for c in COMP)
    mw = [stats.mannwhitneyu(G[("DLPC", k)], G[("DOPC", k)],
                             alternative="two-sided")[1] for k in CONC]
    print(f"  For comparison, per measurement (what the manuscript reports):")
    print(f"       {'KW worst':<12}{kw:>11.1e}")
    for k, p in zip(CONC, mw):
        print(f"       {'MW ' + k:<12}{p:>11.1e}")

    print("\n" + "=" * 78)
    print("3.  WHAT THIS MEANS FOR THE TEXT")
    print("=" * 78)
    print("  Kruskal-Wallis (pitch controls stripe spacing)")
    print("    Significant in 100 % of random partitions at every K. At K = 5 the")
    print("    median p of 0.0025 is close to 0.0019, the smallest value attainable")
    print("    with five samples per group, so the blocks separate almost completely.")
    print("    The conclusion does not depend on the choice of statistical unit.")
    print()
    print("  Mann-Whitney (pure DLPC vs pure DOPC)")
    med = {k: (st.median(G[("DLPC", k)]), st.median(G[("DOPC", k)])) for k in CONC}
    for k in CONC:
        a, b = med[k]
        print(f"    {k:>4} wt.-%   medians {a:.3f} vs {b:.3f} um  "
              f"({100*(b/a-1):+.0f} %)")
    print("    At 2.8 wt.-% the two are so well separated that p < 0.05 in 100 %")
    print("    of partitions; the 'at min' column above gives the share sitting")
    print("    exactly on the smallest value attainable for that K.")
    print("    At 5 wt.-% the difference holds in most partitions. At 10 wt.-%")
    print("    it holds in few, and never in a majority for any setting tested")
    print("    (see part 4). The text therefore reports the difference for the")
    print("    two longer pitches only.")

    print("\n" + "=" * 78)
    print(f"4.  SENSITIVITY TO MIN_BLOCK (currently {MIN_BLOCK})")
    print("=" * 78)
    print("  MIN_BLOCK is a choice, not a measurement, so its effect is shown here.")
    print("  Percentage of 1000 random partitions reaching p < 0.05, at K = 5.\n")
    print(f"  {'MIN':>4}{'KW worst':>11}{'MW 2.8':>10}{'MW 5':>10}{'MW 10':>10}")
    for MB in (3, 5, 8, 12, 15):
        random.seed(SEED)
        hit = {"kw": 0, "0": 0, "1": 0, "2": 0}
        runs = 0
        for _ in range(1000):
            b = {key: random_blocks(v, 5, MB) for key, v in G.items()}
            if any(x is None for x in b.values()):
                break
            runs += 1
            if max(stats.kruskal(*[b[(c, k)] for k in CONC])[1]
                   for c in COMP) < 0.05:
                hit["kw"] += 1
            for idx, k in enumerate(CONC):
                if stats.mannwhitneyu(b[("DLPC", k)], b[("DOPC", k)],
                                      alternative="two-sided")[1] < 0.05:
                    hit[str(idx)] += 1
        if not runs:
            print(f"  {MB:>4}   too few measurements for 5 blocks of {MB}")
            continue
        print(f"  {MB:>4}" + "".join(f"{100*hit[x]/runs:>9.0f} %"
                                     for x in ("kw", "0", "1", "2")))
    print()
    print("  Kruskal-Wallis and the 2.8 wt.-% comparison are unaffected. The")
    print("  5 wt.-% comparison holds throughout. The 10 wt.-% figure does move")
    print("  with MIN_BLOCK, but stays far below a majority at every setting, and")
    print("  falls further as MIN_BLOCK rises towards the realistic value implied")
    print("  by 75-126 measurements spread over 5-6 samples.")


if __name__ == "__main__":
    main()
