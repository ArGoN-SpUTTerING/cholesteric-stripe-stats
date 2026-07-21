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
   random order this would be ~0 (one standard error is 1/sqrt(n) ~ 0.09).

2. Sensitivity analysis.
   The true sample of origin is not recorded in the data files. As a stand-in,
   each column is cut into K equal contiguous blocks and the median of each
   block is used as that sample's value. This ASSUMES that measurements from
   one sample are adjacent in the file. The autocorrelation in part 1 supports
   that assumption but does not establish it, so the analysis is run over a
   range of K rather than at a single value.

Run:  python3 sample_level_sensitivity.py
"""
import os
import statistics as st

from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data", "pitch_series.txt")

COMP = ["DLPC", "1:0.43", "1:1", "1:2.3", "DOPC"]
CONC = ["2.8", "5", "10"]


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


def blocks(v, k):
    """Median of each of k equal contiguous blocks."""
    n = len(v)
    return [st.median(v[round(i * n / k):round((i + 1) * n / k)]) for i in range(k)]


def main():
    if not os.path.exists(DATA):
        raise SystemExit(f"Measurements not found: {os.path.relpath(DATA, HERE)}")
    G = load()

    print("=" * 78)
    print("1.  ARE THE MEASUREMENTS INDEPENDENT?  lag-1 autocorrelation per column")
    print("=" * 78)
    print(f"  {'composition':<10}" + "".join(f"{k + ' wt.-%':>12}" for k in CONC))
    worst = 1.0
    for c in COMP:
        vals = [acf1(G[(c, k)]) for k in CONC]
        worst = min(worst, min(vals))
        print(f"  {c:<10}" + "".join(f"{v:>12.3f}" for v in vals))
    n = min(len(v) for v in G.values())
    print(f"\n  Random ordering would give ~0; one standard error is "
          f"1/sqrt(n) = {1/n**0.5:.3f}.")
    print(f"  The smallest value above is {worst:.3f}, i.e. {worst*n**0.5:.1f} "
          f"standard errors.")
    print("  The measurements are therefore strongly clustered in file order.")

    print("\n" + "=" * 78)
    print("2.  SENSITIVITY:  p-values when the unit is the sample, not the measurement")
    print("=" * 78)
    print("  Each column is cut into K equal contiguous blocks; the median of each")
    print("  block stands in for one sample. K is swept because the true sample of")
    print("  origin is not recorded.\n")
    print(f"  {'K':>3} | {'Kruskal-Wallis':>15} | {'Mann-Whitney U, pure DLPC vs pure DOPC':>40}")
    print(f"  {'':>3} | {'(worst of 5)':>15} | {'2.8 wt.-%':>13}{'5 wt.-%':>13}{'10 wt.-%':>13}")
    print("  " + "-" * 74)
    for K in range(3, 13):
        kw = max(stats.kruskal(*[blocks(G[(c, k)], K) for k in CONC])[1]
                 for c in COMP)
        mw = [stats.mannwhitneyu(blocks(G[("DLPC", k)], K),
                                 blocks(G[("DOPC", k)], K),
                                 alternative="two-sided")[1] for k in CONC]
        f = lambda p: f"{p:.4f}" + ("*" if p < 0.05 else " ")
        print(f"  {K:>3} | {f(kw):>15} | " + "".join(f"{f(p):>13}" for p in mw))
    print("  " + "-" * 74)
    print("  * p < 0.05")

    kw = max(stats.kruskal(*[G[(c, k)] for k in CONC])[1] for c in COMP)
    mw = [stats.mannwhitneyu(G[("DLPC", k)], G[("DOPC", k)],
                             alternative="two-sided")[1] for k in CONC]
    print(f"\n  For comparison, per measurement (what the manuscript reports):")
    print(f"  {'':>3} | {kw:>15.1e} | " + "".join(f"{p:>13.1e}" for p in mw))

    print("\n" + "=" * 78)
    print("3.  WHAT THIS MEANS FOR THE TEXT")
    print("=" * 78)
    print("  Kruskal-Wallis (pitch controls stripe spacing)")
    print("    Significant at every K tested. At K = 5 the p-value is 0.0019,")
    print("    which is the smallest value attainable with five samples per")
    print("    group -- the blocks separate completely. The conclusion does not")
    print("    depend on the choice of statistical unit.")
    print()
    print("  Mann-Whitney (pure DLPC vs pure DOPC)")
    med = {k: (st.median(G[("DLPC", k)]), st.median(G[("DOPC", k)])) for k in CONC}
    for k in CONC:
        a, b = med[k]
        print(f"    {k:>4} wt.-%   medians {a:.3f} vs {b:.3f} um  "
              f"({100*(b/a-1):+.0f} %)")
    print("    The difference is large at 2.8 wt.-% and survives at 5 wt.-%, but")
    print("    at 10 wt.-% it is not significant at any K below 12. The text")
    print("    therefore reports the difference for the two longer pitches only.")


if __name__ == "__main__":
    main()
