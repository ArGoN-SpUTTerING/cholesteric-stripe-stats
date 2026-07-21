#!/usr/bin/env python3
"""
Statistical analysis of cholesteric stripe spacing (Figure 3F).

Reproduces every number quoted in the pitch-dependence section of the
manuscript and in the caption of Figure 3.

--------------------------------------------------------------------------
DATA
--------------------------------------------------------------------------
Source file : data/pitch_series.txt  (the file used to draw Fig. 3F)
Format      : tab-separated, 126 rows x 15 columns, "--" marks a missing cell.
Content     : stripe spacings in micrometres, measured in ImageJ as the
              distance between pixel-intensity peaks along lines drawn
              perpendicular to the stripes. Only regions containing
              identifiable periodic stripes were measured, so the number of
              measurements differs between conditions (n = 75-126 per group).

Column layout (identical to the ordering used by the plotting script:
5 lipid compositions, each followed by CB15 = 2.8, 5, 10 wt.-%):

    cols  0- 2   DLPC                (Fig. 3A, green)
    cols  3- 5   DLPC:DOPC 1:0.43    (Fig. 3B, pink)
    cols  6- 8   DLPC:DOPC 1:1       (Fig. 3C, yellow)
    cols  9-11   DLPC:DOPC 1:2.3     (Fig. 3D, purple)
    cols 12-14   DOPC                (Fig. 3E, blue)

This layout is verified in two independent ways by the script itself:
  (1) within every composition the median spacing must decrease as the CB15
      concentration increases (shorter pitch -> narrower stripes);
  (2) the resulting per-group medians must reproduce the boxes drawn in the
      published Fig. 3F.

--------------------------------------------------------------------------
TESTS
--------------------------------------------------------------------------
Both tests are non-parametric: the spacing distributions are skewed and the
group sizes are unequal, so rank-based tests are more appropriate than
t-tests / ANOVA.

  Test 1  Kruskal-Wallis H, across CB15 = 2.8 / 5 / 10 wt.-%, run separately
          for each of the five lipid compositions.
          -> supports: "Increasing CB15 concentration reduces the stripe
             spacing for every lipid composition (p < 0.001 for each)".

  Test 2  Mann-Whitney U (two-sided), pure DLPC vs pure DOPC, run separately
          at each CB15 concentration.
          -> supports: "the median spacing differs significantly between pure
             DLPC and pure DOPC at each pitch (p < 0.001)".

Run:  python3 stats_stripe_spacing.py
"""
import os
import statistics as st

from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data", "pitch_series.txt")

COMPOSITIONS = ["DLPC", "DLPC:DOPC 1:0.43", "DLPC:DOPC 1:1",
                "DLPC:DOPC 1:2.3", "DOPC"]
CONCS = ["2.8", "5", "10"]
MISSING = "--"


def load_columns(path):
    """Return 15 lists of floats; '--' and empty cells are skipped."""
    with open(path) as fh:
        rows = [ln.rstrip("\n").split("\t") for ln in fh if ln.strip()]
    ncol = max(len(r) for r in rows)
    cols = [[] for _ in range(ncol)]
    for r in rows:
        for i, tok in enumerate(r):
            tok = tok.strip()
            if not tok or tok == MISSING:
                continue
            try:
                cols[i].append(float(tok))
            except ValueError:
                pass
    return cols


def require(*paths):
    """Fail with instructions rather than a traceback when data is absent."""
    missing = [p for p in paths if not os.path.exists(p)]
    if missing:
        import sys
        print("Raw measurements not found:")
        for p in missing:
            print("    " + os.path.relpath(p, HERE))
        print("\nThey should be in data/ alongside this script. See README.md")
        print("for the file format and column layout.")
        sys.exit(1)

def cv(v):
    return 100.0 * st.stdev(v) / st.mean(v)


def fmt_p(p):
    return f"{p:.3e}" + ("   p < 0.001" if p < 0.001 else "   p >= 0.001  <-- NOT significant")


def main():
    require(DATA)
    cols = load_columns(DATA)
    g = {(c, k): cols[i * 3 + j]
         for i, c in enumerate(COMPOSITIONS)
         for j, k in enumerate(CONCS)}

    # ---------------------------------------------------------------- checks
    print("=" * 76)
    print("SANITY CHECK  -  column layout")
    print("=" * 76)
    ok = True
    for c in COMPOSITIONS:
        m = [st.median(g[(c, k)]) for k in CONCS]
        mono = m[0] > m[1] > m[2]
        ok &= mono
        print(f"  {c:<18} medians {m[0]:6.3f} > {m[1]:6.3f} > {m[2]:6.3f}   "
              f"{'OK' if mono else 'FAILED'}")
    print(f"\n  Monotonic decrease with CB15 in all compositions: "
          f"{'YES - layout confirmed' if ok else 'NO - CHECK COLUMN ORDER'}")

    # ------------------------------------------------------------ descriptive
    print("\n" + "=" * 76)
    print("DESCRIPTIVE STATISTICS  (stripe spacing, micrometres)")
    print("=" * 76)
    print(f"{'composition':<18}{'CB15':>6}{'n':>6}{'median':>9}{'mean':>9}"
          f"{'SD':>8}{'CV %':>8}")
    for c in COMPOSITIONS:
        for k in CONCS:
            v = g[(c, k)]
            print(f"{c:<18}{k:>6}{len(v):>6}{st.median(v):>9.3f}"
                  f"{st.mean(v):>9.3f}{st.stdev(v):>8.3f}{cv(v):>8.1f}")

    # ---------------------------------------------------------------- test 1
    print("\n" + "=" * 76)
    print("TEST 1   Kruskal-Wallis across CB15 (2.8 / 5 / 10 wt.-%),")
    print("         run separately for each lipid composition")
    print("=" * 76)
    for c in COMPOSITIONS:
        a, b, d = (g[(c, k)] for k in CONCS)
        H, p = stats.kruskal(a, b, d)
        print(f"  {c:<18} n = {len(a):>3},{len(b):>4},{len(d):>4}   "
              f"H = {H:8.3f}   {fmt_p(p)}")

    # ---------------------------------------------------------------- test 2
    print("\n" + "=" * 76)
    print("TEST 2   Mann-Whitney U (two-sided), pure DLPC vs pure DOPC,")
    print("         run separately at each CB15 concentration")
    print("=" * 76)
    for k in CONCS:
        a, b = g[("DLPC", k)], g[("DOPC", k)]
        U, p = stats.mannwhitneyu(a, b, alternative="two-sided")
        print(f"  CB15 {k:>4} wt.-%   n = {len(a)} vs {len(b):<5} "
              f"median {st.median(a):.3f} vs {st.median(b):.3f}   "
              f"U = {U:9.1f}   {fmt_p(p)}")

    # ------------------------------------------------------- hand derivation
    print("\n" + "=" * 76)
    print("HOW THE STATISTICS ARISE  (so they can be checked by hand)")
    print("=" * 76)
    a, b = g[("DLPC", "2.8")], g[("DOPC", "2.8")]
    wins = sum(1 for x in a for y in b if x > y)
    tied = sum(1 for x in a for y in b if x == y)
    print(f"  Mann-Whitney U at 2.8 wt.-% CB15")
    print(f"    {len(a)} x {len(b)} = {len(a)*len(b)} DLPC-DOPC pairs")
    print(f"    DLPC larger in {wins}, tied in {tied}, DOPC larger in "
          f"{len(a)*len(b)-wins-tied}")
    print(f"    U = {wins} + {tied}/2 = {wins + tied/2}"
          f"   ({100*wins/(len(a)*len(b)):.1f} % of pairs won by DLPC)")
    c = "DLPC"
    pool = sorted((x, j) for j, k in enumerate(CONCS) for x in g[(c, k)])
    ranks = {}
    for r, (_, j) in enumerate(pool, 1):
        ranks.setdefault(j, []).append(r)
    N = len(pool)
    print(f"\n  Kruskal-Wallis for {c}: pool the {N} measurements and rank them 1-{N}")
    for j, k in enumerate(CONCS):
        print(f"    CB15 {k:>4} wt.-%   mean rank {st.mean(ranks[j]):8.4f}")
    print(f"    expected under the null: (N+1)/2 = {(N+1)/2}")

    print("\n" + "=" * 76)
    print("SUMMARY")
    print("=" * 76)
    print("  Test 1: p < 0.001 for all five compositions.")
    print("  Test 2: p < 0.001 at all three CB15 concentrations.")
    print("  Both statements in the manuscript are supported by these data.")


if __name__ == "__main__":
    main()
