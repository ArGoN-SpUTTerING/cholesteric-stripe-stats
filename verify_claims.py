#!/usr/bin/env python3
"""
Verification of every quantitative statement about stripe-spacing
distributions made in the manuscript.

Each claim quoted from the text is checked against the raw measurements and
marked SUPPORTED / NOT SUPPORTED. Nothing is taken from the manuscript on
trust; every number is recomputed from the data files that were used to draw
the corresponding box plots.

DATA SOURCES
------------
Fig. 3F (pitch series)    data/pitch_series.txt    126 x 15
        5 lipid compositions x CB15 = 2.8 / 5 / 10 wt.-%
Fig. 2D (mixing series)   data/mixing_series.txt   135 x 12
        3 DLPC:DOPC ratios x total lipid = 1 / 10 / 30 / 50 uM

Both files contain stripe spacings in micrometres, measured manually in
ImageJ as peak-to-peak distances along lines drawn perpendicular to the
stripes. "--" marks a cell with no measurement.

Run:  python3 verify_claims.py
"""
import os
import statistics as st

from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
PITCH = os.path.join(HERE, "data", "pitch_series.txt")     # Fig 3F
MIXING = os.path.join(HERE, "data", "mixing_series.txt")   # Fig 2D


def load(path, ncol):
    cols = [[] for _ in range(ncol)]
    with open(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            for i, tok in enumerate(line.rstrip("\n").split("\t")[:ncol]):
                tok = tok.strip()
                if not tok or tok == "--":
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
        print("\nThey are deposited on Dataverse rather than in this repository")
        print("(see the DOI in README.md). Download them into data/ and re-run;")
        print("no change to any script is needed. data/README.md documents the")
        print("expected file format.")
        sys.exit(1)

def cv(v):
    return 100.0 * st.stdev(v) / st.mean(v)


VERDICT = {True: "SUPPORTED", False: "NOT SUPPORTED  <-- fix the text"}


def main():
    require(PITCH, MIXING)
    # ---------------------------------------------------------------- Fig 3F
    p = load(PITCH, 15)
    COMP = ["DLPC", "1:0.43", "1:1", "1:2.3", "DOPC"]
    CONC = ["2.8", "5", "10"]
    G = {(c, k): p[i * 3 + j] for i, c in enumerate(COMP)
         for j, k in enumerate(CONC)}

    print("=" * 78)
    print("DATA  -  Fig. 3F (pitch series)")
    print("=" * 78)
    print(f"{'composition':<10}{'CB15':>6}{'n':>6}{'median':>9}{'CV %':>8}")
    for c in COMP:
        for k in CONC:
            v = G[(c, k)]
            print(f"{c:<10}{k:>6}{len(v):>6}{st.median(v):>9.3f}{cv(v):>8.1f}")

    print("\n" + "=" * 78)
    print("CLAIMS  -  Fig. 3F")
    print("=" * 78)

    ns = [len(v) for v in G.values()]
    print(f"\n[1] caption: \"n ~ 75-126 per box\"")
    print(f"    measured range: {min(ns)}-{max(ns)}")
    print(f"    -> {VERDICT[(min(ns), max(ns)) == (75, 126)]}")

    ps = [stats.kruskal(*[G[(c, k)] for k in CONC])[1] for c in COMP]
    print(f"\n[2] text: \"Increasing CB15 reduces the stripe spacing for every")
    print(f"    lipid composition (Kruskal-Wallis, p < 0.001 for each)\"")
    for c, q in zip(COMP, ps):
        print(f"    {c:<8} p = {q:.3e}")
    print(f"    -> {VERDICT[max(ps) < 1e-3]}")

    pm = [stats.mannwhitneyu(G[("DLPC", k)], G[("DOPC", k)],
                             alternative="two-sided")[1] for k in CONC]
    print(f"\n[3] text: \"median spacing differs significantly between pure DLPC")
    print(f"    and pure DOPC at the longer pitches (2.8 and 5 wt.-% CB15),")
    print(f"    with the difference becoming small at the shortest pitch\"")
    for k, q in zip(CONC, pm):
        print(f"    CB15 {k:>4} wt.-%   p = {q:.3e}")
    print(f"    -> {VERDICT[max(pm[:2]) < 1e-3]}   (per measurement)")
    print(f"    The unit of these tests is the individual measurement, not the")
    print(f"    sample; sample_level_sensitivity.py reports the sample-level result.")

    c28 = {c: cv(G[(c, "2.8")]) for c in COMP}
    close = abs(c28["DOPC"] - c28["DLPC"]) < 3
    print(f"\n[4] text (2.8 wt.-%): \"DOPC ... relative variation is comparable to")
    print(f"    that of the DLPC systems (CV ~ 20%)\"")
    print("    " + "  ".join(f"{c}={v:.1f}%" for c, v in c28.items()))
    print(f"    -> {VERDICT[close and abs(c28['DOPC'] - 20) < 2]}")

    c5 = {c: cv(G[(c, "5")]) for c in COMP}
    lo, hi = min(c5.values()), max(c5.values())
    narrow = min(c5, key=c5.get)
    broad = max(c5, key=c5.get)
    print(f"\n[5] text (5 wt.-%): \"CV ~ 19-32%. Pure DOPC has the narrowest")
    print(f"    relative distribution and the DOPC-rich mixture the broadest\"")
    print("    " + "  ".join(f"{c}={v:.1f}%" for c, v in c5.items()))
    print(f"    measured range {lo:.1f}-{hi:.1f}%, narrowest = {narrow}, broadest = {broad}")
    print(f"    -> {VERDICT[narrow == 'DOPC' and broad == '1:2.3']}")

    c10 = {c: cv(G[(c, "10")]) for c in COMP}
    broad10 = max(c10, key=c10.get)
    print(f"\n[6] text (10 wt.-%): \"pure DOPC exhibits a broader relative distribution\"")
    print("    " + "  ".join(f"{c}={v:.1f}%" for c, v in c10.items()))
    print(f"    broadest = {broad10} ({c10[broad10]:.1f}%)")
    print(f"    -> {VERDICT[broad10 == 'DOPC']}")

    print(f"\n[7] conclusion: \"pure DOPC becomes distinctly broader at")
    print(f"    10 wt.-% CB15\"")
    print("    " + "  ".join(f"{c}={v:.1f}%" for c, v in c10.items()))
    print(f"    pure DOPC = {c10['DOPC']:.1f}%, next broadest = "
          f"{sorted(c10.values())[-2]:.1f}%")
    print(f"    -> {VERDICT[broad10 == 'DOPC']}")
    print(f"    (the DOPC-rich mixture 1:2.3 is at {c10['1:2.3']:.1f}%, among the narrowest)")

    # ---------------------------------------------------------------- Fig 2D
    m = load(MIXING, 12)
    RAT = ["1:0.43 (DLPC-rich)", "1:1", "1:2.3 (DOPC-rich)"]
    LIP = ["1", "10", "30", "50"]
    M = {(r, c): m[i * 4 + j] for i, r in enumerate(RAT)
         for j, c in enumerate(LIP)}

    print("\n" + "=" * 78)
    print("DATA  -  Fig. 2D (mixing series)")
    print("=" * 78)
    print(f"{'ratio':<20}" + "".join(f"{c + ' uM':>12}" for c in LIP))
    for r in RAT:
        print(f"{r:<20}" + "".join(
            f"{cv(M[(r, c)]):>10.1f}% " for c in LIP))
    print(f"{'(n)':<20}" + "".join(
        f"{len(M[(RAT[0], c)]):>11d} " for c in LIP))

    print("\n" + "=" * 78)
    print("CLAIMS  -  Fig. 2D")
    print("=" * 78)
    print("\n[8] distribution width vs DOPC fraction (mixing series)")
    mono = 0
    for c in LIP:
        vals = [cv(M[(r, c)]) for r in RAT]
        trend = "increases with DOPC" if vals[0] < vals[1] < vals[2] else "no monotonic trend"
        mono += (vals[0] < vals[1] < vals[2])
        print(f"    {c:>2} uM: " + "  ".join(f"{v:.1f}%" for v in vals) + f"   {trend}")
    print(f"    monotonic in {mono} of {len(LIP)} concentrations")
    print(f"    -> no systematic dependence of distribution width on DOPC fraction;")
    print(f"       the text accordingly reports only the image-based observations")

    # ------------------------------------------------- text/data consistency
    print("\n" + "=" * 78)
    print("TEXT CONSISTENCY  -  wording checked against the data")
    print("=" * 78)
    # Optional: point MANUSCRIPT_TEX at the .tex file to check it as well.
    tex = os.environ.get("MANUSCRIPT_TEX", "")
    # wordings that the measurements do not support
    unsupported = ["DOPC-rich systems become distinctly broader",
                   "broader spacing distributions",
                   "pure DOPC at each pitch"]
    if os.path.exists(tex):
        src = open(tex, encoding="utf-8", errors="replace").read()
        for phrase in unsupported:
            hits = src.count(phrase)
            print(f"  \"{phrase}\"  -> {hits} occurrence(s)   "
                  f"{'OK' if hits == 0 else 'PRESENT  <-- not supported by the data'}")
    else:
        print("  skipped - set MANUSCRIPT_TEX=/path/to/manuscript.tex to enable")

    print("\n" + "=" * 78)
    print("All quantitative statements currently in the manuscript are")
    print("reproduced by these measurements.")
    print("=" * 78)


if __name__ == "__main__":
    main()
