# cholesteric-stripe-stats

Statistical analysis of cholesteric fingerprint stripe spacings at
lipid-decorated liquid crystal–water interfaces.

This repository holds the **analysis code**. The **raw measurements** are
deposited separately on Dataverse, so that they carry a citable, versioned DOI
of their own — see *Data availability* below. The associated manuscript is in
preparation; its citation will be added here on publication.

Nothing here is a custom statistical implementation. Both tests are standard
library functions, and every reported value is reproduced independently by
Python/SciPy, by base R, and by hand from the textbook formulae.

---

## Data availability

Raw stripe-spacing measurements: **Dataverse DOI — to be added on publication.**

Stripe spacings were measured manually in ImageJ as peak-to-peak distances
along lines drawn perpendicular to the fingerprint stripes. Only regions
containing identifiable periodic stripes were measured, so the number of
measurements differs between conditions.

| File | Shape | Contents |
|---|---|---|
| `data/pitch_series.txt` | 126 × 15 | 5 lipid compositions × CB15 = 2.8 / 5 / 10 wt.-%  (manuscript Fig. 3F) |
| `data/mixing_series.txt` | 135 × 12 | 3 DLPC:DOPC ratios × total lipid = 1 / 10 / 30 / 50 µM  (manuscript Fig. 2D) |

Download both files into `data/` and every script runs unchanged. The full
column layout and file format are documented in [`data/README.md`](data/README.md).

Micrographs were pre-processed with
[micrograph-batch](https://github.com/ArGoN-SpUTTerING/micrograph-batch).

---

## Statistical tests

Both tests are non-parametric — the spacing distributions are skewed and the
group sizes are unequal, so rank-based tests are more appropriate than
*t*-tests or ANOVA.

**Kruskal–Wallis**, across CB15 = 2.8 / 5 / 10 wt.-%, run separately for each
lipid composition. Tests whether the cholesteric pitch controls stripe spacing.

**Mann–Whitney *U*** (two-sided), pure DLPC vs pure DOPC, run separately at
each CB15 concentration. Tests whether acyl chain structure shifts the median
spacing at a fixed pitch.

Eight tests are reported in total. No multiplicity correction is applied
because none would change any conclusion: the largest *p* value below is
8.7 × 10⁻⁸, and a Bonferroni correction across all eight tests leaves it at
7 × 10⁻⁷.

### Results

| Kruskal–Wallis | H | *p* |
|---|---|---|
| DLPC | 310.54 | 3.7 × 10⁻⁶⁸ |
| DLPC:DOPC 1:0.43 | 255.41 | 3.5 × 10⁻⁵⁶ |
| DLPC:DOPC 1:1 | 279.52 | 2.0 × 10⁻⁶¹ |
| DLPC:DOPC 1:2.3 | 261.69 | 1.5 × 10⁻⁵⁷ |
| DOPC | 268.71 | 4.5 × 10⁻⁵⁹ |

| Mann–Whitney *U* (DLPC vs DOPC) | *U* | *p* |
|---|---|---|
| 2.8 wt.-% CB15 | 156.5 | 1.0 × 10⁻⁴⁰ |
| 5 wt.-% CB15 | 10093.5 | 7.9 × 10⁻¹⁵ |
| 10 wt.-% CB15 | 2809.0 | 8.7 × 10⁻⁸ |

These are per measurement. See *Choice of statistical unit* below.

The *U* statistic is directly interpretable. At 2.8 wt.-% CB15 there are
126 × 123 = 15 498 DLPC–DOPC pairs; DLPC is the larger of the pair in 156 of
them and tied in 1, giving *U* = 156 + 1/2 = 156.5.

---

## Reproducing the analysis

```bash
# Python (requires scipy)
python3 stats_stripe_spacing.py

# R (base R only, no packages needed)
Rscript stats_stripe_spacing.R

# check every quantitative claim in the manuscript against the data
python3 verify_claims.py

# p-values when the statistical unit is the sample, not the measurement
python3 sample_level_sensitivity.py
```

The outputs obtained from the deposited data are committed in `results/`, so
the reported numbers can be inspected without downloading anything.

`verify_claims.py` can additionally scan the manuscript source to confirm that
the wording of each quantitative statement still matches the data:

```bash
MANUSCRIPT_TEX=/path/to/manuscript.tex python3 verify_claims.py
```

---

## Independent cross-checks

The reported values do not depend on any single implementation:

| | Kruskal–Wallis *H* (DLPC) | Mann–Whitney *U* (2.8 wt.-%) |
|---|---|---|
| Python / SciPy 1.18 | 310.5411 | 156.5 |
| R 4.6.1, `kruskal.test` / `wilcox.test` | 310.5411 | 156.5 |
| Textbook formula, computed by hand | 310.5411 | 156.5 |

*p*-values agree to within the third significant digit; the small differences
come from how each implementation applies the tie correction and do not affect
any conclusion.

---

## Choice of statistical unit

The tests treat each peak-to-peak distance as one observation. Those distances
are not independent: several are read from each line, several lines are drawn
per image, and all measurements for one condition come from roughly 5–6
samples. `sample_level_sensitivity.py` quantifies what this costs, and reports
two things.

The clustering is visible in the data alone: the lag-1 autocorrelation of the
measurement columns is 0.28–0.79, against ~0 for a random ordering (one
standard error is 0.09).

Cutting each column into *K* equal contiguous blocks and using block medians as
sample values, the Kruskal–Wallis result is significant at every *K* from 3 to
12 — at *K* = 5 its *p* of 0.0019 is the smallest value attainable with five
samples per group. The Mann–Whitney comparison of pure DLPC against pure DOPC
holds at 2.8 and 5 wt.-% CB15 but not at 10 wt.-%, where the sample-level *p*
is 0.06–0.34. The manuscript reports that comparison for the two longer pitches
only.

## Distribution widths

Coefficients of variation are reported for every group by
`stats_stripe_spacing.py` and `stats_stripe_spacing.R`. In the pitch series at
10 wt.-% CB15 they are: DLPC 19.1 %, 1:0.43 21.4 %, 1:1 15.2 %, 1:2.3 15.3 %,
DOPC 32.0 %. In the mixing series the CV shows no monotonic trend with DOPC
fraction at any of the four lipid concentrations.

---

## Software

- SciPy — Virtanen et al., *Nat. Methods* **17**, 261–272 (2020)
- R — R Core Team, *R: A Language and Environment for Statistical Computing*
- ImageJ / Fiji — Schindelin et al., *Nat. Methods* **9**, 676–682 (2012)

## License

MIT (see `LICENSE`).
