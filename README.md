# cholesteric-stripe-stats

Statistical analysis of cholesteric fingerprint stripe spacings at
lipid-decorated liquid crystal–water interfaces.

This repository holds the **analysis code**. The **raw measurements** are
deposited separately on Dataverse, so that they carry a citable, versioned DOI
of their own. The associated manuscript is in preparation; its citation will be
added here on publication.

Neither test is a custom implementation, and every reported value is reproduced
independently by Python/SciPy and by base R.

---

## Data availability

Raw stripe-spacing measurements: **Dataverse DOI — to be added on publication.**

Download `pitch_series.txt` into `data/` and the scripts run unchanged. The
column layout and file format are documented in
[`data/README.md`](data/README.md).

Stripe spacings were measured manually in ImageJ as peak-to-peak distances
along lines drawn perpendicular to the fingerprint stripes. Only regions
containing identifiable periodic stripes were measured, so the number of
measurements differs between conditions (75–126 per group). Micrographs were
pre-processed with
[micrograph-batch](https://github.com/ArGoN-SpUTTerING/micrograph-batch).

---

## The tests

Both are non-parametric — the distributions are skewed and the group sizes
unequal, so rank-based tests suit them better than *t*-tests or ANOVA.

**Kruskal–Wallis** across CB15 = 2.8 / 5 / 10 wt.-%, run separately for each
lipid composition. Tests whether the cholesteric pitch controls stripe spacing.

**Mann–Whitney *U*** (two-sided), pure DLPC vs pure DOPC, run separately at
each CB15 concentration. Tests whether acyl chain structure shifts the median
spacing at a fixed pitch.

| Kruskal–Wallis | *H* | *p* |
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

Eight tests in total; a Bonferroni correction across all eight would leave the
largest *p* at about 7 × 10⁻⁷, so none is applied.

*U* is countable by hand, and the scripts print the decomposition: at
2.8 wt.-% CB15 there are 126 × 123 = 15 498 DLPC–DOPC pairs, DLPC is the larger
in 156 of them and tied in 1, giving *U* = 156 + 1/2 = 156.5.

Python and R agree on both test statistics and on seven of the eight
*p*-values. The eighth differs in the third digit, from rounding in R's tail
computation; both are of order 10⁻¹⁵.

---

## Choice of statistical unit

These tests treat each peak-to-peak distance as one observation, but the
distances are not independent: several are read from each line, several lines
from each sample, and each condition draws on roughly 5–6 samples.

`sample_level_sensitivity.py` recomputes both tests with the sample as the
unit. Samples were measured one at a time, so a sample's measurements are
contiguous in the file, though the boundaries are not recorded; the script
therefore splits each column into contiguous blocks of random unequal size and
repeats the analysis over many random partitions, sweeping both the number of
blocks and the minimum block size rather than fixing either.

The Kruskal–Wallis result is unaffected — significant in every partition
tested. The Mann–Whitney comparison holds at 2.8 and 5 wt.-% CB15 but not at
10 wt.-%, where the sample-level *p* is around 0.15. The manuscript reports
that comparison for the two longer pitches only.

---

## Reproducing the analysis

```bash
python3 stats_stripe_spacing.py     # requires scipy
Rscript stats_stripe_spacing.R      # base R only, no packages needed
python3 sample_level_sensitivity.py # p-values if the unit is the sample
```

The output of each is committed in `results/`, so the numbers can be read
without downloading the data.

---

## Software

- SciPy — Virtanen et al., *Nat. Methods* **17**, 261–272 (2020)
- R — R Core Team, *R: A Language and Environment for Statistical Computing*
- ImageJ / Fiji — Schindelin et al., *Nat. Methods* **9**, 676–682 (2012)

## License

MIT (see `LICENSE`).
