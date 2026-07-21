# cholesteric-stripe-stats

Non-parametric statistical analysis of cholesteric fingerprint stripe spacings
measured at lipid-decorated liquid crystal–water interfaces.

No custom statistics: every test is a standard library function, and the Python
and R versions are independent implementations of the same analysis, so running
both cross-checks the result.

## The scripts

**`stats_stripe_spacing.py`** and **`stats_stripe_spacing.R`** — the same two
tests, one in Python/SciPy and one in base R. A Kruskal–Wallis test across
chiral-dopant concentrations, run separately for each lipid composition, and a
two-sided Mann–Whitney *U* test between lipid compositions, run separately at
each concentration. Both are rank-based, because the spacing distributions are
skewed and the group sizes unequal. The scripts also print descriptive
statistics for every group and the pairwise decomposition behind *U*, so the
statistic can be checked by hand.

**`sample_level_sensitivity.py`** — the measurements are clustered: several
peak-to-peak distances are read from each line, several lines from each sample.
This script recomputes both tests with the sample as the statistical unit
instead. Because the sample boundaries are not recorded, it splits each column
into contiguous blocks of random unequal size and repeats the analysis over
many random partitions, sweeping both the number of blocks and the minimum
block size. The random seed is fixed, so the output is reproducible.

## Data

Not in this repository. The raw measurements are deposited on Dataverse so that
they carry a citable DOI of their own — **DOI to be added on publication**.

Download `pitch_series.txt` into `data/` and the scripts run unchanged;
[`data/README.md`](data/README.md) documents the file format and column layout.
The scripts report where to find the data if the directory is empty.

Measurements were made manually in ImageJ as peak-to-peak distances along lines
drawn perpendicular to the fingerprint stripes. Micrographs were pre-processed
with [micrograph-batch](https://github.com/ArGoN-SpUTTerING/micrograph-batch).

## Run

```bash
python3 stats_stripe_spacing.py       # requires scipy
Rscript stats_stripe_spacing.R        # base R only, no packages needed
python3 sample_level_sensitivity.py
```

The output of each is committed in `results/`, so the numbers can be read
without downloading the data.

## Citing

Li, M. *cholesteric-stripe-stats: statistical analysis of cholesteric stripe
spacings* (v1.0.0). Zenodo, 2026. DOI:
[10.5281/zenodo.21477799](https://doi.org/10.5281/zenodo.21477799)

See [`CITATION.cff`](CITATION.cff). The concept DOI
[10.5281/zenodo.21477798](https://doi.org/10.5281/zenodo.21477798) always
resolves to the latest version.

## License

MIT (see `LICENSE`).
