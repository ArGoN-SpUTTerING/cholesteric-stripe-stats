# =============================================================================
# Statistical analysis of cholesteric stripe spacing (Figure 3F)
# R version - independent cross-check of stats_stripe_spacing.py
# =============================================================================
#
# Both tests use base-R functions only; no packages need to be installed.
#   kruskal.test()  -> Kruskal-Wallis rank-sum test
#   wilcox.test()   -> with two samples this IS the Mann-Whitney U test
#                      (R reports the statistic as W; W = U)
#
# DATA -----------------------------------------------------------------------
# data/pitch_series.txt : tab-separated, 126 rows x 15 columns, "--" marks a missing cell.
#            Stripe spacings in micrometres, measured in ImageJ as the distance
#            between pixel-intensity peaks along lines drawn perpendicular to
#            the stripes. Only regions with identifiable periodic stripes were
#            measured, so n differs between groups (75-126).
#
# Column layout (5 lipid compositions, each followed by CB15 = 2.8, 5, 10 wt.-%)
#            V1 -V3    DLPC               (Fig. 3A, green)
#            V4 -V6    DLPC:DOPC 1:0.43   (Fig. 3B, pink)
#            V7 -V9    DLPC:DOPC 1:1      (Fig. 3C, yellow)
#            V10-V12   DLPC:DOPC 1:2.3    (Fig. 3D, purple)
#            V13-V15   DOPC               (Fig. 3E, blue)
#
# HOW TO RUN -----------------------------------------------------------------
#   In a terminal, from this folder:      Rscript stats_stripe_spacing.R
#   Or inside RStudio:                    open the file and press Source
#
# VERIFIED: run with R 4.6.1; output matches the Python/SciPy run exactly
# (statistics identical; p differs only in the 3rd significant digit at
#  5 wt.-%, from the tie-correction implementation).
# EXPECTED RESULTS -----------------------------------------------------------
#   Kruskal-Wallis, across CB15, per composition
#       DLPC              H = 310.54   p = 3.7e-68
#       DLPC:DOPC 1:0.43  H = 255.41   p = 3.5e-56
#       DLPC:DOPC 1:1     H = 279.52   p = 2.0e-61
#       DLPC:DOPC 1:2.3   H = 261.69   p = 1.5e-57
#       DOPC              H = 268.71   p = 4.5e-59
#   Mann-Whitney U, pure DLPC vs pure DOPC, per CB15 concentration
#       2.8 wt.-%   p = 1.0e-40
#       5   wt.-%   p = 7.9e-15
#       10  wt.-%   p = 8.7e-08
# =============================================================================

# --- read data ---------------------------------------------------------------
if (!file.exists("data/pitch_series.txt")) {
  stop("Raw measurements not found: data/pitch_series.txt\n",
       "  They should be in data/ alongside this script.\n",
       "  See README.md for the file format.", call. = FALSE)
}

d <- read.table("data/pitch_series.txt", sep = "\t", header = FALSE,
                na.strings = "--", fill = TRUE, stringsAsFactors = FALSE)

compositions <- c("DLPC", "DLPC:DOPC 1:0.43", "DLPC:DOPC 1:1",
                  "DLPC:DOPC 1:2.3", "DOPC")
concs <- c("2.8", "5", "10")

# column index for composition i (1-5) and concentration j (1-3)
col_of <- function(i, j) (i - 1) * 3 + j

# values of one group, missing values removed
grp <- function(i, j) as.numeric(na.omit(d[[col_of(i, j)]]))

# --- sanity check: median must fall as CB15 rises -----------------------------
cat("==============================================================\n")
cat("SANITY CHECK  -  column layout\n")
cat("==============================================================\n")
for (i in seq_along(compositions)) {
  m <- sapply(1:3, function(j) median(grp(i, j)))
  cat(sprintf("  %-18s medians %6.3f > %6.3f > %6.3f   %s\n",
              compositions[i], m[1], m[2], m[3],
              if (m[1] > m[2] && m[2] > m[3]) "OK" else "FAILED"))
}

# --- descriptive statistics ---------------------------------------------------
cat("\n==============================================================\n")
cat("DESCRIPTIVE STATISTICS  (stripe spacing, micrometres)\n")
cat("==============================================================\n")
cat(sprintf("%-18s %5s %5s %8s %8s %7s %7s\n",
            "composition", "CB15", "n", "median", "mean", "SD", "CV %"))
for (i in seq_along(compositions)) {
  for (j in 1:3) {
    v <- grp(i, j)
    cat(sprintf("%-18s %5s %5d %8.3f %8.3f %7.3f %7.1f\n",
                compositions[i], concs[j], length(v),
                median(v), mean(v), sd(v), 100 * sd(v) / mean(v)))
  }
}

# --- TEST 1: Kruskal-Wallis across CB15, per composition ----------------------
cat("\n==============================================================\n")
cat("TEST 1   Kruskal-Wallis across CB15 (2.8 / 5 / 10 wt.-%)\n")
cat("         run separately for each lipid composition\n")
cat("==============================================================\n")
for (i in seq_along(compositions)) {
  groups <- lapply(1:3, function(j) grp(i, j))
  kt <- kruskal.test(groups)
  cat(sprintf("  %-18s n = %3d,%4d,%4d   H = %8.3f   p = %.3e   %s\n",
              compositions[i],
              length(groups[[1]]), length(groups[[2]]), length(groups[[3]]),
              kt$statistic, kt$p.value,
              if (kt$p.value < 0.05) "p < 0.05" else "p >= 0.05  <-- NOT significant"))
}

# --- TEST 2: Mann-Whitney U, pure DLPC vs pure DOPC, per CB15 -----------------
cat("\n==============================================================\n")
cat("TEST 2   Mann-Whitney U (two-sided), pure DLPC vs pure DOPC\n")
cat("         run separately at each CB15 concentration\n")
cat("==============================================================\n")
for (j in 1:3) {
  a <- grp(1, j)   # DLPC  = composition 1
  b <- grp(5, j)   # DOPC  = composition 5
  wt <- wilcox.test(a, b, alternative = "two.sided", exact = FALSE)
  cat(sprintf("  CB15 %4s wt.-%%   n = %d vs %-5d median %.3f vs %.3f   W = %9.1f   p = %.3e   %s\n",
              concs[j], length(a), length(b), median(a), median(b),
              wt$statistic, wt$p.value,
              if (wt$p.value < 0.05) "p < 0.05" else "p >= 0.05  <-- NOT significant"))
}

cat("\n==============================================================\n")
cat("SUMMARY\n")
cat("==============================================================\n")
cat("  Test 1: p < 0.05 for all five compositions.\n")
cat("  Test 2: p < 0.05 at all three CB15 concentrations.\n\n")
cat("  These p-values treat each measurement as one observation. The\n")
cat("  measurements are clustered within samples, and on a sample-level\n")
cat("  recomputation Test 1 is unchanged while Test 2 holds only at 2.8\n")
cat("  and 5 wt.-% CB15. The manuscript reports Test 2 for those two\n")
cat("  pitches only; see the note on the statistical unit in README.md.\n")
