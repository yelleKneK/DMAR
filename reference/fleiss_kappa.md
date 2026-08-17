# Fleiss's Kappa for Inter-Rater Agreement Among Multiple Raters

Computes Fleiss's (1971) kappa coefficient of agreement among \\m \ge
2\\ raters who classify each of \\N\\ subjects into one of \\k\\ nominal
categories. Returns the point estimate together with the asymptotic
standard error and the Wald confidence interval, plus a test of \\H_0\\:
\kappa = 0\\.

## Usage

``` r
fleiss_kappa(
  ratings,
  conf_level = 0.95,
  ci_method = c("wald", "percentile", "bca"),
  B = 10000L,
  seed = NULL
)
```

## Arguments

- ratings:

  A numeric \\N \times k\\ matrix or `data.frame`: row \\i\\ gives the
  counts of raters who assigned each of the \\k\\ categories to subject
  \\i\\. Each row must sum to the same value \\m\\ (the common number of
  raters per subject).

- conf_level:

  Confidence level for the interval (default `0.95`).

- ci_method:

  Interval method: `"wald"` (the default, the asymptotic interval from
  the Gwet (2008) linearization variance), `"percentile"` (bootstrap
  percentile), or `"bca"` (bootstrap bias-corrected and accelerated).
  The multirater kappa variance literature is unsettled, and Zapf,
  Castell, Morawietz, and Karch (2016) recommend bootstrap intervals in
  this setting: the subjects (rows) are resampled with replacement `B`
  times, kappa is recomputed on each resample, and the interval is read
  off the bootstrap distribution; the BCa variant additionally adjusts
  the quantile positions for median bias and acceleration (Efron &
  Tibshirani, 1993). The `se`, `z_value`, and `p_value` columns keep
  their asymptotic definitions under every `ci_method`; only the
  interval changes.

- B:

  Number of bootstrap replications when `ci_method` is `"percentile"` or
  `"bca"` (default `10000`; ignored for `"wald"`).

- seed:

  Optional integer seed for the bootstrap. The default `NULL` uses the
  current state of the random number generator; a supplied seed is set
  internally and the prior state restored on exit.

## Value

A one-row `data.frame` (class `dmar_tbl`) with columns `kappa`, `se`
(asymptotic standard error of \\\hat\kappa_F\\ used for the interval),
`lower_limit`, `upper_limit`, `z_value`, `p_value` (Wald test of
\\H_0\\: \kappa = 0\\), `n_subjects`, `n_raters` (\\m\\), and
`n_categories` (\\k\\).

## Details

For \\n\_{ij}\\ = the number of raters who assigned subject \\i\\ to
category \\j\\, with \\\sum_j n\_{ij} = m\\ for every \\i\\, define the
marginal proportion of category \\j\\ as \\p_j = \sum_i n\_{ij} /
(Nm)\\, and the per-subject agreement \$\$P_i =
\frac{1}{m(m-1)}\Bigl(\sum_j n\_{ij}^2 - m\Bigr).\$\$ Then Fleiss's
kappa is \$\$\hat\kappa_F = \frac{\bar P - P_e}{1 - P_e}, \qquad \bar P
= \frac{1}{N}\sum_i P_i, \qquad P_e = \sum_j p_j^2.\$\$

**Standard error.** Two variances are involved, because the variance of
\\\hat\kappa_F\\ under \\H_0\\: \kappa = 0\\ is not its variance at a
nonzero value. The test of no agreement uses the null variance of
Fleiss, Nee, and Landis (1979, Equation 12), who corrected the standard
errors given in Fleiss (1971), \$\$\mathrm{Var}\_0(\hat\kappa_F) =
\frac{2\\\bigl(P_e + P_e^2 - 2\sum_j
p_j^3\bigr)}{N\\m\\(m-1)\\(1-P_e)^2},\$\$ and the reported \\z\\
statistic and *p*-value come from it. On the Fleiss (1971) Table 1
example below this gives \\z = 17.65\\, matching `irr::kappam.fleiss`.
The confidence interval instead uses the linearization variance of Gwet
(2008, Section 6), which is consistent at the estimated
\\\hat\kappa_F\\: each subject \\i\\ contributes an influence value
\\\kappa_i^\ast\\ (Gwet's Equations 34 and 35), and
\\\mathrm{Var}(\hat\kappa_F) = \sum_i (\kappa_i^\ast - \hat\kappa_F)^2 /
\\N(N-1)\\\\, Gwet's Equation 33 with the sampling fraction set to zero.
(Gwet derives the variance for the multiple-rater pi statistic, which is
the same estimator as Fleiss's kappa.) The Wald confidence interval is
\\\hat\kappa_F \pm z\_{1-\alpha/2}\\\widehat{\mathrm{SE}}\\, with the
upper limit truncated at 1. Using the null variance for the interval
would understate the standard error and give a spuriously narrow
interval.

**Bootstrap interval.** The variance of multirater kappa is unsettled in
the literature, and Zapf, Castell, Morawietz, and Karch (2016) recommend
a bootstrap interval in this setting. With `ci_method = "percentile"` or
`"bca"` the subjects (the rows of `ratings`) are resampled with
replacement `B` times, kappa is recomputed on each resample, and the
interval is read off the bootstrap distribution: the percentile interval
takes the empirical quantiles, and the BCa interval adjusts the quantile
positions for median bias and for acceleration (Efron & Tibshirani,
1993). Ask for it when \\N\\ is small or \\\hat\kappa_F\\ is near a
boundary, where the Wald interval's coverage is least dependable. The
`se`, `z_value`, and `p_value` columns keep their asymptotic definitions
under every `ci_method`; only the interval changes. Bootstrap results
vary from run to run; supply `seed` for reproducibility.

Fleiss's kappa is purely nominal (no weighting). For ordinal categories
with two raters, use
[`cohen_kappa`](https://yelleknek.github.io/DMAR/reference/cohen_kappa.md)
with quadratic weights; for ordinal categories with three or more
raters, an extension based on the intraclass correlation
([`icc`](https://yelleknek.github.io/DMAR/reference/icc.md)) is more
appropriate.

## References

Fleiss, J. L. (1971). Measuring nominal scale agreement among many
raters. *Psychological Bulletin, 76*(5), 378–382.

Fleiss, J. L., Nee, J. C. M., & Landis, J. R. (1979). Large sample
variance of kappa in the case of different sets of raters.
*Psychological Bulletin, 86*(5), 974–977.

Efron, B., & Tibshirani, R. J. (1993). *An introduction to the
bootstrap*. New York, NY: Chapman & Hall/CRC.

Gwet, K. L. (2008). Computing inter-rater reliability and its variance
in the presence of high agreement. *British Journal of Mathematical and
Statistical Psychology, 61*(1), 29–48.
[doi:10.1348/000711006X126600](https://doi.org/10.1348/000711006X126600)

Zapf, A., Castell, S., Morawietz, L., & Karch, A. (2016). Measuring
inter-rater reliability for nominal data: Which coefficients and
confidence intervals are appropriate? *BMC Medical Research Methodology,
16*, 93.
[doi:10.1186/s12874-016-0200-9](https://doi.org/10.1186/s12874-016-0200-9)

## See also

[`cohen_kappa`](https://yelleknek.github.io/DMAR/reference/cohen_kappa.md),
[`icc`](https://yelleknek.github.io/DMAR/reference/icc.md)

Other reliability:
[`cohen_kappa()`](https://yelleknek.github.io/DMAR/reference/cohen_kappa.md),
[`diagnosis_agreement`](https://yelleknek.github.io/DMAR/reference/diagnosis_agreement.md),
[`icc()`](https://yelleknek.github.io/DMAR/reference/icc.md),
[`reliability()`](https://yelleknek.github.io/DMAR/reference/reliability.md),
[`reliability_H()`](https://yelleknek.github.io/DMAR/reference/reliability_H.md),
[`reliability_alpha()`](https://yelleknek.github.io/DMAR/reference/reliability_alpha.md),
[`reliability_kr20()`](https://yelleknek.github.io/DMAR/reference/reliability_kr20.md),
[`reliability_omega()`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md),
[`reliability_omega_categorical()`](https://yelleknek.github.io/DMAR/reference/reliability_omega_categorical.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Fleiss (1971) Table 1 example: 30 subjects rated by 6 raters into
# 5 diagnostic categories (Depression, Personality Disorder,
# Schizophrenia, Neurosis, Other). Each row of `ratings` gives, for one
# subject, the count of raters who chose each category (rows sum to 6).
# kappa = 0.430, matching Fleiss (1971).
fleiss_1971 <- matrix(c(
  0, 0, 0, 6, 0,
  0, 3, 0, 0, 3,
  0, 1, 4, 0, 1,
  0, 0, 0, 0, 6,
  0, 3, 0, 3, 0,
  2, 0, 4, 0, 0,
  0, 0, 4, 0, 2,
  2, 0, 3, 1, 0,
  2, 0, 0, 4, 0,
  0, 0, 0, 0, 6,
  1, 0, 0, 5, 0,
  1, 1, 0, 4, 0,
  0, 3, 3, 0, 0,
  1, 0, 0, 5, 0,
  0, 2, 0, 3, 1,
  0, 0, 5, 0, 1,
  3, 0, 0, 1, 2,
  5, 1, 0, 0, 0,
  0, 2, 0, 4, 0,
  1, 0, 2, 0, 3,
  0, 0, 0, 0, 6,
  0, 1, 0, 5, 0,
  0, 2, 0, 1, 3,
  2, 0, 0, 4, 0,
  1, 0, 0, 4, 1,
  0, 5, 0, 1, 0,
  4, 0, 0, 0, 2,
  0, 2, 0, 4, 0,
  1, 0, 5, 0, 0,
  0, 0, 0, 0, 6
), nrow = 30, byrow = TRUE)
fleiss_kappa(fleiss_1971)
#>  kappa se     lower_limit upper_limit z_value p_value  n_subjects n_raters
#>  0.43  0.0542 0.324       0.536       17.7    < 0.0001 30         6       
#>  n_categories
#>  5           
#> 
#> Confidence level: 95%

# A bootstrap interval, which resamples the subjects (rows) with
# replacement and recomputes kappa on each resample. Not run here,
# because 2000 refits of kappa is more than a help page should do;
# the call is:
# fleiss_kappa(fleiss_1971, ci_method = "percentile", B = 2000,
#              seed = 113)
```
