# Krippendorff's \\\alpha\\ Inter-Rater Agreement

Computes Krippendorff's (1980, 2004, 2011) \\\alpha\\, a
chance-corrected inter-rater agreement coefficient. Unlike
[`cohen_kappa`](https://yelleknek.github.io/DMAR/reference/cohen_kappa.md)
(two raters; nominal categories, or ordinal ones through weights) or
[`fleiss_kappa`](https://yelleknek.github.io/DMAR/reference/fleiss_kappa.md)
(multiple raters, nominal categories, no missing ratings),
Krippendorff's \\\alpha\\ supports any number of raters, missing values,
and any of four levels of measurement (nominal, ordinal, interval,
ratio) via a user-specified distance metric.

## Usage

``` r
krippendorff_alpha(
  ratings,
  level = c("nominal", "ordinal", "interval", "ratio"),
  conf_level = 0.95,
  boot = FALSE,
  B = 1000L,
  seed = NULL
)
```

## Arguments

- ratings:

  A units \\\times\\ raters matrix (or `data.frame`). Rows = units of
  analysis; columns = raters. `NA` entries are allowed.

- level:

  One of `"nominal"` (default), `"ordinal"`, `"interval"`, or `"ratio"`;
  controls the distance metric used to compute disagreement.

- conf_level:

  Confidence level for the bootstrap CI. Default `0.95`.

- boot:

  Logical. If `TRUE`, fills the interval rows with a bootstrap
  percentile CI. The default `FALSE` runs no bootstrap and leaves those
  rows `NA` (much faster).

- B:

  Number of bootstrap resamples when `boot = TRUE`, a single whole
  number. Default `1000L`, because every resample recomputes the
  coincidences over all pairable values. The fewest accepted is 400 at
  the default `conf_level` (200 at 0.90, and 2000 at 0.99, which the
  default does not meet); a smaller `B` is refused before anything is
  computed (class `dmar_count_invalid`). See *Bootstrap replications*
  under *Details* and
  [`dmar_bootstrap`](https://yelleknek.github.io/DMAR/reference/dmar_bootstrap.md).

- seed:

  Optional integer seed for reproducibility of the bootstrap. Default
  `NULL`, under which the resamples are drawn from, and advance, the
  current random number stream; supply an integer for reproducibility,
  which is set for the call only, the caller's random number state being
  restored when it returns.

## Value

A `data.frame` (class `dmar_tbl`) with the same seven rows whatever
`boot`: the point estimate \\\hat\alpha\\ (`krippendorff_alpha`), the
observed disagreement \\D_o\\ (`D_observed`), the expected disagreement
\\D_e\\ (`D_expected`), the number of pairable values (`n_pairable`),
the lower and upper bootstrap CI limits (`ci_lower`, `ci_upper`), and
`B_used`, the number of resamples that returned a finite value and so
entered the interval. The last three are `NA` when `boot = FALSE`, since
no bootstrap was run. With `boot = TRUE` the table also carries the
`conf_level` attribute (when an interval was formed) and the attributes
`B`, `B_used`, `n_dropped`, `B_required`, and `bootstrap_status` (see
*Bootstrap replications* under Details).

## Details

**Coefficient.** Krippendorff's \\\alpha\\ is \$\$\alpha \\=\\ 1 -
\frac{D_o}{D_e},\$\$ where \\D_o\\ is the observed disagreement (average
squared distance over all within-unit pairs of ratings, scaled by the
number of pairable values), and \\D_e\\ is the expected disagreement
(average squared distance over all between-unit pairs). The metric used
in the squared distance depends on `level`:

- *nominal*: \\d(a, b) = \mathrm{I}(a \ne b)\\

- *ordinal*: distance based on cumulative rank counts

- *interval*: \\d(a, b) = (a - b)^2\\

- *ratio*: \\d(a, b) = ((a - b) / (a + b))^2\\

**CI.** The CI is by case-resampling bootstrap over units (rows): the
rows of `ratings` are resampled with replacement `B` times and
\\\alpha\\ is recomputed on each resample, so units are the sampling
unit and the rater panel is treated as fixed. Only the percentile
interval is offered: the limits are the empirical quantiles of the
bootstrap estimates (Efron & Tibshirani, 1993); there is no
bias-corrected and accelerated (BCa) variant. When the coefficient is
undefined for the data themselves (every value alike), no bootstrap is
run and the interval rows stay `NA`. No closed-form sampling variance is
in general use for Krippendorff's alpha across its measurement levels
and missing data patterns, so the bootstrap is the interval Krippendorff
recommends (Krippendorff, 2011; Hayes & Krippendorff, 2007). `B = 1000L`
typically gives a stable CI to two decimal places. The bootstrap is
opt-in (`boot = FALSE` by default, which returns the point estimate
alone and is much faster); ask for it whenever the coefficient is being
reported rather than explored, since a point estimate on its own says
nothing about how precisely \\\alpha\\ is determined.

**Bootstrap replications.** The two limits are quantiles of the `B`
replicated values of \\\alpha\\, so the interval is run only when `B`
puts at least 10 resamples beyond each limit (400 at the 95 percent
level). A replication in which \\\alpha\\ is undefined (a resample
without enough pairable values, or one whose values are all alike, so
that the expected disagreement is zero) is dropped and counted, never
replaced, and the interval is computed from the rest. The result records
`B` (requested), `B_used` (the number the interval rests on),
`n_dropped`, and `B_required` (the fewest replications a request may ask
for; that count is held to `B`, not to `B_used`) as attributes, and
`B_used` also as a row, with `bootstrap_status` (`"complete"`,
`"dropped"`, `"exploratory"`, `"insufficient_count"`, or
`"insufficient_share"`). One warning of class
`dmar_bootstrap_replications_dropped` states the counts when any
replication was dropped, with a caution when the dropped replications
are at least the 2.5 percent of `B` that lies beyond one limit at the 95
percent level, because they could then decide by themselves where that
limit falls. When more than half were dropped, the two limit rows are
`NA`, the estimate is returned, and the warning has class
`dmar_bootstrap_insufficient`.
`options(dmar.bootstrap_exploratory = TRUE)` lets a smaller `B` (2 or
more) run for a first look or a demonstration; every such call gives a
warning of class `dmar_bootstrap_exploratory` and the result is not one
to report. See
[`dmar_bootstrap`](https://yelleknek.github.io/DMAR/reference/dmar_bootstrap.md)
for the rule and its reasons. Bootstrap results vary from run to run;
supply `seed` for reproducibility.

**Interpretation.** \\\alpha\\ ranges from \\-D_e / D_o\\ (perfect
disagreement) through \\0\\ (chance level) to \\1\\ (perfect agreement).
Report the coefficient with its confidence interval and judge it against
the reliability the application requires; Krippendorff (2004) discusses
how that judgment depends on the cost of acting on unreliable data.

## References

Efron, B., & Tibshirani, R. J. (1993). *An introduction to the
bootstrap*. New York, NY: Chapman & Hall/CRC.

Hayes, A. F., & Krippendorff, K. (2007). Answering the call for a
standard reliability measure for coding data. *Communication Methods and
Measures, 1*(1), 77–89.
[doi:10.1080/19312450709336664](https://doi.org/10.1080/19312450709336664)

Krippendorff, K. (1980). *Content analysis: An introduction to its
methodology*. Sage.

Krippendorff, K. (2004). *Content analysis: An introduction to its
methodology* (2nd ed.). Sage.

Krippendorff, K. (2011). Computing Krippendorff's alpha-reliability.
*Departmental Papers (ASC)*, Annenberg School for Communication,
University of Pennsylvania.

## See also

[`cohen_kappa`](https://yelleknek.github.io/DMAR/reference/cohen_kappa.md),
[`fleiss_kappa`](https://yelleknek.github.io/DMAR/reference/fleiss_kappa.md),
[`icc`](https://yelleknek.github.io/DMAR/reference/icc.md)

Other agreement and measurement:
[`R2_mixed_effects()`](https://yelleknek.github.io/DMAR/reference/R2_mixed_effects.md),
[`content_validity_index()`](https://yelleknek.github.io/DMAR/reference/content_validity_index.md),
[`gwet_ac()`](https://yelleknek.github.io/DMAR/reference/gwet_ac.md),
[`icc_lmer()`](https://yelleknek.github.io/DMAR/reference/icc_lmer.md),
[`limits_of_agreement()`](https://yelleknek.github.io/DMAR/reference/limits_of_agreement.md),
[`lin_ccc()`](https://yelleknek.github.io/DMAR/reference/lin_ccc.md),
[`variance_components_mls()`](https://yelleknek.github.io/DMAR/reference/variance_components_mls.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# 1. Nominal ratings, 4 raters, 12 units (from Krippendorff 2011 Tab. 1):
ratings <- matrix(c(
  1, 2, 3, 3, 2, 1, 4, 1, 2, NA, NA, NA,
  1, 2, 3, 3, 2, 2, 4, 1, 2, 5,  NA, 3,
  NA, 3, 3, 3, 2, 3, 4, 2, 2, 5,  1,  NA,
  1, 2, 3, 3, 2, 4, 4, 1, 2, 5,  1,  NA
), nrow = 12, ncol = 4)
krippendorff_alpha(ratings, level = "nominal")
#>  term               value
#>  krippendorff_alpha 0.743
#>  D_observed         0.200
#>  D_expected         0.779
#>  n_pairable         40   
#>  ci_lower           <NA> 
#>  ci_upper           <NA> 
#>  B_used             <NA> 

# 2. Interval ratings:
set.seed(113)
r1 <- rnorm(30, 0, 1)
r2 <- r1 + rnorm(30, 0, 0.3)
krippendorff_alpha(cbind(r1, r2), level = "interval")
#>  term               value
#>  krippendorff_alpha 0.956
#>  D_observed         0.096
#>  D_expected         2.162
#>  n_pairable         60   
#>  ci_lower           <NA> 
#>  ci_upper           <NA> 
#>  B_used             <NA> 

# The percentile bootstrap interval for the same ratings, which
# recomputes alpha on each of B resamples of the units and fills the
# ci_lower, ci_upper, and B_used rows. B = 400 is the fewest
# resamples a 95 percent interval is reported from and keeps the
# example quick; a reported interval deserves the default B = 1000.
krippendorff_alpha(cbind(r1, r2), level = "interval",
                   boot = TRUE, B = 400L, seed = 113)
#>  term               value
#>  krippendorff_alpha 0.956
#>  D_observed         0.096
#>  D_expected         2.162
#>  n_pairable         60   
#>  ci_lower           0.919
#>  ci_upper           0.975
#>  B_used             400  
#> 
#> Confidence level: 95%
#> 
#> Bootstrap: 400 replications, 400 used (0 dropped; at least 400 must be requested)
```
