# Krippendorff's \\\alpha\\ Inter-Rater Agreement

Computes Krippendorff's (1980, 2004, 2011) \\\alpha\\, the most general
chance-corrected inter-rater agreement coefficient. Unlike
[`cohen_kappa`](https://yelleknek.github.io/DMAR/reference/cohen_kappa.md)
(two raters, nominal data) or
[`fleiss_kappa`](https://yelleknek.github.io/DMAR/reference/fleiss_kappa.md)
(multiple raters, nominal data), Krippendorff's \\\alpha\\ supports any
number of raters, missing values, and any of four levels of measurement
(nominal, ordinal, interval, ratio) via a user-specified distance
metric.

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

  Logical. If `TRUE`, returns a bootstrap percentile CI. Set to `FALSE`
  to return only the point estimate (much faster).

- B:

  Number of bootstrap resamples when `boot = TRUE`. Default `1000L`.

- seed:

  Optional integer seed for reproducibility of the bootstrap. Default
  `NULL`, which leaves the user's current RNG state intact; supply an
  integer for reproducibility.

## Value

A `data.frame` with rows for the point estimate \\\hat\alpha\\, the
observed disagreement \\D_o\\, the expected disagreement \\D_e\\, the
number of pairable values, and, when a bootstrap was run, the lower and
upper bootstrap CI limits and `B_used`, the number of resamples that
returned a finite value and so entered the interval.

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
bias-corrected and accelerated (BCa) variant. Resamples on which the
coefficient cannot be computed (for example, a resample without enough
pairable values) are dropped, and the No closed-form sampling variance
is in general use for Krippendorff's alpha across its measurement levels
and missing-data patterns; the bootstrap is the interval Krippendorff
recommends (Krippendorff, 2011; Hayes & Krippendorff, 2007). The
interval is computed from the resamples that return a finite value; how
many did is reported as the `B_used` row of the result. `B = 1000L`
typically gives a stable CI to two decimal places. Bootstrap results
vary from run to run; supply `seed` for reproducibility.

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
[`lin_ccc()`](https://yelleknek.github.io/DMAR/reference/lin_ccc.md),
[`loa()`](https://yelleknek.github.io/DMAR/reference/loa.md),
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
#>  D_observed         0.2  
#>  D_expected         0.779
#>  n_pairable         40   
#> 
#> Confidence level: 95%

# 2. Interval ratings:
set.seed(113)
r1 <- rnorm(30, 0, 1)
r2 <- r1 + rnorm(30, 0, 0.3)
krippendorff_alpha(cbind(r1, r2), level = "interval",
                   boot = TRUE, B = 500L)
#>  term               value 
#>  krippendorff_alpha 0.956 
#>  D_observed         0.0955
#>  D_expected         2.16  
#>  n_pairable         60    
#>  lower_limit        0.916 
#>  upper_limit        0.976 
#>  B_used             500   
#> 
#> Confidence level: 95%
```
