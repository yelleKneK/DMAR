# Mediation Analysis With Bootstrap Confidence Intervals

Estimates the simple mediation model, a predictor \\X\\ affecting an
outcome \\Y\\ directly and through a mediator \\M\\, and reports the
indirect, direct, and total effects with confidence intervals in one
tidy table. The indirect effect \\a b\\ is the product of the \\X \to
M\\ path and the \\M \to Y\\ path (holding \\X\\), and its sampling
distribution is skewed, which is why the default interval is the
percentile bootstrap rather than a normal approximation; the
bias-corrected and accelerated (BCa) bootstrap, the Monte Carlo
(parametric simulation) interval, and the Sobel normal-theory interval
are available for comparison. This is the analysis counterpart of the
planning functions
[`ss_aipe_indirect_effect`](https://yelleknek.github.io/DMAR/reference/ss_aipe_indirect_effect.md)
and
[`ss_power_indirect_effect`](https://yelleknek.github.io/DMAR/reference/ss_power_indirect_effect.md).

## Usage

``` r
mediate(
  data,
  x,
  m,
  y,
  covariates = NULL,
  ci_method = c("boot_percentile", "boot_bca", "monte_carlo", "sobel"),
  B = 2000,
  conf_level = 0.95,
  seed = NULL
)
```

## Arguments

- data:

  A `data.frame` containing the variables.

- x, m, y:

  Names (single character strings) of the predictor, the mediator, and
  the outcome columns in `data`.

- covariates:

  Optional character vector of covariate column names, entered in both
  the mediator and the outcome models.

- ci_method:

  Confidence interval method for the indirect effect:
  `"boot_percentile"` (default), `"boot_bca"`, `"monte_carlo"` (simulate
  \\a\\ and \\b\\ from their joint normal approximation; MacKinnon,
  Lockwood, & Williams, 2004), or `"sobel"` (the first-order
  normal-theory interval; reported for comparison, not recommended for
  inference). Direct and total effects always carry their ordinary
  *t*-based intervals.

- B:

  Number of bootstrap or Monte Carlo replications. Defaults to 2000;
  published analyses often use 5000 or more, and the BCa interval in
  particular rewards a large `B` (see *Details*).

- conf_level:

  Confidence level. Defaults to 0.95.

- seed:

  Optional integer seed for the resampling, used locally (the caller's
  random number generator state is restored on exit). Default `NULL`
  leaves the random number generator state alone.

## Value

A `data.frame` (class `dmar_tbl`) with rows `indirect_effect` (with its
`ci_method` interval), `direct_effect`, `total_effect`, the paths `a`
and `b`, their standard errors (`se_indirect` per Sobel, `se_a`,
`se_b`), the interval limits (`indirect_lower` / `indirect_upper`,
`direct_lower` / `direct_upper`, `total_lower` / `total_upper`),
`proportion_mediated` (\\ab/c\\; `NA` when the total effect is near
zero, where the ratio is unstable), `N`, and `B`. The interval method is
recorded in the `"ci_method"` attribute and the confidence level in
`"conf_level"`.

## Details

The model is the standard pair of regressions \$\$M = i_M + a X +
\mathbf{g}'\mathbf{C} + e_M, \qquad Y = i_Y + c' X + b M +
\mathbf{h}'\mathbf{C} + e_Y,\$\$ with \\\mathbf{C}\\ the optional
covariates. The indirect effect is \\a b\\, the direct effect \\c'\\,
and the total effect \\c = c' + a b\\ (an identity in linear models with
the same cases, which the implementation exploits as an internal
consistency check).

Bootstrap intervals resample cases (rows) with replacement `B` times,
refitting both regressions in each resample (Efron & Tibshirani, 1993).
`"boot_percentile"` takes the interval limits from the empirical
quantiles of the bootstrapped \\a b\\ estimates; it is not forced to be
symmetric about the estimate, which is the point for a skewed sampling
distribution. `"boot_bca"` (the bias-corrected and accelerated interval)
additionally adjusts the two quantile positions for median bias,
estimated from the bootstrap distribution, and for the rate at which the
variance of the estimator changes with the parameter, the acceleration,
estimated by the jackknife (which adds *N* extra pairs of fits); the
adjustments make it second-order accurate where the percentile interval
is first-order accurate (DiCiccio & Efron, 1996). Because the adjusted
quantile positions sit farther into the tails of the bootstrap
distribution, the BCa interval benefits more than the percentile
interval does from a `B` well above the default. Each resample refits
the two regressions by least squares, a closed-form fit with no
iterative estimation, so in ordinary data all `B` replications enter the
interval. A degenerate resample (one whose refit is rank deficient,
possible with a near-constant predictor) returns no indirect effect;
such replications are dropped with a warning stating how many, and the
interval is computed from the replications that returned a value. The
Sobel standard error is computed by
[`var_indirect_effect`](https://yelleknek.github.io/DMAR/reference/var_indirect_effect.md).
The proportion mediated is one of the effect size measures for mediation
models surveyed by Preacher and Kelley (2011); its instability when the
total effect is small is why it is reported as `NA` near a zero total
effect. Listwise deletion is applied to the analysis variables; for full
information maximum likelihood under missingness, fit the model in
lavaan (see [`mlmr`](https://yelleknek.github.io/DMAR/reference/mlmr.md)
for the package's FIML front end philosophy).

Mediation language implies causal structure: with observational data the
estimates are conditional associations, and the causal reading requires
the usual no-unmeasured-confounding assumptions for both the \\X \to M\\
and \\M \to Y\\ links (MacKinnon, 2008). The function computes; the
design earns the interpretation.

## References

DiCiccio, T. J., & Efron, B. (1996). Bootstrap confidence intervals.
*Statistical Science, 11*(3), 189–228.

Efron, B., & Tibshirani, R. J. (1993). *An introduction to the
bootstrap*. New York, NY: Chapman & Hall/CRC.

MacKinnon, D. P. (2008). *Introduction to statistical mediation
analysis*. Erlbaum.

MacKinnon, D. P., Lockwood, C. M., & Williams, J. (2004). Confidence
limits for the indirect effect: Distribution of the product and
resampling methods. *Multivariate Behavioral Research, 39*(1), 99–128.
[doi:10.1207/s15327906mbr3901_4](https://doi.org/10.1207/s15327906mbr3901_4)

Preacher, K. J., & Hayes, A. F. (2008). Asymptotic and resampling
strategies for assessing and comparing indirect effects in multiple
mediator models. *Behavior Research Methods, 40*(3), 879–891.
[doi:10.3758/BRM.40.3.879](https://doi.org/10.3758/BRM.40.3.879)

Preacher, K. J., & Kelley, K. (2011). Effect size measures for mediation
models: Quantitative strategies for communicating indirect effects.
*Psychological Methods, 16*(2), 93–115.
[doi:10.1037/a0022658](https://doi.org/10.1037/a0022658)

## See also

[`ss_aipe_indirect_effect`](https://yelleknek.github.io/DMAR/reference/ss_aipe_indirect_effect.md)
and
[`ss_power_indirect_effect`](https://yelleknek.github.io/DMAR/reference/ss_power_indirect_effect.md)
for planning the study this function analyzes;
[`var_indirect_effect`](https://yelleknek.github.io/DMAR/reference/var_indirect_effect.md)
for the Sobel variance.

Other mediation:
[`mediation_mbco()`](https://yelleknek.github.io/DMAR/reference/mediation_mbco.md),
[`ss_power_indirect_effect()`](https://yelleknek.github.io/DMAR/reference/ss_power_indirect_effect.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Simulated mediation: X raises M (a = .5), M raises Y (b = .4), and a
# little direct effect remains (c' = .2).
set.seed(113)
n <- 200
x <- rnorm(n)
m <- 0.5 * x + rnorm(n, 0, sqrt(1 - 0.25))
y <- 0.2 * x + 0.4 * m + rnorm(n, 0, 0.8)
d <- data.frame(x = x, m = m, y = y)

# The Sobel interval is closed form, so it is the call that runs here.
# It assumes the product ab is normally distributed, which is why it is
# reported for comparison rather than used for inference. The B row of
# the result is NA because no replications are drawn.
mediate(d, x = "x", m = "m", y = "y", ci_method = "sobel")
#>  term                value 
#>  indirect_effect     0.182 
#>  indirect_lower      0.105 
#>  indirect_upper      0.259 
#>  direct_effect       0.194 
#>  direct_lower        0.0634
#>  direct_upper        0.324 
#>  total_effect        0.375 
#>  total_lower         0.25  
#>  total_upper         0.501 
#>  a                   0.449 
#>  b                   0.404 
#>  se_a                0.0615
#>  se_b                0.0678
#>  se_indirect         0.0393
#>  proportion_mediated 0.484 
#>  N                   200   
#>  B                   <NA>  
#> 
#> Confidence level: 95%

# The percentile bootstrap is the default and is what a reported
# analysis would use. It is not run here because it refits both
# regressions in each of the B resamples; the call is:
# mediate(d, x = "x", m = "m", y = "y", seed = 113)

# The Monte Carlo interval needs no refitting, so it stays quick even
# for a large B, but it still draws B replications. Also not run here:
# mediate(d, x = "x", m = "m", y = "y", ci_method = "monte_carlo",
#         B = 10000, seed = 113)
```
