# Confidence Interval for Generalized Eta Squared (Approximate)

Returns the point estimate of generalized eta squared (\\\eta^2_G\\;
Olejnik & Algina, 2003) along with an *optional* confidence interval
computed by one of two approximate methods. The default is to return
only the point estimate (`method = "none"`), because both available CI
methods are approximations whose coverage properties have not been
broadly validated for this estimand and warrant independent evaluation
before being used in substantive inference.

## Usage

``` r
ci_eta_squared_generalized(
  object = NULL,
  observed = NULL,
  SS_effect = NULL,
  SS_observed = NULL,
  SS_error = NULL,
  F_effect = NULL,
  df_effect = NULL,
  F_observed = NULL,
  df_observed = NULL,
  df_error = NULL,
  N = NULL,
  method = c("none", "parametric", "bootstrap"),
  B = 10000L,
  conf_level = 0.95,
  alpha_lower = NULL,
  alpha_upper = NULL,
  seed = NULL
)
```

## Arguments

- object:

  Optional. A fitted model object of class
  [`aov`](https://rdrr.io/r/stats/aov.html),
  [`lm`](https://rdrr.io/r/stats/lm.html), or `aovlist` (multi-stratum
  aov fit for within-subjects / mixed designs).

- observed:

  Character vector of factor names treated as measured.

- SS_effect, SS_observed, SS_error:

  Sums of squares (option 2 in
  [`eta_squared_generalized`](https://yelleknek.github.io/DMAR/reference/eta_squared_generalized.md)).

- F_effect, df_effect, F_observed, df_observed, df_error:

  *F*-values and degrees of freedom (option 3 in
  [`eta_squared_generalized`](https://yelleknek.github.io/DMAR/reference/eta_squared_generalized.md)).

- N:

  Total sample size. Required when `method = "parametric"` and no fitted
  model is supplied; ignored when `method = "none"` or when a fitted
  model is supplied (derived automatically, via
  [`nobs`](https://rdrr.io/r/stats/nobs.html)`(object)` for
  single-stratum fits, or by summing degrees of freedom across error
  strata for `aovlist`).

- method:

  One of `"none"` (default), `"parametric"`, or `"bootstrap"`. See
  Details.

- B:

  Integer. Number of bootstrap replications when `method = "bootstrap"`.
  Minimum `1000` (enforced); default `10000`. We recommend `10000` or
  more for publication-quality intervals.

- conf_level:

  Desired confidence coverage; default `0.95`.

- alpha_lower, alpha_upper:

  Optional Type I error on the lower and upper side.

- seed:

  Optional integer seed for the bootstrap, for reproducibility. Used
  locally: the caller's random number generator state is restored on
  exit. Default `NULL` leaves the random number generator state alone.

## Value

A `data.frame` with one row per focal effect and the columns `effect`,
`eta_squared_generalized`, `lower_limit`, `upper_limit`, and `method`.
For `aovlist` fits a `stratum` column is also present. When
`method = "none"`, the limit columns contain `NA`.

## Details

**Why CI = "none" is the default.** Confidence interval construction for
\\\eta^2_G\\ is not as settled as for partial \\\eta^2\\ or
\\\omega^2\\, because the denominator mixes sums of squares from
heterogeneous sources (the focal effect, one or more measured factors,
and the error term). No noncentral *F* transformation maps the
population noncentrality parameter directly to \\\eta^2_G\\. Both
methods below are approximations and are exposed for exploration rather
than as defaults.

**`method = "parametric"`.** The function first obtains a confidence
interval for the population noncentrality parameter \\\lambda\\ of the
focal effect's *F*-test via
[`conf_limits_ncf`](https://yelleknek.github.io/DMAR/reference/conf_limits_ncf.md).
The NCP bounds are mapped through the partial-\\\eta^2\\ transformation
\\\eta^2\_{p,\text{bound}} =
\lambda\_{\text{bound}}/(\lambda\_{\text{bound}} + N)\\ (matching the
convention used by
[`ci_pvaf`](https://yelleknek.github.io/DMAR/reference/ci_pvaf.md) and
[`ci_omega_squared`](https://yelleknek.github.io/DMAR/reference/ci_omega_squared.md)),
and then re-expressed as \\\eta^2_G\\ bounds via
\$\$\eta^2\_{G,\text{bound}} =
\frac{r\_{\text{bound}}}{r\_{\text{bound}} + r\_{\text{obs}} + 1},\$\$
where \\r\_{\text{bound}} =
\eta^2\_{p,\text{bound}}/(1-\eta^2\_{p,\text{bound}})\\ and
\\r\_{\text{obs}} = \sum
\mathit{SS}\_{\text{measured}}/\mathit{SS}\_{\text{error}}\\. This
treats the observed-factor sums of squares as fixed at their sample
values, so the interval inherits whatever sampling variability those
contribute. It has not been validated for coverage and should be treated
as preliminary.

**`method = "bootstrap"`.** A residual bootstrap from the fitted model:
the resampling unit is a residual, drawn nonparametrically from the
model's own residuals rather than from a fitted distribution. For each
of the `B` replications, the response is regenerated as \\\hat{y}\_i +
\varepsilon^\*\_i\\ where \\\varepsilon^\*\\ is sampled with replacement
from the model's residuals; the model is refit; \\\eta^2_G\\ is
recomputed; and the percentile interval, the empirical quantiles of the
`B` bootstrap values (Efron & Tibshirani, 1993), is reported. The
percentile interval is the only bootstrap interval offered; there is no
bias-corrected and accelerated (BCa) variant. Replicates whose refit
fails are dropped, and the interval is computed from the replications
that return a value. Bootstrap results vary from run to run; supply
`seed` for reproducibility. Requires a single-stratum `aov`/`lm` fit.
**Not yet supported for `aovlist` (multi-stratum / within-subjects)
fits**, since the bootstrap needs to respect the subject-level
correlation structure, which naive residual resampling does not. Use
`method = "parametric"` for `aovlist` fits. Coverage has not been
broadly validated for this estimand.

**Within-subjects designs (`aovlist`).** For multi-stratum fits the
parametric CI uses each focal effect's stratum-specific *F* test and
degrees of freedom. The denominator ratio is adjusted to reflect the
full set of error strata: the implied \\\mathit{SS}\_{\text{error}}\\
for the focal effect is its own stratum's residual SS, while the
\\r\_{\text{obs}}\\ term includes both measured-factor SS and the
residual SS of all *other* strata. This generalizes the
partial-\\\eta^2\\ CI machinery to the Bakeman (2005) denominator.

## References

Algina, J., Keselman, H. J., & Penfield, R. D. (2005). An alternative to
Cohen's standardized mean difference effect size: A robust parameter and
confidence interval in the two independent groups case. *Psychological
Methods, 10*(3), 317–328.
[doi:10.1037/1082-989X.10.3.317](https://doi.org/10.1037/1082-989X.10.3.317)

Bakeman, R. (2005). Recommended effect size statistics for repeated
measures designs. *Behavior Research Methods, 37*(3), 379–384.
[doi:10.3758/BF03192707](https://doi.org/10.3758/BF03192707)

Efron, B., & Tibshirani, R. J. (1993). *An introduction to the
bootstrap*. New York, NY: Chapman & Hall/CRC.

Kelley, K. (2007). Confidence intervals for standardized effect sizes:
Theory, application, and implementation. *Journal of Statistical
Software, 20*(8), 1–24.
[doi:10.18637/jss.v020.i08](https://doi.org/10.18637/jss.v020.i08)

Kelley, K., & Preacher, K. J. (2012). On effect size. *Psychological
Methods, 17*, 137–152.
[doi:10.1037/a0028086](https://doi.org/10.1037/a0028086)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 3 on \\\eta^2\\, Chapter 7 on factorial
designs, and Chapter 11 on generalized \\\eta^2\\ for within-subjects
designs.)

Olejnik, S., & Algina, J. (2003). Generalized eta and omega squared
statistics: Measures of effect size for some common research designs.
*Psychological Methods, 8*(4), 434–447.
[doi:10.1037/1082-989X.8.4.434](https://doi.org/10.1037/1082-989X.8.4.434)

Smithson, M. (2001). Correct confidence intervals for various regression
effect sizes and parameters: The importance of noncentral distributions
in computing intervals. *Educational and Psychological Measurement, 61*,
605–632.
[doi:10.1177/00131640121971392](https://doi.org/10.1177/00131640121971392)

Steiger, J. H. (2004). Beyond the *F* test: Effect size confidence
intervals and tests of close fit in the analysis of variance and
contrast analysis. *Psychological Methods, 9*(2), 164–182.
[doi:10.1037/1082-989X.9.2.164](https://doi.org/10.1037/1082-989X.9.2.164)

## See also

[`eta_squared_generalized`](https://yelleknek.github.io/DMAR/reference/eta_squared_generalized.md),
[`ci_eta_squared`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared.md)

Other confidence intervals for effect sizes:
[`ci_R2()`](https://yelleknek.github.io/DMAR/reference/ci_R2.md),
[`ci_c()`](https://yelleknek.github.io/DMAR/reference/ci_c.md),
[`ci_c_ancova()`](https://yelleknek.github.io/DMAR/reference/ci_c_ancova.md),
[`ci_c_ancova_bp()`](https://yelleknek.github.io/DMAR/reference/ci_c_ancova_bp.md),
[`ci_correlation`](https://yelleknek.github.io/DMAR/reference/ci_correlation.md),
[`ci_cv()`](https://yelleknek.github.io/DMAR/reference/ci_cv.md),
[`ci_eta_squared()`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared.md),
[`ci_eta_squared_partial()`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared_partial.md),
[`ci_mahalanobis()`](https://yelleknek.github.io/DMAR/reference/ci_mahalanobis.md),
[`ci_omega_squared()`](https://yelleknek.github.io/DMAR/reference/ci_omega_squared.md),
[`ci_pvaf()`](https://yelleknek.github.io/DMAR/reference/ci_pvaf.md),
[`ci_rc()`](https://yelleknek.github.io/DMAR/reference/ci_rc.md),
[`ci_reg_coef()`](https://yelleknek.github.io/DMAR/reference/ci_reg_coef.md),
[`ci_rmsea()`](https://yelleknek.github.io/DMAR/reference/ci_rmsea.md),
[`ci_sc()`](https://yelleknek.github.io/DMAR/reference/ci_sc.md),
[`ci_sc_ancova()`](https://yelleknek.github.io/DMAR/reference/ci_sc_ancova.md),
[`ci_sm()`](https://yelleknek.github.io/DMAR/reference/ci_sm.md),
[`ci_smd()`](https://yelleknek.github.io/DMAR/reference/ci_smd.md),
[`ci_smd_c()`](https://yelleknek.github.io/DMAR/reference/ci_smd_c.md),
[`ci_snr()`](https://yelleknek.github.io/DMAR/reference/ci_snr.md),
[`ci_src()`](https://yelleknek.github.io/DMAR/reference/ci_src.md),
[`ci_srsnr()`](https://yelleknek.github.io/DMAR/reference/ci_srsnr.md),
[`contrast_adjusted()`](https://yelleknek.github.io/DMAR/reference/contrast_adjusted.md),
[`plot_smd()`](https://yelleknek.github.io/DMAR/reference/plot_smd.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# The pygmalion expectancy experiment: treatment is manipulated, while
# grade is a measured classification, so grade belongs in the denominator.
pyg <- pygmalion
pyg$grade <- factor(pyg$grade)
fit <- aov(iq_8 ~ treatment * grade, data = pyg)

# The default returns the point estimate and leaves the limits NA, because
# both interval methods are approximations.
ci_eta_squared_generalized(fit, observed = "grade")
#> No CI computed (method = 'none'). Set method = 'parametric' for an approximate noncentral F transformation, or method = 'bootstrap' for a residual bootstrap (requires a fitted model). Both methods are approximate; see ?ci_eta_squared_generalized.
#>  effect          eta_squared_generalized lower_limit upper_limit method
#>  treatment       0.0202                  <NA>        <NA>        none  
#>  grade           0.044                   <NA>        <NA>        none  
#>  treatment:grade 0.0187                  <NA>        <NA>        none  

# The parametric approximation maps a noncentral F interval for the focal
# effect through the observed sums of squares. Warnings mark the method as
# preliminary and report that the interaction's lower limit is clamped at
# 0; treat the limits accordingly.
ci_eta_squared_generalized(fit, observed = "grade", method = "parametric")
#> Warning: Parametric CI uses an approximate transformation that maps the partial-eta squared CI through the observed sums of squares for measured factors. Coverage properties have not been broadly validated; results should be treated as preliminary and independently evaluated.
#> Warning: The observed F_value is below the alpha_lower critical value of the central F-distribution, so the lower confidence limit on generalized eta squared is 0.
#>  effect          eta_squared_generalized lower_limit upper_limit method    
#>  treatment       0.0202                  0.000974    0.0583      parametric
#>  grade           0.044                   0.00116     0.0807      parametric
#>  treatment:grade 0.0187                  0           0.0413      parametric
#> 
#> Confidence level: 95%

# The third option is a residual bootstrap, which refits the model once per
# replication. It is not run here: B is required to be at least 1000, and
# even that is slower than an example should be, while a reported interval
# deserves B = 10000 or more. The call is
#   ci_eta_squared_generalized(fit, observed = "grade",
#                              method = "bootstrap", B = 10000, seed = 113)

# Within-subjects ANOVA. The parametric CI uses each effect's own stratum.
set.seed(113)
n <- 20
rm_data <- data.frame(
  subject = factor(rep(seq_len(n), each = 3)),
  time    = factor(rep(c("Pre", "Mid", "Post"), n),
                   levels = c("Pre", "Mid", "Post")),
  y       = rnorm(n, sd = 1.5)[rep(seq_len(n), each = 3)] +
            0.7 * rep(1:3, n) + rnorm(n * 3, sd = 1.2)
)
fit_rm <- aov(y ~ time + Error(subject/time), data = rm_data)
ci_eta_squared_generalized(fit_rm, method = "parametric")
#> Warning: Parametric CI uses an approximate transformation that maps the partial-eta squared CI through the observed sums of squares for measured factors. Coverage properties have not been broadly validated; results should be treated as preliminary and independently evaluated.
#>  effect eta_squared_generalized stratum      lower_limit upper_limit method    
#>  time   0.0965                  subject:time 0.00548     0.145       parametric
#> 
#> Confidence level: 95%
```
