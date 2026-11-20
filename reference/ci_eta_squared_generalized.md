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
  [`eta_squared_generalized`](https://yelleknek.github.io/DMAR/reference/eta_squared_generalized.md)):
  `SS_effect` for the focal effect, `SS_error` for its error term, and
  `SS_observed` (a scalar or a vector) for every other source in the
  denominator of \\\eta^2_G\\. The estimate depends on `SS_observed` and
  `SS_error` only through their sum. The parametric interval is built on
  the *F* they imply,
  \\(\mathit{SS}\_{\text{effect}}/df\_{\text{effect}})/(\mathit{SS}\_{\text{error}}/df\_{\text{error}})\\,
  which is the *F* of the effect's test only when `SS_error` is the sum
  of squares that test divides by: the residual, adjusted for any
  covariate, in a design without an `Error()` term, and in a design with
  one the residual of the stratum in which the effect is tested
  (subjects within groups for an effect tested between subjects, and in
  a one-way within-subjects design the interaction of subjects and the
  within-subjects factor). `SS_observed` then holds the other sources
  that Olejnik and Algina (2003, Eq. 5, p. 437) place in the
  denominator: the measured factors, the interactions that contain one,
  the covariates, and the residuals of the other strata, such as that of
  the subjects stratum of a within-subjects design. The interval also
  needs `df_effect`, `df_error`, and `N`.

- F_effect, df_effect, F_observed, df_observed, df_error:

  *F*-values and degrees of freedom (option 3 in
  [`eta_squared_generalized`](https://yelleknek.github.io/DMAR/reference/eta_squared_generalized.md)).
  Each *F* divides a mean square by the focal effect's error mean
  square: `F_effect` is the effect's test, and `F_observed` covers the
  sources that `SS_observed` holds, so for a source in another stratum
  it is not the *F* that
  [`summary()`](https://rdrr.io/r/base/summary.html) of the fit reports.
  `df_effect` and `df_error` go with the sums of squares as well. Given
  with the sums of squares, `F_effect`, `F_observed`, and `df_observed`
  must describe the same test: the *F* the sums of squares imply for the
  effect, and the ratio of `SS_observed` to `SS_error`, must agree with
  the *F* interface's to a relative tolerance of 1e-6, the tolerance
  [`eta_squared_generalized`](https://yelleknek.github.io/DMAR/reference/eta_squared_generalized.md)
  applies to the two estimates, and the interval is then formed from the
  *F*-values. The two can agree on the estimate and still imply
  different *F* statistics for the effect, and so different intervals;
  such a call is refused.

- N:

  Total sample size, with which the confidence limits for the
  noncentrality parameter are converted (see Details). Required when
  `method = "parametric"` and no fitted model is supplied; ignored when
  `method = "none"` or when a fitted model is supplied, in which case a
  single-stratum fit is converted with
  [`nobs`](https://rdrr.io/r/stats/nobs.html)`(object)` and each effect
  of an `aovlist` fit with the size of its own stratum, as in
  [`ci_eta_squared`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared.md).
  For an effect tested in an `Error()` stratum, give the size of that
  stratum instead: the number of subjects for an effect tested between
  subjects, and for an effect tested within subjects the stratum's total
  degrees of freedom, which is `df_effect + df_error` when the effect is
  the only one tested in its stratum. `N` must be finite and at least
  `df_effect + df_error`.

- method:

  One of `"none"` (default), `"parametric"`, or `"bootstrap"`. See
  Details.

- B:

  Number of bootstrap replications when `method = "bootstrap"`, a single
  whole number. Default `10000`. The fewest accepted is 1000 through the
  98 percent level (this function's own floor) and 2000 at the 99
  percent level, where 10 replications beyond each limit take that many;
  a smaller `B` is refused before anything is computed (class
  `dmar_count_invalid`). A reported interval deserves `10000` or more.
  See *Bootstrap replications* under *Details* and
  [`dmar_bootstrap`](https://yelleknek.github.io/DMAR/reference/dmar_bootstrap.md).

- conf_level:

  Desired confidence coverage; default `0.95`.

- alpha_lower, alpha_upper:

  Optional Type I error on the lower and upper side. The table records
  the coverage, `1 - alpha_lower - alpha_upper`, as its level and, when
  the two tail areas differ, the areas themselves, which the line
  beneath the printed table,
  [`results_sentence`](https://yelleknek.github.io/DMAR/reference/dmar_output_helpers.md),
  and the `alpha_lower` and `alpha_upper` columns of `tidy()` state;
  with a tail area of 0 the interval is a one-sided confidence bound.

- seed:

  Optional integer seed for the bootstrap, for reproducibility. Used
  locally: the caller's random number generator state is restored on
  exit. Default `NULL` draws the resamples from the current random
  number stream, which advances that stream.

## Value

A `data.frame` with one row per focal effect and the columns `effect`,
`eta_squared_generalized`, `ci_lower`, `ci_upper`, and `method`. For
`aovlist` fits a `stratum` column is also present. When
`method = "none"`, the limit columns contain `NA`. Otherwise the table
records the confidence level of its intervals in its `conf_level`
attribute: `conf_level`, or `1 - alpha_lower - alpha_upper` when the
tail areas are supplied. With `method = "bootstrap"` the attributes `B`,
`B_used`, `n_dropped`, `B_required`, and `bootstrap_status` record the
replications (see *Bootstrap replications* under Details).

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
[`ci_nc_F`](https://yelleknek.github.io/DMAR/reference/ci_nc_F.md) and
maps its limits through the partial \\\eta^2\\ transformation
\\\eta^2\_{p,\text{bound}} =
\lambda\_{\text{bound}}/(\lambda\_{\text{bound}} + N)\\, as
[`ci_eta_squared`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared.md)
converts them: *N* is the number of observations for a fit without an
`Error()` term (the convention of
[`ci_pvaf`](https://yelleknek.github.io/DMAR/reference/ci_pvaf.md) and
[`ci_omega_squared`](https://yelleknek.github.io/DMAR/reference/ci_omega_squared.md)
as well), the size of the effect's stratum for an effect of an `aovlist`
fit (see below), and `N` as supplied for the raw interfaces. The limits
are then re-expressed as \\\eta^2_G\\ bounds via
\$\$\eta^2\_{G,\text{bound}} =
\frac{r\_{\text{bound}}}{r\_{\text{bound}} + r\_{\text{obs}} + 1},\$\$
where \\r\_{\text{bound}} =
\eta^2\_{p,\text{bound}}/(1-\eta^2\_{p,\text{bound}}) =
\lambda\_{\text{bound}}/N\\ and \\r\_{\text{obs}}\\ is the sum of the
other sums of squares in the denominator of \\\eta^2_G\\ (those of the
measured sources and, for an `aovlist` fit, the residuals of the other
strata) divided by the focal effect's error sum of squares. The point
estimate is the same map applied to
\\\mathit{SS}\_{\text{effect}}/\mathit{SS}\_{\text{error}} = F \cdot
df\_{\text{effect}}/df\_{\text{error}}\\, which for a fitted model
converges to \\\lambda/N\\ as the sample grows. When \\r\_{\text{obs}}\\
is 0, as in a fit without an `Error()` term in which nothing measured
enters the denominator, \\\eta^2_G\\ is the partial \\\eta^2\\, and the
interval is the exact one `ci_eta_squared(estimand = "partial")` gives.
Otherwise the interval is not exact, and no conversion of the limits for
\\\lambda\\ can make it so: \\\eta^2_G\\ then depends on the value that
\\r\_{\text{obs}}\\ converges to as well as on \\\lambda\\, while the
distribution of the focal *F* depends on \\\lambda\\ alone, so the
transformation principle behind the exact intervals (Steiger, 2004,
Proposition 1, p. 167), which maps the limits through a function of that
one parameter, does not apply. The interval takes \\r\_{\text{obs}}\\ as
known and equal to its sample value, which leaves the sampling
variability of \\r\_{\text{obs}}\\ out of it. Its coverage has not been
validated, and it should be treated as preliminary.

**`method = "bootstrap"`.** A residual bootstrap from the fitted model:
the resampling unit is a residual, drawn nonparametrically from the
model's own residuals rather than from a fitted distribution. For each
of the `B` replications, the response is regenerated as \\\hat{y}\_i +
\varepsilon^\*\_i\\ where \\\varepsilon^\*\\ is sampled with replacement
from the model's residuals; the model is refit; \\\eta^2_G\\ is
recomputed; and the percentile interval, the empirical quantiles of the
`B` bootstrap values (Efron & Tibshirani, 1993), is reported. The
percentile interval is the only bootstrap interval offered; there is no
bias-corrected and accelerated (BCa) variant. Bootstrap results vary
from run to run; supply `seed` for reproducibility. Requires a
single-stratum `aov`/`lm` fit. **Not yet supported for `aovlist`
(multi-stratum / within-subjects) fits**, since the bootstrap needs to
respect the subject-level correlation structure, which naive residual
resampling does not. Use `method = "parametric"` for `aovlist` fits.
Coverage has not been broadly validated for this estimand.

**Bootstrap replications.** The limits are quantiles of the `B`
replicated values of each effect, read at `alpha_lower` and
`1 - alpha_upper`. The function requires at least 1000 replications,
more than the 400 that 10 replications beyond each limit take at the 95
percent level; at the 99 percent level that count is 2000 and governs. A
replication in which \\\eta^2_G\\ for some effect is undefined (the
refit fails, or the effect size cannot be computed from it) is dropped
and counted, never replaced, and the interval is computed from the rest.
The result records `B` (requested), `B_used` (the number the interval
rests on), `n_dropped`, and `B_required` (the fewest replications a
request may ask for; that count is held to `B`, not to `B_used`) as
attributes, with `bootstrap_status` (`"complete"`, `"dropped"`,
`"exploratory"`, `"insufficient_count"`, or `"insufficient_share"`). One
warning of class `dmar_bootstrap_replications_dropped` states the counts
when any replication was dropped, with a caution when the dropped
replications are at least the 2.5 percent of `B` that lies beyond one
limit at the 95 percent level, because they could then decide by
themselves where that limit falls. When more than half were dropped, the
limits of every effect are `NA`, the estimate is returned, and the
warning has class `dmar_bootstrap_insufficient`. With several effects,
`B_used` is the fewest replications the limits of any effect rest on.
`options(dmar.bootstrap_exploratory = TRUE)` lets a smaller `B` (2 or
more) run for a first look or a demonstration; every such call gives a
warning of class `dmar_bootstrap_exploratory` and the result is not one
to report. See
[`dmar_bootstrap`](https://yelleknek.github.io/DMAR/reference/dmar_bootstrap.md)
for the rule and its reasons.

**Within-subjects designs (`aovlist`).** Each focal effect's parametric
interval is built from its own stratum's *F* test, and its limits for
\\\lambda\\ are converted with the size of that stratum, as
[`ci_eta_squared`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared.md)
converts them: the stratum's total degrees of freedom for an effect
tested within subjects, and the number of subjects for an effect in the
between-subjects stratum of a mixed design. \\r\_{\text{obs}}\\ includes
the residual sums of squares of the other strata, since every source of
variance that involves subjects belongs to the denominator of
\\\eta^2_G\\ (Olejnik & Algina, 2003, Eq. 5, p. 437; Bakeman, 2005, p.
381). With \\r\_{\text{obs}}\\ at the value it converges to, the
conversion maps each limit for \\\lambda\\ to the corresponding value of
the population \\\eta^2_G\\ of Olejnik and Algina (2003, Eq. 2, p. 436),
the value the estimate converges to as subjects are added. In a one-way
within-subjects design with *n* subjects, *p* occasions, and compound
symmetry, for which Steiger (2004, Eq. 58, p. 177) gives \\\lambda = n
\sum_j (\mu_j - \mu)^2 / \[\sigma^2 (1 - \rho)\]\\, that value is
\\\lambda (1 - \rho)/\[\lambda (1 - \rho) + np\]\\: the variance of the
occasion means, \\\sum_j (\mu_j - \mu)^2 / p\\, over the sum of that
variance and \\\sigma^2\\, which holds the subject, subject by occasion,
and error variances that Olejnik and Algina (2003, p. 443, note 1) place
in the denominator. Given \\\lambda\\, that value still depends on
\\\rho\\, and the distribution of the occasions' *F*, which depends on
\\\lambda\\ alone, carries no further information about \\\rho\\.
Steiger (2004, pp. 176–177) notes that in repeated measures designs
\\\lambda\\ confounds the effect with the correlation among the
observations, and that the root-mean-square standardized effect of such
a design cannot be estimated directly with the exact noncentral methods
unless \\\rho\\ is known. Converted with the number of rows, \\np\\, the
limits would be for a smaller quantity than the one the estimate
estimates: for the three IQ occasions of
[`pygmalion`](https://yelleknek.github.io/DMAR/reference/pygmalion.md)
analyzed in
[`vignette("DMAR", package = "DMAR")`](https://yelleknek.github.io/DMAR/articles/DMAR.md),
930 rows against *F*(2, 618) give the interval \[.019, .037\], below the
estimate of .041; converted with 2 + 618 = 620, the interval is \[.029,
.054\]. As
[`eta_squared_generalized`](https://yelleknek.github.io/DMAR/reference/eta_squared_generalized.md)
does, the function stops when an effect is tested in more than one
stratum, as happens when observations are missing within subjects or a
covariate varies within them, since no one stratum's sum of squares is
then the effect's.

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
[`ci_mahalanobis()`](https://yelleknek.github.io/DMAR/reference/ci_mahalanobis.md),
[`ci_omega_squared()`](https://yelleknek.github.io/DMAR/reference/ci_omega_squared.md),
[`ci_pvaf()`](https://yelleknek.github.io/DMAR/reference/ci_pvaf.md),
[`ci_ratio_of_means()`](https://yelleknek.github.io/DMAR/reference/ci_ratio_of_means.md),
[`ci_ratio_of_means_paired()`](https://yelleknek.github.io/DMAR/reference/ci_ratio_of_means_paired.md),
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
#>  effect          eta_squared_generalized ci_lower ci_upper method
#>  treatment       0.020                   <NA>     <NA>     none  
#>  grade           0.044                   <NA>     <NA>     none  
#>  treatment:grade 0.019                   <NA>     <NA>     none  

# The parametric approximation maps a noncentral F interval for the focal
# effect through the observed sums of squares. Warnings mark the method as
# preliminary and report that the interaction's lower limit is clamped at
# 0; treat the limits accordingly.
ci_eta_squared_generalized(fit, observed = "grade", method = "parametric")
#> Warning: Parametric CI uses an approximate transformation that maps the partial eta squared CI through the observed sums of squares for measured factors. Coverage properties have not been broadly validated; results should be treated as preliminary and independently evaluated.
#> Warning: The observed F_value is below the alpha_lower critical value of the central F-distribution, so the lower confidence limit on generalized eta squared is 0.
#>  effect          eta_squared_generalized ci_lower ci_upper method    
#>  treatment       0.020                   0.001    0.058    parametric
#>  grade           0.044                   0.001    0.081    parametric
#>  treatment:grade 0.019                   0.000    0.041    parametric
#> 
#> Confidence level: 95%

# The third option is a residual bootstrap, which refits the model once per
# replication. B = 100 keeps the example quick and is below the 1000 the
# function requires; the option lets it run and the warning says so. A
# reported interval deserves B = 10000 or more. The seed makes the
# interval reproducible and leaves the random number stream of the
# surrounding session as it was.
old <- options(dmar.bootstrap_exploratory = TRUE)
ci_eta_squared_generalized(fit, observed = "grade",
                           method = "bootstrap", B = 100, seed = 113)
#> Warning: Bootstrap CI resamples the fitted model's residuals nonparametrically (B = 100). Coverage has not been broadly validated for generalized eta squared; results should be treated as preliminary.
#> Warning: ci_eta_squared_generalized(): B = 100 is below the 1000 replications a 95% interval is reported from (10 beyond each limit). The limits are shown because options(dmar.bootstrap_exploratory = TRUE) is set. With 100 replications each limit rests on about 2.5 of them and will move from seed to seed; use the result to explore or to demonstrate, and report an interval only from B of at least 1000. B, B_used, and B_required are recorded in the result.
#>  effect          eta_squared_generalized ci_lower ci_upper method   
#>  treatment       0.020                   0.001    0.059    bootstrap
#>  grade           0.044                   0.022    0.108    bootstrap
#>  treatment:grade 0.019                   0.007    0.083    bootstrap
#> 
#> Confidence level: 95%
#> 
#> Bootstrap: 100 replications, 100 used (0 dropped; at least 1000 must be requested). Exploratory: below the number required for a reported interval.
options(old)

# The parametric treatment row from summary statistics, as an ANOVA table
# reports them: the sums of squares of treatment, of its error term (the
# residual), and of the measured sources in its denominator (grade and
# the treatment by grade interaction), with the degrees of freedom and
# the number of children for the interval.
tab <- anova(fit)
ci_eta_squared_generalized(
  SS_effect   = tab["treatment", "Sum Sq"],
  SS_observed = tab[c("grade", "treatment:grade"), "Sum Sq"],
  SS_error    = tab["Residuals", "Sum Sq"],
  df_effect   = tab["treatment", "Df"], df_error = tab["Residuals", "Df"],
  N = nobs(fit), method = "parametric")
#> Warning: Parametric CI uses an approximate transformation that maps the partial eta squared CI through the observed sums of squares for measured factors. Coverage properties have not been broadly validated; results should be treated as preliminary and independently evaluated.
#>  effect  eta_squared_generalized ci_lower ci_upper method    
#>  overall 0.020                   0.001    0.058    parametric
#> 
#> Confidence level: 95%

# Within-subjects ANOVA. The parametric CI uses each effect's own stratum
# and converts its limits with that stratum's size.
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
#> Warning: Parametric CI uses an approximate transformation that maps the partial eta squared CI through the observed sums of squares for measured factors. Coverage properties have not been broadly validated; results should be treated as preliminary and independently evaluated.
#>  effect eta_squared_generalized stratum      ci_lower ci_upper method    
#>  time   0.097                   subject:time 0.008    0.203    parametric
#> 
#> Confidence level: 95%
```
