# Maximum Likelihood Multiple Regression

Fits a multiple regression model by maximum likelihood, with full
information likelihood handling of missing values by default. The
formula interface and S3 methods mirror
[`lm`](https://rdrr.io/r/stats/lm.html) so that calls such as
[`coef()`](https://rdrr.io/r/stats/coef.html),
[`vcov()`](https://rdrr.io/r/stats/vcov.html),
[`confint()`](https://rdrr.io/r/stats/confint.html),
[`summary()`](https://rdrr.io/r/base/summary.html),
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
[`residuals()`](https://rdrr.io/r/stats/residuals.html), and
[`predict()`](https://rdrr.io/r/stats/predict.html) continue to work.
Confidence intervals default to the likelihood ratio (profile) form;
Wald and bootstrap variants are also available.

## Usage

``` r
mlmr(
  formula,
  data,
  missing = c("fiml", "ml", "listwise", "pairwise", "available.cases"),
  ci_method = c("profile", "wald", "boot"),
  conf_level = 0.95,
  B = 1000L,
  boot_type = c("ordinary", "bollen.stine"),
  boot_seed = NULL,
  estimator = c("ML", "MLR", "MLM", "GLS"),
  se = NULL,
  fixed_x = FALSE,
  auxiliary = NULL,
  effect_sizes = TRUE,
  enforce_es_bounds = FALSE,
  ...,
  wald_df = Inf
)
```

## Arguments

- formula:

  A two-sided [`formula`](https://rdrr.io/r/stats/formula.html) of the
  form `y ~ x1 + x2 + ...`. Factor predictors, interactions (`x1 * x2`),
  polynomial terms (`poly(x, 2)`), and transformations (`I(x^2)`) are
  supported through the usual
  [`model.matrix`](https://rdrr.io/r/stats/model.matrix.html) expansion.
  An intercept-only formula, `y ~ 1`, fits the null model (the mean and
  the residual variance only); it is the natural restricted model in a
  model comparison and pairs with
  [`anova()`](https://rdrr.io/r/stats/anova.html) for a likelihood ratio
  test against a fuller model.

- data:

  A `data.frame` containing the variables in `formula`. Rows with
  missing values on *any* of the modeled variables are retained (under
  `missing = "fiml"`) and contribute to the likelihood through whichever
  components are observed.

- missing:

  Character; how missing values are handled. The default `"fiml"`
  (equivalently `"ml"` in lavaan) uses the full information maximum
  likelihood over all rows that have at least one observed value on the
  modeled variables. Other choices are `"listwise"` (drop rows with any
  missing modeled variable), `"pairwise"` (sample moments computed
  pairwise; not recommended with this model), and `"available.cases"`.

- ci_method:

  Character; the method for confidence intervals on the regression
  coefficients. `"profile"` (default) inverts the likelihood ratio test
  on each parameter through a sequence of constrained refits. `"wald"`
  returns the symmetric estimate \\\pm\\ \\z\_{1 - \alpha/2}\\ standard
  error interval. `"boot"` resamples the rows of `data` `B` times,
  refits the model on each resample, and reports the percentile interval
  of the resampled slopes.

- conf_level:

  Desired level of confidence (the complement of the Type I error rate).
  Defaults to `0.95`.

- B:

  Number of bootstrap resamples when `ci_method = "boot"`, a single
  whole number. Default `1000`, because every resample refits the model.
  The fewest accepted is 400 at the default `conf_level` (200 at 0.90,
  and 2000 at 0.99, which the default does not meet); a smaller `B` is
  refused before anything is computed (class `dmar_count_invalid`). See
  *Bootstrap replications* under *Details* and
  [`dmar_bootstrap`](https://yelleknek.github.io/DMAR/reference/dmar_bootstrap.md).

- boot_type:

  Character; `"ordinary"` (default) for the nonparametric resampling of
  rows, or `"bollen.stine"` for the Bollen and Stine (1992) model-based
  bootstrap implemented in lavaan.

- boot_seed:

  Integer or `NULL`; optional seed for the bootstrap resampling RNG. The
  default `NULL` draws the resamples from, and advances, the current
  random number stream, so successive bootstrap calls draw fresh
  resamples; supply an integer to make the bootstrap reproducible. When
  supplied, the function seeds the RNG locally and restores the prior
  state on exit so the user's global RNG is not polluted.

- estimator:

  Character; the lavaan estimator. One of `"ML"` (default), `"MLR"`,
  `"MLM"`, or `"GLS"`. `ML` is standard maximum likelihood under
  conditional normality of *Y*. `MLR` is robust maximum likelihood with
  Yuan-Bentler scaled standard errors and a Yuan-Bentler scaled test
  statistic, recommended when the conditional distribution of *Y*
  departs from normality and FIML is in use (Yuan & Bentler, 2000).
  `MLM` is Satorra- Bentler scaled \\\chi^2\\ statistics under complete
  data (Satorra & Bentler, 1994). `GLS` is normal theory generalized
  least squares. The ordinal-data estimators (`"DWLS"`, `"WLS"`, `"ULS"`
  and their robust variants) are intentionally not exposed; for ordinal
  outcomes, fit a different model class.

- se:

  Character or `NULL`; the standard error type passed to lavaan. When
  `NULL` (default), the standard error type is chosen automatically from
  `estimator`: `"ML"` and `"GLS"` use `"standard"`, `"MLR"` uses
  `"robust.huber.white"` (Huber-White heteroskedasticity consistent),
  and `"MLM"` uses `"robust.sem"` (Satorra-Bentler). Override only when
  the default does not match the desired analysis. Other values include
  `"robust"`, `"first.order"`, and `"none"`.

- fixed_x:

  Logical; whether to treat predictors as fixed (not modeled jointly) or
  as random (jointly modeled). Defaults to `FALSE`, which is required
  for the full information likelihood to use rows with missing
  predictors. Set to `TRUE` only when the predictors are fully observed
  and the user wants the `lm`-style conditional model.

- auxiliary:

  Character vector of variable names in `data` to include as auxiliary
  variables, or `NULL` (default) for none. Auxiliaries are entered as
  saturated correlates (Graham, 2003): correlated with the outcome
  residual, every predictor, and each other, but not as predictors,
  while the predictors keep their free covariances, so the model is
  saturated, the regression coefficients keep their meaning, and the
  full information maximum likelihood draws on the auxiliaries' observed
  values (the inclusive analysis strategy; Collins, Schafer, & Kam,
  2001). A name in `auxiliary` must be numeric, must be present in
  `data`, must not appear in `formula`, and requires `fixed_x = FALSE`.
  Auxiliaries are refused with `missing = "listwise"`, which would drop
  every row on which an auxiliary is missing and so change the sample
  the coefficients describe, and with `missing = "pairwise"` or
  `"available.cases"`, which compute the moments from each pair of
  variables on its own, so that an auxiliary would leave every estimate
  unchanged.

- effect_sizes:

  Logical; whether to compute regression effect sizes (standardized
  betas, semi-partial \\R^2\\, Cohen's \\f^2\\ per predictor, and the
  overall LR omnibus test). Defaults to `TRUE`. Disabling saves \\p +
  1\\ additional lavaan refits.

- enforce_es_bounds:

  Logical; whether to clamp semi-partial \\R^2\\ and Cohen's \\f^2\\
  estimates to their theoretical lower bound of zero. Defaults to
  `FALSE`: the raw maximum likelihood estimates of
  \\R^2\_{\text{reduced}}\\ and \\R^2\_{\text{full}}\\ are reported
  as-is, and the difference can be slightly negative as a finite-sample
  artifact when the two are nearly equal. Setting to `TRUE` replaces any
  negative value with zero, which yields an estimate that respects the
  parameter space but is no longer the maximum likelihood estimate. When
  the clamp fires, the affected rows of the returned `effect_sizes`
  table carry the attribute `"clamped"` for diagnostics.

- ...:

  Additional arguments forwarded to
  [`lavaan`](https://rdrr.io/pkg/lavaan/man/lavaan.html).

- wald_df:

  The reference distribution of the Wald interval and of the *p*-values:
  `Inf` (the default), the normal distribution; a single positive
  number, the degrees of freedom of a *t* distribution for every
  coefficient; or `"missing_information"`, von Hippel's (2016) degrees
  of freedom \\\max\\3, N(1 - \gamma_j) - k\\\\ for each coefficient,
  from its fraction of missing information \\\gamma_j\\. A finite value
  or the keyword requires `ci_method = "wald"`. The keyword is an error
  where, with values missing, lavaan fit the complete cases in place of
  the FIML fit that `missing = "fiml"` asks for, as it does for
  `se = "robust.sem"` or a test statistic of the Satorra-Bentler kind
  passed through `...`; with no value missing that fit is the FIML fit,
  and the keyword gives \\N - k\\. It follows `...`, so it must be named
  in full. See *Details*.

## Value

An object of class `"mlmr"`, a list with components modeled on the
structure of an `lm` fit:

- `call`:

  The matched call.

- `formula`:

  The model formula.

- `terms`:

  The terms object.

- `model`:

  The model frame (with missing values preserved when
  `missing = "fiml"`).

- `coefficients`:

  Named numeric vector of regression coefficients, with `(Intercept)`
  first when an intercept is in the formula.

- `vcov`:

  The variance-covariance matrix of the regression coefficients,
  returned by [`vcov()`](https://rdrr.io/r/stats/vcov.html).

- `ci`:

  A two-column matrix (`lower`, `upper`) of confidence limits in the
  order of `coefficients`.

- `ci_method`:

  Which method was used to compute `ci`.

- `B`, `B_used`, `n_dropped`, `B_required`, `bootstrap_status`:

  Present only with `ci_method = "boot"`: the resamples requested, the
  number the intervals rest on, their difference, the fewest resamples a
  request may ask for, and the verdict (`"complete"`, `"dropped"`,
  `"exploratory"`, `"insufficient_count"`, or `"insufficient_share"`;
  under the last two the limits are `NA`). See *Bootstrap replications*
  under Details.

- `conf_level`:

  The confidence level used.

- `coef_table`:

  A `data.frame` with columns `term`, `estimate`, `se`,
  `test_statistic`, `p_value`, `ci_lower`, `ci_upper`, carrying the
  confidence level in its `conf_level` attribute, so that
  [`results_sentence`](https://yelleknek.github.io/DMAR/reference/dmar_output_helpers.md)
  on a row reports it.

- `sigma2`:

  Residual variance of *Y*, on the maximum likelihood scale (divisor
  *N*, not \\N - p - 1\\).

- `R2`:

  Model implied squared multiple correlation, \\1 - \hat{\sigma}^2_e /
  \hat{\sigma}^2_Y\\, where both variances come from the FIML estimated
  model implied covariance matrix.

- `adj_R2`:

  Adjusted \\R^2\\ using the number of complete cases (`N_complete`, the
  rows complete on the outcome and every predictor, which are the rows
  that identify the regression) and the number of slopes; the lavaan
  reported `N` can be larger under FIML because it counts rows that
  inform only the predictor distribution.

- `logLik`:

  The log likelihood at the maximum, with attributes `df` and `nobs` for
  compatibility with [`stats::AIC`](https://rdrr.io/r/stats/AIC.html)
  and [`stats::BIC`](https://rdrr.io/r/stats/AIC.html).

- `N`:

  Sample size used by lavaan (rows with at least one observed value when
  `missing = "fiml"`; rows with no missing values when
  `missing = "listwise"`).

- `N_complete`:

  Number of rows that are complete on all modeled variables.

- `fitted.values`:

  Vector of fitted values, length `nrow(data)`, with `NA` for rows
  missing any predictor.

- `residuals`:

  Vector of residuals (`y - fitted`), with `NA` where `y` or any
  predictor was missing.

- `lavaan_fit`:

  The underlying lavaan fit object, returned for advanced users who want
  to apply lavaan accessors directly.

The object carries the attribute `wald_df`, the degrees of freedom of
the *t* reference of the Wald interval and the *p*-values: `Inf` for the
normal reference, the number supplied, or under
`wald_df = "missing_information"` a named vector with one value per
coefficient. In the last case the attribute `wald_df_table` is the
[`df_missing`](https://yelleknek.github.io/DMAR/reference/df_missing.md)
table the values came from, with each coefficient's fraction of missing
information, \\N\\, \\k\\, and whether the floor applied.

## Details

**Why a separate function from `lm`.** `lm` uses ordinary least squares
and listwise deletes any row with a missing value on the outcome or on
any predictor. Two situations motivate a maximum likelihood alternative.

First, when a predictor is missing on some rows, listwise deletion can
be biased if the missingness mechanism depends on other observed
variables (the missing at random or MAR pattern). Full information
maximum likelihood (FIML) jointly models the distribution of \\(Y, X_1,
\ldots, X_p)\\, the outcome and the \\p\\ predictors, and yields
consistent regression estimates under MAR, while listwise estimates can
be biased away from the population values (Enders, 2010; Schafer &
Graham, 2002).

Second, the joint likelihood estimates the predictor distribution as
well, so quantities that depend on the predictor moments (the
standardized coefficients, the model implied \\R^2\\, and the predictor
variances and covariances) draw on every row with an observed predictor,
not only the rows that are complete on the outcome. When the missing
values are confined to the outcome, however, the unstandardized slopes
and their standard errors match listwise deletion up to the maximum
likelihood \\N\\ versus \\N - p - 1\\ variance divisor. The rows with an
observed predictor but a missing outcome inform the marginal
distribution of *X*, not the conditional distribution of *Y* given *X*
that identifies the slopes, so they leave the slope estimates and their
conditional-model standard errors unchanged.

Full information estimation differs from listwise deletion in the slope
estimates themselves when any predictors are missing on some rows, or
when auxiliary variables that correlate with the outcome or with the
missingness mechanism are supplied through `auxiliary` (see below). The
bootstrap gives inference that does not depend on the multivariate
normality assumption.

**Auxiliary variables.** A variable that is not part of the regression
but is correlated with the outcome or with the missingness can be
supplied through `auxiliary`. Auxiliaries are added to the model as
saturated correlates (Graham, 2003): each one is correlated with the
residual of the outcome, with every predictor, and with every other
auxiliary, but is never entered as a predictor, and the predictors keep
their free means, variances, and covariances. The model over the
outcome, the predictors, and the auxiliaries is then saturated (its chi
square is zero on zero degrees of freedom), so the auxiliaries impose no
restriction and the focal regression coefficients keep their meaning: on
complete data they equal those of
[`lm()`](https://rdrr.io/r/stats/lm.html) to working precision, and with
missing values they are the regression implied by the full information
maximum likelihood estimates of the means and covariances of every
variable, auxiliaries included, which
[`cov_fiml`](https://yelleknek.github.io/DMAR/reference/cov_fiml.md)
computes. The full information maximum likelihood uses the auxiliaries'
observed values to make the MAR assumption hold conditional on more of
the observed data and to recover information that listwise deletion
discards. This is the inclusive analysis strategy of Collins, Schafer,
and Kam (2001): a variable that predicts the missingness or the
incomplete outcome belongs in the analysis even when it is of no
substantive interest. Auxiliary variables must be numeric and require
`fixed_x = FALSE` (the default).

**Why likelihood ratio confidence intervals by default.** Wald intervals
(point estimate \\\pm\\ \\z\_{1 - \alpha/2}\\ standard error) are
symmetric by construction and assume the sampling distribution of the
estimator is approximately normal over the relevant range. Likelihood
ratio intervals invert the likelihood ratio test directly: an interval
contains every value of the parameter that would not be rejected at
level \\\alpha\\. Likelihood ratio intervals are invariant under
monotone reparameterizations, often have better coverage in small
samples, and respect parameter boundaries (Pawitan, 2001). The cost is
computational: each parameter requires a sequence of refits with that
parameter constrained. The examples below ask for the Wald and bootstrap
intervals so the help page stays quick; a reported analysis leaves
`ci_method` at its default.

**A *t* reference distribution for the Wald interval.** The Wald
interval and every *p*-value refer the ratio of a coefficient to its
standard error to the normal distribution, which treats the standard
error as known. That holds only in the limit: in small samples, and more
so when values are missing, the normal intervals are too short. In von
Hippel's (2016, pp. 432–433) simulations of bivariate normal data with
half of one variable missing, the normal FIML intervals for the
intercept and slope of a regression covered 5 to 6 points below the
nominal 95% at *n* = 25, and 0 to 3 points below it at *n* = 100. He
kept the maximum likelihood estimate and its standard error and replaced
the normal quantile with a *t* quantile (his Eq. 29), on degrees of
freedom that take the missing information into account. With
`ci_method = "wald"`, `wald_df` does this. A number refers every
coefficient's interval and *p*-value to the *t* distribution on that
many degrees of freedom. The keyword `"missing_information"` applies von
Hippel's rule for regression coefficients (his Eqs. 31 and 34) to each
coefficient separately, \$\$\nu_j = \max\\3, N(1 - \gamma_j) - k\\,\$\$
the complete data degrees of freedom \\N - k\\ of a regression with
\\k\\ parameters (the slopes and the intercept) with the sample size
replaced by the coefficient's effective sample size, where \\N\\ is the
number of cases the fit used and \\\gamma_j\\ the coefficient's fraction
of missing information. The fraction is estimated from the fit as
[`fraction_missing_information`](https://yelleknek.github.io/DMAR/reference/fraction_missing_information.md)
estimates it, by the unstructured information of Chen and Savalei
(2021), which cannot be negative. The model leaves the means and
covariances of the outcome, the predictors, and any auxiliary variables
unrestricted (it is saturated), and the three estimators of Chen and
Savalei then coincide up to numerical error (Chen & Savalei, 2021, pp.
3, 9), so the fraction agrees with lavaan's
`parameterEstimates(fit, fmi = TRUE)` to the precision of lavaan's own
computation. Unlike lavaan's, its computation does not depend on the
units of the variables: an outcome measured in dollars, with a standard
deviation of 100,000,000, gets the degrees of freedom it gets in
millions of dollars (see
[`fraction_missing_information`](https://yelleknek.github.io/DMAR/reference/fraction_missing_information.md)).
It is computed at the fit's estimates, which lavaan can leave short of
the maximum of the likelihood for variables measured in small units
(with every standard deviation near .003 the intercept's degrees of
freedom moved by .15); the call then warns, and `rescale_data = TRUE`,
passed through `...` to lavaan, is the remedy. A fraction that cannot be
computed stops the call with the reason, rather than leave a coefficient
with no degrees of freedom. A coefficient estimated mostly from fully
observed variables keeps nearly all of the complete data degrees of
freedom, and with no missing values every coefficient has \\N - k\\. The
floor of 3 keeps a *t* distribution with finite variance when the
estimated fraction is near one.
[`df_missing`](https://yelleknek.github.io/DMAR/reference/df_missing.md)
computes the same degrees of freedom from any fractions of missing
information, and its help page describes both of von Hippel's rules and
the evidence for them; the table it returns for the fit is the
`wald_df_table` attribute. Two cautions. The *t* interval is an
approximation, since the small-sample distribution of a maximum
likelihood estimate from incomplete data is not exactly *t* (von Hippel,
2016, p. 431). And the standard error is the maximum likelihood one,
whose variance estimate divides by \\N\\ rather than \\N - k\\, so with
complete data the interval on \\N - k\\ degrees of freedom is shorter
than the interval [`lm()`](https://rdrr.io/r/stats/lm.html) reports by
the factor \\\sqrt{(N - k)/N}\\: von Hippel's construction pairs the
maximum likelihood standard error with the *t* quantile, and it is
reproduced here as he published it, without a variance correction. The
default, `Inf`, is the normal reference. The `test_statistic` column
holds the ratio of the estimate to its standard error, the same number
under either reference.
[`summary()`](https://rdrr.io/r/base/summary.html) heads the printed
coefficient table `t value` and `Pr(>|t|)` under a *t* reference, as
[`summary.lm()`](https://rdrr.io/r/stats/summary.lm.html) does, and
names the degrees of freedom directly beneath it.

**The bootstrap interval.** With `ci_method = "boot"` the rows of `data`
are resampled with replacement `B` times (1000 by default) and the model
is refit on each resample; with `boot_type = "bollen.stine"` the
resamples are instead drawn from data transformed to satisfy the fitted
model (Bollen & Stine, 1992), a model-based bootstrap. Only the
percentile interval is offered: each coefficient's limits are the
empirical quantiles of its resampled estimates (Efron & Tibshirani,
1993); there is no bias-corrected and accelerated (BCa) or bootstrap
standard error variant. Raising `B` tightens the Monte Carlo error of
the reported limits. Bootstrap results vary from run to run; supply
`boot_seed` for reproducibility.

**Bootstrap replications.** Each coefficient's limits are quantiles of
its `B` resampled estimates, so the interval is run only when `B` puts
at least 10 resamples beyond each limit (400 at the 95 percent level). A
replication in which some coefficient is undefined (the refit fails or
does not converge, or returns a coefficient that is not a finite number)
is dropped and counted, never replaced, and the interval is computed
from the rest. The result records `B` (requested), `B_used` (the number
the interval rests on), `n_dropped`, and `B_required` (the fewest
replications a request may ask for; that count is held to `B`, not to
`B_used`) as elements of the fitted object, and
[`print()`](https://rdrr.io/r/base/print.html) and
[`summary()`](https://rdrr.io/r/base/summary.html) state them in one
line, with `bootstrap_status` (`"complete"`, `"dropped"`,
`"exploratory"`, `"insufficient_count"`, or `"insufficient_share"`). One
warning of class `dmar_bootstrap_replications_dropped` states the counts
when any replication was dropped, with a caution when the dropped
replications are at least the 2.5 percent of `B` that lies beyond one
limit at the 95 percent level, because they could then decide by
themselves where that limit falls. When more than half were dropped, the
limits of every coefficient are `NA`, the estimate is returned, and the
warning has class `dmar_bootstrap_insufficient`. Such a resample is
dropped for every coefficient, so all intervals rest on the same
resamples. `options(dmar.bootstrap_exploratory = TRUE)` lets a smaller
`B` (2 or more) run for a first look or a demonstration; every such call
gives a warning of class `dmar_bootstrap_exploratory` and the result is
not one to report. See
[`dmar_bootstrap`](https://yelleknek.github.io/DMAR/reference/dmar_bootstrap.md)
for the rule and its reasons.

**Model representation.** Internally the model is fit through lavaan as
a structural equation model in which *Y* is regressed on the predictors
and (when `fixed_x = FALSE`) the predictor distribution is also
estimated. With `fixed_x = FALSE` and complete data, point estimates of
the slopes are identical to `lm` and the residual variance differs only
by the usual \\N\\ versus \\N - p - 1\\ divisor.

**Caveats.** The function assumes that, conditional on the modeled
predictors, the dependent variable is normally distributed with constant
variance. Missingness is assumed to be at most MAR; missing not at
random patterns require selection or pattern mixture models outside the
scope of this function. Factor predictors and interactions are expanded
through [`model.matrix`](https://rdrr.io/r/stats/model.matrix.html) and
entered as numeric covariates, so the same caveats about dummy variable
encoding that apply to `lm` apply here as well.

## References

Bollen, K. A., & Stine, R. A. (1992). Bootstrapping goodness of fit
measures in structural equation models. *Sociological Methods &
Research, 21*, 205–229.
[doi:10.1177/0049124192021002004](https://doi.org/10.1177/0049124192021002004)

Chen, L., & Savalei, V. (2021). Three sample estimates of fraction of
missing information from full information maximum likelihood. *Frontiers
in Psychology, 12*, 667802.
[doi:10.3389/fpsyg.2021.667802](https://doi.org/10.3389/fpsyg.2021.667802)

Collins, L. M., Schafer, J. L., & Kam, C.-M. (2001). A comparison of
inclusive and restrictive strategies in modern missing data procedures.
*Psychological Methods, 6*(4), 330–351.
[doi:10.1037/1082-989X.6.4.330](https://doi.org/10.1037/1082-989X.6.4.330)

Efron, B., & Tibshirani, R. J. (1993). *An introduction to the
bootstrap*. New York, NY: Chapman & Hall/CRC.

Enders, C. K. (2010). *Applied missing data analysis*. New York, NY:
Guilford Press.

Graham, J. W. (2003). Adding missing-data-relevant variables to
FIML-based structural equation models. *Structural Equation Modeling,
10*(1), 80–100.
[doi:10.1207/S15328007SEM1001_4](https://doi.org/10.1207/S15328007SEM1001_4)

Pawitan, Y. (2001). *In all likelihood: Statistical modelling and
inference using likelihood*. Oxford, UK: Oxford University Press.

Rosseel, Y. (2012). lavaan: An R package for structural equation
modeling. *Journal of Statistical Software, 48*(2), 1–36.
[doi:10.18637/jss.v048.i02](https://doi.org/10.18637/jss.v048.i02)

Satorra, A., & Bentler, P. M. (1994). Corrections to test statistics and
standard errors in covariance structure analysis. In A. von Eye & C. C.
Clogg (Eds.), *Latent variables analysis: Applications for developmental
research* (pp. 399–419). Sage.

Schafer, J. L., & Graham, J. W. (2002). Missing data: Our view of the
state of the art. *Psychological Methods, 7*, 147–177.
[doi:10.1037/1082-989X.7.2.147](https://doi.org/10.1037/1082-989X.7.2.147)

von Hippel, P. T. (2016). New confidence intervals and bias comparisons
show that maximum likelihood can beat multiple imputation in small
samples. *Structural Equation Modeling, 23*(3), 422–437.
[doi:10.1080/10705511.2015.1047931](https://doi.org/10.1080/10705511.2015.1047931)

Yuan, K.-H., & Bentler, P. M. (2000). Three likelihood-based methods for
mean and covariance structure analysis with nonnormal missing data.
*Sociological Methodology, 30*(1), 165–200.
[doi:10.1111/0081-1750.00078](https://doi.org/10.1111/0081-1750.00078)

## See also

[`lm`](https://rdrr.io/r/stats/lm.html),
[`sem`](https://rdrr.io/pkg/lavaan/man/sem.html),
[`lavaan`](https://rdrr.io/pkg/lavaan/man/lavaan.html),
[`ci_reg_coef`](https://yelleknek.github.io/DMAR/reference/ci_reg_coef.md),
[`ci_rc`](https://yelleknek.github.io/DMAR/reference/ci_rc.md),
[`ci_src`](https://yelleknek.github.io/DMAR/reference/ci_src.md).

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Complete data: the maximum likelihood estimates agree with lm() to
# working precision. The fit asks for the Wald interval; the choice
# among the intervals is taken up below.
fit_mlmr <- mlmr(t6_paragraph_comprehension ~ t5_general_information +
                   t9_word_meaning,
                 data = holzinger_swineford, ci_method = "wald")
fit_lm   <- lm(t6_paragraph_comprehension ~ t5_general_information +
                 t9_word_meaning,
               data = holzinger_swineford)
cbind(mlmr = coef(fit_mlmr), lm = coef(fit_lm))
#>                              mlmr         lm
#> (Intercept)            2.38038350 2.38038350
#> t5_general_information 0.08482737 0.08482737
#> t9_word_meaning        0.21956217 0.21956217

fit_mlmr
#> 
#> Call:
#> mlmr(formula = t6_paragraph_comprehension ~ t5_general_information + 
#>     t9_word_meaning, data = holzinger_swineford, ci_method = "wald")
#> 
#> Coefficients:
#>            (Intercept)  t5_general_information         t9_word_meaning  
#>                2.38038                 0.08483                 0.21956  
#> 
confint(fit_mlmr)
#>                             2.5 %    97.5 %
#> (Intercept)            1.44486986 3.3158971
#> t5_general_information 0.05263738 0.1170174
#> t9_word_meaning        0.16759660 0.2715277

# summary() adds the intervals, the per-predictor semi-partial R^2
# and Cohen's f^2, and the omnibus likelihood ratio test of all
# slopes equal to zero.
summary(fit_mlmr)
#> 
#> Call:
#> mlmr(formula = t6_paragraph_comprehension ~ t5_general_information + 
#>     t9_word_meaning, data = holzinger_swineford, ci_method = "wald")
#> 
#> Missing data: fiml | Estimator: ML | SE: standard
#> Sample size (used by lavaan): 301   Complete cases: 301
#> 
#> Coefficients:
#>                        Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)             2.38038    0.47731   4.987 6.13e-07 ***
#> t5_general_information  0.08483    0.01642   5.165 2.41e-07 ***
#> t9_word_meaning         0.21956    0.02651   8.281  < 2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Confidence intervals (method: wald, level = 0.95):
#>                          2.5 % 97.5 %
#> (Intercept)            1.44487 3.3159
#> t5_general_information 0.05264 0.1170
#> t9_word_meaning        0.16760 0.2715
#> 
#> Per-predictor effect sizes (semi-partial R^2 and f^2):
#>                           sr^2 Cohen's f^2
#> t5_general_information 0.04101     0.08863
#> t9_word_meaning        0.10542     0.22783
#> 
#> Residual std. error (ML): 2.372 on 298 residual degrees of freedom
#> Model implied R-squared: 0.5373,   Adjusted R-squared: 0.5342,   Cohen's f^2: 1.161
#> Omnibus likelihood ratio test (all slopes = 0): chi square = 232 on 2 df, p = < 2.2e-16
#> Log likelihood: -2792   AIC: 5602   BIC: 5635

# The interval menu is profile, Wald, and bootstrap. The default,
# ci_method = "profile", inverts the likelihood ratio test one
# parameter at a time through a sequence of constrained refits; it
# is what a reported interval deserves, and leaving ci_method at its
# default asks for it. The fits on this page ask for the Wald or the
# bootstrap interval because those refits take longer than a help
# page should. The bootstrap resamples rows and takes percentile
# limits; it is what to ask for when the normality the likelihood
# assumes is doubtful. Every resample refits the model. B = 10 keeps
# the example quick and is below the 400 a reported interval needs;
# the option lets it run and the warning says so. A reported interval
# deserves the default B = 1000. boot_seed fixes the resamples, so
# the limits are reproducible rather than moving from run to run, and
# effect_sizes = FALSE skips the effect size refits, which
# summary() above already showed.
old <- options(dmar.bootstrap_exploratory = TRUE)
fit_boot <- mlmr(t6_paragraph_comprehension ~ t5_general_information +
                   t9_word_meaning, data = holzinger_swineford,
                 ci_method = "boot", B = 10, boot_seed = 113,
                 effect_sizes = FALSE)
#> Warning: mlmr(): B = 10 is below the 400 replications a 95% interval is reported from (10 beyond each limit). The limits are shown because options(dmar.bootstrap_exploratory = TRUE) is set. With 10 replications each limit rests on about 0.3 of them and will move from seed to seed; use the result to explore or to demonstrate, and report an interval only from B of at least 400. B, B_used, and B_required are recorded in the result.
options(old)
confint(fit_boot)
#>                             2.5 %    97.5 %
#> (Intercept)            1.71871728 2.9437862
#> t5_general_information 0.06784019 0.1001241
#> t9_word_meaning        0.20867244 0.2509053

# Missing values on a predictor are where maximum likelihood and
# least squares part company. The full information likelihood keeps
# every row that carries information; listwise deletion keeps only
# the rows that are complete. The Holzinger and Swineford battery
# carries real missingness for this: the revised second-form test
# t26_flags was administered to only 145 of the 301 students, so a
# model using it loses more than half the sample under listwise
# deletion while the full information fit keeps all 301 rows. Notice
# the two sample sizes, and that the two fits do not agree on the
# intercept.
fit_fiml <- mlmr(t6_paragraph_comprehension ~ t7_sentence +
                   t26_flags, data = holzinger_swineford,
                 ci_method = "wald", effect_sizes = FALSE)
fit_lwd  <- mlmr(t6_paragraph_comprehension ~ t7_sentence +
                   t26_flags, data = holzinger_swineford,
                 missing = "listwise",
                 ci_method = "wald", effect_sizes = FALSE)
rbind(FIML = coef(fit_fiml), listwise = coef(fit_lwd))
#>          (Intercept) t7_sentence  t26_flags
#> FIML      -0.4335481   0.4719182 0.04019371
#> listwise  -0.8517917   0.4977943 0.03914096
c(N_fiml = nobs(fit_fiml), N_listwise = nobs(fit_lwd))
#>     N_fiml N_listwise 
#>        301        145 

# An auxiliary variable is not a predictor. The complete speed test
# t13_straight_and_curved_capitals enters as a saturated correlate,
# correlated with the outcome residual and with the predictors, so
# the likelihood can draw on it for the rows where t26_flags is
# missing while the coefficients keep their meaning. Continuing from
# the model above, the coefficients move a little and the slope
# standard errors shrink a little: among the 145 students who took
# t26_flags its correlation with t13 is .34, so the auxiliary
# restores only part of the information the missing scores would
# have carried.
fit_aux <- mlmr(t6_paragraph_comprehension ~ t7_sentence +
                  t26_flags, data = holzinger_swineford,
                ci_method = "wald",
                auxiliary = "t13_straight_and_curved_capitals",
                effect_sizes = FALSE)
cbind(no_aux = coef(fit_fiml), aux = coef(fit_aux))
#>                  no_aux         aux
#> (Intercept) -0.43354809 -0.43606955
#> t7_sentence  0.47191824  0.47561312
#> t26_flags    0.04019371  0.03811271
cbind(se_no_aux = sqrt(diag(vcov(fit_fiml))),
      se_aux = sqrt(diag(vcov(fit_aux))))
#>              se_no_aux     se_aux
#> (Intercept) 0.79537553 0.81583455
#> t7_sentence 0.03046276 0.02967365
#> t26_flags   0.02506817 0.02477070

# With most of a variable missing, the normal reference understates the
# uncertainty in a small sample. Here the outcome t26_flags is observed
# for 145 of the 301 students, and wald_df = "missing_information"
# refers each coefficient to t on von Hippel's degrees of freedom,
# N(1 - FMI) - k, fewer than half of the N - k = 298 that complete
# data would give. The summary lists each coefficient's fraction of
# missing information and degrees of freedom under the intervals.
fit_t <- mlmr(t26_flags ~ t1_visual_perception + t4_lozenges,
              data = holzinger_swineford, ci_method = "wald",
              wald_df = "missing_information", effect_sizes = FALSE)
summary(fit_t)
#> 
#> Call:
#> mlmr(formula = t26_flags ~ t1_visual_perception + t4_lozenges, 
#>     data = holzinger_swineford, ci_method = "wald", effect_sizes = FALSE, 
#>     wald_df = "missing_information")
#> 
#> Missing data: fiml | Estimator: ML | SE: standard
#> Sample size (used by lavaan): 301   Complete cases: 145
#> 
#> Coefficients:
#>                      Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)          19.70907    2.54234   7.752 1.76e-12 ***
#> t1_visual_perception  0.38411    0.09366   4.101 7.03e-05 ***
#> t4_lozenges           0.32774    0.07786   4.209 5.03e-05 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> Wald intervals and p-values: t reference on von Hippel's (2016) missing
#> information degrees of freedom, max(3, N(1 - FMI) - k), with N = 301 and k =
#> 3:
#>    term                 fmi   df   
#>    (Intercept)          0.530 138.4
#>    t1_visual_perception 0.536 136.5
#>    t4_lozenges          0.598 118.0
#> 
#> Confidence intervals (method: wald, level = 0.95):
#>                        2.5 %  97.5 %
#> (Intercept)          14.6822 24.7359
#> t1_visual_perception  0.1989  0.5693
#> t4_lozenges           0.1736  0.4819
#> 
#> Residual std. error (ML): 6.945 on 142 residual degrees of freedom
#> Model implied R-squared: 0.3228,   Adjusted R-squared: 0.3132,   Cohen's f^2: 0.4766
#> Log likelihood: -2556   AIC: 5131   BIC: 5164
```
