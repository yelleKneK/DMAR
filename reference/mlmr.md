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
  ...
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

  Integer; number of bootstrap resamples when `ci_method = "boot"`.
  Defaults to `1000`.

- boot_type:

  Character; `"ordinary"` (default) for the nonparametric resampling of
  rows, or `"bollen.stine"` for the Bollen and Stine (1992) model-based
  bootstrap implemented in lavaan.

- boot_seed:

  Integer or `NULL`; optional seed for the bootstrap resampling RNG. The
  default `NULL` leaves the user's current RNG state untouched, so
  successive bootstrap calls draw fresh resamples; supply an integer to
  make the bootstrap reproducible. When supplied, the function seeds the
  RNG locally and restores the prior state on exit so the user's global
  RNG is not polluted.

- estimator:

  Character; the lavaan estimator. One of `"ML"` (default), `"MLR"`,
  `"MLM"`, or `"GLS"`. `ML` is standard maximum likelihood under
  conditional normality of *Y*. `MLR` is robust maximum likelihood with
  Yuan-Bentler scaled standard errors and a Yuan-Bentler scaled test
  statistic, recommended when the conditional distribution of *Y*
  departs from normality and FIML is in use (Yuan & Bentler, 2000).
  `MLM` is Satorra- Bentler scaled \\\chi^2\\ statistics under complete
  data (Satorra & Bentler, 1994). `GLS` is generalized least squares, an
  alternative ML-family estimator that is less commonly used in modern
  practice. The ordinal-data estimators (`"DWLS"`, `"WLS"`, `"ULS"` and
  their robust variants) are intentionally not exposed; for ordinal
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
  residual, every predictor, and each other, but not as predictors, so
  the regression coefficients keep their meaning while the full
  information maximum likelihood draws on the auxiliaries' observed
  values (the inclusive analysis strategy; Collins, Schafer, & Kam,
  2001). A name in `auxiliary` must be numeric, must be present in
  `data`, must not appear in `formula`, and requires `fixed_x = FALSE`.

- effect_sizes:

  Logical; whether to compute regression effect sizes (standardized
  betas, semi-partial \\R^2\\, Cohen's \\f^2\\ per predictor, and the
  overall LR omnibus test). Defaults to `TRUE`. Disabling saves \\K +
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

- `conf_level`:

  The confidence level used.

- `coef_table`:

  A `data.frame` with columns `term`, `estimate`, `se`, `z_value`,
  `p_value`, `ci_lower`, `ci_upper`.

- `sigma2`:

  Residual variance of *Y*, on the maximum likelihood scale (divisor
  *N*, not \\N - K - 1\\).

- `R2`:

  Model implied squared multiple correlation, \\1 - \hat{\sigma}^2_e /
  \hat{\sigma}^2_Y\\, where both variances come from the FIML estimated
  model implied covariance matrix.

- `adj_R2`:

  Adjusted \\R^2\\ using the lavaan reported sample size and the number
  of slopes.

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

## Details

**Why a separate function from `lm`.** `lm` uses ordinary least squares
and listwise deletes any row with a missing value on the outcome or on
any predictor. Two situations motivate a maximum likelihood alternative.

First, when a predictor is missing on some rows, listwise deletion can
be biased if the missingness mechanism depends on other observed
variables (the missing at random or MAR pattern). Full information
maximum likelihood (FIML) jointly models the distribution of \\(Y, X_1,
\ldots, X_K)\\ and yields consistent regression estimates under MAR,
while listwise estimates can be biased away from the population values
(Enders, 2010; Schafer & Graham, 2002).

Second, the joint likelihood estimates the predictor distribution as
well, so quantities that depend on the predictor moments (the
standardized coefficients, the model implied \\R^2\\, and the predictor
variances and covariances) draw on every row with an observed predictor,
not only the rows that are complete on the outcome. When the missing
values are confined to the outcome, however, the unstandardized slopes
and their standard errors match listwise deletion up to the maximum
likelihood \\N\\ versus \\N - K - 1\\ variance divisor. The rows with an
observed predictor but a missing outcome inform the marginal
distribution of *X*, not the conditional distribution of *Y* given *X*
that identifies the slopes, so they leave the slope estimates and their
conditional-model standard errors unchanged.

The full information advantage is largest when (i) any predictors are
missing on some rows, (ii) auxiliary variables that correlate with the
outcome or with the missingness mechanism are supplied through
`auxiliary` (see below), or (iii) the bootstrap is used to obtain
inference that does not depend on the multivariate normality assumption.

**Auxiliary variables.** A variable that is not part of the regression
but is correlated with the outcome or with the missingness can be
supplied through `auxiliary`. Auxiliaries are added to the model as
saturated correlates (Graham, 2003): each one is correlated with the
residual of the outcome, with every predictor, and with every other
auxiliary, but is never entered as a predictor. The focal regression
coefficients keep their meaning (on complete data they are unchanged to
working precision), while the full information maximum likelihood uses
the auxiliaries' observed values to make the MAR assumption hold
conditional on more of the observed data and to recover information that
listwise deletion discards. This is the inclusive analysis strategy of
Collins, Schafer, and Kam (2001): a variable that predicts the
missingness or the incomplete outcome belongs in the analysis even when
it is of no substantive interest. Auxiliary variables must be numeric
and require `fixed_x = FALSE` (the default).

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
parameter constrained.

**The bootstrap interval.** With `ci_method = "boot"` the rows of `data`
are resampled with replacement `B` times (1000 by default) and the model
is refit on each resample; with `boot_type = "bollen.stine"` the
resamples are instead drawn from data transformed to satisfy the fitted
model (Bollen & Stine, 1992), a model-based bootstrap. Only the
percentile interval is offered: each coefficient's limits are the
empirical quantiles of its resampled estimates (Efron & Tibshirani,
1993); there is no bias-corrected and accelerated (BCa) or bootstrap
standard error variant. Resamples on which the refit does not converge
are dropped, and the interval is computed from the resamples that return
a value. The default `B = 1000` is adequate for the central quantiles a
percentile interval uses; raising it tightens the Monte Carlo error of
the reported limits. Bootstrap results vary from run to run; supply
`boot_seed` for reproducibility.

**Model representation.** Internally the model is fit through lavaan as
a structural equation model in which *Y* is regressed on the predictors
and (when `fixed_x = FALSE`) the predictor distribution is also
estimated. With `fixed_x = FALSE` and complete data, point estimates of
the slopes are identical to `lm` and the residual variance differs only
by the usual \\N\\ versus \\N - K - 1\\ divisor.

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
# Complete data: estimates agree with lm() to working precision.
fit_mlmr <- mlmr(mpg ~ wt + hp, data = mtcars, ci_method = "wald")
fit_lm   <- lm(mpg ~ wt + hp, data = mtcars)
cbind(mlmr = coef(fit_mlmr), lm = coef(fit_lm))
#>                    mlmr          lm
#> (Intercept) 37.22727012 37.22727012
#> wt          -3.87783074 -3.87783074
#> hp          -0.03177295 -0.03177295

# Default print and summary.
fit_mlmr
#> 
#> Call:
#> mlmr(formula = mpg ~ wt + hp, data = mtcars, ci_method = "wald")
#> 
#> Coefficients:
#> (Intercept)           wt           hp  
#>    37.22727     -3.87783     -0.03177  
#> 
summary(fit_mlmr)
#> 
#> Call:
#> mlmr(formula = mpg ~ wt + hp, data = mtcars, ci_method = "wald")
#> 
#> Missing data: fiml | Estimator: ML | SE: standard
#> Sample size (used by lavaan): 32   Complete cases: 32
#> 
#> Coefficients:
#>              Estimate Std. Error z value Pr(>|z|)    
#> (Intercept) 37.227270   1.522000  24.459  < 2e-16 ***
#> wt          -3.877831   0.602344  -6.438 1.21e-10 ***
#> hp          -0.031773   0.008596  -3.696 0.000219 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Confidence intervals (method: wald, level = 0.95):
#>                2.5 %   97.5 %
#> (Intercept) 34.24420 40.21034
#> wt          -5.05840 -2.69726
#> hp          -0.04862 -0.01493
#> 
#> Per-predictor effect sizes (semi-partial R^2 and f^2):
#>       sr^2 Cohen's f^2
#> wt 0.22435      1.2952
#> hp 0.07395      0.4269
#> 
#> Residual std. error (ML): 2.469 on 29 residual degrees of freedom
#> Model implied R-squared: 0.8268,   Adjusted R-squared: 0.8148,   Cohen's f^2: 4.773
#> Omnibus likelihood ratio test (all slopes = 0): chi square = 56.1 on 2 df, p = 6.567e-13
#> Log likelihood: -289.6   AIC: 597.2   BIC: 610.4

# Likelihood ratio CIs (default).
confint(fit_mlmr)
#>                   2.5 %      97.5 %
#> (Intercept) 34.24420415 40.21033609
#> wt          -5.05840396 -2.69725752
#> hp          -0.04862085 -0.01492504

# Auxiliary variable (saturated correlates). qsec is not in the
# model; it is brought in to inform the likelihood. On complete data
# the coefficients are unchanged to working precision, and when the
# outcome is missing as a function of qsec it helps recover them.
# \donttest{
fit_aux <- mlmr(mpg ~ wt + hp, data = mtcars, auxiliary = "qsec")
cbind(no_aux = coef(fit_mlmr), aux = coef(fit_aux))
#>                  no_aux         aux
#> (Intercept) 37.22727012 37.22727014
#> wt          -3.87783074 -3.87783075
#> hp          -0.03177295 -0.03177295

# Demonstrate FIML with missing values on a predictor.
set.seed(113)
d <- mtcars
d$hp[sample.int(nrow(d), 8)] <- NA
fit_fiml <- mlmr(mpg ~ wt + hp, data = d)            # FIML
fit_lwd  <- mlmr(mpg ~ wt + hp, data = d,
                 missing = "listwise")               # listwise
rbind(FIML = coef(fit_fiml), listwise = coef(fit_lwd))
#>          (Intercept)        wt          hp
#> FIML        37.48716 -4.027626 -0.02983439
#> listwise    37.01425 -3.989346 -0.02756755

# Bootstrap confidence intervals; supply boot_seed for reproducibility.
fit_boot <- mlmr(mpg ~ wt + hp, data = mtcars,
                 ci_method = "boot", B = 200, boot_seed = 113)
confint(fit_boot)
#>                   2.5 %     97.5 %
#> (Intercept) 33.27266355 41.2385549
#> wt          -5.06524754 -2.6284807
#> hp          -0.04993342 -0.0194477
# }
```
