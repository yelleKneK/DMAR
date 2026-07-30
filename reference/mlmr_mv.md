# Multivariate Maximum Likelihood Regression With Full Information Missing Data Handling

Fits a multivariate multiple regression model by maximum likelihood,
with full information maximum likelihood handling of missing values by
default. Multiple outcomes are regressed on a shared predictor set
simultaneously, with the residual covariance among outcomes estimated as
part of the model. The formula interface mirrors
[`lm`](https://rdrr.io/r/stats/lm.html)'s multivariate syntax,
`cbind(y1, y2, y3) ~ x1 + x2`, and the returned object supports the same
family of S3 methods as a univariate
[`mlmr`](https://yelleknek.github.io/DMAR/reference/mlmr.md) fit.

## Usage

``` r
mlmr_mv(
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

  A two-sided [`formula`](https://rdrr.io/r/stats/formula.html) whose
  left-hand side is `cbind(y1, y2, ...)` for two or more numeric
  outcomes and whose right-hand side names the shared predictor set.
  Factor predictors, interactions, polynomial terms, and transformations
  are supported as in `mlmr`.

- data:

  A `data.frame` containing the variables in `formula`.

- missing:

  Character; passed to lavaan. Defaults to `"fiml"` (equivalently
  `"ml"`). See
  [`mlmr`](https://yelleknek.github.io/DMAR/reference/mlmr.md) for the
  full list of accepted values.

- ci_method:

  Character; confidence interval method for the regression coefficients.
  `"profile"` (default), `"wald"`, or `"boot"`. The same trade-offs
  apply as in
  [`mlmr`](https://yelleknek.github.io/DMAR/reference/mlmr.md); with *J*
  outcomes and *K* predictors, profile likelihood requires \\O(JK)\\
  constrained refits. `"boot"` resamples the rows of `data` `B` times
  (or draws Bollen-Stine model-based resamples), refits on each, drops
  resamples whose refit does not converge, and reports each
  coefficient's percentile interval, the empirical quantiles of its
  resampled estimates; as in
  [`mlmr`](https://yelleknek.github.io/DMAR/reference/mlmr.md), the
  percentile interval is the only bootstrap interval offered.

- conf_level:

  Desired level of confidence. Defaults to `0.95`.

- B:

  Integer; number of bootstrap resamples when `ci_method = "boot"`.
  Defaults to `1000`.

- boot_type:

  Character; `"ordinary"` (default) or `"bollen.stine"`.

- boot_seed:

  Integer or `NULL`. Defaults to `NULL`, which leaves the user's current
  RNG state intact; supply an integer for reproducible bootstraps. When
  set, the function saves and restores `.Random.seed` on exit.

- estimator:

  Character; one of `"ML"` (default), `"MLR"`, `"MLM"`, `"GLS"`.

- se:

  Character or `NULL`; defaults to `NULL`, meaning "choose automatically
  from `estimator`" (same logic as
  [`mlmr`](https://yelleknek.github.io/DMAR/reference/mlmr.md)).

- fixed_x:

  Logical; defaults to `FALSE` (jointly model the predictor
  distribution, required for FIML to use rows with missing predictors).

- auxiliary:

  Character vector of variable names in `data` to include as auxiliary
  variables (saturated correlates; Graham, 2003), or `NULL` (default)
  for none. Each auxiliary is correlated with every outcome's residual,
  every predictor, and each other auxiliary, but is never entered as a
  predictor, so the per-outcome regression coefficients keep their
  meaning while the full information maximum likelihood draws on the
  auxiliaries' observed values (the inclusive analysis strategy;
  Collins, Schafer, & Kam, 2001). A name in `auxiliary` must be numeric,
  present in `data`, absent from `formula`, and requires
  `fixed_x = FALSE`.

- effect_sizes:

  Logical; whether to compute per-outcome standardized betas,
  semi-partial \\R^2\\, and Cohen's \\f^2\\. Defaults to `TRUE`.

- enforce_es_bounds:

  Logical; if `TRUE`, negative per-outcome \\sr^2\\ or \\f^2\\ estimates
  (finite-sample artifacts) are clamped to zero. Defaults to `FALSE`.

- ...:

  Additional arguments forwarded to
  [`lavaan`](https://rdrr.io/pkg/lavaan/man/lavaan.html).

## Value

An object of class `"mlmr_mv"`, a list with components similar to a
univariate [`mlmr`](https://yelleknek.github.io/DMAR/reference/mlmr.md)
fit but extended for multiple outcomes:

- `call`, `formula`, `terms`, `model`, `xlevels`:

  As in `mlmr`.

- `coefficients`:

  A matrix with predictors (and an intercept row, when present) as rows
  and outcomes as columns, matching `coef.mlm`.

- `coef_table`:

  A long `data.frame` with one row per (outcome, term) combination;
  columns include `outcome`, `term`, `estimate`, `se`, `z_value`,
  `p_value`, `ci_lower`, `ci_upper`, `std_estimate`.

- `vcov`:

  The full joint variance-covariance matrix of all regression
  coefficients across all outcomes.

- `residual_cov`:

  The estimated residual covariance matrix among outcomes (*J* by *J*).

- `R2`:

  Named vector of model implied \\R^2\\ per outcome.

- `adj_R2`:

  Named vector of adjusted \\R^2\\ per outcome.

- `effect_sizes`:

  When `effect_sizes = TRUE`, a long `data.frame` with one row per
  (outcome, predictor) combination giving \\sr^2\\ and Cohen's \\f^2\\.

- `fitted.values`, `residuals`:

  Matrices with rows = observations and columns = outcomes; `NA` in rows
  where any predictor is missing.

- `logLik`, `N`, `N_complete`:

  As in `mlmr`.

- `lavaan_fit`:

  The underlying lavaan fit.

## Details

**Why a separate function from `mlmr`.** A univariate `mlmr` fit handles
the case of one outcome regressed on one or more predictors. `mlmr_mv`
extends to the case of two or more outcomes regressed on the same
predictor set, modeling the residual covariance among outcomes
explicitly. This is the regression problem in which the FIML advantage
over listwise deletion is largest, because rows that are missing on one
outcome still contribute information about the other outcomes (via the
modeled residual covariance) and about the joint distribution of the
predictors. A user who fits separate univariate regressions for each
outcome under listwise deletion can discard a great deal of information
when the outcomes are correlated and missingness patterns differ.

**Same predictor set across outcomes.** The formula
`cbind(y1, y2) ~ x1 + x2` regresses both `y1` and `y2` on `x1` and `x2`.
Per-outcome predictor sets (sometimes called seemingly unrelated
regression with heterogeneous predictors) are not supported; users who
need that can fit separate
[`mlmr`](https://yelleknek.github.io/DMAR/reference/mlmr.md) models or
call [`sem`](https://rdrr.io/pkg/lavaan/man/sem.html) directly with a
custom model string.

**Auxiliary variables.** As in
[`mlmr`](https://yelleknek.github.io/DMAR/reference/mlmr.md), variables
that are not part of the regression but are correlated with an outcome
or with the missingness can be supplied through `auxiliary` and are
entered as saturated correlates (Graham, 2003): each is correlated with
every outcome's residual, every predictor, and each other auxiliary, but
never as a predictor, so the per-outcome coefficients keep their meaning
while the likelihood draws on the auxiliaries' observed values (the
inclusive analysis strategy of Collins, Schafer, & Kam, 2001).

**Caveats.** Same as
[`mlmr`](https://yelleknek.github.io/DMAR/reference/mlmr.md): the
function assumes that, conditional on the predictors, the joint
distribution of the outcomes is multivariate normal with constant
covariance, and that missingness is at most MAR. Factor predictors and
interactions are expanded through
[`model.matrix`](https://rdrr.io/r/stats/model.matrix.html) once and
reused for every outcome.

## See also

[`mlmr`](https://yelleknek.github.io/DMAR/reference/mlmr.md) for the
univariate sibling; [`lm`](https://rdrr.io/r/stats/lm.html) (and the
`"mlm"` object class) for the OLS multivariate analog;
[`sem`](https://rdrr.io/pkg/lavaan/man/sem.html) for the underlying
engine.

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Two outcomes, shared predictor set.
fit <- mlmr_mv(cbind(mpg, disp) ~ wt + hp, data = mtcars,
               ci_method = "wald")
coef(fit)              # matrix: rows = predictors, cols = outcomes
#>                     mpg         disp
#> (Intercept) 37.22727016 -129.9505489
#> wt          -3.87783074   82.1125005
#> hp          -0.03177295    0.6578337
summary(fit)
#> 
#> Call:
#> mlmr_mv(formula = cbind(mpg, disp) ~ wt + hp, data = mtcars, 
#>     ci_method = "wald")
#> 
#> Missing: fiml | Estimator: ML | SE: standard
#> Sample size used: 32   Complete cases: 32
#> 
#> --- Outcome: mpg ---
#>              Estimate Std. Error z value Pr(>|z|)    
#> (Intercept) 37.227270   1.522000  24.459  < 2e-16 ***
#> wt          -3.877831   0.602344  -6.438 1.21e-10 ***
#> hp          -0.031773   0.008596  -3.696 0.000219 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> R^2:        0.8268
#> Adj. R^2:   0.8148
#> 
#> --- Outcome: disp ---
#>              Estimate Std. Error z value Pr(>|z|)    
#> (Intercept) -129.9505    27.7871  -4.677 2.92e-06 ***
#> wt            82.1125    10.9970   7.467 8.21e-14 ***
#> hp             0.6578     0.1569   4.192 2.77e-05 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> R^2:        0.8635
#> Adj. R^2:   0.8541
#> 
#> Residual covariance among outcomes:
#>          mpg      disp
#> mpg   6.0952   -1.9037
#> disp -1.9037 2031.6390
#> 
#> Log likelihood: -456.9   AIC: 941.8   BIC: 962.3
fit$R2                 # per-outcome R^2
#>       mpg      disp 
#> 0.8267854 0.8634722 
fit$residual_cov       # residual covariance among outcomes
#>            mpg        disp
#> mpg   6.095242   -1.903662
#> disp -1.903662 2031.639016

# \donttest{
# FIML versus listwise when one outcome has missing values.
set.seed(113)
d <- mtcars
d$disp[sample.int(nrow(d), 6)] <- NA
fit_fiml <- mlmr_mv(cbind(mpg, disp) ~ wt + hp, data = d,
                    ci_method = "wald")
fit_lwd  <- mlmr_mv(cbind(mpg, disp) ~ wt + hp, data = d,
                    missing = "listwise", ci_method = "wald")
rbind(N_fiml = nobs(fit_fiml), N_listwise = nobs(fit_lwd))
#>            [,1]
#> N_fiml       32
#> N_listwise   26

# Auxiliary variable (saturated correlates): qsec informs the
# likelihood without entering either regression.
fit_aux <- mlmr_mv(cbind(mpg, disp) ~ wt + hp, data = d,
                   ci_method = "wald", auxiliary = "qsec")
coef(fit_aux)
#>                     mpg         disp
#> (Intercept) 37.22726844 -128.9695900
#> wt          -3.87783038   86.2338818
#> hp          -0.03177295    0.5549623
# }
```
