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

  The variance-covariance matrix of the regression coefficients across
  all outcomes, returned by
  [`vcov()`](https://rdrr.io/r/stats/vcov.html). Rows and columns follow
  the outcome-major order of `coef_table` (the column-major flattening
  of `coefficients`) and are named `"outcome:term"`, for example
  `"mpg:wt"`, the naming [`vcov`](https://rdrr.io/r/stats/vcov.html)
  uses for an `"mlm"` fit. The cross-outcome blocks carry the sampling
  covariance between coefficients of different outcomes, so joint Wald
  tests across outcomes compose with
  [`coef()`](https://rdrr.io/r/stats/coef.html).

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

**The bootstrap interval.** `ci_method = "boot"` resamples the rows of
`data` with replacement `B` times (1000 by default), refits the model on
each resample, and reports each coefficient's percentile interval. It is
the interval to ask for when the multivariate normality the likelihood
assumes is doubtful, since its coverage does not rest on that
assumption. The price is `B` refits of a model that already carries *J*
outcomes, so a bootstrap interval is a deliberate request rather than a
default. Bootstrap results vary from run to run; supply `boot_seed` for
reproducibility. The mechanics, including the Bollen-Stine variant, are
given in the `ci_method` argument description and in
[`mlmr`](https://yelleknek.github.io/DMAR/reference/mlmr.md).

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
# Two outcomes on a shared predictor set. The residual covariance
# between the outcomes is estimated as part of the model, which is
# what separates this from two separate regressions. This fit asks
# for the Wald interval and leaves the effect sizes off so that it
# runs at example time. It is the only block here that runs; the
# rest is left as commented code so a reader can see the syntax
# without paying the run time.
fit <- mlmr_mv(cbind(t6_paragraph_comprehension, t9_word_meaning) ~
                 t5_general_information + t7_sentence,
               data = holzinger_swineford,
               ci_method = "wald", effect_sizes = FALSE)
coef(fit)              # matrix: rows = predictors, cols = outcomes
#>                        t6_paragraph_comprehension t9_word_meaning
#> (Intercept)                           -0.25416035      -6.3377552
#> t5_general_information                 0.07653333       0.2844663
#> t7_sentence                            0.36460353       0.5811433
summary(fit)
#> 
#> Call:
#> mlmr_mv(formula = cbind(t6_paragraph_comprehension, t9_word_meaning) ~ 
#>     t5_general_information + t7_sentence, data = holzinger_swineford, 
#>     ci_method = "wald", effect_sizes = FALSE)
#> 
#> Missing: fiml | Estimator: ML | SE: standard
#> Sample size used: 301   Complete cases: 301
#> 
#> --- Outcome: t6_paragraph_comprehension ---
#>                        Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)            -0.25416    0.48954  -0.519    0.604    
#> t5_general_information  0.07653    0.01521   5.031 4.87e-07 ***
#> t7_sentence             0.36460    0.03649   9.993  < 2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> R^2:        0.5734
#> Adj. R^2:   0.5706
#> 
#> --- Outcome: t9_word_meaning ---
#>                        Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)            -6.33776    1.01314  -6.256 3.96e-10 ***
#> t5_general_information  0.28447    0.03148   9.036  < 2e-16 ***
#> t7_sentence             0.58114    0.07551   7.696 1.40e-14 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> R^2:        0.6211
#> Adj. R^2:   0.6186
#> 
#> Residual covariance among outcomes:
#>                            t6_paragraph_comprehension t9_word_meaning
#> t6_paragraph_comprehension                     5.1856          3.0942
#> t9_word_meaning                                3.0942         22.2112
#> 
#> Log likelihood: -3552   AIC: 7132   BIC: 7184
fit$R2                 # per-outcome R^2
#> t6_paragraph_comprehension            t9_word_meaning 
#>                  0.5734131                  0.6211091 
fit$residual_cov       # residual covariance among outcomes
#>                            t6_paragraph_comprehension t9_word_meaning
#> t6_paragraph_comprehension                   5.185583        3.094177
#> t9_word_meaning                              3.094177       22.211178

# The interval menu is profile, Wald, and bootstrap. The default,
# ci_method = "profile", inverts the likelihood ratio test one
# coefficient at a time, and with two outcomes there are twice as
# many coefficients to profile; it is what a reported interval
# deserves. The bootstrap resamples rows and takes percentile
# limits; it is what to ask for when the multivariate normality the
# likelihood assumes is doubtful. Each refits the model many times,
# so neither is run here; the calls are
#   mlmr_mv(cbind(t6_paragraph_comprehension, t9_word_meaning) ~
#             t5_general_information + t7_sentence,
#           data = holzinger_swineford)
#   mlmr_mv(cbind(t6_paragraph_comprehension, t9_word_meaning) ~
#             t5_general_information + t7_sentence,
#           data = holzinger_swineford,
#           ci_method = "boot", B = 1000, boot_seed = 113)
# with boot_seed supplied because bootstrap limits otherwise move
# from run to run.

# The per-outcome effect sizes come back on the fit rather than in
# summary(): one row per outcome and predictor, giving the
# semi-partial R^2 and Cohen's f^2. They cost one constrained refit
# per outcome and predictor, so they are left off above, which also
# leaves the standardized coefficients in coef_table missing. With
# the default effect_sizes = TRUE the table is
#   fit_es <- mlmr_mv(cbind(t6_paragraph_comprehension,
#                           t9_word_meaning) ~
#                       t5_general_information + t7_sentence,
#                     data = holzinger_swineford,
#                     ci_method = "wald")
#   print(fit_es$effect_sizes, row.names = FALSE)

# FIML versus listwise when one outcome has missing values. The
# revised second-form test t26_flags was administered to only 145
# of the 301 students, so it carries real missingness. A row with
# t26_flags missing still informs the likelihood about the other
# outcome, about the predictors, and, through the residual
# covariance, about t26_flags itself, so no row is discarded. Not
# run here because the comparison costs two more fits; the code is:
#   fit_fiml <- mlmr_mv(cbind(t6_paragraph_comprehension,
#                             t26_flags) ~
#                         t7_sentence + t9_word_meaning,
#                       data = holzinger_swineford,
#                       ci_method = "wald", effect_sizes = FALSE)
#   fit_lwd  <- mlmr_mv(cbind(t6_paragraph_comprehension,
#                             t26_flags) ~
#                         t7_sentence + t9_word_meaning,
#                       data = holzinger_swineford,
#                       missing = "listwise", ci_method = "wald",
#                       effect_sizes = FALSE)
#   c(N_fiml = nobs(fit_fiml), N_listwise = nobs(fit_lwd))
#   cbind(FIML = coef(fit_fiml)[, "t6_paragraph_comprehension"],
#         listwise = coef(fit_lwd)[, "t6_paragraph_comprehension"])

# Auxiliary variable (saturated correlates): the complete speed test
# t13_straight_and_curved_capitals informs the likelihood without
# entering either regression. Continuing from the model above, and
# again not run here:
#   fit_aux <- mlmr_mv(cbind(t6_paragraph_comprehension,
#                            t26_flags) ~
#                        t7_sentence + t9_word_meaning,
#                      data = holzinger_swineford,
#                      ci_method = "wald",
#                      auxiliary = "t13_straight_and_curved_capitals",
#                      effect_sizes = FALSE)
#   coef(fit_aux)
```
