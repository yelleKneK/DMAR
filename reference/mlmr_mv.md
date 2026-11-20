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
  ...,
  wald_df = Inf
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
  percentile interval is the only bootstrap interval offered. The
  examples ask for the Wald and bootstrap intervals so the help page
  stays quick; a reported analysis leaves `ci_method` at its default.

- conf_level:

  Desired level of confidence. Defaults to `0.95`.

- B:

  Number of bootstrap resamples when `ci_method = "boot"`, a single
  whole number. Default `1000`, because every resample refits the model
  for every outcome. The fewest accepted is 400 at the default
  `conf_level` (200 at 0.90, and 2000 at 0.99, which the default does
  not meet); a smaller `B` is refused before anything is computed (class
  `dmar_count_invalid`). See *Bootstrap replications* under *Details*
  and
  [`dmar_bootstrap`](https://yelleknek.github.io/DMAR/reference/dmar_bootstrap.md).

- boot_type:

  Character; `"ordinary"` (default) or `"bollen.stine"`.

- boot_seed:

  Integer or `NULL`. Defaults to `NULL`, under which the resamples are
  drawn from, and advance, the current random number stream; supply an
  integer for reproducible bootstraps. When set, the seed applies for
  the duration of the call and the caller's generator state is restored
  on exit.

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
  predictor, while the predictors keep their free covariances, so the
  model is saturated, the per-outcome regression coefficients keep their
  meaning, and the full information maximum likelihood draws on the
  auxiliaries' observed values (the inclusive analysis strategy;
  Collins, Schafer, & Kam, 2001). A name in `auxiliary` must be numeric,
  present in `data`, absent from `formula`, and requires
  `fixed_x = FALSE`. Auxiliaries are refused with
  `missing = "listwise"`, which would drop every row on which an
  auxiliary is missing and so change the sample the coefficients
  describe, and with `missing = "pairwise"` or `"available.cases"`,
  which compute the moments from each pair of variables on its own, so
  that an auxiliary would leave every estimate unchanged.

- effect_sizes:

  Logical; whether to compute per-outcome standardized betas,
  semi-partial \\R^2\\, and Cohen's \\f^2\\. Defaults to `TRUE`.

- enforce_es_bounds:

  Logical; if `TRUE`, negative per-outcome \\sr^2\\ or \\f^2\\ estimates
  (finite-sample artifacts) are clamped to zero. Defaults to `FALSE`.

- ...:

  Additional arguments forwarded to
  [`lavaan`](https://rdrr.io/pkg/lavaan/man/lavaan.html).

- wald_df:

  The reference distribution of the Wald interval and of the *p*-values:
  `Inf` (the default), the normal distribution; a single positive
  number, the degrees of freedom of a *t* distribution for every
  coefficient; or `"missing_information"`, von Hippel's (2016) degrees
  of freedom for each coefficient from its fraction of missing
  information (see
  [`mlmr`](https://yelleknek.github.io/DMAR/reference/mlmr.md)). A
  finite value or the keyword requires `ci_method = "wald"`; the keyword
  is an error where, with values missing, lavaan fit the complete cases
  in place of the FIML fit asked for, as in
  [`mlmr`](https://yelleknek.github.io/DMAR/reference/mlmr.md). It
  follows `...`, so it must be named in full.

## Value

An object of class `"mlmr_mv"`, a list with components similar to a
univariate [`mlmr`](https://yelleknek.github.io/DMAR/reference/mlmr.md)
fit but extended for multiple outcomes:

- `call`, `formula`, `terms`, `model`, `xlevels`:

  As in `mlmr`.

- `coefficients`:

  A matrix with predictors (and an intercept row, when present) as rows
  and outcomes as columns, matching `coef.mlm`.

- `B`, `B_used`, `n_dropped`, `B_required`, `bootstrap_status`:

  Present only with `ci_method = "boot"`: the resamples requested, the
  number the intervals rest on, their difference, the fewest resamples a
  request may ask for, and the verdict (`"complete"`, `"dropped"`,
  `"exploratory"`, `"insufficient_count"`, or `"insufficient_share"`;
  under the last two the limits are `NA`). See *Bootstrap replications*
  under Details.

- `coef_table`:

  A long `data.frame` with one row per (outcome, term) combination;
  columns include `outcome`, `term`, `estimate`, `se`, `test_statistic`,
  `p_value`, `ci_lower`, `ci_upper`, `std_estimate`, carrying the
  confidence level in its `conf_level` attribute, so that
  [`results_sentence`](https://yelleknek.github.io/DMAR/reference/dmar_output_helpers.md)
  on a row reports it.

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

The object carries the attributes `wald_df` and, under
`wald_df = "missing_information"`, `wald_df_table`, as an
[`mlmr`](https://yelleknek.github.io/DMAR/reference/mlmr.md) fit does,
with the coefficients named `"outcome:term"` in the order of
`coef_table`.

## Details

**Why a separate function from `mlmr`.** A univariate `mlmr` fit handles
the case of one outcome regressed on one or more predictors. `mlmr_mv`
extends to the case of two or more outcomes regressed on the same
predictor set, modeling the residual covariance among outcomes
explicitly. Under FIML a row that is missing on one outcome still
contributes information about the other outcomes (through the modeled
residual covariance) and about the joint distribution of the predictors,
where listwise deletion drops the row. A user who fits separate
univariate regressions for each outcome under listwise deletion can
discard a great deal of information when the outcomes are correlated and
missingness patterns differ.

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
never as a predictor, and the predictors keep their free means,
variances, and covariances. The model over the outcomes, the predictors,
and the auxiliaries is then saturated, so the per-outcome coefficients
keep their meaning (on complete data they equal those of
[`lm()`](https://rdrr.io/r/stats/lm.html)) while the likelihood draws on
the auxiliaries' observed values (the inclusive analysis strategy of
Collins, Schafer, & Kam, 2001).

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

**A *t* reference distribution for the Wald interval.** With
`ci_method = "wald"`, `wald_df` refers every coefficient's interval and
*p*-value to a *t* distribution rather than the normal one, keeping the
estimate and its standard error (von Hippel, 2016, Eq. 29), exactly as
in [`mlmr`](https://yelleknek.github.io/DMAR/reference/mlmr.md), whose
help page gives the reasoning. A number applies to every coefficient of
every outcome. Under `"missing_information"` each coefficient gets von
Hippel's \\\max\\3, N(1 - \gamma) - k\\\\ from its own fraction of
missing information \\\gamma\\, with \\k\\ the number of regression
parameters of its outcome (the slopes and the intercept), so the
coefficients of a fully observed outcome keep nearly \\N - k\\ degrees
of freedom while those of an outcome observed for few cases get fewer.

**Caveats.** Same as
[`mlmr`](https://yelleknek.github.io/DMAR/reference/mlmr.md): the
function assumes that, conditional on the predictors, the joint
distribution of the outcomes is multivariate normal with constant
covariance, and that missingness is at most MAR. Factor predictors and
interactions are expanded through
[`model.matrix`](https://rdrr.io/r/stats/model.matrix.html) once and
reused for every outcome.

## References

Collins, L. M., Schafer, J. L., & Kam, C.-M. (2001). A comparison of
inclusive and restrictive strategies in modern missing data procedures.
*Psychological Methods, 6*(4), 330–351.
[doi:10.1037/1082-989X.6.4.330](https://doi.org/10.1037/1082-989X.6.4.330)

Graham, J. W. (2003). Adding missing-data-relevant variables to
FIML-based structural equation models. *Structural Equation Modeling,
10*(1), 80–100.
[doi:10.1207/S15328007SEM1001_4](https://doi.org/10.1207/S15328007SEM1001_4)

von Hippel, P. T. (2016). New confidence intervals and bias comparisons
show that maximum likelihood can beat multiple imputation in small
samples. *Structural Equation Modeling, 23*(3), 422–437.
[doi:10.1080/10705511.2015.1047931](https://doi.org/10.1080/10705511.2015.1047931)

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
# what separates this from two separate regressions. The fit asks
# for the Wald interval, a choice taken up below, and keeps the
# default effect_sizes = TRUE: the per-outcome effect sizes come
# back on the fit rather than in summary(), one row per outcome and
# predictor, giving the semi-partial R^2 and Cohen's f^2, and they
# are what fills the standardized coefficients in coef_table.
fit <- mlmr_mv(cbind(t6_paragraph_comprehension, t9_word_meaning) ~
                 t5_general_information + t7_sentence,
               data = holzinger_swineford,
               ci_method = "wald")
coef(fit)              # matrix: rows = predictors, cols = outcomes
#>                        t6_paragraph_comprehension t9_word_meaning
#> (Intercept)                           -0.25416032      -6.3377552
#> t5_general_information                 0.07653333       0.2844663
#> t7_sentence                            0.36460353       0.5811433
summary(fit)
#> 
#> Call:
#> mlmr_mv(formula = cbind(t6_paragraph_comprehension, t9_word_meaning) ~ 
#>     t5_general_information + t7_sentence, data = holzinger_swineford, 
#>     ci_method = "wald")
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
print(fit$effect_sizes, row.names = FALSE)
#>                     outcome                   term        sr2        f2
#>  t6_paragraph_comprehension t5_general_information 0.03587460 0.0840968
#>  t6_paragraph_comprehension            t7_sentence 0.14153107 0.3317754
#>             t9_word_meaning t5_general_information 0.10277307 0.2712471
#>             t9_word_meaning            t7_sentence 0.07456057 0.1967864

# The interval menu is profile, Wald, and bootstrap. The default,
# ci_method = "profile", inverts the likelihood ratio test one
# coefficient at a time, and with two outcomes there are twice as
# many coefficients to profile; it is what a reported interval
# deserves, and leaving ci_method at its default asks for it. The
# fits on this page ask for the Wald or the bootstrap interval
# because those refits take longer than a help page should. The
# bootstrap resamples rows and takes percentile limits; it is what
# to ask for when the multivariate normality the likelihood assumes
# is doubtful. Every resample refits the two-outcome model. B = 10
# keeps the example quick and is below the 400 a reported interval
# needs; the option lets it run and the warning says so. A reported
# interval deserves the default B = 1000. boot_seed fixes the
# resamples, so the limits are reproducible rather than moving from
# run to run. The effect sizes cost one constrained refit per outcome
# and predictor, so this fit and the ones after it leave them off.
old <- options(dmar.bootstrap_exploratory = TRUE)
fit_boot <- mlmr_mv(cbind(t6_paragraph_comprehension, t9_word_meaning) ~
                      t5_general_information + t7_sentence,
                    data = holzinger_swineford,
                    ci_method = "boot", B = 10, boot_seed = 113,
                    effect_sizes = FALSE)
#> Warning: mlmr_mv(): B = 10 is below the 400 replications a 95% interval is reported from (10 beyond each limit). The limits are shown because options(dmar.bootstrap_exploratory = TRUE) is set. With 10 replications each limit rests on about 0.3 of them and will move from seed to seed; use the result to explore or to demonstrate, and report an interval only from B of at least 400. B, B_used, and B_required are recorded in the result.
options(old)
confint(fit_boot)
#>                      outcome                   term       2.5 %      97.5 %
#> 1 t6_paragraph_comprehension            (Intercept) -1.00690701  0.27654636
#> 2 t6_paragraph_comprehension t5_general_information  0.07182434  0.09680779
#> 3 t6_paragraph_comprehension            t7_sentence  0.29716900  0.40923300
#> 4            t9_word_meaning            (Intercept) -7.78549657 -4.94677137
#> 5            t9_word_meaning t5_general_information  0.26989256  0.33955085
#> 6            t9_word_meaning            t7_sentence  0.45180978  0.67755562

# FIML versus listwise when one outcome has missing values. The
# revised second-form test t26_flags was administered to only 145
# of the 301 students, so it carries real missingness. A row with
# t26_flags missing still informs the likelihood about the other
# outcome, about the predictors, and, through the residual
# covariance, about t26_flags itself, so no row is discarded. Notice
# the two sample sizes, and that the coefficients of the complete
# outcome differ between the fits: listwise deletion drops 156 of
# its observed rows along with the missing t26_flags values.
fit_fiml <- mlmr_mv(cbind(t6_paragraph_comprehension,
                          t26_flags) ~
                      t7_sentence + t9_word_meaning,
                    data = holzinger_swineford,
                    ci_method = "wald", effect_sizes = FALSE)
fit_lwd  <- mlmr_mv(cbind(t6_paragraph_comprehension,
                          t26_flags) ~
                      t7_sentence + t9_word_meaning,
                    data = holzinger_swineford,
                    missing = "listwise", ci_method = "wald",
                    effect_sizes = FALSE)
c(N_fiml = nobs(fit_fiml), N_listwise = nobs(fit_lwd))
#>     N_fiml N_listwise 
#>        301        145 
cbind(FIML = coef(fit_fiml)[, "t6_paragraph_comprehension"],
      listwise = coef(fit_lwd)[, "t6_paragraph_comprehension"])
#>                      FIML  listwise
#> (Intercept)     1.1169741 0.9709397
#> t7_sentence     0.3174143 0.3134578
#> t9_word_meaning 0.1669888 0.1777868

# Auxiliary variable (saturated correlates): the complete speed test
# t13_straight_and_curved_capitals informs the likelihood without
# entering either regression. Continuing from the model above, the
# coefficients of the complete outcome are unchanged to working
# precision, while those of t26_flags move, since the auxiliary
# carries information about the rows where t26_flags is missing.
fit_aux <- mlmr_mv(cbind(t6_paragraph_comprehension,
                         t26_flags) ~
                     t7_sentence + t9_word_meaning,
                   data = holzinger_swineford,
                   ci_method = "wald",
                   auxiliary = "t13_straight_and_curved_capitals",
                   effect_sizes = FALSE)
coef(fit_aux)
#>                 t6_paragraph_comprehension  t26_flags
#> (Intercept)                      1.1169741 27.0919099
#> t7_sentence                      0.3174143  0.3408605
#> t9_word_meaning                  0.1669888  0.1749609
```
