# Reliability Coefficient With a Confidence Interval (General Dispatch)

General-purpose entry point for the reliability family. Dispatches to
[`reliability_alpha`](https://yelleknek.github.io/DMAR/reference/reliability_alpha.md),
[`reliability_kr20`](https://yelleknek.github.io/DMAR/reference/reliability_kr20.md),
[`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md),
[`reliability_omega_categorical`](https://yelleknek.github.io/DMAR/reference/reliability_omega_categorical.md),
or
[`reliability_H`](https://yelleknek.github.io/DMAR/reference/reliability_H.md)
according to the requested `type`. When `type` is not specified, the
function picks a reasonable default from the supplied input following
the recommendations of Kelley and Pornprasertmanit (2016).

## Usage

``` r
reliability(
  data = NULL,
  S = NULL,
  N = NULL,
  type = NULL,
  estimator = c("analytic", "model_implied"),
  denominator = c("observed", "model_implied"),
  missing = NULL,
  aux = NULL,
  ci_method = NULL,
  conf_level = 0.95,
  B = NULL,
  seed = NULL,
  ncpus = 1L,
  improper = c("keep", "drop"),
  impute = NULL,
  imputations = NULL,
  wald_df = Inf,
  progress = NULL
)
```

## Arguments

- data:

  A numeric matrix or data frame of item scores, or `NULL`.

- S:

  A symmetric covariance matrix among the items, or `NULL`. If supplied,
  `N` must also be supplied; methods that require raw data are then
  unavailable.

- N:

  Total sample size; required when `S` is supplied.

- type:

  Character; one of `"alpha"`, `"kr20"`, `"omega"`,
  `"omega_categorical"` (`"omega_c"` is accepted as a shorthand), `"H"`
  (coefficient *H*, the reliability of the optimally weighted composite,
  from raw data or a covariance matrix;
  [`reliability_H`](https://yelleknek.github.io/DMAR/reference/reliability_H.md)
  also takes a fitted model or standardized loadings), or `NULL` for
  auto-detection. See *Details*.

- estimator:

  For `type = "alpha"` only: how coefficient \\\alpha\\ is estimated,
  `"analytic"` (default; the closed-form equation applied to the
  observed covariance matrix) or `"model_implied"` (the reliability
  implied by the \\\tau\\-equivalent single-factor model fit by maximum
  likelihood); forwarded to
  [`reliability_alpha`](https://yelleknek.github.io/DMAR/reference/reliability_alpha.md),
  whose help page discusses the choice and which interval methods each
  estimator supports. Supplying it with any other `type` is an error, as
  supplying `denominator` is for any type but `"omega"` and
  `"omega_categorical"`.

- denominator:

  For `type = "omega"` and `type = "omega_categorical"` only: how the
  total variance in the denominator is estimated, `"observed"` (default;
  robust omega, and for categorical omega the unrestricted polychoric
  correlations) or `"model_implied"` (for categorical omega, the
  coefficient as Green and Yang, 2009, defined it); forwarded to
  [`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md)
  or
  [`reliability_omega_categorical`](https://yelleknek.github.io/DMAR/reference/reliability_omega_categorical.md),
  whose help pages discuss the choice.

- missing:

  How incomplete rows of `data` are handled, forwarded to the family
  function when supplied, which checks it against its own choices and
  whose help page discusses them: for `type = "alpha"` and
  `type = "omega"`, `"listwise"`, `"fiml"` (full information maximum
  likelihood), or `"two_stage"` (two-stage maximum likelihood); for
  `type = "H"`, `"listwise"` or `"fiml"`; for
  `type = "omega_categorical"`, `"listwise"`, `"pairwise"`, or `"mi"`
  (multiple imputation). The default, `NULL`, leaves the family
  function's default in place: for alpha, omega, and H, full information
  maximum likelihood when an item value is missing and the complete data
  analysis when none is, since the complete-data intervals of alpha and
  omega need that analysis
  ([`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md)
  fits complete data by FIML too; a row with every item missing does not
  count, since both treatments drop it), and `"listwise"` for
  categorical omega. Supplying it with `type = "kr20"` is an error.

- aux:

  Optional character vector naming auxiliary variable columns of `data`,
  forwarded to the family function: saturated correlates under full
  information maximum likelihood, and variables entered into the first
  stage under two-stage estimation, for `type = "alpha"`,
  `type = "omega"`, and `type = "H"` (supplying `aux` without `missing`
  implies `missing = "fiml"`); predictors in the imputation model for
  `type = "omega_categorical"` (it implies `missing = "mi"`). Supplying
  it with `type = "kr20"` is an error.

- ci_method:

  Method for constructing the confidence interval, or `NULL` to use the
  chosen family function's default. See the help page for the chosen
  `reliability_*` function for the full list of accepted values.

- conf_level:

  Confidence level. Defaults to `0.95`.

- B:

  Number of bootstrap replications when a bootstrap method is selected.
  The default, `NULL`, uses the family function's default: 10000, or
  1000 for categorical omega under `missing = "mi"`. The family function
  checks `B` against the number its interval method and level require
  (400 for `"percentile"`, 2000 for `"bca"`, 100 for the two bootstrap
  standard error methods, at the default `conf_level`) and refuses a
  smaller one (class `dmar_count_invalid`); see
  [`dmar_bootstrap`](https://yelleknek.github.io/DMAR/reference/dmar_bootstrap.md).

- seed:

  Random number seed used for bootstrap reproducibility. Defaults to
  `NULL`, under which the resamples are drawn from, and advance, the
  current random number stream; supply an integer for a reproducible
  interval, which is set for the call only, the caller's random number
  state being restored when it returns.

- ncpus:

  The number of processes over which the bootstrap replications are
  spread; forwarded to every family function. The default `1` runs them
  in this R session, and a given `seed` returns the same interval
  whatever `ncpus` is (see
  [`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md)).

- improper:

  For `type = "alpha"` and `type = "omega"` only: how a bootstrap
  replication whose single-factor solution is improper but still yields
  a coefficient is treated, `"keep"` (the default) or `"drop"`;
  forwarded to
  [`reliability_alpha`](https://yelleknek.github.io/DMAR/reference/reliability_alpha.md)
  or
  [`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md),
  whose help pages describe it. Categorical omega is not defined for an
  improper solution, and KR-20 fits no model, so supplying it with
  another `type` is an error.

- impute, imputations:

  For `type = "omega_categorical"` with `missing = "mi"` only: the
  imputation model and the number of imputations of each bootstrap
  replication, forwarded to
  [`reliability_omega_categorical`](https://yelleknek.github.io/DMAR/reference/reliability_omega_categorical.md),
  whose help page describes them. The default, `NULL`, leaves that
  function's defaults (`"pmm"` and 5) in place; supplying either with
  another `type` is an error.

- wald_df:

  For `type = "alpha"` and `type = "omega"` only: the reference
  distribution of a Wald interval, `Inf` (the default, the normal
  reference), a single positive number of degrees of freedom of a *t*
  reference, or, for `type = "omega"` with
  `denominator = "model_implied"` only, `"missing_information"`, the
  degrees of freedom of a *t* reference from omega's fraction of missing
  information; forwarded to
  [`reliability_alpha`](https://yelleknek.github.io/DMAR/reference/reliability_alpha.md)
  or
  [`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md),
  whose help pages describe it (the keyword is an error for alpha, for
  which no rule passed the calibration study's decision rule). A finite
  value or the keyword with another `type` is an error.

- progress:

  `NULL` (the default), or a function of two arguments that the
  bootstrap calls as it runs, so that an application can show how far a
  long bootstrap has gone and let its user stop it (a Shiny app, for
  example). It is called as `progress(done, total)`, with the number of
  replications attempted so far and `B`: after each replication when
  `ncpus = 1`, and with more processes after each block of replications
  they complete together, about a fiftieth of `B`
  (`ncpus * ceiling(B / (50 * ncpus))` replications, the last block
  possibly fewer), so that `done` advances by the block. A BCa interval
  takes its acceleration from a regression on the resampling frequencies
  when at least as many replications return a value as there are cases,
  and from the jackknife, one leave-one-out refit per case, when fewer
  do or when that regression fails or gives influence values that are
  not all finite; the jackknife runs at most once in a call, and not at
  all when no interval is formed or every returned replicate has the
  same value. While it runs, the callback is called again after each
  refit (with more processes, after each block of
  `ncpus * ceiling(n / (50 * ncpus))` refits, `n` the number of cases
  resampled), with `done` and `total` both `B`, so it can be called more
  than `B` times, and `done` equal to `total` does not mean that the
  bootstrap has finished. If the callback returns `FALSE`, the bootstrap
  stops before the next replication or block with an error of class
  `dmar_bootstrap_cancelled`, whose elements `done` and `total` say how
  far it got (see
  [`dmar_metadata`](https://yelleknek.github.io/DMAR/reference/dmar_metadata.md));
  any other value is ignored. An error inside the callback stops the
  call with an error of class `dmar_progress_error`, whose message names
  the callback and which carries `done`, `total`, and, as `parent`, the
  callback's own condition with its class intact, so that an application
  can tell a failing callback from a failing analysis. The callback runs
  outside the random number stream, so a given `seed` returns the same
  result with or without it, and a callback that draws random numbers
  leaves the caller's stream as it found it. With `NULL` nothing is
  called. Ignored unless a bootstrap method is selected.

## Value

The `data.frame` returned by the dispatched `reliability_*` function
(rows: `estimate`, `se`, `ci_lower`, `ci_upper`, `conf_level`, `N`,
`N_complete`, `n_items`). The `se` row is on the coefficient scale; the
transformation-based intervals (`"fisher"`, `"bonett"`,
`"hakstian_whalen"`) add an `se_transformed` row carrying the
transformation-scale standard error, with the scale named in the
`se_transform_scale` attribute. The `coefficient` attribute identifies
which coefficient was computed, coefficients alpha and omega record
`wald_df`, and a bootstrap call also carries the `bootstrap` attribute
with its replication counts (see
[`reliability_alpha`](https://yelleknek.github.io/DMAR/reference/reliability_alpha.md)).

## Details

Auto-detection rules (used only when `type = NULL`):

- If raw items are integer-valued and every item has at most 10 distinct
  values, `type = "omega_categorical"` (categorical omega; appropriate
  when items are ordered categorical and the relationship between the
  underlying factor and the observed items is non-linear).

- Otherwise, `type = "omega"` (McDonald's coefficient \\\omega\\ from a
  single-factor CFA). This includes the case where only a covariance
  matrix `S` and sample size `N` are supplied; lavaan fits the CFA on
  the covariance matrix.

The rule reads the items only. The columns of `data` named in `aux` are
auxiliary variables, which never enter the composite, so their scale
says nothing about the items' response format and they take no part.
Auto-detection emits a single
[`message()`](https://rdrr.io/r/base/message.html) indicating which
`type` was chosen, so the choice is stated rather than made silently. It
never chooses `type = "H"`, which describes a different composite (the
optimally weighted one) rather than a different set of assumptions about
the summed one, and is computed only when asked for.

The selected family function determines which `ci_method` values are
accepted; see the help page for the chosen function for the full list.
When `ci_method` is left at its default (`NULL`), the family function's
own default is used:

- `reliability_alpha`: `"bonett"`, or `"ml"` when incomplete items are
  analyzed by full information maximum likelihood (the default of
  `missing` for them).

- `reliability_alpha(estimator = "model_implied")`: `"mlr"` with raw
  data; `"ml"` with covariance input.

- `reliability_kr20`: `"feldt"`.

- `reliability_omega`: for `denominator = "observed"` (robust omega, the
  default), the point estimate with no interval, since its interval is
  bootstrap based and no bootstrap runs unless requested; for
  `denominator = "model_implied"`, `"mlr_logistic"` with raw data and
  `"ml_logistic"` with covariance input, the limits formed on the logit
  scale (`"mlr"` and `"ml"` under `wald_df = "missing_information"`).

- `reliability_omega_categorical`: the point estimate with no interval,
  for the same reason; request `"bca"`.

- `reliability_H`: `"mlr_logistic"` with raw data; `"ml_logistic"` with
  covariance input.

A bootstrap is never run by default anywhere in the family. When a
bootstrap method is requested, `B = 10000` replications is the default,
except for categorical omega under `missing = "mi"`, whose bootstrap
then impute replications each impute the data and default to `B = 1000`.

Several reliability coefficients exist because their assumptions differ.
Coefficient \\\alpha\\ (and its dichotomous specialization KR-20) equals
the population reliability under essential \\\tau\\-equivalence (equal
loadings); McDonald's \\\omega\\ relaxes that assumption to a congeneric
single-factor model; and `reliability_omega(denominator = "observed")`
further relaxes the requirement that the single-factor model fit
perfectly by using the observed composite variance in the denominator
(see the
[`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md)
help page for the properties of that choice). For well-behaved
homogeneous measurement instruments these coefficients typically yield
very similar values. For ordered categorical items the relationship
between the latent factor and the observed responses is non-linear, and
`reliability_omega_categorical` handles that case explicitly via a
probit-link single-factor model.

## References

Green, S. B., & Yang, Y. (2009). Reliability of summed item scores using
structural equation modeling: An alternative to coefficient alpha.
*Psychometrika, 74*, 155–167.
[doi:10.1007/s11336-008-9099-3](https://doi.org/10.1007/s11336-008-9099-3)

Kelley, K., & Cheng, Y. (2012). Estimation of and confidence interval
formation for reliability coefficients of homogeneous measurement
instruments. *Methodology, 8*, 39–50.
[doi:10.1027/1614-2241/a000036](https://doi.org/10.1027/1614-2241/a000036)

Kelley, K., & Pornprasertmanit, S. (2016). Confidence intervals for
population reliability coefficients: Evaluation of methods,
recommendations, and software for composite measures. *Psychological
Methods, 21*, 69–92.
[doi:10.1037/a0040086](https://doi.org/10.1037/a0040086)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.

Terry, L. J., & Kelley, K. (2012). Sample size planning for composite
reliability coefficients: Accuracy in parameter estimation via narrow
confidence intervals. *British Journal of Mathematical and Statistical
Psychology, 65*, 371–401.
[doi:10.1111/j.2044-8317.2011.02030.x](https://doi.org/10.1111/j.2044-8317.2011.02030.x)

## See also

[`reliability_alpha`](https://yelleknek.github.io/DMAR/reference/reliability_alpha.md),
[`reliability_kr20`](https://yelleknek.github.io/DMAR/reference/reliability_kr20.md),
[`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md),
[`reliability_omega_categorical`](https://yelleknek.github.io/DMAR/reference/reliability_omega_categorical.md),
[`ss_aipe_reliability`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reliability.md),
[`cfa_1`](https://yelleknek.github.io/DMAR/reference/cfa_1.md).

Other reliability:
[`cohen_kappa()`](https://yelleknek.github.io/DMAR/reference/cohen_kappa.md),
[`diagnosis_agreement`](https://yelleknek.github.io/DMAR/reference/diagnosis_agreement.md),
[`fleiss_kappa()`](https://yelleknek.github.io/DMAR/reference/fleiss_kappa.md),
[`icc()`](https://yelleknek.github.io/DMAR/reference/icc.md),
[`reliability_H()`](https://yelleknek.github.io/DMAR/reference/reliability_H.md),
[`reliability_alpha()`](https://yelleknek.github.io/DMAR/reference/reliability_alpha.md),
[`reliability_kr20()`](https://yelleknek.github.io/DMAR/reference/reliability_kr20.md),
[`reliability_omega()`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md),
[`reliability_omega_categorical()`](https://yelleknek.github.io/DMAR/reference/reliability_omega_categorical.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
set.seed(113)
J <- 6
loadings <- seq(0.4, 0.8, length.out = J)
eta <- rnorm(200)
errors <- matrix(rnorm(200 * J), 200, J) %*% diag(sqrt(1 - loadings^2))
items <- sweep(matrix(rep(eta, J), 200, J), 2, loadings, `*`) + errors
colnames(items) <- paste0("y", seq_len(J))

# Auto-detection picks coefficient omega for continuous data.
reliability(data = items)
#> Auto-detected type = "omega" (McDonald's coefficient omega from a single-factor CFA).
#> Robust omega is reported without a confidence interval by default because its interval is bootstrap based. Request it with ci_method = "percentile" (or "bca"); B = 10000 replications is the default when you do.
#>  term       value
#>  estimate   0.781
#>  se         <NA> 
#>  ci_lower   <NA> 
#>  ci_upper   <NA> 
#>  conf_level 0.950
#>  N          200  
#>  N_complete 200  
#>  n_items    6    
#> 
#> Denominator: observed composite variance (robust omega)
#> Largest residual correlation under the single-factor model: -0.048 (y2 with y5)

# Explicit type.
reliability(data = items, type = "alpha")
#>  term           value
#>  estimate       0.767
#>  se             0.026
#>  se_transformed 0.110
#>  ci_lower       0.710
#>  ci_upper       0.812
#>  conf_level     0.950
#>  N              200  
#>  N_complete     200  
#>  n_items        6    

# A covariance matrix and its sample size stand in for raw data, and
# auto-detection again picks coefficient omega. The point estimate is
# the same as from the raw data; the message explains why no interval
# accompanies it, since the bootstrap behind robust omega resamples
# rows that a covariance matrix does not carry.
reliability(S = cov(items), N = 200)
#> Auto-detected type = "omega" (McDonald's coefficient omega from a single-factor CFA).
#> Robust omega is reported without a confidence interval: its interval is bootstrap based, which requires raw data rather than a covariance matrix.
#>  term       value
#>  estimate   0.781
#>  se         <NA> 
#>  ci_lower   <NA> 
#>  ci_upper   <NA> 
#>  conf_level 0.950
#>  N          200  
#>  N_complete 200  
#>  n_items    6    
#> 
#> Denominator: observed composite variance (robust omega)
#> Largest residual correlation under the single-factor model: -0.048 (y2 with y5)
```
