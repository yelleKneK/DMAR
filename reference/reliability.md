# Reliability Coefficient With a Confidence Interval (General Dispatch)

General-purpose entry point for the reliability family. Dispatches to
[`reliability_alpha`](https://yelleknek.github.io/DMAR/reference/reliability_alpha.md),
[`reliability_kr20`](https://yelleknek.github.io/DMAR/reference/reliability_kr20.md),
[`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md),
or
[`reliability_omega_categorical`](https://yelleknek.github.io/DMAR/reference/reliability_omega_categorical.md)
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
  missing = c("listwise", "fiml"),
  aux = NULL,
  ci_method = NULL,
  conf_level = 0.95,
  B = 10000,
  seed = NULL
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
  `"omega_categorical"` (`"omega_c"` is accepted as a shorthand), or
  `NULL` for auto-detection. See *Details*.

- estimator:

  For `type = "alpha"` only: how coefficient \\\alpha\\ is estimated,
  `"analytic"` (default; the closed-form equation applied to the
  observed covariance matrix) or `"model_implied"` (the reliability
  implied by the \\\tau\\-equivalent single-factor model fit by maximum
  likelihood); forwarded to
  [`reliability_alpha`](https://yelleknek.github.io/DMAR/reference/reliability_alpha.md),
  whose help page discusses the choice and which interval methods each
  estimator supports. Supplying it with any other `type` is an error, as
  `denominator` is for any type but `"omega"`.

- denominator:

  For `type = "omega"` only: how the total variance in the denominator
  of \\\omega\\ is estimated, `"observed"` (default; robust omega) or
  `"model_implied"`; forwarded to
  [`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md),
  whose help page discusses the choice.

- missing:

  For `type = "alpha"` and `type = "omega"` only: how incomplete rows of
  `data` are handled, `"listwise"` (the default) or `"fiml"` (full
  information maximum likelihood); forwarded to the family function,
  whose help page discusses the choice. Supplying it with any other
  `type` is an error.

- aux:

  For `type = "alpha"` and `type = "omega"` only: optional character
  vector naming auxiliary variable columns of `data`, entered as
  saturated correlates under full information maximum likelihood;
  forwarded to the family function. Supplying `aux` implies
  `missing = "fiml"`.

- ci_method:

  Method for constructing the confidence interval, or `NULL` to use the
  chosen family function's default. See the help page for the chosen
  `reliability_*` function for the full list of accepted values.

- conf_level:

  Confidence level. Defaults to `0.95`.

- B:

  Number of bootstrap replications when a bootstrap method is selected.
  Defaults to `10000`.

- seed:

  Random number seed used for bootstrap reproducibility. Defaults to
  `NULL`, which leaves the user's current RNG state intact; supply an
  integer for reproducibility.

## Value

The `data.frame` returned by the dispatched `reliability_*` function
(rows: `estimate`, `se`, `lower_limit`, `upper_limit`, `conf_level`,
`N`, `N_complete`, `J`). The `coefficient` attribute identifies which
coefficient was computed.

## Details

Auto-detection rules (used only when `type = NULL`):

- If raw items are integer-valued and every column has at most 10
  distinct values, `type = "omega_categorical"` (categorical omega;
  appropriate when items are ordered-categorical and the relationship
  between the underlying factor and the observed items is non-linear).

- Otherwise, `type = "omega"` (McDonald's coefficient \\\omega\\ from a
  single-factor CFA). This includes the case where only a covariance
  matrix `S` and sample size `N` are supplied; lavaan fits the CFA on
  the covariance matrix.

Auto-detection emits a single
[`message()`](https://rdrr.io/r/base/message.html) indicating which
`type` was chosen, so it never surprises the user silently.

The selected family function determines which `ci_method` values are
accepted; see the help page for the chosen function for the full list.
When `ci_method` is left at its default (`NULL`), the family function's
own default is used:

- `reliability_alpha`: `"bonett"`.

- `reliability_alpha(estimator = "model_implied")`: `"mlr"` with raw
  data; `"ml"` with covariance input.

- `reliability_kr20`: `"feldt"`.

- `reliability_omega`: for `denominator = "observed"` (robust omega, the
  default), the point estimate with no interval, since its interval is
  bootstrap based and no bootstrap runs unless requested; for
  `denominator = "model_implied"`, `"mlr"`.

- `reliability_omega_categorical`: the point estimate with no interval,
  for the same reason; request `"bca"`.

A bootstrap is never run by default anywhere in the family. When a
bootstrap method is requested, `B = 10000` replications is the default.

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
very similar values. For ordered-categorical items the relationship
between the latent factor and the observed responses is non-linear, and
`reliability_omega_categorical` handles that case explicitly via a
probit-link single-factor model.

## References

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
# \donttest{
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
#>  term        value
#>  estimate    0.781
#>  se          <NA> 
#>  lower_limit <NA> 
#>  upper_limit <NA> 
#>  conf_level  0.95 
#>  N           200  
#>  N_complete  200  
#>  J           6    

# Explicit type.
reliability(data = items, type = "alpha")
#>  term        value
#>  estimate    0.767
#>  se          0.11 
#>  lower_limit 0.71 
#>  upper_limit 0.812
#>  conf_level  0.95 
#>  N           200  
#>  N_complete  200  
#>  J           6    

# From a covariance matrix (also picks coefficient omega).
reliability(S = cov(items), N = 200)
#> Auto-detected type = "omega" (McDonald's coefficient omega from a single-factor CFA).
#> Robust omega is reported without a confidence interval: its interval is bootstrap based, which requires raw data rather than a covariance matrix.
#>  term        value
#>  estimate    0.781
#>  se          <NA> 
#>  lower_limit <NA> 
#>  upper_limit <NA> 
#>  conf_level  0.95 
#>  N           200  
#>  N_complete  200  
#>  J           6    
# }
```
