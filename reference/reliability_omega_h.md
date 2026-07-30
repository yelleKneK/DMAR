# Hierarchical omega with a bootstrap confidence interval

Estimates hierarchical omega (\\\omega_H\\) for a homogeneous composite
score and returns a bootstrap confidence interval for the population
coefficient. Hierarchical omega differs from
[`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md)
(McDonald's coefficient \\\omega\\) in that it uses the *observed*
variance of the composite in the denominator rather than the model
implied total variance, which makes \\\omega_H\\ robust to minor-factor
misspecification of the single-factor model.

## Usage

``` r
reliability_omega_h(
  data = NULL,
  S = NULL,
  N = NULL,
  ci_method = c("percentile", "bca", "bootstrap_se", "bootstrap_se_logistic", "none"),
  conf_level = 0.95,
  B = 10000,
  seed = NULL
)
```

## Arguments

- data:

  A numeric matrix or data frame of item scores. Required for confidence
  interval estimation (the bootstrap resamples raw data). Rows with any
  missing values are listwise-deleted.

- S:

  A symmetric covariance matrix among the items. Used only for the point
  estimate when raw data are not available; in that case `N` must be
  supplied and `ci_method` must be `"none"`.

- N:

  Total sample size; required when `S` (rather than `data`) is supplied.

- ci_method:

  Method for constructing the confidence interval. See *Details*.
  Defaults to `"percentile"`.

- conf_level:

  Confidence level for the interval. Defaults to `0.95`.

- B:

  Number of bootstrap replications. Defaults to `10000`.

- seed:

  Random number seed used for bootstrap reproducibility. Defaults to
  `NULL`, which leaves the user's current RNG state intact; supply an
  integer for reproducibility.

## Value

A `data.frame` with columns `term` and `value` and rows `"estimate"`
(sample \\\omega_H\\), `"se"` (bootstrap standard deviation across
replications, `NA` for `ci_method = "none"`), `"lower_limit"` and
`"upper_limit"` (clamped to \[0, 1\]), `"conf_level"`, `"N"`, and `"J"`.
Attributes `coefficient` (`"omega_h"`), `ci_method`, and `B` record the
computation.

## Details

Hierarchical omega was introduced by Kelley and Pornprasertmanit (2016,
Eq. 15) as a generalization of McDonald's (1999) coefficient \\\omega\\.
For a \\J\\-item composite \\Y = \sum_j X_j\\ with population factor
loadings \\\lambda_j\\ from a single-factor model, the population
coefficient is \$\$\omega_H = \frac{\left(\sum\_{j}
\lambda\_{j}\right)^{2}}{\sigma\_{Y}^{2}},\$\$ where \\\sigma\_{Y}^{2}\\
is the variance of the composite. In sample form, the numerator is the
squared sum of loadings from a single-factor CFA fit, and the
denominator is the observed variance of the composite, equivalently
`sum(stats::cov(data))` on raw items.

Using the observed (rather than model implied) total variance is the key
distinction from coefficient \\\omega\\. When the single-factor model
fits the items perfectly, \\\omega_H\\ and \\\omega\\ are equal; for
well-behaved homogeneous measurement instruments they are typically very
similar. When there are minor factors or correlated errors among the
items, the model implied total variance and the observed total variance
diverge; \\\omega_H\\ absorbs that divergence into the denominator and
so is unaffected by it. In their Monte Carlo studies, Kelley and
Pornprasertmanit (2016) found that bootstrap confidence intervals for
\\\omega_H\\ achieved acceptable coverage across the conditions
evaluated, including conditions with non-trivial minor-factor
misspecification.

Available confidence interval methods (set via `ci_method`):

- `"percentile"`:

  Percentile bootstrap (Efron & Tibshirani, 1993). Default; performed
  best in Kelley and Pornprasertmanit's (2016) Monte Carlo studies for
  continuous items.

- `"bca"`:

  Bias-corrected and accelerated bootstrap.

- `"bootstrap_se"`, `"bootstrap_se_logistic"`:

  Wald intervals using the bootstrap standard deviation as a standard
  error; the `_logistic` variant applies Browne's (1982) logit
  transformation.

- `"none"`:

  Return only the point estimate.

No closed-form standard error has been derived for \\\omega_H\\ because
the denominator is the observed composite variance rather than a model
parameter; Wald / delta method intervals from the single-factor fit
(e.g., `ml`, `mlr`, `adf`) describe coefficient \\\omega\\, not
\\\omega_H\\. Use
[`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md)
when a closed-form CI is needed and the single-factor model is assumed
to fit perfectly.

**Comparison with other packages.** The psych package's
[`omega`](https://rdrr.io/pkg/psych/man/omega.html) function reports
\\\omega_h\\ (lowercase *h*) defined relative to a Schmid-Leiman general
factor extracted from a hierarchical factor model with specific group
factors. That quantity is conceptually distinct from the \\\omega_H\\
estimated here, which is derived from a single common factor and uses
the observed composite variance to absorb minor-factor variability. The
two coincide in the special case of a homogeneous scale with no specific
factors.

## References

Browne, M. W. (1982). Covariance structures. In D. M. Hawkins (Ed.),
*Topics in applied multivariate analysis* (pp. 72–141). Cambridge, UK:
Cambridge University Press.

Efron, B., & Tibshirani, R. J. (1993). *An introduction to the
bootstrap*. New York, NY: Chapman & Hall/CRC.

Kelley, K., & Cheng, Y. (2012). Estimation of and confidence interval
formation for reliability coefficients of homogeneous measurement
instruments. *Methodology, 8*, 39–50.

Kelley, K., & Pornprasertmanit, S. (2016). Confidence intervals for
population reliability coefficients: Evaluation of methods,
recommendations, and software for composite measures. *Psychological
Methods, 21*, 69–92.

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.

McDonald, R. P. (1999). *Test theory: A unified treatment*. Mahwah, NJ:
Lawrence Erlbaum Associates.

Terry, L. J., & Kelley, K. (2012). Sample size planning for composite
reliability coefficients: Accuracy in parameter estimation via narrow
confidence intervals. *British Journal of Mathematical and Statistical
Psychology, 65*, 371–401.

Zinbarg, R. E., Yovel, I., Revelle, W., & McDonald, R. P. (2006).
Estimating generalizability to a latent variable common to all of a
scale's indicators: A comparison of estimators for \\\omega_h\\.
*Applied Psychological Measurement, 30*, 121–144.

## See also

[`reliability`](https://yelleknek.github.io/DMAR/reference/reliability.md)
(general wrapper),
[`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md),
[`reliability_omega_c`](https://yelleknek.github.io/DMAR/reference/reliability_omega_c.md),
[`reliability_alpha`](https://yelleknek.github.io/DMAR/reference/reliability_alpha.md),
[`cfa_1`](https://yelleknek.github.io/DMAR/reference/cfa_1.md)
(single-factor CFA used internally),
[`omega`](https://rdrr.io/pkg/psych/man/omega.html).

Other reliability:
[`reliability`](https://yelleknek.github.io/DMAR/reference/reliability.md)`()`,
[`reliability_H`](https://yelleknek.github.io/DMAR/reference/reliability_H.md)`()`,
[`reliability_alpha`](https://yelleknek.github.io/DMAR/reference/reliability_alpha.md)`()`,
[`reliability_kr20`](https://yelleknek.github.io/DMAR/reference/reliability_kr20.md)`()`,
[`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md)`()`,
[`reliability_omega_c`](https://yelleknek.github.io/DMAR/reference/reliability_omega_c.md)`()`

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# \donttest{
set.seed(113)
J <- 6
loadings <- seq(0.5, 0.8, length.out = J)
eta <- rnorm(200)
errors <- matrix(rnorm(200 * J), 200, J) %*% diag(sqrt(1 - loadings^2))
items <- sweep(matrix(rep(eta, J), 200, J), 2, loadings, `*`) + errors
colnames(items) <- paste0("y", seq_len(J))

# Default percentile bootstrap (paper's recommendation) with reduced B.
reliability_omega_h(data = items, B = 200)
#>  term        value 
#>  estimate    0.816 
#>  se          0.0203
#>  lower_limit 0.769 
#>  upper_limit 0.849 
#>  conf_level  0.95  
#>  N           200   
#>  J           6     
# }
```
