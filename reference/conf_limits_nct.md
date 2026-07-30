# Confidence Limits for a Noncentrality Parameter From a *t*-distribution

Finds the noncentrality parameters of a noncentral *t*-distribution that
bracket an observed *t*-value with the requested tail probabilities,
giving a confidence interval on the population noncentrality parameter.
Together with
[`conf_limits_ncf`](https://yelleknek.github.io/DMAR/reference/conf_limits_ncf.md)
and
[`conf_limits_nc_chisq`](https://yelleknek.github.io/DMAR/reference/conf_limits_nc_chisq.md),
this is one of the low-level noncentral distribution workhorses on which
the `ci_*` confidence interval functions (e.g.,
[`ci_smd`](https://yelleknek.github.io/DMAR/reference/ci_smd.md),
[`ci_smd_c`](https://yelleknek.github.io/DMAR/reference/ci_smd_c.md),
[`ci_cv`](https://yelleknek.github.io/DMAR/reference/ci_cv.md)) are
built; most analyses reach it through those functions rather than
calling it directly.

## Usage

``` r
conf_limits_nct(
  ncp,
  df,
  conf_level = 0.95,
  alpha_lower = NULL,
  alpha_upper = NULL,
  t_value,
  tol = 1e-09,
  verbose = TRUE,
  ...
)
```

## Arguments

- ncp:

  The noncentrality parameter (e.g., observed *t*-value) of interest

- df:

  The degrees of freedom

- conf_level:

  The level of confidence for a symmetric confidence interval

- alpha_lower:

  The proportion of values beyond the lower limit of the confidence
  interval (cannot be used with `conf_level`)

- alpha_upper:

  The proportion of values beyond the upper limit of the confidence
  interval (cannot be used with `conf_level`)

- t_value:

  Alias for `ncp`

- tol:

  The convergence tolerance passed to
  [`uniroot`](https://rdrr.io/r/stats/uniroot.html) when locating each
  limit

- verbose:

  If `TRUE` (the default), the returned data frame additionally reports
  the achieved tail probabilities at each limit; if `FALSE`, only `term`
  and `value` are returned

- ...:

  Additional arguments forwarded to
  [`uniroot`](https://rdrr.io/r/stats/uniroot.html)

## Value

A `data.frame` with one row per confidence limit and the columns:

- term:

  Either `"lower_limit"` or `"upper_limit"`.

- value:

  The noncentrality parameter at that limit. `-Inf` when
  `alpha_lower = 0`; `Inf` when `alpha_upper = 0`.

- prob_less:

  (`verbose = TRUE`) The probability \\P(T \le \mathrm{ncp})\\ that a
  *t*-statistic from the noncentral *t*-distribution centered at the
  row's limit falls at or below the observed `ncp`. By construction this
  equals `alpha_upper` on the `upper_limit` row and \\1 -
  \mathtt{alpha\\lower}\\ on the `lower_limit` row.

- prob_greater:

  (`verbose = TRUE`) The complementary probability \\P(T \ge
  \mathrm{ncp})\\. By construction this equals `alpha_lower` on the
  `lower_limit` row and \\1 - \mathtt{alpha\\upper}\\ on the
  `upper_limit` row.

## Details

Each confidence limit is the noncentrality parameter of a noncentral
*t*-distribution with `df` degrees of freedom whose appropriate tail at
the observed `ncp` contains the requested probability:

- the lower limit satisfies \\P(T \ge \mathrm{ncp}) =
  \mathtt{alpha\\lower}\\;

- the upper limit satisfies \\P(T \le \mathrm{ncp}) =
  \mathtt{alpha\\upper}\\.

Each tail probability is continuous and strictly monotone in the
noncentrality parameter, so each limit is the unique root of a
one-dimensional equation. The roots are located with
[`uniroot`](https://rdrr.io/r/stats/uniroot.html) starting from a
bracket centered on `ncp` with half-width scaled by the asymptotic
standard error of the noncentrality estimator; `extendInt` is used to
widen the bracket if needed.

This function is especially useful for forming confidence intervals
around standardized mean differences (Cohen's *d*, Glass's *g*, Hedges'
*g*), standardized regression coefficients, and coefficients of
variation.

## Warning

As of R 4.0.0, the largest `ncp` that R can accurately handle is 37.62.

## References

Cumming, G., & Finch, S. (2001). A primer on the understanding, use, and
calculation of confidence intervals that are based on central and
noncentral distributions. *Educational and Psychological Measurement,
61*(4), 532–574.
[doi:10.1177/0013164401614002](https://doi.org/10.1177/0013164401614002)

Kelley, K. (2005). The effects of nonnormal distributions on confidence
intervals around the standardized mean difference: Bootstrap and
parametric confidence intervals, *Educational and Psychological
Measurement, 65*, 51–69.
[doi:10.1177/0013164404264850](https://doi.org/10.1177/0013164404264850)

Kelley, K. (2007). Confidence intervals for standardized effect sizes:
Theory, application, and implementation. *Journal of Statistical
Software, 20*(8), 1–24.
[doi:10.18637/jss.v020.i08](https://doi.org/10.18637/jss.v020.i08)

Steiger, J. H., & Fouladi, R. T. (1997). Noncentrality interval
estimation and the evaluation of statistical methods. In L. L. Harlow,
S. A. Mulaik, & J. H. Steiger (Eds.), *What if there were no
significance tests?* (pp. 221–257). Mahwah, NJ: Lawrence Erlbaum.

## See also

[`stats::pt()`](https://rdrr.io/r/stats/TDist.html),
[`stats::qt()`](https://rdrr.io/r/stats/TDist.html),
[`uniroot`](https://rdrr.io/r/stats/uniroot.html),
[`ci_smd`](https://yelleknek.github.io/DMAR/reference/ci_smd.md),
[`ci_smd_c`](https://yelleknek.github.io/DMAR/reference/ci_smd_c.md),
[`conf_limits_ncf`](https://yelleknek.github.io/DMAR/reference/conf_limits_ncf.md),
[`conf_limits_nc_chisq`](https://yelleknek.github.io/DMAR/reference/conf_limits_nc_chisq.md)

Other noncentral distribution confidence limits:
[`conf_limits_nc_chisq()`](https://yelleknek.github.io/DMAR/reference/conf_limits_nc_chisq.md),
[`conf_limits_ncf()`](https://yelleknek.github.io/DMAR/reference/conf_limits_ncf.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Suppose observed t-value based on 'df'=126 is 2.83. Finding the lower
# and upper critical values for the population noncentrality parameter
# with a symmetric confidence interval with 95\% confidence is given as:
conf_limits_nct(ncp = 2.83, df = 126, conf_level = .95)
#>          term     value prob_less prob_greater
#> 1 lower_limit 0.8337503     0.975        0.025
#> 2 upper_limit 4.8153591     0.025        0.975

# Modifying the above example so that a nonsymmetric 95% confidence interval
# can be formed:
conf_limits_nct(ncp = 2.83, df = 126, alpha_lower = .01, alpha_upper = .04, conf_level = NULL)
#>          term    value prob_less prob_greater
#> 1 lower_limit 0.461692      0.99         0.01
#> 2 upper_limit 4.602743      0.04         0.96

# Modifying the above example so that a single-sided 95% confidence interval
# can be formed:
conf_limits_nct(ncp = 2.83, df = 126, alpha_lower = 0, alpha_upper = .05, conf_level = NULL)
#>          term    value prob_less prob_greater
#> 1 lower_limit     -Inf      1.00         0.00
#> 2 upper_limit 4.495225      0.05         0.95
```
