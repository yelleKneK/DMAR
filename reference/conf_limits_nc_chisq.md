# Confidence Limits for the Noncentrality Parameter of a Noncentral Chi Square Distribution

Finds the noncentrality parameters of a noncentral chi square
distribution that bracket an observed chi square value with the
requested tail probabilities, giving a confidence interval on the
population noncentrality parameter. Together with
[`conf_limits_nct`](https://yelleknek.github.io/DMAR/reference/conf_limits_nct.md)
and
[`conf_limits_ncf`](https://yelleknek.github.io/DMAR/reference/conf_limits_ncf.md),
this is one of the low-level noncentral distribution workhorses on which
the `ci_*` confidence interval functions are built; most analyses reach
it through those functions rather than calling it directly.

## Usage

``` r
conf_limits_nc_chisq(
  chi_square = NULL,
  conf_level = 0.95,
  df = NULL,
  alpha_lower = NULL,
  alpha_upper = NULL,
  tol = 1e-09,
  verbose = TRUE,
  ...
)
```

## Arguments

- chi_square:

  The observed chi square value

- conf_level:

  The desired degree of confidence for a symmetric interval

- df:

  The degrees of freedom

- alpha_lower:

  The proportion of values beyond the lower limit (cannot be used with
  `conf_level`)

- alpha_upper:

  The proportion of values beyond the upper limit (cannot be used with
  `conf_level`)

- tol:

  The convergence tolerance passed to
  [`uniroot`](https://rdrr.io/r/stats/uniroot.html)

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

  The noncentrality parameter at that limit. `0` when `alpha_lower = 0`
  or the lower limit is unattainable; `Inf` when `alpha_upper = 0`; `NA`
  when the observed `chi_square` is so small that even at \\\lambda =
  0\\ the lower-tail probability is already at or below `alpha_upper`,
  so the upper limit is undefined (a warning is issued).

- prob_less:

  (`verbose = TRUE`) The probability \\P(X \le \mathtt{chi\\square})\\
  that an observation from the noncentral chi square distribution
  centered at the row's limit falls at or below the observed
  `chi_square`.

- prob_greater:

  (`verbose = TRUE`) The complementary probability \\P(X \ge
  \mathtt{chi\\square})\\. By construction this equals `alpha_lower` on
  the `lower_limit` row and \\1 - \mathtt{alpha\\upper}\\ on the
  `upper_limit` row.

## Details

Each confidence limit is the noncentrality parameter \\\lambda \ge 0\\
of a noncentral chi square distribution with `df` degrees of freedom
whose appropriate tail at the observed `chi_square` contains the
requested probability:

- the lower limit satisfies \\P(X \ge \mathtt{chi\\square}) =
  \mathtt{alpha\\lower}\\;

- the upper limit satisfies \\P(X \le \mathtt{chi\\square}) =
  \mathtt{alpha\\upper}\\.

The two conditions run in opposite directions in \\\lambda\\: the
lower-tail probability \\P(X \le \mathtt{chi\\square})\\ is continuous
and strictly decreasing in the noncentrality parameter, so the
upper-tail probability \\P(X \ge \mathtt{chi\\square})\\ is continuous
and strictly increasing in it. The lower limit is the \\\lambda\\ at
which the upper tail has grown to `alpha_lower`, and the upper limit is
the \\\lambda\\ at which the lower tail has shrunk to `alpha_upper`.
Each is therefore the unique non-negative root of a one-dimensional
equation, and both are located with
[`uniroot`](https://rdrr.io/r/stats/uniroot.html) on the decreasing
lower-tail scale; `extendInt` is used to widen the search bracket if
needed.

Because the noncentrality parameter is bounded below by zero, the lower
limit is set to zero whenever the observed `chi_square` is smaller than
the `alpha_lower` critical value of the central chi square distribution
(i.e., the data is consistent with \\\lambda = 0\\ at the requested
confidence level). A warning is issued in that case, and the achieved
probabilities reported in the output reflect the actual values at
\\\lambda = 0\\.

Symmetrically, when the observed `chi_square` is so small that even at
\\\lambda = 0\\ the lower-tail probability is already at or below
`alpha_upper`, no \\\lambda \ge 0\\ places as much as `alpha_upper` mass
at or below `chi_square`; the upper limit is undefined and is returned
as `NA`, with a warning.

## References

Cumming, G., & Finch, S. (2001). A primer on the understanding, use, and
calculation of confidence intervals that are based on central and
noncentral distributions. *Educational and Psychological Measurement,
61*(4), 532–574.
[doi:10.1177/0013164401614002](https://doi.org/10.1177/0013164401614002)

Kelley, K. (2007). Confidence intervals for standardized effect sizes:
Theory, application, and implementation. *Journal of Statistical
Software, 20*(8), 1–24.
[doi:10.18637/jss.v020.i08](https://doi.org/10.18637/jss.v020.i08)

## See also

[`conf_limits_nct`](https://yelleknek.github.io/DMAR/reference/conf_limits_nct.md),
[`conf_limits_ncf`](https://yelleknek.github.io/DMAR/reference/conf_limits_ncf.md),
[`stats::pchisq()`](https://rdrr.io/r/stats/Chisquare.html),
[`stats::qchisq()`](https://rdrr.io/r/stats/Chisquare.html),
[`uniroot`](https://rdrr.io/r/stats/uniroot.html)

Other noncentral distribution confidence limits:
[`conf_limits_ncf()`](https://yelleknek.github.io/DMAR/reference/conf_limits_ncf.md),
[`conf_limits_nct()`](https://yelleknek.github.io/DMAR/reference/conf_limits_nct.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# A typical call to the function.
conf_limits_nc_chisq(chi_square = 30, conf_level = .95, df = 15)
#>          term     value prob_less prob_greater
#> 1 lower_limit  1.407074     0.975        0.025
#> 2 upper_limit 38.876511     0.025        0.975

# A one-sided (upper) confidence interval.
conf_limits_nc_chisq(chi_square = 30, alpha_lower = 0, alpha_upper = .05,
                     conf_level = NULL, df = 15)
#>          term   value prob_less prob_greater
#> 1 lower_limit  0.0000      1.00         0.00
#> 2 upper_limit 34.6284      0.05         0.95
```
