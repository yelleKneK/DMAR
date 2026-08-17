# Confidence Limits for the Noncentrality Parameter of a Noncentral *F*-distribution

Finds the noncentrality parameters of a noncentral *F*-distribution that
bracket an observed *F*-value with the requested tail probabilities,
giving a confidence interval on the population noncentrality parameter.
Together with
[`conf_limits_nct`](https://yelleknek.github.io/DMAR/reference/conf_limits_nct.md)
and
[`conf_limits_nc_chisq`](https://yelleknek.github.io/DMAR/reference/conf_limits_nc_chisq.md),
this is one of the low-level noncentral distribution workhorses on which
the `ci_*` confidence interval functions (e.g.,
[`ci_pvaf`](https://yelleknek.github.io/DMAR/reference/ci_pvaf.md),
[`ci_snr`](https://yelleknek.github.io/DMAR/reference/ci_snr.md),
[`ci_R2`](https://yelleknek.github.io/DMAR/reference/ci_R2.md)) are
built; most analyses reach it through those functions rather than
calling it directly.

## Usage

``` r
conf_limits_ncf(
  F_value = NULL,
  conf_level = 0.95,
  df_1 = NULL,
  df_2 = NULL,
  alpha_lower = NULL,
  alpha_upper = NULL,
  tol = 1e-09,
  verbose = TRUE,
  ...
)
```

## Arguments

- F_value:

  The observed *F*-value

- conf_level:

  The desired degree of confidence for a symmetric interval

- df_1:

  The numerator degrees of freedom

- df_2:

  The denominator degrees of freedom

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
  on the `upper_limit` row when the observed `F_value` is so small that
  even at \\\lambda = 0\\ the lower-tail probability is already at or
  below `alpha_upper`, leaving the upper noncentrality limit undefined
  (a warning is issued).

- prob_less:

  (`verbose = TRUE`) The probability \\P(F \le \mathtt{F\\value})\\ that
  an *F*-statistic from the noncentral *F*-distribution centered at the
  row's limit falls at or below the observed `F_value`.

- prob_greater:

  (`verbose = TRUE`) The complementary probability \\P(F \ge
  \mathtt{F\\value})\\. By construction this equals `alpha_lower` on the
  `lower_limit` row and \\1 - \mathtt{alpha\\upper}\\ on the
  `upper_limit` row.

## Details

Each confidence limit is the noncentrality parameter \\\lambda \ge 0\\
of a noncentral *F*-distribution with `df_1` and `df_2` degrees of
freedom whose appropriate tail at the observed `F_value` contains the
requested probability:

- the lower limit satisfies \\P(F \ge \mathtt{F\\value}) =
  \mathtt{alpha\\lower}\\;

- the upper limit satisfies \\P(F \le \mathtt{F\\value}) =
  \mathtt{alpha\\upper}\\.

The two conditions run in opposite directions in \\\lambda\\: the
lower-tail probability \\P(F \le \mathtt{F\\value})\\ is continuous and
strictly decreasing in the noncentrality parameter, so the upper-tail
probability \\P(F \ge \mathtt{F\\value})\\ is continuous and strictly
increasing in it. The lower limit is the \\\lambda\\ at which the upper
tail has grown to `alpha_lower`, and the upper limit is the \\\lambda\\
at which the lower tail has shrunk to `alpha_upper`. Each is therefore
the unique non-negative root of a one-dimensional equation, and both are
located with [`uniroot`](https://rdrr.io/r/stats/uniroot.html) on the
decreasing lower-tail scale; `extendInt` is used to widen the search
bracket if needed.

Because the noncentrality parameter is bounded below by zero, the lower
limit is set to zero whenever the observed `F_value` is smaller than the
`alpha_lower` critical value of the central *F*-distribution (i.e., the
data is consistent with \\\lambda = 0\\ at the requested confidence
level). A warning is issued in that case, and the achieved probabilities
reported in the output reflect the actual values at \\\lambda = 0\\
rather than the requested `alpha_lower`. The warning carries the
condition class `dmar_ncf_clamp`, so a caller that inverts the
noncentral *F* repeatedly can muffle or deduplicate it by class.

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

[`ss_aipe_R2`](https://yelleknek.github.io/DMAR/reference/ss_aipe_R2.md),
[`ci_R2`](https://yelleknek.github.io/DMAR/reference/ci_R2.md),
[`conf_limits_nct`](https://yelleknek.github.io/DMAR/reference/conf_limits_nct.md),
[`conf_limits_nc_chisq`](https://yelleknek.github.io/DMAR/reference/conf_limits_nc_chisq.md),
[`stats::pf()`](https://rdrr.io/r/stats/Fdist.html),
[`stats::qf()`](https://rdrr.io/r/stats/Fdist.html),
[`uniroot`](https://rdrr.io/r/stats/uniroot.html)

Other noncentral distribution confidence limits:
[`conf_limits_nc_chisq()`](https://yelleknek.github.io/DMAR/reference/conf_limits_nc_chisq.md),
[`conf_limits_nct()`](https://yelleknek.github.io/DMAR/reference/conf_limits_nct.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
conf_limits_ncf(F_value = 5, conf_level = .95, df_1 = 5, df_2 = 100)
#>          term     value prob_less prob_greater
#> 1 lower_limit  5.353713     0.975        0.025
#> 2 upper_limit 45.276111     0.025        0.975

# A one-sided (upper) confidence interval.
conf_limits_ncf(F_value = 5, conf_level = NULL, df_1 = 5, df_2 = 100,
                alpha_lower = 0, alpha_upper = .05)
#>          term    value prob_less prob_greater
#> 1 lower_limit  0.00000      1.00         0.00
#> 2 upper_limit 40.74672      0.05         0.95
```
