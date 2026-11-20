# Confidence Interval for the Noncentrality Parameter of a Noncentral *t* Distribution

Finds the noncentrality parameters of a noncentral *t*-distribution that
bracket an observed *t*-value with the requested tail probabilities,
giving a confidence interval on the population noncentrality parameter.
Together with
[`ci_nc_F`](https://yelleknek.github.io/DMAR/reference/ci_nc_F.md) and
[`ci_nc_chisq`](https://yelleknek.github.io/DMAR/reference/ci_nc_chisq.md),
this is one of the low-level noncentral distribution workhorses on which
the `ci_*` confidence interval functions (e.g.,
[`ci_smd`](https://yelleknek.github.io/DMAR/reference/ci_smd.md),
[`ci_smd_c`](https://yelleknek.github.io/DMAR/reference/ci_smd_c.md),
[`ci_cv`](https://yelleknek.github.io/DMAR/reference/ci_cv.md)) are
built; most analyses reach it through those functions rather than
calling it directly. The function corresponds to `conf.limits.nct()` in
MBESS; its name places it in the `ci_*` family because a confidence
interval is what it computes.

## Usage

``` r
ci_nc_t(
  ncp,
  df,
  conf_level = 0.95,
  alpha_lower = NULL,
  alpha_upper = NULL,
  tol = 1e-09,
  verbose = TRUE,
  ...
)
```

## Arguments

- ncp:

  The observed *t*-value, which is the estimate of the noncentrality
  parameter the interval is formed for

- df:

  The degrees of freedom

- conf_level:

  The level of confidence for a symmetric confidence interval. A level
  above .999999998, whose tail areas are below 1e-9, is refused (see
  "Small tail areas" in Details)

- alpha_lower:

  The proportion of values beyond the lower limit of the confidence
  interval (cannot be used with `conf_level`): 0, for an interval with
  no lower limit, or at least 1e-9 (see "Small tail areas" in Details)

- alpha_upper:

  The proportion of values beyond the upper limit of the confidence
  interval (cannot be used with `conf_level`): 0, for an interval with
  no upper limit, or at least 1e-9

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

  Either `"ci_lower"` or `"ci_upper"`.

- value:

  The noncentrality parameter at that limit. `-Inf` when
  `alpha_lower = 0`; `Inf` when `alpha_upper = 0`.

- prob_less:

  (`verbose = TRUE`) The probability \\P(T \le \mathrm{ncp})\\ that a
  *t*-statistic from the noncentral *t*-distribution centered at the
  row's limit falls at or below the observed `ncp`. By construction this
  equals `alpha_upper` on the `ci_upper` row and \\1 -
  \mathtt{alpha\\lower}\\ on the `ci_lower` row.

- prob_greater:

  (`verbose = TRUE`) The complementary probability \\P(T \ge
  \mathrm{ncp})\\. By construction this equals `alpha_lower` on the
  `ci_lower` row and \\1 - \mathtt{alpha\\upper}\\ on the `ci_upper`
  row.

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

**Computation.** R's [`pt()`](https://rdrr.io/r/stats/TDist.html) sums a
series for the noncentral *t* distribution up to a noncentrality of
37.62 and with at most 400,000 degrees of freedom, and returns a normal
approximation otherwise. Compared with the exact probability over a grid
of degrees of freedom, noncentralities, and quantiles (from \\10^{-6}\\
to \\1 - 10^{-6}\\), it was within \\10^{-9}\\ up to a noncentrality of
37.62 with at most 2,000 degrees of freedom and up to 33 with more (the
series stops short in the upper tail near the top of its range: the
error was \\10^{-6}\\ at 34 and .001 at 36 with 100,000 degrees of
freedom), while the normal approximation was off by as much as .14 with
one degree of freedom, .008 with 38, and .001 with 1,000; below one
degree of freedom it is not accurate at any noncentrality other than 0.
Each limit is therefore found with
[`pt()`](https://rdrr.io/r/stats/TDist.html) where the limit lies in
that range (with at least one degree of freedom), and otherwise with the
exact probability \$\$P(T \le t) = \int_0^\infty \Phi(t s - \lambda)\\
f(s)\\ ds,\$\$ where \\\Phi\\ is the standard normal distribution
function, \\f\\ the density of \\S = \sqrt{V/\nu}\\ for \\V\\ chi square
on \\\nu\\ degrees of freedom, and \\T = (Z + \lambda)/S\\ with \\Z\\
standard normal, independent of \\V\\; the integral is computed
numerically. The limits agreed with bisection on a second exact form,
which conditions on the normal variate instead, to within \\3.5 \times
10^{-9}\\ of the width of the interval in 400 random cases with 1 to
1,000,000 degrees of freedom. The tail probabilities in `prob_less` and
`prob_greater` are computed the same way as the limit in their row. A
limit outside the range takes about a millisecond; one inside it, a
tenth of that.

**Small tail areas.** A nonzero tail area below 1e-9, or a `conf_level`
above .999999998, is refused with an error of class
`dmar_argument_invalid` whose `reason` is `"tail_below_floor"`. A limit
is the root of a tail probability, and
[`pt()`](https://rdrr.io/r/stats/TDist.html) sums its series to an
absolute tolerance of \\10^{-12}\\ and forms one tail as 1 minus the
other, so the error of the probability becomes a larger part of the tail
as the tail shrinks. Against the exact probability, the limits with
equal tails at \\t\\ = 5 on 47 degrees of freedom, .5 on 10, 2 on 200,
10 on 20, -3 on 30, and 1 on 1,000, at 66 tail areas from \\10^{-6}\\
down to \\10^{-9}\\, were off by at most \\1.7 \times 10^{-5}\\ at tail
areas of \\10^{-8}\\ and above, \\7.3 \times 10^{-5}\\ at \\2 \times
10^{-9}\\ and above, and \\1.4 \times 10^{-4}\\ at \\10^{-9}\\ (the
upper limit at \\t\\ = 10 on 20 degrees of freedom, \\6.0 \times
10^{-6}\\ of the width of the interval); at \\10^{-10}\\ the largest
error was \\1.1 \times 10^{-3}\\, and at \\10^{-20}\\ the upper limit at
\\t\\ = 5 on 47 degrees of freedom was 7.6 percent too small. Every
interval and plan whose limits come from this function keeps the floor,
among them
[`ci_smd`](https://yelleknek.github.io/DMAR/reference/ci_smd.md),
[`ci_smd_c`](https://yelleknek.github.io/DMAR/reference/ci_smd_c.md),
[`ci_sm`](https://yelleknek.github.io/DMAR/reference/ci_sm.md),
[`ci_cv`](https://yelleknek.github.io/DMAR/reference/ci_cv.md),
[`ci_sc`](https://yelleknek.github.io/DMAR/reference/ci_sc.md),
[`ci_sc_ancova`](https://yelleknek.github.io/DMAR/reference/ci_sc_ancova.md),
[`ci_src`](https://yelleknek.github.io/DMAR/reference/ci_src.md), and
[`ci_reg_coef`](https://yelleknek.github.io/DMAR/reference/ci_reg_coef.md)
with `noncentral = TRUE`, and the planners of those intervals. The
central intervals of a contrast, a correlation, and a regression
coefficient
([`ci_c`](https://yelleknek.github.io/DMAR/reference/ci_c.md),
[`ci_r`](https://yelleknek.github.io/DMAR/reference/ci_correlation.md),
[`ci_rc`](https://yelleknek.github.io/DMAR/reference/ci_rc.md), and
[`ci_reg_coef`](https://yelleknek.github.io/DMAR/reference/ci_reg_coef.md)
with `noncentral = FALSE`) and their planners
([`ss_aipe_c`](https://yelleknek.github.io/DMAR/reference/ss_aipe_c.md),
[`ss_aipe_r`](https://yelleknek.github.io/DMAR/reference/ss_aipe_r.md),
and
[`ss_aipe_rc`](https://yelleknek.github.io/DMAR/reference/ss_aipe_rc.md))
keep it too, so that an interval and the plan for it accept the same
tail areas.

**A limit past 37.62.** Each limit is found by its own rule, so an
interval can have its lower limit from
[`pt()`](https://rdrr.io/r/stats/TDist.html) and its upper limit from
the exact probability. The interval for the standardized mean difference
of two groups of 20 at \\d = 10\\ with 99 percent confidence
([`ci_smd`](https://yelleknek.github.io/DMAR/reference/ci_smd.md)) is
such an interval: it is \[7.003, 13.099\], 6.096 wide. In 100,000
samples of two groups of 20 with \\\delta = 13\\, this 99 percent
interval fell below \\\delta\\ in .52 percent and above it in .48
percent (.50 percent each by design).

This function is especially useful for forming confidence intervals
around standardized mean differences (Cohen's *d*, Glass's *g*, Hedges'
*g*), standardized regression coefficients, and coefficients of
variation.

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
estimation and the evaluation of statistical models. In L. L. Harlow, S.
A. Mulaik, & J. H. Steiger (Eds.), *What if there were no significance
tests?* (pp. 221–257). Mahwah, NJ: Lawrence Erlbaum.

## See also

[`stats::pt()`](https://rdrr.io/r/stats/TDist.html),
[`stats::qt()`](https://rdrr.io/r/stats/TDist.html),
[`uniroot`](https://rdrr.io/r/stats/uniroot.html),
[`ci_smd`](https://yelleknek.github.io/DMAR/reference/ci_smd.md),
[`ci_smd_c`](https://yelleknek.github.io/DMAR/reference/ci_smd_c.md),
[`ci_nc_F`](https://yelleknek.github.io/DMAR/reference/ci_nc_F.md),
[`ci_nc_chisq`](https://yelleknek.github.io/DMAR/reference/ci_nc_chisq.md)

Other noncentral distribution confidence intervals:
[`ci_nc_F()`](https://yelleknek.github.io/DMAR/reference/ci_nc_F.md),
[`ci_nc_chisq()`](https://yelleknek.github.io/DMAR/reference/ci_nc_chisq.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Suppose observed t-value based on 'df'=126 is 2.83. Finding the lower
# and upper critical values for the population noncentrality parameter
# with a symmetric confidence interval with 95\% confidence is given as:
ci_nc_t(ncp = 2.83, df = 126, conf_level = .95)
#>       term     value prob_less prob_greater
#> 1 ci_lower 0.8337503     0.975        0.025
#> 2 ci_upper 4.8153591     0.025        0.975

# Modifying the above example so that a nonsymmetric 95% confidence interval
# can be formed:
ci_nc_t(ncp = 2.83, df = 126, alpha_lower = .01, alpha_upper = .04, conf_level = NULL)
#>       term    value prob_less prob_greater
#> 1 ci_lower 0.461692      0.99         0.01
#> 2 ci_upper 4.602743      0.04         0.96

# Modifying the above example so that a single-sided 95% confidence interval
# can be formed:
ci_nc_t(ncp = 2.83, df = 126, alpha_lower = 0, alpha_upper = .05, conf_level = NULL)
#>       term    value prob_less prob_greater
#> 1 ci_lower     -Inf      1.00         0.00
#> 2 ci_upper 4.495225      0.05         0.95
```
