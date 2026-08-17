# Power of the TOST or Noninferiority Test for a Linear Contrast

Computes the exact power of the Schuirmann (1987) two one-sided tests
procedure, or of the one-sided noninferiority test, for a linear
contrast of group means \\\psi = \sum_j c_j \mu_j\\ with one pooled
error term. For equivalence, the power is the probability that the
\\(1 - 2\alpha)\\ confidence interval for \\\psi\\ lies entirely inside
\\(-\delta_L, \delta_U)\\, computed by numerical integration over the
chi distribution of the estimated error standard deviation; for
noninferiority, the power is a noncentral *t* probability in closed
form. This is the contrast generalization of
[`power_equivalence_md`](https://yelleknek.github.io/DMAR/reference/power_equivalence_md.md).

## Usage

``` r
power_equivalence_c(
  c_weights,
  n,
  sigma,
  delta_lower = NULL,
  delta_upper = NULL,
  true_psi = 0,
  alpha_level = 0.05,
  side = c("equivalence", "noninferiority"),
  df_error = NULL
)
```

## Arguments

- c_weights:

  The contrast weights. The weights must sum to zero with the positive
  weights summing to 1 and the negative weights to -1, so that the
  bounds are on the raw scale of the response.

- n:

  Sample sizes per group (if length 1, equal group sizes are assumed).
  Together with `c_weights`, `n` determines the standard error factor
  \\\sqrt{\sum_j c_j^2 / n_j}\\.

- sigma:

  The error standard deviation (the square root of the mean square
  error).

- delta_lower, delta_upper:

  Equivalence bounds on the raw scale of the response. Both must be
  positive; the equivalence region is \\(-\delta_L, +\delta_U)\\. If
  only `delta_upper` is supplied, the bounds are symmetric.
  Noninferiority uses \\-\delta_L\\ alone.

- true_psi:

  The population value of the contrast at which the power is evaluated.
  Default `0`, the most favorable point for an equivalence declaration.

- alpha_level:

  One-sided significance level for each test. Default `0.05`.

- side:

  `"equivalence"` (default) for the TOST power, or `"noninferiority"`
  for the one-sided test against \\-\delta_L\\.

- df_error:

  The error degrees of freedom. Defaults to \\N - J\\; supply it
  directly when the error term comes from a model with more groups or
  additional predictors than the contrast involves (for example, a
  five-group model supplying the pooled error for a two-group contrast).

## Value

A one-row `data.frame` with columns `term` (`"power"`) and `value` (the
computed power, in \\\[0, 1\]\\).

## Details

**Equivalence power.** Conditional on the estimated error standard
deviation \\S\\, the \\(1 - 2\alpha)\\ CI fits inside the bounds on a
computable event, and the unconditional power integrates that event over
the scaled chi distribution of \\S\\ on `df_error` degrees of freedom.
With `c_weights = c(1, -1)` and equal `n`, the result reproduces
[`power_equivalence_md`](https://yelleknek.github.io/DMAR/reference/power_equivalence_md.md)
exactly.

**Noninferiority power.** The one-sided test rejects when \\t =
(\hat\psi + \delta_L)/\mathrm{SE}(\hat\psi)\\ exceeds
\\t\_{1-\alpha,\nu}\\, so the power is \\\Pr(T'\_{\nu}(\lambda) \>
t\_{1-\alpha,\nu})\\ with noncentrality \\\lambda = (\psi +
\delta_L)/(\sigma \sqrt{\sum_j c_j^2/n_j})\\.

**The feasibility condition.** If the expected half-width of the CI is
not smaller than the bounds allow, the equivalence power is zero or near
zero regardless of `true_psi`: an imprecise design cannot declare
equivalence even when the arms are truly identical. Planning should
target a half-width of about half the bound; see
[`ss_power_equivalence_c`](https://yelleknek.github.io/DMAR/reference/ss_power_equivalence_c.md)
and
[`ss_aipe_c`](https://yelleknek.github.io/DMAR/reference/ss_aipe_c.md).

## References

Chattopadhyay, B., Bandyopadhyay, T., Kelley, K., & Padalunkal, J. J.
(2025). A sequential approach for noninferiority or equivalence of a
linear contrast under cost constraints. *Psychological Methods, 30*(2),
425–439. [doi:10.1037/met0000570](https://doi.org/10.1037/met0000570)

Phillips, K. F. (1990). Power of the two one-sided tests procedure in
bioequivalence. *Journal of Pharmacokinetics and Biopharmaceutics,
18*(2), 139–144.
[doi:10.1007/BF01063556](https://doi.org/10.1007/BF01063556)

Schuirmann, D. J. (1987). A comparison of the two one-sided tests
procedure and the power approach for assessing the equivalence of
average bioavailability. *Journal of Pharmacokinetics and
Biopharmaceutics, 15*(6), 657–680.

## See also

[`power_equivalence_md`](https://yelleknek.github.io/DMAR/reference/power_equivalence_md.md),
[`ss_power_equivalence_c`](https://yelleknek.github.io/DMAR/reference/ss_power_equivalence_c.md),
[`equivalence_c`](https://yelleknek.github.io/DMAR/reference/equivalence_c.md),
[`ss_aipe_c`](https://yelleknek.github.io/DMAR/reference/ss_aipe_c.md)

Other equivalence testing:
[`equivalence_c()`](https://yelleknek.github.io/DMAR/reference/equivalence_c.md),
[`equivalence_r()`](https://yelleknek.github.io/DMAR/reference/equivalence_r.md),
[`equivalence_smd()`](https://yelleknek.github.io/DMAR/reference/equivalence_smd.md),
[`plot_equivalence()`](https://yelleknek.github.io/DMAR/reference/plot_equivalence.md),
[`power_density_equivalence_md()`](https://yelleknek.github.io/DMAR/reference/power_density_equivalence_md.md),
[`power_equivalence_md()`](https://yelleknek.github.io/DMAR/reference/power_equivalence_md.md),
[`power_equivalence_md_plot()`](https://yelleknek.github.io/DMAR/reference/power_equivalence_md_plot.md),
[`ss_power_equivalence_c()`](https://yelleknek.github.io/DMAR/reference/ss_power_equivalence_c.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# 1. Two groups of 61 and 113 sharing a five-group pooled error term
#    (so df_error = 404 - 5 = 399), bounds of 5 raw-scale points:
#    the design's probability of declaring equivalence when the
#    groups are truly identical.
power_equivalence_c(c_weights = c(1, -1), n = c(61, 113),
                    sigma = 15.67, delta_upper = 5,
                    true_psi = 0, df_error = 399)
#>  term  value
#>  power 0.281

# 2. The same design's noninferiority power at the same point.
power_equivalence_c(c_weights = c(1, -1), n = c(61, 113),
                    sigma = 15.67, delta_upper = 5,
                    true_psi = 0, df_error = 399,
                    side = "noninferiority")
#>  term  value
#>  power 0.641

# 3. Agreement with power_equivalence_md() in the two-group case
#    (Phillips, 1990, Table 1: expected 0.8029678).
power_equivalence_c(c_weights = c(1, -1), n = 24, sigma = 0.20,
                    delta_lower = 0.2, delta_upper = 0.2,
                    true_psi = 0.05, df_error = 22)
#>  term  value
#>  power 0.803
```
