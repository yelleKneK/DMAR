# Power of the TOST or Noninferiority Test for a Linear Contrast

Computes the exact power of the Schuirmann (1987) two one-sided tests
procedure, or of the one-sided noninferiority test, for a linear
contrast of group means \\\psi = \sum_j c_j \mu_j\\ with one pooled
error term, or for a contrast measured within each unit of a paired
design. For equivalence, the power is the probability that the \\(1 -
2\alpha)\\ confidence interval for \\\psi\\ lies entirely inside
\\(-\delta_L, \delta_U)\\, computed by numerical integration over the
chi distribution of the estimated error standard deviation; for
noninferiority, the power is a noncentral *t* probability in closed
form. This is the contrast generalization of
[`power_equivalence_md`](https://yelleknek.github.io/DMAR/reference/power_equivalence_md.md).

## Usage

``` r
power_equivalence_c(
  c_weights = NULL,
  n_per_group = NULL,
  N = NULL,
  true_sd,
  contrast_lower = NULL,
  contrast_upper = NULL,
  true_contrast = 0,
  alpha_level = 0.05,
  side = c("equivalence", "noninferiority"),
  df_error = NULL,
  design = c("independent", "paired")
)
```

## Arguments

- c_weights:

  The contrast weights, required for `design = "independent"`. The
  weights must sum to zero with the positive weights summing to 1 and
  the negative weights to -1, so that the bounds are on the raw scale of
  the response. Leave `NULL` for `design = "paired"`, where the weights
  are already in the per-unit contrast scores.

- n_per_group:

  For `design = "independent"`, the sample sizes per group (if length 1,
  equal group sizes are assumed); together with `c_weights`,
  `n_per_group` determines the standard error factor \\\sqrt{\sum_j
  c_j^2 / n_j}\\. Leave `NULL` for `design = "paired"`, which has no
  groups of different units and takes `N`; a paired call that gives
  `n_per_group` is refused.

- N:

  For `design = "paired"`, the number of units, each measured under
  every condition the contrast involves. Every unit contributes one
  contrast score, so this is the total sample size. Leave `NULL` for
  `design = "independent"`, where the sizes are given per group in
  `n_per_group`; an independent call that gives `N` is refused.

- true_sd:

  For `design = "independent"`, the error standard deviation (the square
  root of the mean square error). For `design = "paired"`, the standard
  deviation of the per-unit contrast scores (for two conditions, of the
  difference scores). It is \\\sigma\\ in the formulas below.

- contrast_lower, contrast_upper:

  Bounds on the raw scale of the response, each positive; the
  equivalence region is \\(-\delta_L, +\delta_U)\\, where \\\delta_L\\
  is `contrast_lower` and \\\delta_U\\ is `contrast_upper`. For
  `side = "equivalence"`, `contrast_upper` is required and
  `contrast_lower` defaults to it (symmetric bounds). For
  `side = "noninferiority"`, the test is against \\-\delta_L\\ alone, so
  `contrast_lower` is the argument to give; a call that gives only the
  upper bound uses it as the margin.

- true_contrast:

  The population value of the contrast at which the power is evaluated
  (\\\psi\\ in the formulas). Default `0`, the most favorable point for
  an equivalence declaration.

- alpha_level:

  One-sided significance level for each test. Default `0.05`.

- side:

  `"equivalence"` (default) for the TOST power, or `"noninferiority"`
  for the one-sided test against \\-\delta_L\\.

- df_error:

  The error degrees of freedom, \\\nu\\ in the formulas below. Defaults
  to \\\sum_j n_j - J\\ for `design = "independent"`, the sum of
  `n_per_group` less the number of groups, and to \\N - 1\\ for
  `design = "paired"`; supply it directly when the error term comes from
  a model with more groups or additional predictors than the contrast
  involves (for example, a five-group model supplying the pooled error
  for a two-group contrast).

- design:

  `"independent"` (default) for groups of different units sharing one
  pooled error term, or `"paired"` for a contrast measured within each
  unit (two conditions, or a pre and a post measurement, on the same
  units), analyzed as the mean of the \\N\\ per-unit contrast scores.
  The design decides which size argument is read: `n_per_group` for
  independent groups, `N` for a paired design. See Details.

## Value

A one-row `data.frame` with columns `term` (`"power"`) and `value` (the
computed power, in \\\[0, 1\]\\).

## Details

**Equivalence power.** Conditional on the estimated error standard
deviation \\S\\, the \\(1 - 2\alpha)\\ CI fits inside the bounds on a
computable event, and the unconditional power integrates that event over
the scaled chi distribution of \\S\\ on `df_error` degrees of freedom.
With `c_weights = c(1, -1)` and equal `n_per_group`, the result
reproduces
[`power_equivalence_md`](https://yelleknek.github.io/DMAR/reference/power_equivalence_md.md)
exactly.

**Noninferiority power.** The one-sided test rejects when \\t =
(\hat\psi + \delta_L)/\mathrm{SE}(\hat\psi)\\ exceeds
\\t\_{1-\alpha,\nu}\\, so the power is \\\Pr(T'\_{\nu}(\lambda) \>
t\_{1-\alpha,\nu})\\ with noncentrality \\\lambda = (\psi +
\delta_L)/(\sigma \sqrt{\sum_j c_j^2/n_j})\\.

**Paired designs.** With `design = "paired"`, each unit contributes one
contrast score, \\d_i = \sum_j c_j y\_{ij}\\, and the analysis is of
their mean: the standard error is \\\sigma / \sqrt{N}\\ on \\N - 1\\
degrees of freedom, where \\N\\ is the number of units (the argument
`N`), so \\\sqrt{\sum_j c_j^2 / n_j}\\ above becomes \\1 / \sqrt{N}\\
and both powers are computed as for the independent design. The weights
that form the scores should follow the same rule as `c_weights`
(positive weights summing to 1, negative to -1), so that the bounds are
on the raw scale; for two conditions the score is the difference. Its
standard deviation `true_sd` carries the correlation between the
measures: for two measures with a common standard deviation \\\sigma_y\\
correlated \\\rho\\, it is \\\sigma_y \sqrt{2(1 - \rho)}\\, smaller than
\\\sigma_y \sqrt{2}\\, the standard deviation of a difference between
independent observations, whenever \\\rho \> 0\\; that is the precision
a paired design buys.

**The feasibility condition.** If the expected half-width of the CI is
not smaller than the bounds allow, the equivalence power is zero or near
zero regardless of `true_contrast`: an imprecise design cannot declare
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
18*(2), 137–144.
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
[`equivalence_proportion_paired()`](https://yelleknek.github.io/DMAR/reference/equivalence_proportion_paired.md),
[`equivalence_r()`](https://yelleknek.github.io/DMAR/reference/equivalence_r.md),
[`equivalence_rmsea()`](https://yelleknek.github.io/DMAR/reference/equivalence_rmsea.md),
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
power_equivalence_c(c_weights = c(1, -1), n_per_group = c(61, 113),
                    true_sd = 15.67, contrast_upper = 5,
                    true_contrast = 0, df_error = 399)
#>  term  value
#>  power 0.281

# 2. The same design's noninferiority power at the same point.
power_equivalence_c(c_weights = c(1, -1), n_per_group = c(61, 113),
                    true_sd = 15.67, contrast_upper = 5,
                    true_contrast = 0, df_error = 399,
                    side = "noninferiority")
#>  term  value
#>  power 0.641

# 3. Agreement with power_equivalence_md() in the two-group case
#    (Phillips, 1990, Table 1: expected 0.8029678).
power_equivalence_c(c_weights = c(1, -1), n_per_group = 24, true_sd = 0.20,
                    contrast_lower = 0.2, contrast_upper = 0.2,
                    true_contrast = 0.05, df_error = 22)
#>  term  value
#>  power 0.803

# 4. A paired design: 50 units measured under both conditions, with
#    difference scores whose SD is 1, and a noninferiority margin of
#    0.5. No upper bound is needed for a noninferiority question.
power_equivalence_c(N = 50, true_sd = 1, contrast_lower = 0.5,
                    side = "noninferiority", design = "paired")
#>  term  value
#>  power 0.967
```
