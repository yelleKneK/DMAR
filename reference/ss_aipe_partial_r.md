# Sample Size for AIPE on a Partial Correlation

Takes a planning value of the population partial correlation \\\rho\_{XY
\cdot Z_1 \cdots Z\_{p-1}}\\, the number of predictors (*X* and the
variables partialed out), and a desired confidence interval width, and
returns the smallest sample size at which the interval is expected to be
that narrow (accuracy in parameter estimation, AIPE; Kelley, 2008). With
`assurance` it returns the smallest sample size at which the realized
interval is that narrow with at least the stated probability.

## Usage

``` r
ss_aipe_partial_r(
  true_cor,
  p,
  width,
  which_width = c("Full", "Lower", "Upper"),
  conf_level = 0.95,
  fisher_z = FALSE,
  assurance = NULL
)
```

## Arguments

- true_cor:

  Population partial correlation \\\rho\_{XY \cdot Z_1 \cdots
  Z\_{p-1}}\\ (written \\\rho\\ in Details) the sample size is planned
  for, in \\(-1, 1)\\. The planner treats this value as the truth;
  [`ss_aipe_partial_r_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_partial_r_sensitivity.md)
  shows what the plan delivers when the truth (`true_cor` there) differs
  from the value planned with (`supposed_cor` there).

- p:

  Number of predictors in the regression of *Y* on *X* and the variables
  partialed out, so `p - 1` variables are partialed out and `p` is at
  least 2. A categorical variable with *a* levels counts as *a* - 1
  predictors. It is the `p` of
  [`ss_power_R2`](https://yelleknek.github.io/DMAR/reference/ss_power_R2.md),
  [`ci_R2`](https://yelleknek.github.io/DMAR/reference/ci_R2.md), and
  [`ci_reg_coef`](https://yelleknek.github.io/DMAR/reference/ci_reg_coef.md):
  the correlation concerns one predictor of a regression with `p`
  predictors. In the notation of a model comparison, \\q\\ is the number
  of predictors removed from the full model to give the restricted model
  (\\q = p\_{full} - p\_{restricted}\\, the predictors being tested),
  and the number of variables partialed out is \\p - q\\; this
  correlation concerns one predictor, so \\q = 1\\ and \\p - 1\\
  variables are partialed out.

- width:

  Desired full width of the confidence interval on the partial
  correlation.

- which_width:

  Whether `width` is the `"Full"` width of the interval (default) or a
  half-width: `"Lower"` and `"Upper"` both interpret `width` as half the
  full width, so they plan for a full width of twice `width` and return
  the same sample size. Because the interval is not generally symmetric
  about the estimate (with `fisher_z = TRUE` in particular), its
  realized lower and upper half-widths can differ from each other and
  from half the full width; the planner does not target them separately.
  A genuinely one-sided width target is not currently offered.

- conf_level:

  Desired confidence level (default `0.95`).

- fisher_z:

  Logical. If `TRUE`, the plan is for the interval formed on Fisher's
  *Z* scale (variance \\1/(N - p - 2)\\) and back-transformed with
  \\\tanh\\. If `FALSE` (default), the plan is for the interval formed
  on the raw scale from the first-order variance.

- assurance:

  Optional. The probability that the realized interval is no wider than
  `width`. The function returns the smallest sample size, at or above
  the sample size planned without it, at which that probability is at
  least `assurance`. A single number above 0.5 and below 1; when `NULL`
  (the default) the plan is for the expected width only.

## Value

A `data.frame` with the columns `term` and `value` and one row per term:

- `necessary_N`:

  the planned total sample size.

- `expected_width`:

  the full width of the interval at `necessary_N` when the sample
  partial correlation equals `true_cor`.

- `true_cor`, `p`, `width`, `conf_level`:

  the inputs, with `width` the `width` as given.

- `assurance`:

  the requested probability; present only when `assurance` was given.

- `assurance_attained`:

  the probability, at `necessary_N` and `true_cor`, that the realized
  interval is no wider than the target; present only when `assurance`
  was given. It is at least `assurance` and exceeds it because the
  sample size is a whole number.

## Details

**Estimand and interval.** The quantity planned for is the population
partial correlation between \\X\\ and \\Y\\ given \\Z_1, \ldots,
Z\_{p-1}\\, for independent observations from a multivariate normal
population. The interval is two-sided at `conf_level` and is one of two
Wald intervals, chosen by `fisher_z`; the plan is for the one chosen,
and a study that reports the other is not the study that was planned.
Below, \\\rho\\ is that partial correlation (`true_cor`) and \\1 -
\alpha\\ is `conf_level`.

**Raw scale, first-order variance (`fisher_z = FALSE`).** The
first-order variance of a simple correlation is \\(1 - \rho^2)^2 / N\\
(Olkin & Finn, 1995, Equation 3). The sample partial correlation given
\\p - 1\\ variables is distributed as a simple correlation from \\N -
(p - 1)\\ observations (Fisher, 1924, p. 330, for normally distributed
variables), and the functions here use \\(1 - \rho^2)^2 / (N - p)\\,
which differs from that at order \\1/N^2\\ (see
[`var_partial_r`](https://yelleknek.github.io/DMAR/reference/var_partial_r.md)).
The interval is \\r \pm z\_{1 - \alpha/2} (1 - r^2) / \sqrt{N - p}\\,
with half-width \$\$w\_{1/2}(r, N) \\=\\ z\_{1 - \alpha/2}\\ \frac{1 -
r^2}{\sqrt{N - p}}.\$\$ Setting \\r = \rho\\ and solving for the
smallest whole \\N\\ with \\w\_{1/2} \le\\ the target half-width gives
\$\$N \\=\\ p + \Big\lceil z\_{1 - \alpha/2}^2 \\ (1 - \rho^2)^2 /
w\_{1/2}^{2} \Big\rceil.\$\$

**Fisher's *Z* scale (`fisher_z = TRUE`).** Bonett and Wright (2000,
Section 3, pp. 24–25) plan the sample size for a simple correlation so
that the interval formed on Fisher's *Z* scale (Fisher, 1921, 1924) and
back-transformed has the desired width. The same plan is applied here to
the partial correlation, with the standard error \\1/\sqrt{N - p - 2}\\
in place of \\1/\sqrt{N - 3}\\. With \\Z = \mathrm{arctanh}(r)\\ and
standard error \\1/\sqrt{N - p - 2}\\, the interval is \\\tanh(Z \pm
z\_{1 - \alpha/2} / \sqrt{N - p - 2})\\, and half its width is
\$\$w\_{1/2}(r, N) \\=\\ \tfrac{1}{2}\Big\\
\tanh\\\Big(\mathrm{arctanh}(r) + \frac{z\_{1 - \alpha/2}}{\sqrt{N - p -
2}}\Big) - \tanh\\\Big(\mathrm{arctanh}(r) - \frac{z\_{1 -
\alpha/2}}{\sqrt{N - p - 2}}\Big) \Big\\.\$\$ The planner finds the
smallest \\N\\ of at least \\p + 3\\ at which this is no larger than the
target at \\r = \rho\\, by a one-dimensional search (the half-width
decreases in \\N\\). This interval stays inside \\(-1, 1)\\ and is
generally the more accurate of the two when \\N - p\\ is small or
\\\|\rho\|\\ is large.

**The assurance step.** A plan without `assurance` makes the interval as
narrow as the target when the sample partial correlation equals
`true_cor`; in about half of all samples it is wider. With `assurance`
the function finds the sample size at which the realized interval is no
wider than the target with at least that probability, from the sampling
distribution of the sample partial correlation:

1.  At a given \\N\\ the width of either interval depends on the sample
    only through \\r\\, and both half-widths above decrease as \\\|r\|\\
    grows. The realized interval is therefore no wider than the target
    exactly when \\\|r\| \ge r^\*(N)\\, where \\r^\*(N)\\ is the
    \\\|r\|\\ at which the half-width equals the target half-width
    \\h\\: on the raw scale \\r^\*(N) = \sqrt{1 - h \sqrt{N - p} /
    z\_{1 - \alpha/2}}\\, and on Fisher's scale the root of
    \\w\_{1/2}(r, N) = h\\. When the half-width at \\r = 0\\, the widest
    the interval can be, already meets the target, every sample does and
    the probability is 1.

2.  Under multivariate normality the sample partial correlation given
    \\p - 1\\ variables from \\N\\ observations is distributed as a
    simple correlation from \\m = N - (p - 1)\\ observations. Its
    distribution function is computed as \$\$F(c) \\=\\ \int_0^1
    T\_{2m - 3}\\\Big(\sqrt{2m - 3}\\ \Big\\\frac{c\\\sqrt{1 -
    b}}{\sqrt{1 - c^2}} - \frac{\rho\\\sqrt{b}}{\sqrt{1 -
    \rho^2}}\Big\\\Big)\\ g(b)\\ db,\$\$ where \\T\_{2m - 3}\\ is the
    distribution function of Student's *t* on \\2m - 3\\ degrees of
    freedom and \\g\\ is the density of the beta distribution with
    parameters \\(m - 1)/2\\ and \\(m - 2)/2\\. The form follows from
    writing \\r / \sqrt{1 - r^2}\\ as \\(Z + \theta \sqrt{U}) /
    \sqrt{V}\\, with \\Z\\ standard normal, \\U\\ and \\V\\ chi square
    on \\m - 1\\ and \\m - 2\\ degrees of freedom, all independent, and
    \\\theta = \rho / \sqrt{1 - \rho^2}\\, and then conditioning on \\U
    / (U + V)\\. It uses the central *t* only, because R's noncentral
    *t* is an approximation at the noncentrality values that plans for
    narrow intervals reach.

3.  The probability at \\N\\ is \\P_N = 1 - F(r^\*(N)) + F(-r^\*(N))\\.

4.  The function returns the smallest \\N\\, at or above the sample size
    planned without `assurance`, with \\P_N \ge\\ `assurance`, found by
    doubling a step upward and then bisecting. `assurance_attained` is
    \\P_N\\ at that \\N\\.

The probability is exact, up to numerical integration, under
multivariate normality for the interval that is planned. It was checked
by simulation in two ways: in 112 planning sets of 200,000 simulated
correlations each, the proportion of intervals no wider than the target
at the planned sample size was within .0023 of \\P_N\\; and in 16
further planning sets simulated from regression data (20,000 data sets
each, with the variables partialed out by least squares) it was within
.004 of \\P_N\\, about two Monte Carlo standard errors, and between .002
and .012 above the request.

Four properties of the step are worth knowing.

- The answer is never below the sample size planned without `assurance`.
  Below that sample size \\P_N\\ is not monotone in \\N\\ (at a very
  small \\N\\ the sample correlation is variable enough to exceed
  \\r^\*(N)\\ by chance, and does so less often as \\N\\ grows, before
  the narrowing of the interval takes over), so a smaller \\N\\ that
  happens to meet the request is not a plan. At or above it \\P_N\\ did
  not decrease in any of 72 planning sets scanned over every whole
  \\N\\.

- When `true_cor` is 0 the plan without `assurance` already makes the
  interval as narrow as the target in every sample, since the interval
  is widest at \\r = 0\\; `assurance` then adds nothing and
  `assurance_attained` is 1. Near 0 the same holds approximately and the
  step adds little.

- As \\\|\rho\|\\ grows the width becomes more variable from sample to
  sample and the step adds more: for a width of .10 with one variable
  partialed out and an assurance of .80, the plan moves from 867 to 915
  at \\\rho = .50\\ and from 121 to 151 at \\\rho = .85\\.

- `assurance` must be above 0.5 and below 1; other values are refused.

**Why the step is not a chi square factor.** Multiplying the sample size
planned without `assurance` by a chi square quantile over its degrees of
freedom is the assurance step for an interval whose width is
proportional to a sample standard deviation (Kelley & Maxwell, 2003);
the width of an interval for a correlation depends on the sample through
\\r\\ instead, and the factor gives too few cases when \\\|\rho\|\\ is
large and too many when it is small. For one variable partialed out (\\p
= 2\\), a width of .10, and a requested .80 it gives 426 at \\\rho =
.70\\, where the probability is .675 (448 are needed), and 1318 at
\\\rho = .30\\, where 1310 suffice; in the 16 planning sets simulated
from regression data its sample sizes delivered from .14 below the
request to .15 above it.

**When to use partial or simple correlation planning.** Use this
function when the target is the population correlation between \\X\\ and
\\Y\\ after statistically controlling for \\Z_1, \ldots, Z\_{p-1}\\. For
the simple Pearson correlation, see
[`ss_aipe_r`](https://yelleknek.github.io/DMAR/reference/ss_aipe_r.md).
[`ss_aipe_partial_r_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_partial_r_sensitivity.md)
simulates the plan, including under a population value other than the
planning value.

## References

Algina, J., & Olejnik, S. (2003). Sample size tables for correlation
analysis with applications in partial correlation and multiple
regression analysis. *Multivariate Behavioral Research, 38*(3), 309–323.
[doi:10.1207/s15327906mbr3803_02](https://doi.org/10.1207/s15327906mbr3803_02)

Bonett, D. G., & Wright, T. A. (2000). Sample size requirements for
estimating Pearson, Kendall and Spearman correlations. *Psychometrika,
65*(1), 23–28.
[doi:10.1007/BF02294183](https://doi.org/10.1007/BF02294183)

Fisher, R. A. (1921). On the "probable error" of a coefficient of
correlation deduced from a small sample. *Metron, 1*, 3–32.

Fisher, R. A. (1924). The distribution of the partial correlation
coefficient. *Metron, 3*, 329–332.

Kelley, K. (2008). Sample size planning for the squared multiple
correlation coefficient: Accuracy in parameter estimation via narrow
confidence intervals. *Multivariate Behavioral Research, 43*(4),
524–555.
[doi:10.1080/00273170802490632](https://doi.org/10.1080/00273170802490632)

Kelley, K., & Maxwell, S. E. (2003). Sample size for multiple
regression: Obtaining regression coefficients that are accurate, not
simply significant. *Psychological Methods, 8*(3), 305–321.
[doi:10.1037/1082-989X.8.3.305](https://doi.org/10.1037/1082-989X.8.3.305)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 3 on the one-way ANOVA and Chapter 4 on
contrasts.)

Olkin, I., & Finn, J. D. (1995). Correlations redux. *Psychological
Bulletin, 118*(1), 155–164.
[doi:10.1037/0033-2909.118.1.155](https://doi.org/10.1037/0033-2909.118.1.155)

## See also

[`var_partial_r`](https://yelleknek.github.io/DMAR/reference/var_partial_r.md),
[`expected_partial_r`](https://yelleknek.github.io/DMAR/reference/expected_partial_r.md),
[`ss_aipe_partial_r_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_partial_r_sensitivity.md),
[`ss_aipe_semipartial_r`](https://yelleknek.github.io/DMAR/reference/ss_aipe_semipartial_r.md),
[`ss_aipe_r`](https://yelleknek.github.io/DMAR/reference/ss_aipe_r.md)

[`design_consequences`](https://yelleknek.github.io/DMAR/reference/design_consequences.md)
for what a chosen design delivers: power, the Type S (sign) and Type M
(exaggeration) errors of the significance filter, and the expected
confidence interval width.

Other AIPE sample size planning:
[`ss_aipe_c_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_c_sensitivity.md),
[`ss_aipe_cliff_delta()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cliff_delta.md),
[`ss_aipe_cliff_delta_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cliff_delta_sensitivity.md),
[`ss_aipe_composite_sem()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_composite_sem.md),
[`ss_aipe_crossed()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_crossed.md),
[`ss_aipe_crossed_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_crossed_sensitivity.md),
[`ss_aipe_equivalence_r()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_equivalence_r.md),
[`ss_aipe_equivalence_r_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_equivalence_r_sensitivity.md),
[`ss_aipe_equivalence_smd()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_equivalence_smd.md),
[`ss_aipe_equivalence_smd_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_equivalence_smd_sensitivity.md),
[`ss_aipe_icc()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_icc.md),
[`ss_aipe_icc_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_icc_sensitivity.md),
[`ss_aipe_indirect_effect()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_indirect_effect.md),
[`ss_aipe_indirect_effect_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_indirect_effect_sensitivity.md),
[`ss_aipe_mean()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mean.md),
[`ss_aipe_mean_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mean_sensitivity.md),
[`ss_aipe_mixed_effects_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mixed_effects_sensitivity.md),
[`ss_aipe_omega_squared()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_omega_squared.md),
[`ss_aipe_omega_squared_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_omega_squared_sensitivity.md),
[`ss_aipe_partial_r_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_partial_r_sensitivity.md),
[`ss_aipe_pcm_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_pcm_sensitivity.md),
[`ss_aipe_r()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_r.md),
[`ss_aipe_r_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_r_sensitivity.md),
[`ss_aipe_reliability_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reliability_sensitivity.md),
[`ss_aipe_semipartial_r()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_semipartial_r.md),
[`ss_aipe_semipartial_r_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_semipartial_r_sensitivity.md),
[`ss_aipe_two_phase()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_two_phase.md),
[`ss_aipe_two_phase_cost()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_two_phase_cost.md),
[`ss_aipe_two_phase_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_two_phase_sensitivity.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# 1. Plan N so the 95% CI on rho_XY.Z (two variables partialed out, so
#        p = 3 predictors) has full width <= 0.20, when the planning
#        value of the partial correlation is 0.30.
ss_aipe_partial_r(true_cor = 0.30, p = 3, width = 0.20)
#>  term           value
#>  necessary_N    322  
#>  expected_width 0.200
#>  true_cor       0.300
#>  p              3    
#>  width          0.200
#>  conf_level     0.950
#> 
#> Confidence level: 95%

# 2. The same problem for the interval formed on Fisher's Z scale
#        (Bonett & Wright, 2000).
ss_aipe_partial_r(true_cor = 0.30, p = 3, width = 0.20, fisher_z = TRUE)
#>  term           value
#>  necessary_N    322  
#>  expected_width 0.200
#>  true_cor       0.300
#>  p              3    
#>  width          0.200
#>  conf_level     0.950
#> 
#> Confidence level: 95%

# 3. With an assurance of .80: the smallest N at which the realized
#        interval is no wider than 0.20 with probability at least .80.
#        The last row is the probability at the returned N.
ss_aipe_partial_r(true_cor = 0.30, p = 3, width = 0.20, assurance = 0.80)
#>  term               value
#>  necessary_N        338  
#>  expected_width     0.195
#>  true_cor           0.300
#>  p                  3    
#>  width              0.200
#>  conf_level         0.950
#>  assurance          0.800
#>  assurance_attained 0.804
#> 
#> Confidence level: 95%

# 4. A larger correlation needs fewer cases for the same width, and
#        the assurance step adds proportionally more of them.
ss_aipe_partial_r(true_cor = 0.70, p = 3, width = 0.20)
#>  term           value
#>  necessary_N    103  
#>  expected_width 0.200
#>  true_cor       0.700
#>  p              3    
#>  width          0.200
#>  conf_level     0.950
#> 
#> Confidence level: 95%
ss_aipe_partial_r(true_cor = 0.70, p = 3, width = 0.20, assurance = 0.80)
#>  term               value
#>  necessary_N        126  
#>  expected_width     0.180
#>  true_cor           0.700
#>  p                  3    
#>  width              0.200
#>  conf_level         0.950
#>  assurance          0.800
#>  assurance_attained 0.809
#> 
#> Confidence level: 95%
```
