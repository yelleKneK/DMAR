# Sample Size Planning From the AIPE Perspective for a Mean Under Two-Phase Sampling

Plans the number of labeled units, \\n\\, on which a costly target
measurement must be taken so that the confidence interval for the target
mean from
[`ci_two_phase`](https://yelleknek.github.io/DMAR/reference/ci_two_phase.md)
is sufficiently narrow, when an inexpensive proxy measurement is
available on every unit of a pool of \\M\\ units. The plan follows the
accuracy in parameter estimation (AIPE) approach: the goal is an
interval no wider than `width`, either at the planning values (the
expected width) or with a stated assurance that the realized interval
will be that narrow (Kelley & Rausch, 2006). Beside the answer it
reports the number of labeled units the same goal needs without the
proxy, so the saving is in view.

## Usage

``` r
ss_aipe_two_phase(
  true_sd_Y,
  true_cor,
  width,
  n_pool = Inf,
  conf_level = 0.95,
  assurance = NULL,
  true_kurtosis = 0
)
```

## Arguments

- true_sd_Y:

  Planning value of the standard deviation of the target measurement,
  \\\sigma_Y\\.

- true_cor:

  Planning value of the correlation between the target and the proxy,
  \\\rho\\.

- width:

  The desired full width of the confidence interval (the upper limit
  minus the lower limit), in the units of the target.

- n_pool:

  The number of units in the pool, \\M\\, every one of which carries the
  proxy. The default `Inf` plans as though the pool were unlimited.

- conf_level:

  Confidence level (default `0.95`).

- assurance:

  `NULL` (default) to plan for the expected width, or a probability in
  (0, 1) that the realized width will be no wider than `width`, as in
  [`ss_aipe_c`](https://yelleknek.github.io/DMAR/reference/ss_aipe_c.md).

- true_kurtosis:

  Excess kurtosis (the scale of
  [`kurtosis`](https://yelleknek.github.io/DMAR/reference/kurtosis.md))
  of the residual of the target on the proxy; default `0`, its normal
  value. Used only with `assurance`, where it adjusts the degrees of
  freedom of the variance estimate (see Details).

## Value

A `data.frame` (classes `dmar_ss_aipe` and `dmar_tbl`) with rows

- necessary_n_labeled:

  the planned number of labeled units.

- necessary_n_target_only:

  the number the same goal needs without the proxy.

- ratio:

  the first divided by the second.

- expected_width:

  the full width at `necessary_n_labeled` and the planning values.

- minimum_width:

  the narrowest width the pool allows, the normal-quantile width with
  all but one unit of the pool labeled (0 for an unlimited pool; see
  Details).

- width, true_sd_Y, true_cor, n_pool, conf_level:

  the inputs.

- assurance, true_kurtosis:

  the inputs, present only when an assurance was supplied.

## Details

**The interval planned for.** The regression estimator of double
sampling (Tamhane, 1978, Section 2.1, p. 420; the design is that of
Neyman, 1938, who formed strata from the inexpensive variable) that
[`ci_two_phase`](https://yelleknek.github.io/DMAR/reference/ci_two_phase.md)
reports has, to first order, variance \\\sigma_Y^2(1 - \rho^2)/n +
\rho^2\sigma_Y^2/M\\, and its interval uses the *t* distribution on
\\n - 2\\ degrees of freedom. With \\\omega\\ the desired width
(`width`), \\1 - \alpha\\ the `conf_level`, and \\t_n =
t\_{1-\alpha/2,\\ n-2}\\, the plan for the expected width is the
smallest \\n\\ with \$\$2\\ t_n \left\\\frac{\sigma_Y^2(1 -
\rho^2)}{n} + \frac{\rho^2\sigma_Y^2}{M}\right\\^{1/2} \le \omega.\$\$
In an evaluation of a system whose outputs an automated judge labels,
for example, \\n\\ is the number of outputs a human assessor must label
when the judge has labeled \\M\\ of them.

**Assurance.** The realized width replaces \\\sigma_Y^2(1 - \rho^2)\\ by
the residual mean square \\S_e^2\\ and \\\rho^2 \sigma_Y^2\\ by \\\hat
b^2 S_X^2\\. Treating the second, which is divided by the large \\M\\,
as fixed at its planning value, and \\\nu S_e^2 / \sigma_e^2\\ as chi
square on \\\nu\\ degrees of freedom, the probability that the realized
width is at most \\\omega\\ is \$\$F\_{\chi^2\_\nu}\left(\frac{\nu\\ n\\
\\\omega^2 / (4 t_n^2) - \rho^2 \sigma_Y^2 / M\\}{\sigma_Y^2 (1 -
\rho^2)}\right),\$\$ and the plan is the smallest \\n\\ from which it is
at least `assurance` at every larger \\n\\. Under normal residuals \\\nu
= n - 2\\. For residuals with excess kurtosis \\\kappa\\
(`true_kurtosis`) the variance of \\S_e^2\\ is approximately
\\\sigma_e^4\\2/(n - 2) + \kappa/n\\\\, and matching it to a scaled chi
square gives \\\nu = 2 / \\2/(n - 2) + \kappa/n\\\\: for an assurance of
.80 or more, a lighter-tailed residual (\\\kappa \< 0\\) needs fewer
units and a heavier-tailed one more. A binary target is not light-tailed
as a rule: a 0/1 variable with proportion \\p\\ has excess kurtosis
\\\\1 - 6p(1 - p)\\/\\p(1 - p)\\\\, which is negative only for \\p\\
between about .21 and .79 and is about 5.1 at \\p = .90\\. The planning
value is best taken from the residuals of a pilot regression of the
target on the proxy, through
[`kurtosis`](https://yelleknek.github.io/DMAR/reference/kurtosis.md).
With a large \\\kappa\\ and a small \\n\\, \\\nu\\ is small, and the
two-moment approximation can reach the assurance at a few small \\n\\,
fall below it over a stretch of larger ones, and reach it again; the
plan is the start of the stretch over which it holds for good, which the
search finds exactly for any `assurance` in (0, 1). Under normal
residuals and with an assurance of one half or more the two readings
agree.

**The pool sets a floor.** The labeled units are drawn from the pool,
and
[`ci_two_phase`](https://yelleknek.github.io/DMAR/reference/ci_two_phase.md)
needs at least one unit that carries the proxy alone, so \\n \le M -
1\\. No number of labeled units can therefore bring the width below \\2
z\_{1-\alpha/2} \\\sigma_Y^2(1 - \rho^2)/(M - 1) +
\rho^2\sigma_Y^2/M\\^{1/2}\\ (the `minimum_width` row), which is close
to \\2 z\_{1-\alpha/2}\\ \sigma_Y / \sqrt{M}\\, the width of the target
mean were every unit of the pool labeled; a target at or below it is an
error. A target just above it can still be out of reach once the *t*
quantile and the assurance are accounted for; when the plan needs more
than \\M - 1\\ labeled units, that too is an error, naming the pool as
the constraint, rather than a missing value.

**The target-only comparison.** `necessary_n_target_only` is the planned
size of the usual *t* interval for the target mean from labeled units
alone, the plan
[`ss_aipe_mean()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mean.md)
makes for a single mean, so the two functions agree. Without an
assurance it is the smallest \\n\\ with \\2\\ t\_{1-\alpha/2,\\ n-1}\\
c_4(n)\\ \sigma_Y / \sqrt{n} \le \omega\\, where \\c_4(n)\sigma_Y\\ is
the expected sample standard deviation under normality, so that the left
side is the exact expected width; with `assurance` it is the same chi
square argument on \\n - 1\\ degrees of freedom. It is not bounded by
`n_pool`, because a design without the proxy has no use for a pool. The
`true_kurtosis` supplied is the residual's, and the comparison uses it
for the target as well, an approximation: the target's own kurtosis can
differ, and the comparison for it is
[`ss_aipe_mean()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mean.md)
with that value. The labeled-unit plan, by contrast, sets the two
estimates in its standard error at their planning values, since the
expected value of the two-phase standard error has no closed form once
the pool term enters; this overstates its expected width slightly, by a
factor of at most about \\1 + 1/\\4(n - 2)\\\\, and can add a unit to
`necessary_n_labeled`.

## References

Kelley, K., & Rausch, J. R. (2006). Sample size planning for the
standardized mean difference: Accuracy in parameter estimation via
narrow confidence intervals. *Psychological Methods, 11*(4), 363–385.
[doi:10.1037/1082-989X.11.4.363](https://doi.org/10.1037/1082-989X.11.4.363)

Neyman, J. (1938). Contribution to the theory of sampling human
populations. *Journal of the American Statistical Association, 33*(201),
101–116.
[doi:10.1080/01621459.1938.10503378](https://doi.org/10.1080/01621459.1938.10503378)

Tamhane, A. C. (1978). Inference based on regression estimator in double
sampling. *Biometrika, 65*(2), 419–427.
[doi:10.1093/biomet/65.2.419](https://doi.org/10.1093/biomet/65.2.419)

## See also

[`ci_two_phase`](https://yelleknek.github.io/DMAR/reference/ci_two_phase.md)
for the interval,
[`ss_aipe_two_phase_cost`](https://yelleknek.github.io/DMAR/reference/ss_aipe_two_phase_cost.md)
for the least costly numbers of labeled and pool units,
[`ss_aipe_two_phase_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_two_phase_sensitivity.md)
to evaluate a plan by simulation, and
[`ss_seq_two_phase`](https://yelleknek.github.io/DMAR/reference/ss_seq_two_phase.md)
for the sequential rule.

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
[`ss_aipe_partial_r()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_partial_r.md),
[`ss_aipe_partial_r_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_partial_r_sensitivity.md),
[`ss_aipe_pcm_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_pcm_sensitivity.md),
[`ss_aipe_r()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_r.md),
[`ss_aipe_r_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_r_sensitivity.md),
[`ss_aipe_reliability_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reliability_sensitivity.md),
[`ss_aipe_semipartial_r()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_semipartial_r.md),
[`ss_aipe_semipartial_r_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_semipartial_r_sensitivity.md),
[`ss_aipe_two_phase_cost()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_two_phase_cost.md),
[`ss_aipe_two_phase_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_two_phase_sensitivity.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Planning values from a pilot: a target standard deviation of 0.43 and
# a correlation of 0.55 with the proxy, which is available on a pool of
# 20,000 units. A 90% interval no wider than 0.10:
ss_aipe_two_phase(true_sd_Y = 0.43, true_cor = 0.55, width = 0.10,
                  n_pool = 20000,
                  conf_level = 0.90)
#>  term                    value
#>  necessary_n_labeled     142  
#>  necessary_n_target_only 202  
#>  ratio                   0.703
#>  expected_width          0.100
#>  minimum_width           0.010
#>  width                   0.100
#>  true_sd_Y               0.430
#>  true_cor                0.550
#>  n_pool                  20000
#>  conf_level              0.900
#> 
#> Confidence level: 90%

# The same goal with 80% assurance that the realized interval is no
# wider than 0.10:
ss_aipe_two_phase(true_sd_Y = 0.43, true_cor = 0.55, width = 0.10,
                  n_pool = 20000,
                  conf_level = 0.90, assurance = 0.80)
#>  term                    value
#>  necessary_n_labeled     156  
#>  necessary_n_target_only 218  
#>  ratio                   0.716
#>  expected_width          0.095
#>  minimum_width           0.010
#>  width                   0.100
#>  true_sd_Y               0.430
#>  true_cor                0.550
#>  n_pool                  20000
#>  conf_level              0.900
#>  assurance               0.800
#>  true_kurtosis           0.000
#> 
#> Confidence level: 90%

# A smaller pool raises the floor and the number of labeled units.
ss_aipe_two_phase(true_sd_Y = 0.43, true_cor = 0.55, width = 0.10,
                  n_pool = 1000,
                  conf_level = 0.90, assurance = 0.80)
#>  term                    value
#>  necessary_n_labeled     165  
#>  necessary_n_target_only 218  
#>  ratio                   0.757
#>  expected_width          0.096
#>  minimum_width           0.045
#>  width                   0.100
#>  true_sd_Y               0.430
#>  true_cor                0.550
#>  n_pool                  1000 
#>  conf_level              0.900
#>  assurance               0.800
#>  true_kurtosis           0.000
#> 
#> Confidence level: 90%
```
