# Sample Size Planning for a Mean From the Accuracy in Parameter Estimation (AIPE) Perspective

Plans the sample size *N* needed for the confidence interval for a
population mean to be sufficiently narrow, from the accuracy in
parameter estimation (AIPE) perspective. The mean may be that of a
single sample or the mean of paired differences: in a pretest-posttest
or matched pairs design each unit contributes one difference score, the
interval for the mean difference is the single sample interval computed
on those scores, and the plan is the same with `true_sd` the standard
deviation of the differences. Without an `assurance` the plan makes the
expected width of the interval no wider than `width`; with one, it makes
the realized width no wider than `width` with that probability.

## Usage

``` r
ss_aipe_mean(
  true_sd,
  width,
  conf_level = 0.95,
  assurance = NULL,
  true_kurtosis = 0
)
```

## Arguments

- true_sd:

  The planning value of the population standard deviation of the scores
  whose mean is estimated. For a paired design it is the standard
  deviation of the differences, which for two measures with standard
  deviations \\\sigma_1\\ and \\\sigma_2\\ and correlation \\\rho\\ is
  \\\sqrt{\sigma_1^2 + \sigma_2^2 - 2\rho\\\sigma_1\sigma_2}\\. Must be
  positive.

- width:

  The desired full width of the confidence interval (twice its
  half-width), in the units of the scores.

- conf_level:

  The confidence level of the interval (default `0.95`).

- assurance:

  `NULL` (the default) to plan the expected width, or the desired
  probability, strictly between 0.5 and 1, that the realized width will
  be no wider than `width`.

- true_kurtosis:

  The planning value of the population excess kurtosis of the scores, on
  the scale
  [`kurtosis`](https://yelleknek.github.io/DMAR/reference/kurtosis.md)
  reports (0 for a normal population). It enters only the assurance
  adjustment (see Details). The default 0 is normal theory. The value
  must be greater than -2: no distribution has an excess kurtosis below
  -2, and only a population of two equally likely values reaches it. A
  pilot's
  [`kurtosis()`](https://yelleknek.github.io/DMAR/reference/kurtosis.md)
  can usually be passed as it is, but
  [`kurtosis()`](https://yelleknek.github.io/DMAR/reference/kurtosis.md)
  corrects its estimate for small-sample bias, and for a pilot whose
  scores take nearly two values it can return -2 or less
  (`kurtosis(rep(c(0, 1), 123))` is about -2.016), which is refused.
  Such a pilot says the population is close to that two-valued limit.
  Pass instead a value a little above -2, such as -1.9, or 0, the normal
  theory plan, which for an assurance of .80 or more is never the
  smaller of the two.

## Value

A `data.frame` of class `dmar_ss_aipe` with columns `term` and `value`
and the rows `necessary_N` (the planned sample size), `expected_width`
(the expected width of the interval at that size under normality), the
input echoes `true_sd`, `width`, and `conf_level`, and, when an
assurance was supplied, `assurance` and `true_kurtosis`.

## Details

**The interval.** The \\100(1 - \alpha)\\% interval for the mean is
\\\bar{y} \pm t\_{1-\alpha/2,\\N-1}\\ S/\sqrt{N}\\, the interval
[`t.test()`](https://rdrr.io/r/stats/t.test.html) reports for a single
sample and for paired data. Its full width \\W = 2\\
t\_{1-\alpha/2,\\N-1}\\ S/\sqrt{N}\\ is random because the sample
standard deviation \\S\\ is, so a plan is a statement either about the
expected width or about the probability that \\W\\ is no larger than the
desired width \\\omega\\. These are the two criteria Kelley and Rausch
(2006) developed for the standardized mean difference and that Kelley,
Maxwell, and Rausch (2003) developed for the unstandardized mean
difference.

**Expected width.** Without an assurance, `necessary_N` is the smallest
*N* at which \$\$E(W) = 2\\ t\_{1-\alpha/2,\\N-1}\\ c_4(N)\\ \sigma /
\sqrt{N} \le \omega,\$\$ where, under normality, \$\$c_4(N) =
\frac{E(S)}{\sigma} = \sqrt{\frac{2}{N - 1}}\\
\frac{\Gamma(N/2)}{\Gamma((N - 1)/2)}.\$\$ The expected width is
therefore exact under normality. It falls as *N* rises, so the
criterion, once met, holds at every larger *N*. The `expected_width` row
reports it at `necessary_N`.

**How the plan differs from
[`ss_aipe_c()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_c.md).**
[`ss_aipe_c()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_c.md),
which plans a contrast among independent groups, sets \\S\\ equal to
\\\sigma\\ in the width. Because \\E(S) = c_4(N)\\\sigma \< \sigma\\,
the width it plans is larger than the expected width, and the same rule
applied to a single mean gives either the *N* returned here or one more.
With an assurance the difference disappears: the criterion is then the
probability that \\S\\ is small enough, which does not involve \\E(S)\\,
and under normality it is the criterion
[`ss_aipe_c()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_c.md)
uses for a contrast, here applied to a single mean.
([`ss_aipe_c()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_c.md)
cannot plan a single mean, since it requires contrast weights that sum
to zero.)

**Assurance.** With assurance \\\gamma\\, `necessary_N` is the smallest
*N* from which \\P(W \le \omega) \ge \gamma\\ holds at every larger *N*.
Because \\W \le \omega\\ exactly when \$\$\frac{S^2}{\sigma^2} \le
\frac{N \omega^2}{4\\ t^2\_{1-\alpha/2,\\N-1}\\ \sigma^2},\$\$ the
probability is that of \\S^2\\ falling below a bound. Under normality
\\(N - 1) S^2/\sigma^2\\ is chi square on \\N - 1\\ degrees of freedom,
so the probability is exact, and the criterion, once met, holds at every
larger *N*; `necessary_N` is then simply the smallest *N* at which \\P(W
\le \omega) \ge \gamma\\. The two readings part only under the kurtosis
adjustment (see below).

**Kurtosis.** For *N* independent scores whose excess kurtosis is
\\\kappa\\, the variance of the sample variance is \$\$\mathrm{Var}(S^2)
= \sigma^4 \left\\\frac{2}{N - 1} + \frac{\kappa}{N}\right\\,\$\$
against \\2\sigma^4/(N - 1)\\ under normality. With a nonzero
`true_kurtosis` the chi square law is replaced by a scaled chi square,
\\S^2/\sigma^2 \approx \chi^2\_\nu/\nu\\, with \$\$\nu =
\frac{2}{2/(N - 1) + \kappa/N}\$\$ chosen so that the approximation has
the exact mean and variance of \\S^2/\sigma^2\\; \\\nu = N - 1\\ when
\\\kappa = 0\\. Tails heavier than the normal's (\\\kappa \> 0\\) make
\\S\\ more variable, and for an assurance of .80 or more they make it
harder to reach, while lighter tails (\\\kappa \< 0\\, as for scores on
a bounded scale) make it easier. Precisely, with an assurance of .80 or
more, when the plan made with a given \\\kappa\\ has at least \\\kappa +
3\\ observations, no smaller planning value of \\\kappa\\ gives a larger
plan, so an overstated \\\kappa\\ errs on the side of a larger *N*. With
an assurance close to one half the criterion turns instead on the median
of \\S^2\\, which heavier tails pull further below \\\sigma^2\\, and a
larger \\\kappa\\ can give a slightly smaller plan. The argument is on
the excess scale of
[`kurtosis`](https://yelleknek.github.io/DMAR/reference/kurtosis.md), so
`kurtosis(pilot)` can usually be passed as it is (see the
`true_kurtosis` argument for a pilot that gives -2 or less). A pilot's
kurtosis is itself imprecise. The kurtosis does not enter the expected
width plan: there it would lower \\E(S)\\ by a relative amount of about
\\\kappa/(8N)\\ and *N* by about \\\kappa/4\\, less than one when
\\\|\kappa\| \< 4\\, so the plan keeps normal theory's \\c_4(N)\\ and
the function warns when a nonzero `true_kurtosis` is supplied without an
`assurance`.

The two-moment approximation is coarse when \\\kappa\\ is large and *N*
is small. Once *N* is below about \\\kappa/2\\, \\\nu\\ is below one,
and the smaller \\\nu\\ is, the more of the scaled chi square's
probability sits near zero, much as most small samples from a population
with heavy tails miss the tails. The probability can then meet the
assurance at a very small *N*, fall below it over a stretch of larger
*N*, and meet it again. This is where the definition above matters:
`necessary_N` is the start of the stretch over which the criterion holds
at every larger *N*, not the first *N* at which it holds, and in this
region a larger \\\kappa\\ can give a smaller plan. The function finds
that *N* exactly: it bounds the criterion over whole ranges of *N*,
using only that \\\nu\\ and the bound on \\S^2/\sigma^2\\ both rise with
*N* and that a chi square quantile rises with its degrees of freedom,
and evaluates it directly where the bounds cannot decide.
[`ss_aipe_mean_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mean_sensitivity.md),
with a pilot's scores as the population, checks a plan without the
approximation, and it also shows the coverage of the *t* interval under
the shape of those scores, which the plan does not address.

**When no planning value is trustworthy.** A fixed *N* is only as good
as the planning value of \\\sigma\\. Sequential fixed-width confidence
interval procedures for a mean keep sampling until the interval reaches
the target width and need no planning value of \\\sigma\\.

## References

Kelley, K., Maxwell, S. E., & Rausch, J. R. (2003). Obtaining power or
obtaining precision: Delineating methods of sample-size planning.
*Evaluation & the Health Professions, 26*(3), 258–287.
[doi:10.1177/0163278703255242](https://doi.org/10.1177/0163278703255242)

Kelley, K., & Rausch, J. R. (2006). Sample size planning for the
standardized mean difference: Accuracy in parameter estimation via
narrow confidence intervals. *Psychological Methods, 11*(4), 363–385.
[doi:10.1037/1082-989X.11.4.363](https://doi.org/10.1037/1082-989X.11.4.363)

## See also

[`ss_aipe_mean_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mean_sensitivity.md)
for a Monte Carlo check of a plan;
[`ss_aipe_c`](https://yelleknek.github.io/DMAR/reference/ss_aipe_c.md)
for a contrast among independent groups;
[`ss_aipe_sm`](https://yelleknek.github.io/DMAR/reference/ss_aipe_sm.md)
for the standardized mean;
[`kurtosis`](https://yelleknek.github.io/DMAR/reference/kurtosis.md);
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
[`ss_aipe_two_phase()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_two_phase.md),
[`ss_aipe_two_phase_cost()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_two_phase_cost.md),
[`ss_aipe_two_phase_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_two_phase_sensitivity.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# A pilot: the IQ gains (end of year minus pretest) of the 246 children
# in the control group of the pygmalion data. A gain is a paired
# difference, so its mean is planned as a single mean.
gain <- pygmalion$iq_gain[pygmalion$treatment == "Control"]
sd(gain)
#> [1] 13.29708
kurtosis(gain)
#> [1] 1.196691

# The N at which the expected full width of the 95 percent interval for
# the mean gain is at most 4 IQ points.
ss_aipe_mean(true_sd = sd(gain), width = 4)
#>  term           value 
#>  necessary_N    172   
#>  expected_width 3.997 
#>  true_sd        13.297
#>  width          4     
#>  conf_level     0.950 
#> 
#> Confidence level: 95%

# With 90 percent assurance that the realized width is at most 4 points,
# first under normal theory and then with the heavier than normal tails
# the pilot shows.
ss_aipe_mean(true_sd = sd(gain), width = 4, assurance = 0.90)
#>  term           value 
#>  necessary_N    195   
#>  expected_width 3.751 
#>  true_sd        13.297
#>  width          4     
#>  conf_level     0.950 
#>  assurance      0.900 
#>  true_kurtosis  0.000 
#> 
#> Confidence level: 95%
ss_aipe_mean(true_sd = sd(gain), width = 4, assurance = 0.90,
             true_kurtosis = kurtosis(gain))
#>  term           value 
#>  necessary_N    201   
#>  expected_width 3.694 
#>  true_sd        13.297
#>  width          4     
#>  conf_level     0.950 
#>  assurance      0.900 
#>  true_kurtosis  1.197 
#> 
#> Confidence level: 95%
```
