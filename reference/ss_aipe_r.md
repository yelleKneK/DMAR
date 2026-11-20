# Sample Size for AIPE on a Pearson Correlation

Determines the sample size needed for a confidence interval on a
population Pearson correlation \\\rho\\ to have a desired width
(accuracy in parameter estimation; Kelley & Maxwell, 2003). The interval
planned for is the Fisher's *Z* interval that
[`correlations_test`](https://yelleknek.github.io/DMAR/reference/correlations_test.md)
reports for `method = "pearson"`: the correlation is transformed as
\\Z(r) = \mathrm{atanh}(r)\\, an interval with standard error
\\1/\sqrt{N - 3}\\ is formed on the *Z* scale, and the limits are
back-transformed through \\\tanh(\cdot)\\ (Fisher, 1921; Bonett &
Wright, 2000, Equation 2). Because the plan targets the same interval
the analysis will report, the planned width and the analyzed width
agree.

## Usage

``` r
ss_aipe_r(true_cor, width, conf_level = 0.95, assurance = NULL)
```

## Arguments

- true_cor:

  Population Pearson correlation \\\rho\\ the sample size is planned
  for, in \\(-1, 1)\\. The planner treats this value as the truth;
  [`ss_aipe_r_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_r_sensitivity.md)
  shows what the plan delivers when the truth (`true_cor` there) differs
  from the value planned with (`supposed_cor` there).

- width:

  Desired full width of the confidence interval on the correlation.

- conf_level:

  Desired confidence level (default `0.95`). A level above .999999998,
  whose tail areas are below 1e-9, is refused (see "Small tail areas" in
  [`ci_nc_t`](https://yelleknek.github.io/DMAR/reference/ci_nc_t.md)).

- assurance:

  Optional. The probability that the realized interval is no wider than
  `width`. The function returns the smallest sample size, at or above
  the sample size planned without it, at which that probability is at
  least `assurance`. A single number above 0.5 and below 1; when `NULL`
  (the default) the plan is for the width at \\r = \rho\\ only.

## Value

A `data.frame` with the columns `term` and `value` and one row per term:

- `necessary_N`:

  the planned total sample size.

- `expected_width`:

  the full width of the interval at `necessary_N` when the sample
  correlation equals `true_cor`.

- `true_cor`, `width`, `conf_level`:

  the inputs.

- `assurance`:

  the requested probability; present only when `assurance` was given.

- `assurance_attained`:

  the probability, at `necessary_N` and `true_cor`, that the realized
  interval is no wider than `width`; present only when `assurance` was
  given. It is at least `assurance`, and above it because the sample
  size is a whole number or because the plan without `assurance` already
  exceeds the request.

## Details

**Closed-form first pass.** On the Fisher's *Z* scale the interval has
half-width \\z\_{1 - \alpha/2} / \sqrt{N - 3}\\, and the delta method
maps it back to the correlation scale as approximately \$\$w \\\approx\\
2\\ z\_{1 - \alpha/2} \\ \frac{1 - \rho^2}{\sqrt{N - 3}},\$\$ which
solves to the first-stage approximation of Bonett and Wright (2000),
\$\$n_0 \\=\\ 3 + \Big\lceil 4\\ (z\_{1 - \alpha/2})^2 (1 - \rho^2)^2 /
w^2 \Big\rceil.\$\$ Here and below \\\rho\\ is `true_cor`, \\w\\ is
`width`, and \\1 - \alpha\\ is `conf_level`.

**Exact iteration.** The back-transformed width depends on \\\rho\\
through \\\tanh(\cdot)\\, so the delta method approximation can land a
few observations off in either direction. Starting from \\n_0\\, the
function evaluates the exact back-transformed width \\\tanh(z\_\rho +
z\_{1 - \alpha/2}/\sqrt{N - 3}) - \tanh(z\_\rho - z\_{1 -
\alpha/2}/\sqrt{N - 3})\\ and steps the integer \\N\\ until it is the
smallest sample size whose width is at or below `width`. Where Bonett
and Wright (2000) stop after a single second-stage adjustment, this
search is exact.

**The assurance step.** A plan without `assurance` makes the interval as
narrow as the target when the sample correlation equals `true_cor`; in
many samples it is wider. With `assurance` the function finds the sample
size at which the realized interval is no wider than the target with at
least that probability, from the sampling distribution of the sample
correlation *r* of *N* independent observations from a bivariate normal
population:

1.  At a given *N* the width of the interval depends on the sample only
    through *r*, and it decreases as \\\|r\|\\ grows. The realized
    interval is therefore no wider than the target exactly when \\\|r\|
    \ge r^\*(N)\\, where \\r^\*(N)\\ is the root of
    \\\tanh(\mathrm{atanh}(r) + z\_{1 - \alpha/2}/\sqrt{N - 3}) -
    \tanh(\mathrm{atanh}(r) - z\_{1 - \alpha/2}/\sqrt{N - 3}) = w\\.
    When the width at \\r = 0\\, the widest the interval can be, already
    meets the target, every sample does and the probability is 1.

2.  The distribution function of *r* is computed as \$\$F(c) \\=\\
    \int_0^1 T\_{2N - 3}\\\Big(\sqrt{2N - 3}\\ \Big\\\frac{c\\\sqrt{1 -
    b}}{\sqrt{1 - c^2}} - \frac{\rho\\\sqrt{b}}{\sqrt{1 -
    \rho^2}}\Big\\\Big)\\ g(b)\\ db,\$\$ where \\T\_{2N - 3}\\ is the
    distribution function of Student's *t* on \\2N - 3\\ degrees of
    freedom and \\g\\ is the density of the beta distribution with
    parameters \\(N - 1)/2\\ and \\(N - 2)/2\\. The form follows from
    writing \\r / \sqrt{1 - r^2}\\ as \\(Z + \theta \sqrt{U}) /
    \sqrt{V}\\, with \\Z\\ standard normal, \\U\\ and \\V\\ chi square
    on \\N - 1\\ and \\N - 2\\ degrees of freedom, all independent, and
    \\\theta = \rho / \sqrt{1 - \rho^2}\\, and then conditioning on \\U
    / (U + V)\\. It is the computation
    [`ss_aipe_partial_r`](https://yelleknek.github.io/DMAR/reference/ss_aipe_partial_r.md)
    uses with no variables partialed out.

3.  The probability at *N* is \\P_N = 1 - F(r^\*(N)) + F(-r^\*(N))\\,
    and the function returns the smallest *N*, at or above the sample
    size planned without `assurance`, with \\P_N \ge\\ `assurance`,
    found by doubling a step upward and then bisecting.
    `assurance_attained` is \\P_N\\ at that *N*.

The probability is exact, up to numerical integration, under bivariate
normality. It was checked against the density of *r* computed from its
hypergeometric series, and against the proportion of intervals no wider
than the target among samples of bivariate normal data. The search does
not go below the plan without `assurance`, so a plan with `assurance` is
never smaller than the plan without it. Near \\\rho = 0\\ a slightly
smaller *N* can attain the probability, because there the interval at
\\r = 0\\ is barely wider than the target (at \\\rho = 0\\ and a width
of .10 the probability is .917 at 1537, one below the plan without
`assurance`). At and above the plan without `assurance`, \\P_N\\ did not
decrease from one *N* to the next in any of the 15 planning sets scanned
(\\\|\rho\|\\ of 0, .1, .5, .9, and .99 with widths of .05, .1, and .2).
The result depends on `true_cor` only through \\\|\rho\|\\. At
`true_cor = 0` the plan without `assurance` already makes every interval
as narrow as the target, since the interval is widest at \\r = 0\\;
`assurance` then adds nothing and `assurance_attained` is 1. As
\\\|\rho\|\\ grows the width becomes more variable from sample to sample
and the step adds more: for a width of .10 and an assurance of .80 the
plan moves from 62 to 83 at \\\rho = .90\\.

**Why the step is not a chi square factor.** Multiplying the sample size
planned without `assurance` by \\\chi^2\_{\gamma, N - 2} / (N - 2)\\,
the \\\gamma\\ quantile of a chi square on \\N - 2\\ degrees of freedom
over its degrees of freedom, with \\\gamma\\ the `assurance`, does not
serve here. That factor is the \\\gamma\\ quantile of the ratio of a
sample variance to its population value, so it applies to an interval
whose width is proportional to a sample standard deviation; the width of
the interval for a correlation depends on the sample through *r*
instead. At `true_cor = .90`, a width of .10, and an assurance of .80
the factor gives 72, where the probability is .672 (83 are needed); with
an assurance of .99 it gives 92, where the probability is .887 (122 are
needed). At `true_cor = 0` it gives more than the plan without
`assurance`, which already attains a probability of 1.

**The planning value matters least near zero.** At a fixed sample size
the back-transformed width is largest at \\\rho = 0\\ and shrinks as
\\\|\rho\|\\ grows, so a planning value closer to zero yields a larger,
more conservative sample size. When little is known about the population
correlation, `true_cor = 0` gives the sample size that suffices for any
population value.

**When to use simple vs. partial correlation planning.** Use this
function when the inferential target is the correlation between two
variables with nothing partialed out. When the target is the correlation
after statistically controlling for other variables, see
[`ss_aipe_partial_r`](https://yelleknek.github.io/DMAR/reference/ss_aipe_partial_r.md).

The Monte Carlo companion
[`ss_aipe_r_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_r_sensitivity.md)
evaluates how the plan behaves when the population correlation differs
from the planning value.

## References

Bonett, D. G., & Wright, T. A. (2000). Sample size requirements for
estimating Pearson, Kendall and Spearman correlations. *Psychometrika,
65*(1), 23–28.
[doi:10.1007/BF02294183](https://doi.org/10.1007/BF02294183)

Fisher, R. A. (1921). On the "probable error" of a coefficient of
correlation deduced from a small sample. *Metron, 1*, 3–32.

Kelley, K., & Maxwell, S. E. (2003). Sample size for multiple
regression: Obtaining regression coefficients that are accurate, not
simply significant. *Psychological Methods, 8*(3), 305–321.
[doi:10.1037/1082-989X.8.3.305](https://doi.org/10.1037/1082-989X.8.3.305)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 9 on correlations.)

## See also

[`ss_aipe_r_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_r_sensitivity.md),
[`correlations_test`](https://yelleknek.github.io/DMAR/reference/correlations_test.md),
[`ss_aipe_partial_r`](https://yelleknek.github.io/DMAR/reference/ss_aipe_partial_r.md),
[`ss_power_r`](https://yelleknek.github.io/DMAR/reference/ss_power_r.md),
[`convert_r_Z`](https://yelleknek.github.io/DMAR/reference/convert_r_Z.md)

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
[`ss_aipe_partial_r()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_partial_r.md),
[`ss_aipe_partial_r_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_partial_r_sensitivity.md),
[`ss_aipe_pcm_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_pcm_sensitivity.md),
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
# Plan N so the 95% CI on the Pearson correlation has full width
# at most 0.20, when the anticipated correlation is 0.30.
ss_aipe_r(true_cor = 0.30, width = 0.20)
#>  term           value
#>  necessary_N    320  
#>  expected_width 0.200
#>  true_cor       0.300
#>  width          0.200
#>  conf_level     0.950
#> 
#> Confidence level: 95%

# A narrower target width requires a larger sample size.
ss_aipe_r(true_cor = 0.30, width = 0.10)
#>  term           value
#>  necessary_N    1274 
#>  expected_width 0.100
#>  true_cor       0.300
#>  width          0.100
#>  conf_level     0.950
#> 
#> Confidence level: 95%

# The smallest N at which the realized interval is no wider than
# 0.20 with probability at least .80. The last row is the probability
# at the returned N.
ss_aipe_r(true_cor = 0.30, width = 0.20, assurance = 0.80)
#>  term               value
#>  necessary_N        336  
#>  expected_width     0.195
#>  true_cor           0.300
#>  width              0.200
#>  conf_level         0.950
#>  assurance          0.800
#>  assurance_attained 0.805
#> 
#> Confidence level: 95%

# At a large correlation the width varies more from sample to sample,
# and the assurance step adds more cases.
ss_aipe_r(true_cor = 0.90, width = 0.10)
#>  term           value
#>  necessary_N    62   
#>  expected_width 0.100
#>  true_cor       0.900
#>  width          0.100
#>  conf_level     0.950
#> 
#> Confidence level: 95%
ss_aipe_r(true_cor = 0.90, width = 0.10, assurance = 0.80)
#>  term               value
#>  necessary_N        83   
#>  expected_width     0.085
#>  true_cor           0.900
#>  width              0.100
#>  conf_level         0.950
#>  assurance          0.800
#>  assurance_attained 0.807
#> 
#> Confidence level: 95%

# Planning at true_cor = 0 gives the sample size that suffices for any
# population correlation, since the interval is widest there.
ss_aipe_r(true_cor = 0, width = 0.20)
#>  term           value
#>  necessary_N    385  
#>  expected_width 0.200
#>  true_cor       0    
#>  width          0.200
#>  conf_level     0.950
#> 
#> Confidence level: 95%
```
