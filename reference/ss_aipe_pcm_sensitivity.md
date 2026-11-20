# Sensitivity Analysis for Sample Size Planning From the AIPE Perspective for a Polynomial Change Parameter

Quantifies how much misspecification of the population between-subject
slope variance and within-subject error variance distorts an AIPE-based
sample size plan for the group-by-time polynomial change parameter. On
each replication the function simulates two independent groups of *n*
subjects each, measured at \\M = f \times D + 1\\ timepoints, where
every subject has a true linear slope drawn from a normal distribution
with mean 0 and variance `true_var_trend`, and within-subject
observations have residual variance `true_var_error`. Subject-level OLS
slopes are computed in each group, the between-group difference in mean
slopes (the change parameter \\\beta\_{m1}\\ that
[`ss_aipe_pcm`](https://yelleknek.github.io/DMAR/reference/ss_aipe_pcm.md)
plans for) is estimated, and a two-group *t*-confidence interval on that
difference (pooled standard error, \\2n - 2\\ degrees of freedom) is
recorded. The function only handles `trend = "linear"` in the simulator;
for quadratic / cubic trends, the planner's closed-form solution is
still available via
[`ss_aipe_pcm`](https://yelleknek.github.io/DMAR/reference/ss_aipe_pcm.md).

## Usage

``` r
ss_aipe_pcm_sensitivity(
  true_var_trend = NULL,
  true_var_error = NULL,
  supposed_var_trend = NULL,
  supposed_var_error = NULL,
  duration,
  frequency,
  width,
  specified_n_per_group = NULL,
  conf_level = 0.95,
  assurance = NULL,
  G = 1000,
  print_iter = FALSE,
  filename = NULL,
  seed = NULL
)
```

## Arguments

- true_var_trend:

  Population between-subject variance of the polynomial change
  coefficient (the data generating \\\sigma^2\_{\upsilon_m}\\ of Kelley
  & Rausch, 2011).

- true_var_error:

  Population within-subject error variance (\\\sigma^2\_\epsilon\\).

- supposed_var_trend:

  Planning value of `true_var_trend` passed to
  [`ss_aipe_pcm`](https://yelleknek.github.io/DMAR/reference/ss_aipe_pcm.md);
  supply this and `supposed_var_error`, or supply
  `specified_n_per_group`.

- supposed_var_error:

  Planning value of `true_var_error` passed to
  [`ss_aipe_pcm`](https://yelleknek.github.io/DMAR/reference/ss_aipe_pcm.md).

- duration:

  Study duration (in time units).

- frequency:

  Number of measurements per unit time. The number of measurement
  occasions is \\f \times D + 1\\ (the row `n_occasions`).

- width:

  Desired full width of the CI on the between-group difference in change
  parameters (\\\beta\_{m1}\\).

- specified_n_per_group:

  Number of subjects per group to evaluate (incompatible with the
  supposed-variance arguments).

- conf_level:

  Confidence level (default `0.95`).

- assurance:

  Optional assurance probability passed to the planner.

- G:

  Number of Monte Carlo replications.

- print_iter:

  Logical.

- filename:

  Optional path to a CSV file; when supplied, the per-replication
  results (the slope difference, its interval limits and width, and the
  two tail misses) are written there, appended when the file already
  exists, and a throwaway run should point it at
  `tempfile(fileext = ".csv")`; the default `NULL` writes nothing.

- seed:

  `NULL` (the default) to draw from the current random number stream, or
  a single number that makes the call reproducible; a supplied seed is
  set for this call only, and the caller's random number state is
  restored when the function returns.

## Value

A `data.frame` with rows for mean / median / SD of the realized
estimated slope difference and CI width, the proportion of intervals at
or below `width`, tail-specific and overall non-coverage of the
population slope difference (0 by construction in this simulator), and
the input echoes, including `assurance` (present only when an assurance
was supplied).

## The record of the simulation

Every sensitivity function of the package summarizes its replications
over one denominator: the replications that produced a confidence
interval (for a power sensitivity function, a test). The proportion of
intervals that met the width target (`pct_ci_less_width`), the
proportions of intervals that missed the true value, and the summaries
of the widths are all taken over those replications. A replication
without an interval is in none of these rows. It is counted, and the
returned table carries the counts in its attribute `simulation`, a list
read with `attr(x, "simulation")`:

- `attempted`:

  the number of replications asked for (`G`).

- `returned`:

  the number that produced an interval (a test).

- `failed`:

  the number that did not (`attempted - returned`).

- `denominator`:

  the number of replications the rows named above are taken over, which
  is `returned`.

- `width_target_probability`:

  the proportion of the replications with an interval whose interval met
  the width target, the value of the `pct_ci_less_width` row; `NA` for a
  function or a call without a width target, and when no replication
  produced an interval.

- `width_target_mc_ci`:

  the exact binomial (Clopper-Pearson) limits for that proportion, at
  the confidence level of the table. They describe the error of the
  simulation, which shrinks as `G` grows, and nothing else.

A function that records more about its run adds elements to the list,
and its help page names them. The line "Replications" beneath the
printed table gives `G`, and the three counts when a replication was
without an interval. With many replications without an interval the rows
describe only the studies in which the interval could be formed, so
`failed` is to be reported beside them.

## References

Kelley, K., & Rausch, J. R. (2011). Sample size planning for
longitudinal models: Accuracy in parameter estimation for polynomial
change parameters. *Psychological Methods, 16*(4), 391–405.
[doi:10.1037/a0023352](https://doi.org/10.1037/a0023352)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapters 11 and 15.)

## See also

[`ss_aipe_pcm`](https://yelleknek.github.io/DMAR/reference/ss_aipe_pcm.md),
[`ss_aipe_mixed_effects_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mixed_effects_sensitivity.md)

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
# Every replication simulates two full groups of subjects, fits a slope
# for each subject, and forms a confidence interval on the difference in
# mean slopes. G = 10 keeps the example quick; a reported sensitivity
# study deserves the default G = 1000. With the planning values equal to
# the population values, the realized mean width should sit close to
# the target.
set.seed(113)
ss_aipe_pcm_sensitivity(
  true_var_trend = 0.003, true_var_error = 0.0262,
  supposed_var_trend = 0.003, supposed_var_error = 0.0262,
  duration = 4, frequency = 1, width = 0.05,
  G = 10, print_iter = FALSE
)
#>  term               value
#>  mean_slope_diff    0.004
#>  median_slope_diff  0.002
#>  sd_slope_diff      0.012
#>  mean_ci_width      0.051
#>  median_ci_width    0.051
#>  sd_ci_width        0.003
#>  pct_ci_less_width  0.500
#>  pct_ci_miss_low    0.100
#>  pct_ci_miss_high   0.000
#>  total_type_I_error 0.100
#>  n_per_group        71   
#>  n_occasions        5    
#>  true_var_trend     0.003
#>  true_var_error     0.026
#>  supposed_var_trend 0.003
#>  supposed_var_error 0.026
#>  duration           4    
#>  frequency          1    
#>  width              0.050
#>  conf_level         0.950
#> 
#> Replications: 10
#> 
#> Confidence level: 95%
```
