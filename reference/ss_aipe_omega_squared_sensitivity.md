# Sensitivity Analysis for Sample Size Planning From the AIPE Perspective for Omega Squared

Quantifies how much misspecification of the population \\\omega^2\\
distorts an AIPE-based sample size plan. The planner
[`ss_aipe_omega_squared`](https://yelleknek.github.io/DMAR/reference/ss_aipe_omega_squared.md)
solves for the smallest *N* that yields an expected CI width below the
target at the planning value. Here we generate `G` datasets from a
balanced one-way ANOVA with population \\\omega^2 =
\\`true_omega_squared` and `df_effect + 1` groups at the
planner-recommended *N*, compute the noncentral *F* confidence interval
on each replication via
[`ci_omega_squared`](https://yelleknek.github.io/DMAR/reference/ci_omega_squared.md),
and summarize the realized widths and coverage of `true_omega_squared`.

## Usage

``` r
ss_aipe_omega_squared_sensitivity(
  true_omega_squared = NULL,
  supposed_omega_squared = NULL,
  df_effect,
  width,
  specified_N = NULL,
  conf_level = 0.95,
  assurance = NULL,
  G = 1000,
  print_iter = FALSE,
  filename = NULL,
  seed = NULL
)
```

## Arguments

- true_omega_squared:

  Population \\\omega^2\\ (the data generating value); in \\\[0, 1)\\.

- supposed_omega_squared:

  The supposed \\\omega^2\\, the planning value from which
  [`ss_aipe_omega_squared`](https://yelleknek.github.io/DMAR/reference/ss_aipe_omega_squared.md)
  plans the study (not an estimate computed from data); supply this or
  `specified_N` but not both.

- df_effect:

  Numerator degrees of freedom for the omnibus *F*, equal to the number
  of groups minus 1.

- width:

  Desired full width of the confidence interval on \\\omega^2\\.

- specified_N:

  Total sample size to evaluate (incompatible with
  `supposed_omega_squared`).

- conf_level:

  Confidence level (default `0.95`).

- assurance:

  Optional assurance probability passed to
  [`ss_aipe_omega_squared`](https://yelleknek.github.io/DMAR/reference/ss_aipe_omega_squared.md)
  when resolving the planned sample size.

- G:

  Number of Monte Carlo replications (default 1000).

- print_iter:

  Logical. Print iteration index per replication.

- filename:

  Optional path for a comma separated file recording every replication
  (the sample \\\hat\omega^2\\, the two confidence limits, the interval
  width, and two indicators of whether the interval missed
  `true_omega_squared` below or above): nothing is written when
  `filename` is `NULL` (the default), a new file with a header row is
  created otherwise, an existing file at that path is appended to, and a
  throwaway run should point it at `tempfile(fileext = ".csv")`.

- seed:

  `NULL` (the default) to draw from the current random number stream, or
  a single number that makes the call reproducible; a supplied seed is
  set for this call only, and the caller's random number state is
  restored when the function returns.

## Value

A `data.frame` with rows for mean / median / SD of the realized
\\\hat\omega^2\\ and interval width, the proportion of intervals at or
below `width`, tail-specific and overall empirical non-coverage of
`true_omega_squared`, and the input echoes, including `assurance`
(present only when an assurance was supplied).

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

Kelley, K. (2007). Confidence intervals for standardized effect sizes:
Theory, application, and implementation. *Journal of Statistical
Software, 20*(8), 1–24.
[doi:10.18637/jss.v020.i08](https://doi.org/10.18637/jss.v020.i08)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 3 on effect size measures.)

## See also

[`ss_aipe_omega_squared`](https://yelleknek.github.io/DMAR/reference/ss_aipe_omega_squared.md),
[`ci_omega_squared`](https://yelleknek.github.io/DMAR/reference/ci_omega_squared.md)

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
# Well-specified: planner used omega^2 = 0.10, truth is 0.10.
# G is kept small here so the example runs quickly; raise it for a
# stable sensitivity estimate.
set.seed(113)
ss_aipe_omega_squared_sensitivity(
  true_omega_squared     = 0.10,
  supposed_omega_squared = 0.10,
  df_effect = 2, width = 0.10,
  G = 25, print_iter = FALSE
)
#> During the iterative sample size search, the noncentral F lower-limit clamp in ci_nc_F() fired in 10 intermediate evaluations.
#>  term                   value
#>  mean_omega_squared     0.067
#>  median_omega_squared   0.070
#>  sd_omega_squared       0.020
#>  mean_ci_width          0.085
#>  median_ci_width        0.088
#>  sd_ci_width            0.010
#>  pct_ci_less_width      0.920
#>  pct_ci_miss_low        0.000
#>  pct_ci_miss_high       0.320
#>  total_type_I_error     0.320
#>  total_N                471  
#>  n_per_group            157  
#>  true_omega_squared     0.100
#>  supposed_omega_squared 0.100
#>  width                  0.100
#>  conf_level             0.950
#> 
#> Replications: 25
#> 
#> Confidence level: 95%
```
