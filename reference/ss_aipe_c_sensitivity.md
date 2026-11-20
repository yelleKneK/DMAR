# Sensitivity Analysis for Sample Size Planning From the AIPE Perspective for an Unstandardized Contrast

Quantifies how much misspecification of the population error variance
distorts an AIPE-based sample size plan for an unstandardized contrast
of means. Because the half-width of the confidence interval on \\\psi =
\sum_j c_j \mu_j\\ depends on the error variance, the contrast weights,
and the per-group sample size, but not on the value of \\\psi\\ itself,
this sensitivity analysis varies the planning value of the error
variance. On each replication the function simulates \\n\\ observations
per group from a normal population with variance `true_var_error`,
builds the confidence interval via
[`ci_c`](https://yelleknek.github.io/DMAR/reference/ci_c.md), and
summarizes the realized widths and coverage of `true_contrast`.

## Usage

``` r
ss_aipe_c_sensitivity(
  true_var_error = NULL,
  supposed_var_error = NULL,
  c_weights,
  width,
  true_contrast = 0,
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

- true_var_error:

  Population error variance (the data generating value). Must be
  positive.

- supposed_var_error:

  Error variance used to plan the study (the value passed to
  [`ss_aipe_c`](https://yelleknek.github.io/DMAR/reference/ss_aipe_c.md)).
  Supply this or `specified_n_per_group` but not both.

- c_weights:

  Contrast weight vector. Must sum to zero.

- width:

  Desired full width of the confidence interval on the unstandardized
  contrast.

- true_contrast:

  Population value of the contrast; the simulator places group means
  such that \\\sum_j c_j \mu_j = \\`true_contrast`. The width of the
  interval does not depend on this value but the realized coverage of
  `true_contrast` does. Default `0`.

- specified_n_per_group:

  Per-group sample size to evaluate (incompatible with
  `supposed_var_error`); when used, the planner is bypassed.

- conf_level:

  Confidence level (default `0.95`). A level above .999999998, whose
  tail areas are below 1e-9, is refused (see "Small tail areas" in
  [`ci_nc_t`](https://yelleknek.github.io/DMAR/reference/ci_nc_t.md)).

- assurance:

  Optional probability that the realized interval is no wider than
  `width`; passed to
  [`ss_aipe_c`](https://yelleknek.github.io/DMAR/reference/ss_aipe_c.md)
  when resolving the planned sample size.

- G:

  Number of Monte Carlo replications (default 1000).

- print_iter:

  Logical. Print the iteration index after each replication (helpful for
  long runs); default `FALSE`.

- filename:

  Optional path of a CSV file to receive the per-replication results
  (the contrast estimate, the confidence limits, the interval width, and
  the two tail misses), appended when the file already exists and
  created otherwise; the default `NULL` writes nothing, and a throwaway
  run that wants the file should point it at
  `tempfile(fileext = ".csv")`.

- seed:

  `NULL` (the default) to draw from the current random number stream, or
  a single number that makes the call reproducible; a supplied seed is
  set for this call only, and the caller's random number state is
  restored when the function returns.

## Value

A `data.frame` with rows for mean / median / SD of the realized contrast
(`mean_contrast`, `median_contrast`, `sd_contrast`) and of the interval
width, the proportion of intervals at or below `width`, the
tail-specific and overall empirical non-coverage of `true_contrast`, and
the input echoes (per-group sample size, total sample size, the true and
supposed error variances `true_var_error` and `supposed_var_error`,
`true_contrast`, width, confidence level, and, when one was supplied,
assurance).

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

Lai, K., & Kelley, K. (2012). Accuracy in parameter estimation for
ANCOVA and ANOVA contrasts: Sample size planning via narrow confidence
intervals. *British Journal of Mathematical and Statistical Psychology,
65*, 350–370.
[doi:10.1111/j.2044-8317.2011.02029.x](https://doi.org/10.1111/j.2044-8317.2011.02029.x)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 4 on individual comparisons.)

## See also

[`ss_aipe_c`](https://yelleknek.github.io/DMAR/reference/ss_aipe_c.md),
[`ci_c`](https://yelleknek.github.io/DMAR/reference/ci_c.md),
[`ss_aipe_sc_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_sc_sensitivity.md)

[`design_consequences`](https://yelleknek.github.io/DMAR/reference/design_consequences.md)
for what a chosen design delivers: power, the Type S (sign) and Type M
(exaggeration) errors of the significance filter, and the expected
confidence interval width.

Other AIPE sample size planning:
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
[`ss_aipe_two_phase()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_two_phase.md),
[`ss_aipe_two_phase_cost()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_two_phase_cost.md),
[`ss_aipe_two_phase_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_two_phase_sensitivity.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# G = 50 keeps the example quick; a reported analysis deserves the
# default G of 1000. First the well-specified case: the planner used an
# error variance of 4, and the truth is 4.
set.seed(113)
ss_aipe_c_sensitivity(
  true_var_error     = 4,
  supposed_var_error = 4,
  c_weights = c(-1, 0, 1),
  width = 1, G = 50, print_iter = FALSE
)
#>  term               value 
#>  mean_contrast      -0.010
#>  median_contrast    -0.006
#>  sd_contrast        0.213 
#>  mean_ci_width      0.995 
#>  median_ci_width    1.000 
#>  sd_ci_width        0.032 
#>  pct_ci_less_width  0.520 
#>  pct_ci_miss_low    0.000 
#>  pct_ci_miss_high   0.020 
#>  total_type_I_error 0.020 
#>  n_per_group        124   
#>  total_N            372   
#>  true_var_error     4     
#>  supposed_var_error 4     
#>  true_contrast      0     
#>  width              1     
#>  conf_level         0.950 
#> 
#> Replications: 50
#> 
#> Confidence level: 95%

# Misspecified: planner used 4, truth is 9. Realized widths inflate.
set.seed(113)
ss_aipe_c_sensitivity(
  true_var_error     = 9,
  supposed_var_error = 4,
  c_weights = c(-1, 0, 1),
  width = 1, G = 50, print_iter = FALSE
)
#>  term               value 
#>  mean_contrast      -0.015
#>  median_contrast    -0.009
#>  sd_contrast        0.320 
#>  mean_ci_width      1.493 
#>  median_ci_width    1.499 
#>  sd_ci_width        0.048 
#>  pct_ci_less_width  0.000 
#>  pct_ci_miss_low    0.000 
#>  pct_ci_miss_high   0.020 
#>  total_type_I_error 0.020 
#>  n_per_group        124   
#>  total_N            372   
#>  true_var_error     9     
#>  supposed_var_error 4     
#>  true_contrast      0     
#>  width              1     
#>  conf_level         0.950 
#> 
#> Replications: 50
#> 
#> Confidence level: 95%
```
