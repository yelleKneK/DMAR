# Sensitivity Analysis for Sample Size Planning From the AIPE Perspective for a Polynomial Change Parameter

Quantifies how much misspecification of the population between-subject
slope variance and within-subject error variance distorts an AIPE-based
sample size plan for the group-by-time polynomial change parameter. On
each replication the function simulates two independent groups of *n*
subjects each, measured at \\M = f \times D + 1\\ timepoints, where
every subject has a true linear slope drawn from \\N(0,
\mathrm{true\\variance\\trend})\\ and within-subject observations have
residual variance `true_error_variance`. Subject-level OLS slopes are
computed in each group, the between-group difference in mean slopes (the
change parameter \\\beta\_{m1}\\ that
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
  true_variance_trend = NULL,
  true_error_variance = NULL,
  estimated_variance_trend = NULL,
  estimated_error_variance = NULL,
  duration,
  frequency,
  width,
  n_per_group = NULL,
  conf_level = 0.95,
  assurance = NULL,
  G = 1000,
  print_iter = FALSE,
  save = FALSE,
  filename = "ss_aipe_pcm_sensitivity_result.csv"
)
```

## Arguments

- true_variance_trend:

  Population between-subject variance of the polynomial change
  coefficient (the data generating \\\sigma^2\_{\upsilon_m}\\ of Kelley
  & Rausch, 2011).

- true_error_variance:

  Population within-subject error variance (\\\sigma^2\_\epsilon\\).

- estimated_variance_trend:

  Planning value of `variance_trend` passed to
  [`ss_aipe_pcm`](https://yelleknek.github.io/DMAR/reference/ss_aipe_pcm.md);
  supply this and `estimated_error_variance`, or supply `n_per_group`.

- estimated_error_variance:

  Planning value of `error_variance` passed to
  [`ss_aipe_pcm`](https://yelleknek.github.io/DMAR/reference/ss_aipe_pcm.md).

- duration:

  Study duration (in time units).

- frequency:

  Number of measurements per unit time. Total timepoints = \\f \times
  D + 1\\.

- width:

  Desired full width of the CI on the between-group difference in change
  parameters (\\\beta\_{m1}\\).

- n_per_group:

  Number of subjects to evaluate (incompatible with the
  estimated-variance arguments).

- conf_level:

  Confidence level (default `0.95`).

- assurance:

  Optional assurance probability passed to the planner.

- G:

  Number of Monte Carlo replications.

- print_iter:

  Logical.

- save:

  Logical. Save per-replication CSV.

- filename:

  Path used when `save = TRUE`.

## Value

A `data.frame` with rows for mean / median / SD of the realized
estimated slope difference and CI width, the proportion of intervals at
or below `width`, tail-specific and overall non-coverage of the
population slope difference (0 by construction in this simulator), and
the input echoes, including `assurance` (present only when an assurance
was supplied).

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
[`ss_aipe_equivalence_r()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_equivalence_r.md),
[`ss_aipe_equivalence_r_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_equivalence_r_sensitivity.md),
[`ss_aipe_equivalence_smd()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_equivalence_smd.md),
[`ss_aipe_equivalence_smd_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_equivalence_smd_sensitivity.md),
[`ss_aipe_icc()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_icc.md),
[`ss_aipe_icc_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_icc_sensitivity.md),
[`ss_aipe_indirect_effect()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_indirect_effect.md),
[`ss_aipe_indirect_effect_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_indirect_effect_sensitivity.md),
[`ss_aipe_mixed_effects_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mixed_effects_sensitivity.md),
[`ss_aipe_omega_squared()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_omega_squared.md),
[`ss_aipe_omega_squared_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_omega_squared_sensitivity.md),
[`ss_aipe_partial_r()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_partial_r.md),
[`ss_aipe_partial_r_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_partial_r_sensitivity.md),
[`ss_aipe_r()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_r.md),
[`ss_aipe_r_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_r_sensitivity.md),
[`ss_aipe_reliability_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reliability_sensitivity.md),
[`ss_aipe_semipartial_r()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_semipartial_r.md),
[`ss_aipe_semipartial_r_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_semipartial_r_sensitivity.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Every replication simulates two full groups of subjects, fits a
# slope for each subject, and forms a confidence interval on the
# difference in mean slopes, so the sweep is not run at example time.
# The G below is far smaller than a reported sensitivity study would
# use; the default of 1000 is the realistic setting. The call is:
# set.seed(113)
# ss_aipe_pcm_sensitivity(
#   true_variance_trend       = 0.003,
#   true_error_variance       = 0.0262,
#   estimated_variance_trend  = 0.003,
#   estimated_error_variance  = 0.0262,
#   duration  = 4, frequency = 1,
#   width     = 0.05,
#   G = 20, print_iter = FALSE
# )
```
