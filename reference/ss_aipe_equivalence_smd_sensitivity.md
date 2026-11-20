# Sensitivity Analysis for Sample Size Planning From the AIPE Perspective for an Equivalence-Test SMD

Quantifies how much misspecification of the population standardized mean
difference distorts an AIPE-based sample size plan for the
two-one-sided-tests (TOST) confidence interval on the standardized mean
difference. On each replication the function simulates two normal groups
of size *n* per group with population standardized mean difference
`true_smd`, computes the sample standardized mean difference and the
interval
[`equivalence_smd`](https://yelleknek.github.io/DMAR/reference/equivalence_smd.md)
reports for it, and summarizes the realized widths, the proportion of
replications in which the interval falls entirely inside the equivalence
bounds (`equivalent`), and the proportions that miss `true_smd`.

## Usage

``` r
ss_aipe_equivalence_smd_sensitivity(
  true_smd = 0,
  supposed_smd = NULL,
  width,
  smd_lower = NULL,
  smd_upper = NULL,
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

- true_smd:

  Population standardized mean difference (the data generating value), a
  single finite number. Defaults to `0` (perfect equivalence).

- supposed_smd:

  The supposed population SMD, the planning value passed to
  [`ss_aipe_equivalence_smd`](https://yelleknek.github.io/DMAR/reference/ss_aipe_equivalence_smd.md)
  (not an estimate computed from data), a single finite number; supply
  this or `specified_n_per_group` but not both.

- width:

  Desired full width of the two-sided CI on the SMD.

- smd_lower, smd_upper:

  Equivalence bounds on the SMD, as positive magnitudes with the same
  meaning as in
  [`equivalence_smd`](https://yelleknek.github.io/DMAR/reference/equivalence_smd.md):
  the region is \\(-\delta_L, +\delta_U)\\, with `smd_lower` as
  \\\delta_L\\ and `smd_upper` as \\\delta_U\\. `smd_upper` is required;
  `smd_lower` defaults to `smd_upper` (a symmetric region). The
  simulator records whether the realized CI falls entirely inside the
  region.

- specified_n_per_group:

  Per-group sample size to evaluate, a whole number of at least 2.

- conf_level:

  Confidence level of the interval (default `0.95`); the two one-sided
  tests it decides are at level (1 - `conf_level`) / 2.

- assurance:

  Optional assurance probability, passed to
  [`ss_aipe_equivalence_smd`](https://yelleknek.github.io/DMAR/reference/ss_aipe_equivalence_smd.md)
  when the sample size is planned from `supposed_smd`.

- G:

  Number of Monte Carlo replications.

- print_iter:

  Logical.

- filename:

  Optional path for a comma separated file recording every replication
  (the sample standardized mean difference, the two confidence limits,
  the interval width, whether the interval fell inside the equivalence
  region, and two indicators of whether the interval missed `true_smd`
  below or above): nothing is written when `filename` is `NULL` (the
  default), a new file with a header row is created otherwise, an
  existing file at that path is appended to, and a throwaway run should
  point it at `tempfile(fileext = ".csv")`.

- seed:

  `NULL` (the default) to draw from the current random number stream, or
  a single number that makes the call reproducible; a supplied seed is
  set for this call only, and the caller's random number state is
  restored when the function returns.

## Value

A `data.frame` with rows for mean / median / SD of the realized SMD and
CI width, the proportion of intervals at or below `width`, tail-specific
and overall non-coverage of `true_smd`, the proportion of intervals
classified as `equivalent` (CI fully inside the bounds), and the input
echoes, including `assurance` (present only when an assurance was
supplied).

## Details

**The interval simulated.** The interval is the one
[`equivalence_smd`](https://yelleknek.github.io/DMAR/reference/equivalence_smd.md)
reports and
[`ss_aipe_equivalence_smd`](https://yelleknek.github.io/DMAR/reference/ss_aipe_equivalence_smd.md)
plans for, \$\$\hat d \pm t\_{1 - \alpha,\\ 2n - 2} \sqrt{(8 + \hat
d^{\\2}) / (4n)},\$\$ a large-sample standard error of \\\hat d\\ with a
Student *t* critical value, at confidence level `conf_level`, so that
\\\alpha\\ is (1 - `conf_level`) / 2: the default 0.95 is the interval
of the two one-sided tests at \\\alpha = .025\\, and `conf_level = 0.90`
is the interval of the tests at .05. A sample size from `supposed_smd`
is planned by
`ss_aipe_equivalence_smd(alpha_level = (1 - conf_level) / 2)`, so the
plan, the analysis, and the simulation describe one interval at one
level: a simulated sample with standardized mean difference \\\hat d\\
gives the limits
`equivalence_smd(smd = d, n_1 = n, n_2 = n, alpha_level = (1 - conf_level) / 2)`
reports. When `true_smd` equals `supposed_smd`, `pct_ci_less_width`
estimates the probability that the interval at the planned size is no
wider than `width`, which the planner computes and reports as
`assurance_attained` when `assurance` is given; the `simulation`
attribute (below) gives the Monte Carlo limits of the estimate.

The interval is a large-sample approximation, so its coverage of
`true_smd` (the `pct_ci_miss_low`, `pct_ci_miss_high`, and
`total_type_I_error` rows) can differ from 1 - `conf_level`, most in
small samples; the noncentral *t* interval of
[`ci_smd`](https://yelleknek.github.io/DMAR/reference/ci_smd.md), which
[`ss_aipe_smd_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_smd_sensitivity.md)
simulates, has exact coverage for the normal model.

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

Kelley, K., & Rausch, J. R. (2006). Sample size planning for the
standardized mean difference: Accuracy in parameter estimation via
narrow confidence intervals. *Psychological Methods, 11*, 363–385.
[doi:10.1037/1082-989X.11.4.363](https://doi.org/10.1037/1082-989X.11.4.363)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.

## See also

[`ss_aipe_equivalence_smd`](https://yelleknek.github.io/DMAR/reference/ss_aipe_equivalence_smd.md),
[`equivalence_smd`](https://yelleknek.github.io/DMAR/reference/equivalence_smd.md),
[`ss_aipe_smd_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_smd_sensitivity.md)

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
# Reduced Monte Carlo sweep (small G) for a fast, illustrative run.
set.seed(113)
ss_aipe_equivalence_smd_sensitivity(
  true_smd      = 0.0,
  supposed_smd  = 0.0,
  width         = 0.30,
  smd_upper     = 0.20,
  G = 50, print_iter = FALSE
)
#>  term               value 
#>  mean_smd           0.000 
#>  median_smd         0.010 
#>  sd_smd             0.087 
#>  mean_ci_width      0.300 
#>  median_ci_width    0.300 
#>  sd_ci_width        0.000 
#>  pct_ci_less_width  0.620 
#>  pct_equivalent     0.460 
#>  pct_ci_miss_low    0.020 
#>  pct_ci_miss_high   0.040 
#>  total_type_I_error 0.060 
#>  n_per_group        343   
#>  total_N            686   
#>  true_smd           0     
#>  supposed_smd       0     
#>  width              0.300 
#>  conf_level         0.950 
#>  smd_lower          -0.200
#>  smd_upper          0.200 
#> 
#> Replications: 50
#> 
#> Confidence level: 95%
```
