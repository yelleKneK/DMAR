# Sensitivity Analysis for Sample Size Planning From the AIPE Perspective for an Indirect Effect

Quantifies how much misspecification of the population mediation path
coefficients \\a\\ and \\b\\ distorts an AIPE-based sample size plan for
the indirect effect \\ab\\. On each replication the function simulates a
three-variable mediation system \\X \to M \to Y\\ of size *n* with
population standardized path coefficients `true_std_reg_coef_a` and
`true_std_reg_coef_b`, fits the two regressions of *M* on *X* and *Y* on
*M* and *X*, computes the sample indirect effect \\\hat a \hat b\\, and
forms the interval the plan targeted: the symmetric Wald interval from
the delta method standard error (`method = "closed_form"`) or the Monte
Carlo interval (`method = "monte_carlo"`), matching
[`ss_aipe_indirect_effect`](https://yelleknek.github.io/DMAR/reference/ss_aipe_indirect_effect.md).

## Usage

``` r
ss_aipe_indirect_effect_sensitivity(
  true_std_reg_coef_a = NULL,
  true_std_reg_coef_b = NULL,
  supposed_std_reg_coef_a = NULL,
  supposed_std_reg_coef_b = NULL,
  width,
  specified_N = NULL,
  method = c("closed_form", "monte_carlo"),
  conf_level = 0.95,
  B = 5000L,
  G = 1000,
  print_iter = FALSE,
  filename = NULL,
  seed = NULL
)
```

## Arguments

- true_std_reg_coef_a:

  Population standardized regression coefficient of path *a* (from *X*
  to *M*); the data generating value.

- true_std_reg_coef_b:

  Population standardized regression coefficient of path *b* (from *M*
  to *Y* after controlling for *X*); the data generating value.

- supposed_std_reg_coef_a, supposed_std_reg_coef_b:

  The standardized path coefficients supposed when the study is planned
  (passed to
  [`ss_aipe_indirect_effect`](https://yelleknek.github.io/DMAR/reference/ss_aipe_indirect_effect.md)
  as its `true_std_reg_coef_a` and `true_std_reg_coef_b`). They are
  planning assumptions, not estimates from data, and may differ from the
  data generating values. Supply both or neither (if neither, supply
  `specified_N`).

- width:

  Desired full width of the CI on \\ab\\.

- specified_N:

  Sample size to evaluate (incompatible with `supposed_std_reg_coef_a`
  and `supposed_std_reg_coef_b`).

- method:

  One of `"closed_form"` (default) or `"monte_carlo"`; the interval
  computed on each replication, also forwarded to the planner when the
  sample size is planned from `supposed_std_reg_coef_a` and
  `supposed_std_reg_coef_b`. A planning call with
  `method = "monte_carlo"` runs the planner's a priori Monte Carlo
  search at its default `G`, so it takes a few seconds.

- conf_level:

  Confidence level (default `0.95`).

- B:

  Number of Monte Carlo draws used for the indirect-effect CI when
  `method = "monte_carlo"` (default 5000). The interval is read from the
  quantiles of the draws, so `B` follows the rule of
  [`dmar_bootstrap`](https://yelleknek.github.io/DMAR/reference/dmar_bootstrap.md):
  at least 400 at the 95 percent level, unless
  `options(dmar.bootstrap_exploratory = TRUE)` is set.

- G:

  Number of outer simulation replications (default 1000).

- print_iter:

  Logical.

- filename:

  Optional path for a comma separated file recording every replication
  (the sample indirect effect \\\hat a \hat b\\ in the column
  `std_indirect_effect`, the two confidence limits, the interval width,
  and two indicators of whether the interval missed
  `true_std_reg_coef_a * true_std_reg_coef_b` below or above): nothing
  is written when `filename` is `NULL` (the default), a new file with a
  header row is created otherwise, an existing file at that path is
  appended to, and a throwaway run should point it at
  `tempfile(fileext = ".csv")`.

- seed:

  `NULL` (the default) to draw from the current random number stream, or
  a single number that makes the call reproducible; a supplied seed is
  set for this call only, and the caller's random number state is
  restored when the function returns.

## Value

A `data.frame` with rows for the mean, median, and SD of the sample
indirect effect \\\hat a \hat b\\ (`mean_std_indirect_effect`,
`median_std_indirect_effect`, `sd_std_indirect_effect`) and of the CI
width, the proportion of intervals at or below `width`, tail-specific
and overall non-coverage of the population value
`true_std_indirect_effect`
(`true_std_reg_coef_a * true_std_reg_coef_b`), and the input echoes
`total_N`, `true_std_reg_coef_a`, `true_std_reg_coef_b`,
`supposed_std_reg_coef_a`, `supposed_std_reg_coef_b` (`NA` when
`specified_N` was given), `width`, and `conf_level`. Under
`method = "monte_carlo"` the result carries the attributes `B`,
`B_used`, `n_dropped`, `B_required`, and `bootstrap_status` of
[`dmar_bootstrap`](https://yelleknek.github.io/DMAR/reference/dmar_bootstrap.md);
they record the request each replication ran (`B_used` equals `B`), and
`bootstrap_status` is `"exploratory"` when `B` is below `B_required`.

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

Preacher, K. J., & Kelley, K. (2011). Effect size measures for mediation
models: Quantitative strategies for communicating indirect effects.
*Psychological Methods, 16*(2), 93–115.
[doi:10.1037/a0022658](https://doi.org/10.1037/a0022658)

Tofighi, D., & Kelley, K. (2020). Improved inference in mediation
analysis: Introducing the model-based constrained optimization
procedure. *Psychological Methods, 25*(4), 496–515.
[doi:10.1037/met0000259](https://doi.org/10.1037/met0000259)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.

## See also

[`ss_aipe_indirect_effect`](https://yelleknek.github.io/DMAR/reference/ss_aipe_indirect_effect.md),
[`var_indirect_effect`](https://yelleknek.github.io/DMAR/reference/var_indirect_effect.md)

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
# Reduced replications and a wide target width keep this fast.
set.seed(113)
ss_aipe_indirect_effect_sensitivity(
  true_std_reg_coef_a = 0.4, true_std_reg_coef_b = 0.3,
  supposed_std_reg_coef_a = 0.4, supposed_std_reg_coef_b = 0.3,
  width = 0.40, method = "closed_form",
  G = 50, print_iter = FALSE
)
#>  term                       value
#>  mean_std_indirect_effect   0.128
#>  median_std_indirect_effect 0.134
#>  sd_std_indirect_effect     0.107
#>  mean_ci_width              0.435
#>  median_ci_width            0.434
#>  sd_ci_width                0.149
#>  pct_ci_less_width          0.420
#>  pct_ci_miss_low            0.000
#>  pct_ci_miss_high           0.100
#>  total_type_I_error         0.100
#>  total_N                    27   
#>  true_std_reg_coef_a        0.400
#>  true_std_reg_coef_b        0.300
#>  true_std_indirect_effect   0.120
#>  supposed_std_reg_coef_a    0.400
#>  supposed_std_reg_coef_b    0.300
#>  width                      0.400
#>  conf_level                 0.950
#> 
#> Replications: 50
#> 
#> Confidence level: 95%
```
