# Sensitivity Analysis for Sample Size Planning From the AIPE Perspective for a Semipartial Correlation

Quantifies how much misspecification of the population semipartial
correlation, or of the full-model squared multiple correlation, distorts
an AIPE-based sample size plan. The function constructs a \\(p + 1)
\times (p + 1)\\ population covariance matrix for *Y*, *X*, and \\Z_1,
\ldots, Z\_{p-1}\\ whose semipartial correlation between *Y* and *X*
(partialing \\Z_1, \ldots, Z\_{p-1}\\ out of *X* only, not out of *Y*)
equals `true_r_sp` and whose full-model \\R^2\\ equals `true_R2_full`.
On each replication it draws an *N*-row sample, computes the sample
semipartial correlation and the sample full-model \\R^2\\, and forms the
Wald interval that
[`ss_aipe_semipartial_r`](https://yelleknek.github.io/DMAR/reference/ss_aipe_semipartial_r.md)
plans for: the estimate plus and minus the normal critical value times
the square root of the first-order variance
([`var_semipartial_r`](https://yelleknek.github.io/DMAR/reference/var_semipartial_r.md)
with `under_null = FALSE`) at the two sample values.

## Usage

``` r
ss_aipe_semipartial_r_sensitivity(
  true_r_sp = NULL,
  supposed_r_sp = NULL,
  p,
  width,
  specified_N = NULL,
  conf_level = 0.95,
  assurance = NULL,
  G = 1000,
  print_iter = FALSE,
  filename = NULL,
  true_R2_full = NULL,
  supposed_R2_full = NULL,
  seed = NULL
)
```

## Arguments

- true_r_sp:

  Population semipartial correlation \\\rho\_{sp}\\; must lie in \\(-1,
  1)\\.

- supposed_r_sp:

  Planning value passed to
  [`ss_aipe_semipartial_r`](https://yelleknek.github.io/DMAR/reference/ss_aipe_semipartial_r.md);
  supply this or `specified_N` but not both.

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
  variables are partialed out. The simulated data have \\p + 1\\
  columns.

- width:

  Desired full width of the CI on the semipartial correlation.

- specified_N:

  Sample size to evaluate; a whole number of at least \\p + 2\\.

- conf_level:

  Confidence level (default `0.95`).

- assurance:

  Optional assurance probability passed to
  [`ss_aipe_semipartial_r`](https://yelleknek.github.io/DMAR/reference/ss_aipe_semipartial_r.md);
  used only with `supposed_r_sp`.

- G:

  Number of Monte Carlo replications (default 1000).

- print_iter:

  Logical. Print the iteration index of each replication.

- filename:

  Optional path for a comma separated file recording every replication
  (the sample semipartial correlation, the two confidence limits, the
  interval width, and two indicators of whether the interval missed
  `true_r_sp` below or above): nothing is written when `filename` is
  `NULL` (the default), a new file with a header row is created
  otherwise, an existing file at that path is appended to, and a
  throwaway run should point it at `tempfile(fileext = ".csv")`.

- true_R2_full:

  Population squared multiple correlation of the full model (*Y* on *X*
  and \\Z_1, \ldots, Z\_{p-1}\\); at least `true_r_sp` squared and
  below 1. When `NULL` (the default) it is `true_r_sp` squared, a
  population in which the variables partialed out explain none of the
  outcome.

- supposed_R2_full:

  Optional planning value of the full-model \\R^2\\, passed to
  [`ss_aipe_semipartial_r`](https://yelleknek.github.io/DMAR/reference/ss_aipe_semipartial_r.md)
  as `true_R2_full`; when `NULL` that function plans for the largest
  variance over every possible value. It cannot be given with
  `specified_N`.

- seed:

  `NULL` (the default) to draw from the current random number stream, or
  a single number that makes the run reproducible; a supplied seed is
  set for this call only, and the caller's random number state is
  restored when the function returns. With `assurance`, the seed also
  makes the simulated plan of
  [`ss_aipe_semipartial_r`](https://yelleknek.github.io/DMAR/reference/ss_aipe_semipartial_r.md)
  reproducible.

## Value

A `data.frame` with rows for the realized semipartial correlation, the
interval width, the proportion of intervals at or below `width`,
tail-specific and overall non-coverage of `true_r_sp`, and the input
echoes, including `true_R2_full` (the value simulated),
`supposed_R2_full` (`NA` when none was supplied), and `assurance`
(present only when an assurance was supplied).

## Details

Under multivariate normality the sampling distribution of the sample
semipartial correlation depends on the population only through the
semipartial correlation, the full-model \\R^2\\, \\p\\, and the sample
size, because the statistic is unchanged by any nonsingular linear
transformation of \\Z_1, \ldots, Z\_{p-1}\\, by adding any linear
function of them to *X*, and by positive rescaling of *X* and *Y*. The
population simulated here (uncorrelated standardized *X* and \\Z_1,
\ldots, Z\_{p-1}\\, with *Y* related to *X* and to \\Z_1\\) therefore
stands for every normal population with the same two values.

The confidence limits are not truncated at \\-1\\ and \\1\\, so that the
widths are those the planner plans; with a large \\\|\rho\_{sp}\|\\ and
a small sample a limit can pass 1 in absolute value.

The function simulates \\p - 1\\ variables partialed out, the count
[`ss_aipe_semipartial_r`](https://yelleknek.github.io/DMAR/reference/ss_aipe_semipartial_r.md)
plans for.

With `assurance`, the sample size comes from the simulated assurance
step of
[`ss_aipe_semipartial_r`](https://yelleknek.github.io/DMAR/reference/ss_aipe_semipartial_r.md)
(at its default number of draws), so two unseeded calls can simulate
sample sizes that differ by a case or two.

A value of `true_R2_full` or `supposed_R2_full` that falls below the
squared semipartial correlation by no more than rounding to two decimals
allows is taken to mean that the variables partialed out explain none of
the outcome, with a warning of class `dmar_semipartial_r_rounded`; a
smaller value is refused.

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

Kelley, K., & Maxwell, S. E. (2003). Sample size for multiple
regression: Obtaining regression coefficients that are accurate, not
simply significant. *Psychological Methods, 8*(3), 305–321.
[doi:10.1037/1082-989X.8.3.305](https://doi.org/10.1037/1082-989X.8.3.305)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.

## See also

[`ss_aipe_semipartial_r`](https://yelleknek.github.io/DMAR/reference/ss_aipe_semipartial_r.md),
[`ss_aipe_partial_r_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_partial_r_sensitivity.md)

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
[`ss_aipe_r()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_r.md),
[`ss_aipe_r_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_r_sensitivity.md),
[`ss_aipe_reliability_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reliability_sensitivity.md),
[`ss_aipe_semipartial_r()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_semipartial_r.md),
[`ss_aipe_two_phase()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_two_phase.md),
[`ss_aipe_two_phase_cost()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_two_phase_cost.md),
[`ss_aipe_two_phase_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_two_phase_sensitivity.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# G = 50 keeps the example quick; a reported check deserves the
# default G = 1000 or more.
ss_aipe_semipartial_r_sensitivity(
  true_r_sp = 0.30, supposed_r_sp = 0.30, p = 4, width = 0.20,
  G = 50, print_iter = FALSE,
  true_R2_full = 0.50, supposed_R2_full = 0.50, seed = 113
)
#>  term               value
#>  mean_r_sp          0.297
#>  median_r_sp        0.293
#>  sd_r_sp            0.045
#>  mean_ci_width      0.197
#>  median_ci_width    0.196
#>  sd_ci_width        0.007
#>  pct_ci_less_width  0.680
#>  pct_ci_miss_low    0.000
#>  pct_ci_miss_high   0.020
#>  total_type_I_error 0.020
#>  total_N            193  
#>  p                  4    
#>  true_r_sp          0.300
#>  supposed_r_sp      0.300
#>  width              0.200
#>  conf_level         0.950
#>  true_R2_full       0.500
#>  supposed_R2_full   0.500
#> 
#> Replications: 50
#> 
#> Confidence level: 95%
```
