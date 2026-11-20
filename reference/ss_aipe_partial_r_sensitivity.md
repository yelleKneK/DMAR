# Sensitivity Analysis for Sample Size Planning From the AIPE Perspective for a Partial Correlation

Simulates a sample size plan for a partial correlation to show what it
delivers when the population value differs from the planning value, or
at a sample size chosen in advance. It takes the population partial
correlation, either a planning value or a sample size, the number of
predictors (*X* and the variables partialed out), and a target width,
and returns summaries of the estimate, the interval width, the
proportion of intervals no wider than the target, and coverage over `G`
simulated samples.

## Usage

``` r
ss_aipe_partial_r_sensitivity(
  true_cor = NULL,
  supposed_cor = NULL,
  p,
  width,
  specified_N = NULL,
  conf_level = 0.95,
  assurance = NULL,
  G = 1000,
  print_iter = FALSE,
  filename = NULL,
  fisher_z = FALSE,
  seed = NULL
)
```

## Arguments

- true_cor:

  Population partial correlation \\\rho\\ between *Y* and *X*
  controlling for \\Z_1, \ldots, Z\_{p-1}\\; must lie in \\(-1, 1)\\.

- supposed_cor:

  Planning value of the partial correlation passed to
  [`ss_aipe_partial_r`](https://yelleknek.github.io/DMAR/reference/ss_aipe_partial_r.md);
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

  Desired full width of the CI on the partial correlation.

- specified_N:

  Sample size to evaluate (incompatible with `supposed_cor`); a whole
  number of at least \\p + 2\\, or \\p + 3\\ with `fisher_z = TRUE`.

- conf_level:

  Confidence level (default `0.95`).

- assurance:

  Optional assurance probability passed to
  [`ss_aipe_partial_r`](https://yelleknek.github.io/DMAR/reference/ss_aipe_partial_r.md);
  used only with `supposed_cor`.

- G:

  Number of Monte Carlo replications (default 1000).

- print_iter:

  Logical. Print iteration index per replication.

- filename:

  Optional path for a comma separated file recording every replication
  (the sample partial correlation, the two confidence limits, the
  interval width, and two indicators of whether the interval missed
  `true_cor` below or above): nothing is written when `filename` is
  `NULL` (the default), a new file with a header row is created
  otherwise, an existing file at that path is appended to, and a
  throwaway run should point it at `tempfile(fileext = ".csv")`.

- fisher_z:

  Logical. The interval that is planned and simulated. `FALSE` (the
  default) is the raw-scale interval, \\r \pm z\_{1 - \alpha/2} (1 -
  r^2) / \sqrt{N - p}\\; `TRUE` is the interval formed on Fisher's *Z*
  scale with standard error \\1/\sqrt{N - p - 2}\\ and back-transformed.
  The value is passed to
  [`ss_aipe_partial_r`](https://yelleknek.github.io/DMAR/reference/ss_aipe_partial_r.md),
  so the interval simulated is always the one planned.

- seed:

  `NULL` (the default) to draw from the current random number stream, or
  a single number that makes the run reproducible; a supplied seed is
  set for this call only, and the caller's random number state is
  restored when the function returns.

## Value

A `data.frame` with the columns `term` and `value`: `mean_partial_r`,
`median_partial_r`, and `sd_partial_r` (the sample partial
correlations); `mean_ci_width`, `median_ci_width`, and `sd_ci_width`
(the full widths of the intervals); `pct_ci_less_width` (the proportion
of intervals no wider than `width`); `pct_ci_miss_low`,
`pct_ci_miss_high`, and `total_type_I_error` (non-coverage of
`true_cor`, by side and in total); and the echoes `total_N` (the sample
size simulated), `p`, `true_cor`, `supposed_cor` (`NA` when
`specified_N` was given), `width`, `conf_level`, and `assurance`
(present only when an assurance was supplied). The proportions are on
the 0 to 1 scale. The attribute `interval` is `"raw"` or `"fisher_z"`,
and the printed table names the interval beneath it.

## Details

**Estimand and interval.** The population partial correlation between
*Y* and *X* given \\Z_1, \ldots, Z\_{p-1}\\, estimated from *N*
independent multivariate normal observations by the correlation of the
two sets of least squares residuals. The interval is two-sided at
`conf_level` and is the one `fisher_z` names. The raw-scale limits are
not truncated at \\-1\\ and \\1\\, so that the widths are those the
planner plans; with a large \\\|\rho\|\\ and a small sample a limit can
pass 1 in absolute value.

**What is simulated.** The function constructs a \\(p + 1) \times (p +
1)\\ population covariance matrix for *Y*, *X*, and \\Z_1, \ldots,
Z\_{p-1}\\, all standardized and uncorrelated but for a correlation of
`true_cor` between *Y* and *X*, so that the partial correlation is
`true_cor`. Under normality the distribution of the sample partial
correlation depends on the population only through that value, \\p\\,
and *N*, so this population stands for every normal population with the
same partial correlation. With `supposed_cor`, *N* is the plan of
[`ss_aipe_partial_r`](https://yelleknek.github.io/DMAR/reference/ss_aipe_partial_r.md)
for that value (with `assurance` and `fisher_z` passed on); with
`specified_N` it is that number. Each of `G` replications draws *N*
rows, computes the sample partial correlation, and forms the interval.
The \\p - 1\\ variables partialed out are the count
[`ss_aipe_partial_r`](https://yelleknek.github.io/DMAR/reference/ss_aipe_partial_r.md)
plans for and
[`var_partial_r`](https://yelleknek.github.io/DMAR/reference/var_partial_r.md)
uses. A sample with nothing partialed out is the case of
[`ss_aipe_r_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_r_sensitivity.md).

**Reading the result.** `pct_ci_less_width` is the proportion of the `G`
intervals no wider than `width`: about one half for a plan without
`assurance` whose planning value is correct, and at least `assurance`,
within Monte Carlo error, for a plan with it. The Monte Carlo standard
error of a proportion \\A\\ is \\\sqrt{A (1 - A) / G}\\.
`pct_ci_miss_low` counts intervals that lay entirely above `true_cor`
and `pct_ci_miss_high` those entirely below; their sum is
`total_type_I_error`. Every replication has an interval, so all rows are
over `G` replications.

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

[`ss_aipe_partial_r`](https://yelleknek.github.io/DMAR/reference/ss_aipe_partial_r.md),
[`ss_aipe_semipartial_r_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_semipartial_r_sensitivity.md)

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
# A plan for a full width of .40 at a planning value of .40 that is
# correct. G = 50 keeps the example quick; a reported check deserves
# the default G = 1000 or more.
ss_aipe_partial_r_sensitivity(
  true_cor = 0.40, supposed_cor = 0.40, p = 4, width = 0.40,
  G = 50, seed = 113
)
#>  term               value
#>  mean_partial_r     0.400
#>  median_partial_r   0.409
#>  sd_partial_r       0.106
#>  mean_ci_width      0.394
#>  median_ci_width    0.396
#>  sd_ci_width        0.040
#>  pct_ci_less_width  0.540
#>  pct_ci_miss_low    0.020
#>  pct_ci_miss_high   0.020
#>  total_type_I_error 0.040
#>  total_N            72   
#>  p                  4    
#>  true_cor           0.400
#>  supposed_cor       0.400
#>  width              0.400
#>  conf_level         0.950
#> 
#> Replications: 50
#> 
#> Confidence level: 95%
#> 
#> Interval: raw scale, r plus or minus z (1 - r^2) / sqrt(N - p)

# The same plan when the population value is smaller than planned: the
# intervals are wider and fewer meet the target.
ss_aipe_partial_r_sensitivity(
  true_cor = 0.20, supposed_cor = 0.40, p = 4, width = 0.40,
  G = 50, seed = 113
)
#>  term               value
#>  mean_partial_r     0.201
#>  median_partial_r   0.209
#>  sd_partial_r       0.121
#>  mean_ci_width      0.449
#>  median_ci_width    0.455
#>  sd_ci_width        0.024
#>  pct_ci_less_width  0.020
#>  pct_ci_miss_low    0.020
#>  pct_ci_miss_high   0.040
#>  total_type_I_error 0.060
#>  total_N            72   
#>  p                  4    
#>  true_cor           0.200
#>  supposed_cor       0.400
#>  width              0.400
#>  conf_level         0.950
#> 
#> Replications: 50
#> 
#> Confidence level: 95%
#> 
#> Interval: raw scale, r plus or minus z (1 - r^2) / sqrt(N - p)

# Fisher's Z interval, planned and simulated.
ss_aipe_partial_r_sensitivity(
  true_cor = 0.40, supposed_cor = 0.40, p = 4, width = 0.40,
  G = 50, fisher_z = TRUE, seed = 113
)
#>  term               value
#>  mean_partial_r     0.398
#>  median_partial_r   0.408
#>  sd_partial_r       0.110
#>  mean_ci_width      0.394
#>  median_ci_width    0.395
#>  sd_ci_width        0.039
#>  pct_ci_less_width  0.540
#>  pct_ci_miss_low    0.020
#>  pct_ci_miss_high   0.020
#>  total_type_I_error 0.040
#>  total_N            73   
#>  p                  4    
#>  true_cor           0.400
#>  supposed_cor       0.400
#>  width              0.400
#>  conf_level         0.950
#> 
#> Replications: 50
#> 
#> Confidence level: 95%
#> 
#> Interval: Fisher's Z scale, back-transformed (standard error 1 / sqrt(N - p - 2))
```
