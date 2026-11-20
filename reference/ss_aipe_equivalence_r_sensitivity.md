# Sensitivity Analysis for Sample Size Planning From the AIPE Perspective for an Equivalence-Test Correlation

Quantifies how much misspecification of the population correlation
distorts an AIPE-based sample size plan for the two-one-sided-tests
(TOST) confidence interval on the Pearson correlation. On each
replication the function simulates *N* bivariate normal pairs with
population correlation `true_cor`, computes the sample correlation and
its Fisher's \\Z\\ confidence interval via
[`ci_r`](https://yelleknek.github.io/DMAR/reference/ci_correlation.md),
and summarizes the realized widths and the proportion of replications in
which the computed interval falls entirely inside (`equivalent`) the
specified equivalence bounds.

## Usage

``` r
ss_aipe_equivalence_r_sensitivity(
  true_cor = 0,
  supposed_cor = NULL,
  width,
  cor_lower = NULL,
  cor_upper = NULL,
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

- true_cor:

  Population correlation (the data generating value). Defaults to `0`
  (no association, the exact equivalence case).

- supposed_cor:

  Planning value of the population correlation passed to
  [`ss_aipe_equivalence_r`](https://yelleknek.github.io/DMAR/reference/ss_aipe_equivalence_r.md);
  supply this or `specified_N` but not both.

- width:

  Desired full width of the two-sided CI on the correlation.

- cor_lower, cor_upper:

  Equivalence bounds on the correlation, as positive magnitudes with the
  same meaning as in
  [`equivalence_r`](https://yelleknek.github.io/DMAR/reference/equivalence_r.md):
  the region is \\(-\rho_L, +\rho_U)\\, where \\\rho_L\\ is `cor_lower`
  and \\\rho_U\\ is `cor_upper`. `cor_upper` is required; `cor_lower`
  defaults to `cor_upper` (a symmetric region). The simulator records
  whether the realized CI falls entirely inside the region.

- specified_N:

  Sample size to evaluate.

- conf_level:

  Confidence level (default `0.95`). A level above .999999998, whose
  tail areas are below 1e-9, is refused (see "Small tail areas" in
  [`ci_nc_t`](https://yelleknek.github.io/DMAR/reference/ci_nc_t.md)).

- assurance:

  Optional assurance probability.

- G:

  Number of Monte Carlo replications.

- print_iter:

  Logical.

- filename:

  Optional path for a comma separated file recording every replication
  (the sample correlation, the two confidence limits, the interval
  width, whether the interval fell inside the equivalence region, and
  two indicators of whether the interval missed `true_cor` below or
  above): nothing is written when `filename` is `NULL` (the default), a
  new file with a header row is created otherwise, an existing file at
  that path is appended to, and a throwaway run should point it at
  `tempfile(fileext = ".csv")`.

- seed:

  `NULL` (the default) to draw from the current random number stream, or
  a single number that makes the call reproducible; a supplied seed is
  set for this call only, and the caller's random number state is
  restored when the function returns.

## Value

A `data.frame` with rows for mean / median / SD of the realized
correlation and CI width, the proportion of intervals at or below
`width`, tail-specific and overall non-coverage of `true_cor`, the
proportion of intervals classified as `equivalent` (CI fully inside the
bounds), and the input echoes, including `assurance` (present only when
an assurance was supplied).

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

Counsell, A., & Cribbie, R. A. (2015). Equivalence tests for comparing
correlation and regression coefficients. *British Journal of
Mathematical and Statistical Psychology, 68*(2), 292–309.
[doi:10.1111/bmsp.12045](https://doi.org/10.1111/bmsp.12045)

Goertzen, J. R., & Cribbie, R. A. (2010). Detecting a lack of
association: An equivalence testing approach. *British Journal of
Mathematical and Statistical Psychology, 63*(3), 527–537.
[doi:10.1348/000711009X475853](https://doi.org/10.1348/000711009X475853)

## See also

[`ss_aipe_equivalence_r`](https://yelleknek.github.io/DMAR/reference/ss_aipe_equivalence_r.md),
[`equivalence_r`](https://yelleknek.github.io/DMAR/reference/equivalence_r.md),
[`ss_aipe_r_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_r_sensitivity.md),
[`ss_aipe_equivalence_smd_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_equivalence_smd_sensitivity.md)

Other AIPE sample size planning:
[`ss_aipe_c_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_c_sensitivity.md),
[`ss_aipe_cliff_delta()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cliff_delta.md),
[`ss_aipe_cliff_delta_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cliff_delta_sensitivity.md),
[`ss_aipe_composite_sem()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_composite_sem.md),
[`ss_aipe_crossed()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_crossed.md),
[`ss_aipe_crossed_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_crossed_sensitivity.md),
[`ss_aipe_equivalence_r()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_equivalence_r.md),
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
# Reduced Monte Carlo sweep (small G) for a fast, illustrative run.
set.seed(113)
ss_aipe_equivalence_r_sensitivity(
  true_cor     = 0.0,
  supposed_cor = 0.0,
  width        = 0.30,
  cor_upper    = 0.20,
  G = 50, print_iter = FALSE
)
#>  term               value 
#>  mean_r             0.009 
#>  median_r           0.004 
#>  sd_r               0.076 
#>  mean_ci_width      0.298 
#>  median_ci_width    0.298 
#>  sd_ci_width        0.002 
#>  pct_ci_less_width  1.000 
#>  pct_equivalent     0.500 
#>  pct_ci_miss_low    0.020 
#>  pct_ci_miss_high   0.000 
#>  total_type_I_error 0.020 
#>  total_N            172   
#>  true_cor           0     
#>  supposed_cor       0     
#>  width              0.300 
#>  conf_level         0.950 
#>  cor_lower          -0.200
#>  cor_upper          0.200 
#> 
#> Replications: 50
#> 
#> Confidence level: 95%
```
