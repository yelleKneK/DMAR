# Sensitivity Analysis for Sample Size Planning From the AIPE Perspective for a Pearson Correlation

Quantifies how much misspecification of the population Pearson
correlation distorts an AIPE-based sample size plan. On each replication
the function draws an *n*-row sample from a bivariate normal
distribution with correlation `true_cor` and computes the sample
correlation and its Fisher's \\Z\\ CI, the interval
[`ss_aipe_r`](https://yelleknek.github.io/DMAR/reference/ss_aipe_r.md)
plans for and
[`correlations_test`](https://yelleknek.github.io/DMAR/reference/correlations_test.md)
reports. Because the back-transformed width is largest at \\\rho = 0\\
and shrinks as \\\|\rho\|\\ grows, a planning value whose magnitude
overstates the population correlation yields realized intervals wider
than planned, and the summary rows report by how much.

## Usage

``` r
ss_aipe_r_sensitivity(
  true_cor = NULL,
  supposed_cor = NULL,
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

- true_cor:

  Population Pearson correlation \\\rho\\; must lie in \\(-1, 1)\\.

- supposed_cor:

  Planning value of the correlation passed to
  [`ss_aipe_r`](https://yelleknek.github.io/DMAR/reference/ss_aipe_r.md);
  supply this or `specified_N` but not both.

- width:

  Desired full width of the CI on the correlation.

- specified_N:

  Sample size to evaluate (incompatible with `supposed_cor`).

- conf_level:

  Confidence level (default `0.95`). A level above .999999998, whose
  tail areas are below 1e-9, is refused (see "Small tail areas" in
  [`ci_nc_t`](https://yelleknek.github.io/DMAR/reference/ci_nc_t.md)).

- assurance:

  Optional assurance probability passed to
  [`ss_aipe_r`](https://yelleknek.github.io/DMAR/reference/ss_aipe_r.md).

- G:

  Number of Monte Carlo replications (default 1000).

- print_iter:

  Logical. Print iteration index per replication.

- filename:

  Optional path for a comma separated file recording every replication
  (the sample correlation, the two confidence limits, the interval
  width, and two indicators of whether the interval missed `true_cor`
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

A `data.frame` with rows for the realized correlation, the interval
width, the proportion of intervals at or below `width`, tail-specific
and overall non-coverage of `true_cor`, and the input echoes, including
`assurance` (present only when an assurance was supplied).

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

Bonett, D. G., & Wright, T. A. (2000). Sample size requirements for
estimating Pearson, Kendall and Spearman correlations. *Psychometrika,
65*(1), 23–28.
[doi:10.1007/BF02294183](https://doi.org/10.1007/BF02294183)

Kelley, K., & Maxwell, S. E. (2003). Sample size for multiple
regression: Obtaining regression coefficients that are accurate, not
simply significant. *Psychological Methods, 8*(3), 305–321.
[doi:10.1037/1082-989X.8.3.305](https://doi.org/10.1037/1082-989X.8.3.305)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.

## See also

[`ss_aipe_r`](https://yelleknek.github.io/DMAR/reference/ss_aipe_r.md),
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
# Reduced replications and a wide target interval keep this fast.
set.seed(113)
ss_aipe_r_sensitivity(
  true_cor = 0.30, supposed_cor = 0.30, width = 0.40,
  G = 50, print_iter = FALSE
)
#>  term               value
#>  mean_r             0.289
#>  median_r           0.278
#>  sd_r               0.087
#>  mean_ci_width      0.399
#>  median_ci_width    0.404
#>  sd_ci_width        0.022
#>  pct_ci_less_width  0.440
#>  pct_ci_miss_low    0.000
#>  pct_ci_miss_high   0.000
#>  total_type_I_error 0.000
#>  total_N            81   
#>  true_cor           0.300
#>  supposed_cor       0.300
#>  width              0.400
#>  conf_level         0.950
#> 
#> Replications: 50
#> 
#> Confidence level: 95%
```
