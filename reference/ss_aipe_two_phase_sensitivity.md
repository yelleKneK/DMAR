# Sensitivity Analysis for Sample Size Planning From the AIPE Perspective for a Mean Under Two-Phase Sampling

Evaluates a two-phase design by simulation: how wide the confidence
intervals from
[`ci_two_phase`](https://yelleknek.github.io/DMAR/reference/ci_two_phase.md)
turn out to be, how often they are no wider than planned, and how often
they miss the population mean of the target, when the design is carried
out in a population that may differ from the one the plan assumed. The
design is either planned by
[`ss_aipe_two_phase`](https://yelleknek.github.io/DMAR/reference/ss_aipe_two_phase.md)
from planning values, so that the distance between planning values and
population values is what is studied, or fixed by the user.

## Usage

``` r
ss_aipe_two_phase_sensitivity(
  true_sd_Y = NULL,
  true_cor = NULL,
  supposed_sd_Y = NULL,
  supposed_cor = NULL,
  specified_n_labeled = NULL,
  n_pool,
  width,
  population = NULL,
  conf_level = 0.95,
  assurance = NULL,
  true_kurtosis = 0,
  G = 1000,
  print_iter = FALSE,
  filename = NULL,
  seed = NULL
)
```

## Arguments

- true_sd_Y, true_cor:

  The population standard deviation of the target measurement and its
  correlation with the proxy, from which the data are generated when
  `population` is `NULL`.

- supposed_sd_Y, supposed_cor:

  Planning values passed to
  [`ss_aipe_two_phase`](https://yelleknek.github.io/DMAR/reference/ss_aipe_two_phase.md)
  to choose the number of labeled units; supply both, or
  `specified_n_labeled` instead.

- specified_n_labeled:

  The number of labeled units to evaluate, in place of a planned one.

- n_pool:

  The number of units in the pool, \\M\\, a whole number larger than the
  number of labeled units.

- width:

  The desired full width of the confidence interval, the planning target
  and the yardstick for `pct_ci_less_width`.

- population:

  `NULL` (default) for a bivariate normal population, or a data frame of
  observed (target, proxy) pairs with numeric columns `y` and `x`, to be
  resampled with replacement as the population; `true_sd_Y` and
  `true_cor` are then computed from it and must not be supplied.

- conf_level:

  Confidence level (default `0.95`).

- assurance, true_kurtosis:

  Passed to
  [`ss_aipe_two_phase`](https://yelleknek.github.io/DMAR/reference/ss_aipe_two_phase.md)
  when the design is planned; see that page.

- G:

  The number of Monte Carlo replications (default 1000).

- print_iter:

  Logical; print the replication number as the loop runs.

- filename:

  Optional path for a comma separated file recording every replication
  (the estimate `mean_Y`, the two confidence limits, the interval width,
  and two indicators of whether the interval missed the population mean
  below or above): nothing is written when `filename` is `NULL` (the
  default), a new file with a header row is created otherwise, an
  existing file at that path is appended to, and a throwaway run should
  point it at `tempfile(fileext = ".csv")`.

- seed:

  `NULL` (default) to draw from the current random number stream, or a
  number that makes the run reproducible without disturbing the stream
  of the calling script.

## Value

A `data.frame` (class `dmar_tbl`) with `term` and `value` columns:
`mean_mean_Y`, `median_mean_Y`, and `sd_mean_Y` for the estimates of the
target mean; the realized widths (`mean_ci_width`, `median_ci_width`,
`sd_ci_width`); the proportion of intervals no wider than `width`
(`pct_ci_less_width`); the proportions that missed the population mean
from above (`pct_ci_miss_low`, the mean below the lower limit) and from
below (`pct_ci_miss_high`) and their sum (`total_type_I_error`), all on
the 0 to 1 scale; then `n_labeled`, `total_N` (the pool size \\M\\, the
design's total number of units, which the rest of the two-phase family
calls `n_pool`; the sensitivity family's shared schema names the size of
a whole design `total_N`), `true_mean_Y` (the population target mean),
`true_sd_Y`, `true_cor`, `supposed_sd_Y` and `supposed_cor` (`NA` when
`specified_n_labeled` was given), `width`, `conf_level`, and, only when
one was supplied, `assurance` and `true_kurtosis`. The number of
replications used is in the `simulation` attribute (its element
`returned`), described below.

## Details

**The data generating mechanism.** With `population = NULL`, each
replication draws the \\n\\ labeled units as (target, proxy) pairs from
a bivariate normal distribution with standard deviations `true_sd_Y` and
1, correlation `true_cor`, and means 0, and the \\M - n\\ unlabeled
units' proxies from the same proxy distribution. The interval's width
and coverage do not depend on the two means or on the proxy's scale, so
these choices lose no generality; the population target mean the
intervals are scored against is then 0. Because the regression estimator
uses the unlabeled proxies only through their mean and their sum of
squares, the simulation draws those two from their exact joint
distribution (independent normal and scaled chi square) rather than
drawing every unit of the pool, which keeps a large pool cheap. With a
`population` data frame, each replication resamples the labeled pairs
and the unlabeled proxies with replacement from its rows, a population
that keeps whatever skew, bounds, or nonlinearity the observed pairs
have; the target mean scored against is then the mean of its `y` column.
In the evaluation of a system whose outputs an automated judge labels,
for example, a pilot set of outputs carrying both a human assessor's
label and the judge's makes such a population.

**What is evaluated.** Each replication computes the regression interval
of
[`ci_two_phase`](https://yelleknek.github.io/DMAR/reference/ci_two_phase.md)
(the interval the planner plans for) and records its width and whether
it missed the population target mean below or above. A replication whose
labeled units all share one proxy value has no regression interval, and
one whose labeled units all share one target value would have an
interval of zero width, which
[`ci_two_phase`](https://yelleknek.github.io/DMAR/reference/ci_two_phase.md)
refuses. Either can happen only with a `population` data frame, chiefly
one whose proxy takes few distinct values or whose target is binary with
a proportion near 0 or 1. Such replications are dropped, with a single
warning that counts each kind, and the number kept is
`attr(x, "simulation")$returned` (see the section on the record of the
simulation). Many dropped replications say that the design has too few
labeled units for that population, and the summaries, which describe the
replications that produced an interval, then flatter the design.

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
narrow confidence intervals. *Psychological Methods, 11*(4), 363–385.
[doi:10.1037/1082-989X.11.4.363](https://doi.org/10.1037/1082-989X.11.4.363)

Neyman, J. (1938). Contribution to the theory of sampling human
populations. *Journal of the American Statistical Association, 33*(201),
101–116.
[doi:10.1080/01621459.1938.10503378](https://doi.org/10.1080/01621459.1938.10503378)

## See also

[`ss_aipe_two_phase`](https://yelleknek.github.io/DMAR/reference/ss_aipe_two_phase.md),
[`ci_two_phase`](https://yelleknek.github.io/DMAR/reference/ci_two_phase.md),
[`ss_aipe_two_phase_cost`](https://yelleknek.github.io/DMAR/reference/ss_aipe_two_phase_cost.md)

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
[`ss_aipe_semipartial_r_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_semipartial_r_sensitivity.md),
[`ss_aipe_two_phase()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_two_phase.md),
[`ss_aipe_two_phase_cost()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_two_phase_cost.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# A plan made with a correlation of 0.55, carried out where the
# correlation is only 0.45. G = 200 keeps the example quick; a reported
# evaluation deserves the default G = 1000 or more.
ss_aipe_two_phase_sensitivity(true_sd_Y = 0.43, true_cor = 0.45,
                              supposed_sd_Y = 0.43, supposed_cor = 0.55,
                              n_pool = 20000, width = 0.10,
                              conf_level = 0.90, G = 200, seed = 113)
#>  term               value
#>  mean_mean_Y        0.000
#>  median_mean_Y      0.001
#>  sd_mean_Y          0.032
#>  mean_ci_width      0.107
#>  median_ci_width    0.106
#>  sd_ci_width        0.007
#>  pct_ci_less_width  0.180
#>  pct_ci_miss_low    0.045
#>  pct_ci_miss_high   0.060
#>  total_type_I_error 0.105
#>  n_labeled          142  
#>  total_N            20000
#>  true_mean_Y        0    
#>  true_sd_Y          0.430
#>  true_cor           0.450
#>  supposed_sd_Y      0.430
#>  supposed_cor       0.550
#>  width              0.100
#>  conf_level         0.900
#> 
#> Replications: 200
#> 
#> Confidence level: 90%

# A fixed design, 150 labeled units in a pool of 5,000, evaluated in a
# population of binary (target, proxy) pairs resampled from a pilot of
# 300 units: how often is the 90% interval narrower than 0.085?
set.seed(113)
proxy  <- rbinom(300, 1, 0.6)
target <- ifelse(runif(300) < 0.8, proxy, rbinom(300, 1, 0.6))
ss_aipe_two_phase_sensitivity(specified_n_labeled = 150, n_pool = 5000,
                              width = 0.085, conf_level = 0.90,
                              population = data.frame(y = target, x = proxy),
                              G = 200, seed = 113)
#>  term               value
#>  mean_mean_Y        0.565
#>  median_mean_Y      0.565
#>  sd_mean_Y          0.024
#>  mean_ci_width      0.084
#>  median_ci_width    0.084
#>  sd_ci_width        0.009
#>  pct_ci_less_width  0.520
#>  pct_ci_miss_low    0.035
#>  pct_ci_miss_high   0.035
#>  total_type_I_error 0.070
#>  n_labeled          150  
#>  total_N            5000 
#>  true_mean_Y        0.567
#>  true_sd_Y          0.496
#>  true_cor           0.789
#>  supposed_sd_Y      <NA> 
#>  supposed_cor       <NA> 
#>  width              0.085
#>  conf_level         0.900
#> 
#> Replications: 200
#> 
#> Confidence level: 90%
```
