# Sensitivity Analysis for Sample Size Planning From the AIPE Perspective for Cliff's Delta

Quantifies how much misspecification of the population Cliff's delta
(\\\delta = \Pr(X \> Y) - \Pr(X \< Y)\\) distorts an AIPE-based sample
size plan. On each replication the function simulates two independent
samples whose population Cliff's delta equals `true_delta`, computes the
sample
[`cliff_delta`](https://yelleknek.github.io/DMAR/reference/cliff_delta.md)
and its CI, and summarizes the realized widths and coverage.

**Replications without an interval.** When the two simulated groups are
completely separated (every score of one above every score of the
other), or every score is equal, the estimated variance of Cliff's delta
is 0 and
[`cliff_delta`](https://yelleknek.github.io/DMAR/reference/cliff_delta.md)
reports both limits equal to the estimate, with a warning that they are
not a confidence interval. Such a replication has no interval here. Its
estimate still enters `mean_cliff_delta`, `median_cliff_delta`, and
`sd_cliff_delta`. The width, width-target, and coverage rows
(`mean_ci_width`, `median_ci_width`, `sd_ci_width`, `pct_ci_less_width`,
`pct_ci_miss_low`, `pct_ci_miss_high`, `total_type_I_error`) describe
the replications that have an interval, and are `NA` when none has;
`pct_ci_less_width` is the number of intervals no wider than `width`
divided by the number of replications with an interval. The warning of
each replication is not passed on. One warning of class
`dmar_replications_without_interval`, with the fields
`n_without_interval` and `G`, states how many there were, and the
`simulation` attribute records the counts. Separation is common when
`true_delta` is large in absolute value and the groups are small: with
`true_delta = 0.9` and six per group, about four replications in ten are
separated.

**Data generating mechanism.** The simulator draws each sample from a
normal distribution and chooses the mean shift so that the implied
Cliff's delta equals `true_delta`. For normal samples \\\delta = 2
\Phi(\Delta/\sqrt{2}) - 1\\ where \\\Delta\\ is the standardized mean
difference, so the simulator sets \\\Delta = \sqrt{2} \cdot
\Phi^{-1}((1 + \delta)/2)\\.

## Usage

``` r
ss_aipe_cliff_delta_sensitivity(
  true_delta = NULL,
  supposed_delta = NULL,
  ratio = 1,
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

- true_delta:

  Population Cliff's delta (the data generating value); in \\(-1, 1)\\.

- supposed_delta:

  The supposed population Cliff's delta, the planning value passed to
  [`ss_aipe_cliff_delta`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cliff_delta.md)
  (not an estimate computed from data); supply this or `specified_N` but
  not both.

- ratio:

  Allocation ratio \\n_1 / n_2\\ (default 1).

- width:

  Desired full width of the CI on Cliff's delta.

- specified_N:

  Total sample size to evaluate (split per `ratio`).

- conf_level:

  Confidence level (default `0.95`).

- assurance:

  Optional assurance probability.

- G:

  Number of Monte Carlo replications.

- print_iter:

  Logical.

- filename:

  Optional path for a comma separated file recording every replication
  (the sample Cliff's delta, the two confidence limits, the interval
  width, and two indicators of whether the interval missed `true_delta`
  below or above): nothing is written when `filename` is `NULL` (the
  default), a new file with a header row is created otherwise, an
  existing file at that path is appended to, and a throwaway run should
  point it at `tempfile(fileext = ".csv")`. A replication without an
  interval is recorded with its estimate and with `NA` in the other
  columns.

- seed:

  Optional integer seed. The default `NULL` draws from, and advances,
  the current random number stream; an integer makes the call
  reproducible, and the caller's random number state is restored when
  the call returns.

## Value

A `data.frame` (class `dmar_tbl`) with columns `term` and `value` and
these rows: `mean_cliff_delta`, `median_cliff_delta`, and
`sd_cliff_delta` (summaries of the `G` estimates); `mean_ci_width`,
`median_ci_width`, and `sd_ci_width` (summaries of the widths of the
replications with an interval); `pct_ci_less_width` (the proportion of
all `G` replications whose interval is no wider than `width`);
`pct_ci_miss_low`, `pct_ci_miss_high`, and `total_type_I_error` (among
the replications with an interval, the proportions whose interval lies
wholly above `true_delta`, wholly below it, and either); `n_1`, `n_2`,
and `total_N` (the sample sizes simulated); and the inputs `true_delta`,
`supposed_delta` (`NA` when `specified_N` was given), `ratio`, `width`,
`conf_level`, and `assurance` (present only when an assurance was
supplied). The `simulation` attribute records the numbers of
replications attempted, with an interval, and without one (see the
section on the record of the simulation).

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

Cliff, N. (1993). Dominance statistics: Ordinal analyses to answer
ordinal questions. *Psychological Bulletin, 114*(3), 494–509.
[doi:10.1037/0033-2909.114.3.494](https://doi.org/10.1037/0033-2909.114.3.494)

Feng, D., & Cliff, N. (2004). Monte Carlo evaluation of ordinal d with
improved confidence interval. *Journal of Modern Applied Statistical
Methods, 3*(2), 322–332.
[doi:10.22237/jmasm/1099267560](https://doi.org/10.22237/jmasm/1099267560)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.

## See also

[`ss_aipe_cliff_delta`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cliff_delta.md),
[`cliff_delta`](https://yelleknek.github.io/DMAR/reference/cliff_delta.md)

[`design_consequences`](https://yelleknek.github.io/DMAR/reference/design_consequences.md)
for what a chosen design delivers: power, the Type S (sign) and Type M
(exaggeration) errors of the significance filter, and the expected
confidence interval width.

Other AIPE sample size planning:
[`ss_aipe_c_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_c_sensitivity.md),
[`ss_aipe_cliff_delta()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cliff_delta.md),
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
# Small G keeps the Monte Carlo sweep fast; raise G for a real plan.
ss_aipe_cliff_delta_sensitivity(
  true_delta = 0.30, supposed_delta = 0.30,
  width = 0.30, G = 25, print_iter = FALSE, seed = 113
)
#>  term               value
#>  mean_cliff_delta   0.310
#>  median_cliff_delta 0.303
#>  sd_cliff_delta     0.068
#>  mean_ci_width      0.240
#>  median_ci_width    0.241
#>  sd_ci_width        0.008
#>  pct_ci_less_width  1.000
#>  pct_ci_miss_low    0.080
#>  pct_ci_miss_high   0.000
#>  total_type_I_error 0.080
#>  n_1                156  
#>  n_2                156  
#>  total_N            312  
#>  true_delta         0.300
#>  supposed_delta     0.300
#>  ratio              1    
#>  width              0.300
#>  conf_level         0.950
#> 
#> Replications: 25
#> 
#> Confidence level: 95%

# A large effect with six per group: some replications are completely
# separated and have no interval. The width and coverage rows describe
# the others, and one warning says how many were without an interval.
ss_aipe_cliff_delta_sensitivity(
  true_delta = 0.90, specified_N = 12, width = 0.40, G = 100,
  seed = 113
)
#> Warning: 43 of 100 replications had completely separated groups (or no variation), where the estimated variance of Cliff's delta is 0 and no confidence interval exists. The width, width-target, and coverage rows describe the other 57 replications.
#>  term               value
#>  mean_cliff_delta   0.914
#>  median_cliff_delta 0.944
#>  sd_cliff_delta     0.103
#>  mean_ci_width      0.959
#>  median_ci_width    0.909
#>  sd_ci_width        0.100
#>  pct_ci_less_width  0.000
#>  pct_ci_miss_low    0.000
#>  pct_ci_miss_high   0.000
#>  total_type_I_error 0.000
#>  n_1                6    
#>  n_2                6    
#>  total_N            12   
#>  true_delta         0.900
#>  supposed_delta     <NA> 
#>  ratio              1    
#>  width              0.400
#>  conf_level         0.950
#> 
#> Replications: 100 attempted, 57 with an interval, 43 without
#> 
#> Confidence level: 95%
```
