# Sensitivity Analysis for Sample Size Planning for a Mean From the Accuracy in Parameter Estimation (AIPE) Perspective

Evaluates a sample size plan for the confidence interval for a mean by
Monte Carlo simulation. Each of `G` replications draws a sample of the
planned (or a specified) size from the data generating population, forms
the *t* interval, and records its width and whether it missed the true
mean. The summary shows what a plan from
[`ss_aipe_mean`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mean.md)
delivers when the planning value of the standard deviation is wrong,
when the population is not normal, or both. The population is either
normal, with mean `true_mean` and standard deviation `true_sd`, or the
empirical distribution of a vector of observed scores, resampled with
replacement, which is how a plan is checked against the skewness, tails,
and discreteness of real data such as a pilot's difference scores.

## Usage

``` r
ss_aipe_mean_sensitivity(
  true_sd = NULL,
  supposed_sd = NULL,
  specified_N = NULL,
  width,
  true_mean = 0,
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

- true_sd:

  The data generating standard deviation of a normal population. Supply
  it, or `population`, but not both.

- supposed_sd:

  The supposed standard deviation, the planning value passed to
  [`ss_aipe_mean`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mean.md)
  to plan the sample size (not an estimate computed from data). Supply
  this or `specified_N`, but not both.

- specified_N:

  A sample size to evaluate in place of a planned one; the planner is
  then bypassed. A single whole number of at least 2.

- width:

  The desired full width of the confidence interval: the planning
  target, and the bound against which `pct_ci_less_width` counts the
  realized widths.

- true_mean:

  The data generating mean of a normal population (default 0). The width
  does not depend on it; the misses are counted against it. Not supplied
  with `population`.

- population:

  `NULL` (the default) for a normal population, or a numeric vector of
  scores whose empirical distribution is the population: each
  replication draws its sample from these scores with replacement. The
  true mean is then their mean and the true standard deviation their
  standard deviation with divisor equal to their number, the standard
  deviation of the distribution resampled. Missing values are refused
  rather than dropped.

- conf_level:

  The confidence level of the interval (default `0.95`).

- assurance:

  `NULL` (the default) or the assurance passed to
  [`ss_aipe_mean`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mean.md)
  when it plans the sample size; when one is supplied it is also echoed,
  as the value `pct_ci_less_width` is read against.

- true_kurtosis:

  The planning value of the excess kurtosis passed to
  [`ss_aipe_mean`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mean.md)
  (default 0), on the scale
  [`kurtosis`](https://yelleknek.github.io/DMAR/reference/kurtosis.md)
  reports. It changes the plan only together with `supposed_sd` and an
  `assurance`; the function warns when a nonzero value would go unused.
  As there, it must be greater than -2; the `true_kurtosis` argument of
  [`ss_aipe_mean`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mean.md)
  says what to pass for a pilot whose
  [`kurtosis()`](https://yelleknek.github.io/DMAR/reference/kurtosis.md)
  is -2 or less.

- G:

  The number of Monte Carlo replications (default 1000).

- print_iter:

  Logical. Print the index of each replication as it runs, which is
  helpful in a long run; default `FALSE`.

- filename:

  Optional path of a CSV file to receive the per-replication results
  (the sample mean and standard deviation, the interval limits and
  width, and the two tail misses), appended when the file already exists
  and created otherwise. The default `NULL` writes nothing; a throwaway
  run that wants the file should point it at
  `tempfile(fileext = ".csv")`.

- seed:

  `NULL` (the default) to draw from the current random number stream, or
  a single number that makes the run reproducible; a supplied seed is
  set for this call only, and the caller's random number state is
  restored when the function returns.

## Value

A `data.frame` with columns `term` and `value`. The rows `mean_mean`,
`median_mean`, and `sd_mean` summarize the `G` sample means;
`mean_ci_width`, `median_ci_width`, and `sd_ci_width` summarize the
interval widths; `pct_ci_less_width` is the proportion of intervals no
wider than `width`, to be read against the assurance of the plan (or
against a value near one half for an expected width plan);
`pct_ci_miss_low` is the proportion of intervals lying entirely above
the true mean (the true mean below the lower limit) and
`pct_ci_miss_high` the proportion lying entirely below it, with
`total_type_I_error` their sum. These proportions are on the 0 to 1
scale. The input echoes are `total_N` (the sample size evaluated),
`true_sd` and `true_mean` (computed from `population` when one is
supplied), `supposed_sd` (NA when `specified_N` was supplied instead),
`width`, `conf_level`, and, when an assurance was supplied, `assurance`.
The attribute `population` records whether the population was `"normal"`
or `"resampled"`.

## Details

**What each replication does.** A sample of `total_N` scores is drawn,
and the interval \\\bar{y} \pm t\_{1-\alpha/2,\\N-1}\\ S/\sqrt{N}\\ is
formed, the one [`t.test()`](https://rdrr.io/r/stats/t.test.html)
reports and the one the plan was made for. Every replication returns an
interval, so none is dropped. A resample whose scores are all equal,
which a discrete population allows at a small *N*, gives an interval of
zero width, and it is counted as it falls.

**What to read.** The width rows against `width` show the effect of a
planning value of the standard deviation that differs from the true one:
the width is proportional to \\S\\, so a population standard deviation
larger than the planning value widens every interval.
`pct_ci_less_width` against the plan's assurance shows whether the plan
delivered the probability it promised; with a resampled population it is
the check, free of the two-moment approximation, on the planner's
`true_kurtosis` adjustment. The two miss rows, each against
\\\alpha/2\\, show the coverage of the *t* interval in the population at
hand: under normality both tails miss at the nominal rate, and a skewed
population misses more often in one tail than in the other, which the
plan does not address.

**A pilot as the population.** Passing a pilot sample's scores as
`population` and its standard deviation as `supposed_sd` evaluates the
plan in a population with exactly the pilot's shape, so the planning
value is correct (up to the divisor) and what the simulation isolates is
the effect of the shape. The pilot is itself a sample, so its shape, and
above all its tails, carries sampling error of its own.

Sensitivity analysis by simulation, in which the population value that
generates the data and the value the plan assumed are set separately,
was introduced for the standardized mean difference by Kelley and Rausch
(2006); the `ss_aipe_*_sensitivity()` family carries it to the other
AIPE planners, with a common set of rows.

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

## See also

[`ss_aipe_mean`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mean.md),
[`kurtosis`](https://yelleknek.github.io/DMAR/reference/kurtosis.md),
[`ss_aipe_c_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_c_sensitivity.md)

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
# The pilot of the ss_aipe_mean() page: the IQ gains of the 246 children
# in the control group of the pygmalion data.
gain <- pygmalion$iq_gain[pygmalion$treatment == "Control"]

# A plan for 90 percent assurance made under normal theory, evaluated by
# resampling the pilot's own gains, whose tails are heavier than the
# normal's: the share of intervals no wider than 4 points (the row
# pct_ci_less_width) falls short of .90.
ss_aipe_mean_sensitivity(supposed_sd = sd(gain), width = 4,
                         population = gain, assurance = 0.90, seed = 113)
#>  term               value 
#>  mean_mean          8.660 
#>  median_mean        8.687 
#>  sd_mean            0.928 
#>  mean_ci_width      3.736 
#>  median_ci_width    3.733 
#>  sd_ci_width        0.234 
#>  pct_ci_less_width  0.870 
#>  pct_ci_miss_low    0.021 
#>  pct_ci_miss_high   0.034 
#>  total_type_I_error 0.055 
#>  total_N            195   
#>  true_sd            13.270
#>  true_mean          8.663 
#>  supposed_sd        13.297
#>  width              4     
#>  conf_level         0.950 
#>  assurance          0.900 
#> 
#> Replications: 1000
#> 
#> Confidence level: 95%

# The plan made with the pilot's excess kurtosis is six children larger,
# and the share comes close to .90.
ss_aipe_mean_sensitivity(supposed_sd = sd(gain), width = 4,
                         population = gain, assurance = 0.90,
                         true_kurtosis = kurtosis(gain), seed = 113)
#>  term               value 
#>  mean_mean          8.657 
#>  median_mean        8.669 
#>  sd_mean            0.941 
#>  mean_ci_width      3.680 
#>  median_ci_width    3.684 
#>  sd_ci_width        0.230 
#>  pct_ci_less_width  0.913 
#>  pct_ci_miss_low    0.025 
#>  pct_ci_miss_high   0.026 
#>  total_type_I_error 0.051 
#>  total_N            201   
#>  true_sd            13.270
#>  true_mean          8.663 
#>  supposed_sd        13.297
#>  width              4     
#>  conf_level         0.950 
#>  assurance          0.900 
#> 
#> Replications: 1000
#> 
#> Confidence level: 95%

# A normal population whose standard deviation is 15 when the plan
# assumed 13.3: the realized intervals run wider than the target.
ss_aipe_mean_sensitivity(true_sd = 15, supposed_sd = 13.3, width = 4,
                         seed = 113)
#>  term               value 
#>  mean_mean          -0.003
#>  median_mean        0.023 
#>  sd_mean            1.137 
#>  mean_ci_width      4.512 
#>  median_ci_width    4.511 
#>  sd_ci_width        0.241 
#>  pct_ci_less_width  0.013 
#>  pct_ci_miss_low    0.031 
#>  pct_ci_miss_high   0.021 
#>  total_type_I_error 0.052 
#>  total_N            172   
#>  true_sd            15    
#>  true_mean          0     
#>  supposed_sd        13.300
#>  width              4     
#>  conf_level         0.950 
#> 
#> Replications: 1000
#> 
#> Confidence level: 95%
```
