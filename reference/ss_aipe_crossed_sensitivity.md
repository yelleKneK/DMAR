# Sensitivity Analysis for Sample Size Planning From the AIPE Perspective for the Mean of a Crossed Design of Units and Raters

Shows by Monte Carlo simulation what a plan from
[`ss_aipe_crossed`](https://yelleknek.github.io/DMAR/reference/ss_aipe_crossed.md)
delivers when the true variance components differ from the planning
values. The design is either planned from the supposed components, as a
researcher would plan it, or fixed by the user; complete arrays of that
size are then drawn from the population with the true components, the
interval of
[`ci_crossed`](https://yelleknek.github.io/DMAR/reference/ci_crossed.md)
is formed from each, and the realized widths and the coverage of the
true mean are summarized. A planning value of the rater component that
is too small is the costly mistake in this design, since that component
sets the floor below which no number of units can bring the interval.

## Usage

``` r
ss_aipe_crossed_sensitivity(
  true_var_unit,
  true_var_rater,
  true_var_residual,
  supposed_var_unit = NULL,
  supposed_var_rater = NULL,
  supposed_var_residual = NULL,
  specified_n_units = NULL,
  specified_n_raters = NULL,
  width,
  conf_level = 0.95,
  assurance = NULL,
  n_raters = NULL,
  unit_cost = 0,
  rater_cost = 0,
  judgment_cost = 1,
  quantile = c("t", "normal"),
  max_units = 1e+05,
  max_raters = 500,
  true_grand_mean = 0,
  G = 1000,
  print_iter = FALSE,
  filename = NULL,
  seed = NULL
)
```

## Arguments

- true_var_unit, true_var_rater, true_var_residual:

  The data generating variance components of the unit effects, the rater
  effects, and the residuals (see
  [`ss_aipe_crossed`](https://yelleknek.github.io/DMAR/reference/ss_aipe_crossed.md)).
  The first two must be nonnegative and the third positive.

- supposed_var_unit, supposed_var_rater, supposed_var_residual:

  The supposed variance components, the planning values the design is
  planned from (not estimates computed from data), all three together.
  Supply these or `specified_n_units` and `specified_n_raters`, but not
  both.

- specified_n_units, specified_n_raters:

  A design to evaluate (whole numbers of at least 2, both together) in
  place of a planned one; the planner is then bypassed.

- width:

  Desired full width of the confidence interval for the mean, against
  which the realized widths are judged.

- conf_level:

  Confidence level of the interval (default `0.95`).

- assurance:

  Optional assurance, strictly between 0 and 1, passed to
  [`ss_aipe_crossed`](https://yelleknek.github.io/DMAR/reference/ss_aipe_crossed.md)
  when the design is planned; echoed in the output when supplied, with a
  specified design as well.

- n_raters, unit_cost, rater_cost, judgment_cost, quantile, max_units,
  max_raters:

  Passed to
  [`ss_aipe_crossed`](https://yelleknek.github.io/DMAR/reference/ss_aipe_crossed.md)
  when the design is planned, and unused with a specified design.

- true_grand_mean:

  The data generating mean (default `0`). The widths and the coverage do
  not depend on it; it sets the scale of the reported estimates.

- G:

  Number of Monte Carlo replications (default `1000`).

- print_iter:

  Logical. Print the replication index as the simulation runs; default
  `FALSE`.

- filename:

  Optional path of a CSV file to receive the per-replication results
  (the grand mean, the two limits, the width, the degrees of freedom,
  and the two tail misses), appended when the file already holds a table
  and created otherwise; the default `NULL` writes nothing, and a
  throwaway run that wants the file should point it at
  `tempfile(fileext = ".csv")`.

- seed:

  Optional integer seed. The default `NULL` draws from, and advances,
  the current random number stream; an integer makes the call
  reproducible, and the caller's random number state is restored when
  the call returns.

## Value

A `data.frame` (class `dmar_tbl`) with rows `mean_grand_mean`,
`median_grand_mean`, and `sd_grand_mean` (summaries of the \\G\\
estimates); `mean_ci_width`, `median_ci_width`, and `sd_ci_width` (of
the realized widths); `pct_ci_less_width` (the proportion of intervals
no wider than `width`); `pct_ci_miss_low` (the proportion lying entirely
above the true mean), `pct_ci_miss_high` (entirely below it), and
`total_type_I_error` (their sum); `mean_df` (the mean Satterthwaite
degrees of freedom); and the echoes `n_units` and `n_raters` (the design
evaluated), the three `true_var_*` and three `supposed_var_*` components
(the latter `NA` for a specified design), `true_grand_mean`, `width`,
`conf_level`, and `assurance` when one was supplied. The proportions are
on the 0 to 1 scale.

## Details

Each replication draws \\n_u\\ unit effects, \\n_r\\ rater effects, and
\\n_u n_r\\ residuals from normal distributions with the true variances,
forms the complete array, and computes the interval exactly as
[`ci_crossed`](https://yelleknek.github.io/DMAR/reference/ci_crossed.md)
does: the ANOVA estimates of the components, negative ones truncated at
zero, plugged into the variance of the grand mean, with Satterthwaite
(1946) degrees of freedom. The summary separates two sources of error:
the planning values (compare a well specified run with a misspecified
one) and the approximations of the plan itself (a well specified run
whose `pct_ci_less_width` is well below one half, or below the
assurance, or whose coverage is far from `conf_level`; a plan without an
assurance meets its width somewhat more than half the time, as
[`ss_aipe_crossed`](https://yelleknek.github.io/DMAR/reference/ss_aipe_crossed.md)
explains). A plan with `quantile = "normal"` shows the second kind,
since its target ignores the *t* quantile the interval uses; with a
rater component at or near zero it can fall far short.

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

Satterthwaite, F. E. (1946). An approximate distribution of estimates of
variance components. *Biometrics Bulletin, 2*(6), 110–114.
[doi:10.2307/3002019](https://doi.org/10.2307/3002019)

## See also

[`ss_aipe_crossed`](https://yelleknek.github.io/DMAR/reference/ss_aipe_crossed.md),
[`ci_crossed`](https://yelleknek.github.io/DMAR/reference/ci_crossed.md)

Other AIPE sample size planning:
[`ss_aipe_c_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_c_sensitivity.md),
[`ss_aipe_cliff_delta()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cliff_delta.md),
[`ss_aipe_cliff_delta_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cliff_delta_sensitivity.md),
[`ss_aipe_composite_sem()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_composite_sem.md),
[`ss_aipe_crossed()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_crossed.md),
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
# G = 200 keeps the example quick; a reported study deserves the
# default G = 1000 or more. The first plan on the help page of
# ss_aipe_crossed(), evaluated when its planning values are right:
ss_aipe_crossed_sensitivity(
  true_var_unit = 0.50, true_var_rater = 0.10, true_var_residual = 0.40,
  supposed_var_unit = 0.50, supposed_var_rater = 0.10,
  supposed_var_residual = 0.40,
  width = 0.40, G = 200, seed = 113)
#>  term                  value 
#>  mean_grand_mean       -0.015
#>  median_grand_mean     -0.003
#>  sd_grand_mean         0.103 
#>  mean_ci_width         0.396 
#>  median_ci_width       0.395 
#>  sd_ci_width           0.033 
#>  pct_ci_less_width     0.555 
#>  pct_ci_miss_low       0.010 
#>  pct_ci_miss_high      0.040 
#>  total_type_I_error    0.050 
#>  mean_df               75.376
#>  n_units               94    
#>  n_raters              22    
#>  true_var_unit         0.500 
#>  true_var_rater        0.100 
#>  true_var_residual     0.400 
#>  supposed_var_unit     0.500 
#>  supposed_var_rater    0.100 
#>  supposed_var_residual 0.400 
#>  true_grand_mean       0     
#>  width                 0.400 
#>  conf_level            0.950 
#> 
#> Replications: 200
#> 
#> Confidence level: 95%

# The same plan when the raters in fact differ twice as much as
# planned: nearly every interval is now wider than the target, since
# the design has too few raters for the larger rater component.
ss_aipe_crossed_sensitivity(
  true_var_unit = 0.50, true_var_rater = 0.20, true_var_residual = 0.40,
  supposed_var_unit = 0.50, supposed_var_rater = 0.10,
  supposed_var_residual = 0.40,
  width = 0.40, G = 200, seed = 113)
#>  term                  value 
#>  mean_grand_mean       -0.019
#>  median_grand_mean     -0.009
#>  sd_grand_mean         0.121 
#>  mean_ci_width         0.479 
#>  median_ci_width       0.477 
#>  sd_ci_width           0.049 
#>  pct_ci_less_width     0.035 
#>  pct_ci_miss_low       0.015 
#>  pct_ci_miss_high      0.045 
#>  total_type_I_error    0.060 
#>  mean_df               51.404
#>  n_units               94    
#>  n_raters              22    
#>  true_var_unit         0.500 
#>  true_var_rater        0.200 
#>  true_var_residual     0.400 
#>  supposed_var_unit     0.500 
#>  supposed_var_rater    0.100 
#>  supposed_var_residual 0.400 
#>  true_grand_mean       0     
#>  width                 0.400 
#>  conf_level            0.950 
#> 
#> Replications: 200
#> 
#> Confidence level: 95%
```
