# Sensitivity Analysis for Sample Size Planning to Make the Omnibus Test of \\R^2\\ Sufficiently Powerful

Monte Carlo sensitivity analysis for the power of the omnibus *F*-test
of the squared multiple correlation coefficient (\\R^2\\). Given a
`supposed_R2` used for sample size planning and a `true_R2` that
actually obtains in the population (the two need not agree), the
function draws `G` replications, fits the regression, compares \\F\\ to
its critical value, and reports the realized empirical power and a
summary of the realized \\R^2\\ and \\F\\ distributions. The simulation
honors the same `random_predictors` / `generate_random_predictors`
crossing as
[`ss_aipe_R2_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_R2_sensitivity.md),
so the user can examine the effect of planning under one regression
model (fixed or random predictors) but actually realizing the other.

## Usage

``` r
ss_power_R2_sensitivity(
  true_R2 = NULL,
  supposed_R2 = NULL,
  desired_power = 0.85,
  p = NULL,
  alpha_level = 0.05,
  random_predictors = TRUE,
  specified_N = NULL,
  generate_random_predictors = TRUE,
  true_cor_YX = 0.3,
  true_cor_XX = 0.3,
  G = 10000,
  print_iter = TRUE,
  filename = NULL,
  seed = NULL
)
```

## Arguments

- true_R2:

  Value of the population squared multiple correlation coefficient

- supposed_R2:

  Value of the squared multiple correlation coefficient used for sample
  size planning. Either `supposed_R2` or `specified_N` must be supplied
  (not both).

- desired_power:

  Desired degree of statistical power used for planning

- p:

  Number of predictors

- alpha_level:

  Type I error rate

- random_predictors:

  Whether the sample size planning step treats predictors as random
  (`TRUE`, the default) or fixed (`FALSE`)

- specified_N:

  Sample size at which the realized power should be computed;
  alternative to specifying `supposed_R2`

- generate_random_predictors:

  Whether the internal simulation should generate predictors as random
  (`TRUE`, the default) or fixed (`FALSE`)

- true_cor_YX:

  Correlation between the dependent variable (*Y*) and each of the *X*
  variables

- true_cor_XX:

  Correlation among the *X* variables (off-diagonal of the predictor
  correlation matrix)

- G:

  Number of Monte Carlo replications

- print_iter:

  Whether to print the iteration number during the simulation

- filename:

  Optional path of a CSV file to receive the per-replication results
  (the observed \\R^2\\ `R2` and \\F\\ statistic `F_value`), overwriting
  any file already at that path; the default `NULL` writes nothing, and
  a throwaway run that wants the file should point it at
  `tempfile(fileext = ".csv")`.

- seed:

  `NULL` (the default) to draw from the current random number stream, or
  a single number that makes the call reproducible; a supplied seed is
  set for this call only, and the caller's random number state is
  restored when the function returns.

## Value

A `data.frame` with columns `term` and `value` summarizing the Monte
Carlo sensitivity analysis across `G` replications. The `term` entries
are: `total_N` (the sample size evaluated), `empirical_power` (the
proportion of replications on which \\F\\ exceeded the critical value),
`analytic_power` (computed from
[`ss_power_R2`](https://yelleknek.github.io/DMAR/reference/ss_power_R2.md)
under the same model as planning), `mean_R2` / `median_R2` / `sd_R2` and
`mean_F_value` / `median_F_value` / `sd_F_value` (summaries of the
realized \\R^2\\ and \\F\\), `F_crit` (the critical value), and the
input echoes `p`, `true_R2`, `supposed_R2` and `desired_power` (both NA
when `specified_N` was supplied instead), and `alpha_level`. The result
carries the `dmar_ss_power_sensitivity` class, so
[`tidy`](https://generics.r-lib.org/reference/tidy.html) reports the
planned sample size beside the empirical and analytic power, and
[`glance`](https://generics.r-lib.org/reference/glance.html) adds the
simulated \\R^2\\ and \\F\\ distribution beside the echoed inputs.

## Details

When `supposed_R2` equals `true_R2`, the function performs a straight
Monte Carlo evaluation of the planning procedure (no misspecification).
Pass `specified_N` to evaluate realized power at a specified sample
size; in that case `supposed_R2` must not be supplied. The crossing of
`random_predictors` (used in planning) with `generate_random_predictors`
(used in the simulation) lets the user inspect the consequences of
planning under one regression model but realizing the other. See
Gatsonis and Sampson (1989) for the comparison of fixed and random
predictor power for the omnibus test.

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

Cohen, J. (1988). *Statistical power analysis for the behavioral
sciences* (2nd ed.). Hillsdale, NJ: Lawrence Erlbaum.

Gatsonis, C., & Sampson, A. R. (1989). Multiple correlation: Exact power
and sample size calculations. *Psychological Bulletin, 106*(3), 516–524.

Kelley, K., & Maxwell, S. E. (2008). Sample size planning with
applications to multiple regression: Power and accuracy for omnibus and
targeted effects. In P. Alasuutari, L. Bickman, & J. Brannen (Eds.),
*The Sage handbook of social research methods* (pp. 166–192). Sage.

Lee, Y. S. (1971). Some results on the sampling distribution of the
multiple correlation coefficient. *Journal of the Royal Statistical
Society, Series B, 33*(1), 117–130.

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 3 on \\R^2\\ as a model comparison effect
size.)

## See also

[`ss_power_R2`](https://yelleknek.github.io/DMAR/reference/ss_power_R2.md),
[`ss_aipe_R2_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_R2_sensitivity.md),
[`ci_R2`](https://yelleknek.github.io/DMAR/reference/ci_R2.md)

[`design_consequences`](https://yelleknek.github.io/DMAR/reference/design_consequences.md)
for what a chosen design delivers: power, the Type S (sign) and Type M
(exaggeration) errors of the significance filter, and the expected
confidence interval width.

Other sample size for power:
[`power_fisher_exact()`](https://yelleknek.github.io/DMAR/reference/power_fisher_exact.md),
[`ss_aipe_mixed_effects()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mixed_effects.md),
[`ss_power_R2()`](https://yelleknek.github.io/DMAR/reference/ss_power_R2.md),
[`ss_power_R2_change()`](https://yelleknek.github.io/DMAR/reference/ss_power_R2_change.md),
[`ss_power_c()`](https://yelleknek.github.io/DMAR/reference/ss_power_c.md),
[`ss_power_c_ancova()`](https://yelleknek.github.io/DMAR/reference/ss_power_c_ancova.md),
[`ss_power_composite_ancova()`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_ancova.md),
[`ss_power_composite_ancova_2group()`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_ancova_2group.md),
[`ss_power_composite_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_anova.md),
[`ss_power_composite_factorial_ancova()`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_factorial_ancova.md),
[`ss_power_composite_factorial_ancova_het()`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_factorial_ancova_het.md),
[`ss_power_composite_factorial_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_factorial_anova.md),
[`ss_power_composite_sem()`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_sem.md),
[`ss_power_contrast()`](https://yelleknek.github.io/DMAR/reference/ss_power_contrast.md),
[`ss_power_equivalence_c()`](https://yelleknek.github.io/DMAR/reference/ss_power_equivalence_c.md),
[`ss_power_factorial_ancova()`](https://yelleknek.github.io/DMAR/reference/ss_power_factorial_ancova.md),
[`ss_power_factorial_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_factorial_anova.md),
[`ss_power_indirect_effect()`](https://yelleknek.github.io/DMAR/reference/ss_power_indirect_effect.md),
[`ss_power_mixed_effects()`](https://yelleknek.github.io/DMAR/reference/ss_power_mixed_effects.md),
[`ss_power_one_way_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_one_way_anova.md),
[`ss_power_pcm()`](https://yelleknek.github.io/DMAR/reference/ss_power_pcm.md),
[`ss_power_r()`](https://yelleknek.github.io/DMAR/reference/ss_power_r.md),
[`ss_power_rc()`](https://yelleknek.github.io/DMAR/reference/ss_power_rc.md),
[`ss_power_reg_coef()`](https://yelleknek.github.io/DMAR/reference/ss_power_reg_coef.md),
[`ss_power_reg_coef_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_power_reg_coef_sensitivity.md),
[`ss_power_rm_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_rm_anova.md),
[`ss_power_sc()`](https://yelleknek.github.io/DMAR/reference/ss_power_sc.md),
[`ss_power_sem()`](https://yelleknek.github.io/DMAR/reference/ss_power_sem.md),
[`ss_power_smd()`](https://yelleknek.github.io/DMAR/reference/ss_power_smd.md),
[`ss_power_split_plot_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_split_plot_anova.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
set.seed(113)
# Realized power when planning under the fixed-predictor model but the
# data are actually generated with random predictors. G = 200 keeps the
# example quick; a reported analysis deserves the default G of 10000.
ss_power_R2_sensitivity(true_R2 = 0.30, supposed_R2 = 0.30,
                        desired_power = 0.80, p = 5,
                        random_predictors = FALSE,
                        generate_random_predictors = TRUE,
                        G = 200, print_iter = FALSE)
#>  term            value
#>  total_N         36   
#>  empirical_power 0.795
#>  analytic_power  0.802
#>  mean_R2         0.400
#>  median_R2       0.397
#>  sd_R2           0.118
#>  mean_F_value    4.430
#>  median_F_value  3.953
#>  sd_F_value      2.255
#>  F_crit          2.534
#>  p               5    
#>  true_R2         0.300
#>  supposed_R2     0.300
#>  desired_power   0.800
#>  alpha_level     0.050
#> 
#> Replications: 200
```
