# Sensitivity Analysis for Sample Size Planning to Make the Omnibus Test of \\R^2\\ Sufficiently Powerful

Monte Carlo sensitivity analysis for the power of the omnibus *F*-test
of the squared multiple correlation coefficient (\\R^2\\). Given an
`estimated_R2` used for sample size planning and a `true_R2` that
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
  estimated_R2 = NULL,
  desired_power = 0.85,
  p = NULL,
  alpha_level = 0.05,
  random_predictors = TRUE,
  specified_N = NULL,
  generate_random_predictors = TRUE,
  rho_yx = 0.3,
  rho_xx = 0.3,
  G = 10000,
  print_iter = TRUE,
  save = FALSE,
  filename = "ss_power_r2_sensitivity_result.csv",
  ...
)
```

## Arguments

- true_R2:

  Value of the population squared multiple correlation coefficient

- estimated_R2:

  Value of the squared multiple correlation coefficient used for sample
  size planning. Either `estimated_R2` or `specified_N` must be supplied
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
  alternative to specifying `estimated_R2`

- generate_random_predictors:

  Whether the internal simulation should generate predictors as random
  (`TRUE`, the default) or fixed (`FALSE`)

- rho_yx:

  Correlation between the dependent variable (*Y*) and each of the *X*
  variables

- rho_xx:

  Correlation among the *X* variables (off-diagonal of the predictor
  correlation matrix)

- G:

  Number of Monte Carlo replications

- print_iter:

  Whether to print the iteration number during the simulation

- save:

  Whether to write the per-replication results to a CSV file

- filename:

  Name of the CSV file written when `save = TRUE`

- ...:

  Additional arguments forwarded to internal helpers

## Value

A `data.frame` with columns `term` and `value` summarizing the Monte
Carlo sensitivity analysis across `G` replications. Rows include the
planning sample size used (`total_N`), the empirical power
(`empirical_power`, the proportion of replications on which \\F\\
exceeded the critical value), the analytic power computed from
[`ss_power_R2`](https://yelleknek.github.io/DMAR/reference/ss_power_R2.md)
under the same model as planning (`analytic_power`), the mean / median /
SD of the realized \\R^2\\ and \\F\\, and the critical value (`f_crit`).
The result carries the `dmar_ss_power_sensitivity` class, so
[`tidy`](https://generics.r-lib.org/reference/tidy.html) reports the
planned sample size beside the empirical and analytic power, and
[`glance`](https://generics.r-lib.org/reference/glance.html) adds the
simulated \\R^2\\ and \\F\\ distribution.

## Details

When `estimated_R2` equals `true_R2`, the function performs a straight
Monte Carlo evaluation of the planning procedure (no misspecification).
Pass `specified_N` to evaluate realized power at a specified sample
size; in that case `estimated_R2` must not be supplied. The crossing of
`random_predictors` (used in planning) with `generate_random_predictors`
(used in the simulation) lets the user inspect the consequences of
planning under one regression model but realizing the other. See
Gatsonis and Sampson (1989) for the comparison of fixed and random
predictor power for the omnibus test.

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
[`ss_aipe_tost_smd()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_tost_smd.md),
[`ss_power_R2()`](https://yelleknek.github.io/DMAR/reference/ss_power_R2.md),
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
# data are actually generated with random predictors. G is small here
# for illustration; use G = 10,000 in practice.
ss_power_R2_sensitivity(true_R2 = 0.30, estimated_R2 = 0.30,
                        desired_power = 0.80, p = 5,
                        random_predictors = FALSE,
                        generate_random_predictors = TRUE,
                        G = 200, print_iter = FALSE)
#>  term            value
#>  total_N         36   
#>  empirical_power 0.77 
#>  analytic_power  0.802
#>  mean_r2         0.396
#>  median_r2       0.386
#>  sd_r2           0.122
#>  mean_f          4.41 
#>  median_f        3.78 
#>  sd_f            2.51 
#>  f_crit          2.53 
```
