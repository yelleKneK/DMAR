# Sensitivity Analysis for the Power of a Targeted Regression Coefficient

Monte Carlo sensitivity analysis for the statistical power of the
*t*-test of a targeted regression coefficient. Given a planned
(`supposed_*`) covariance structure and a true (`true_*`) covariance
structure, the function draws `G` replications, fits the multiple
regression, and reports the empirical proportion of replications on
which the *t*-test of the targeted coefficient rejects, together with
the realized distribution of the estimated unstandardized coefficient
\\\hat b_j\\, its standard error, and the test statistic.
`ss_power_reg_coef_sensitivity()` is the power-oriented sibling of
[`ss_aipe_reg_coef_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reg_coef_sensitivity.md)
(which is CI-width oriented).

## Usage

``` r
ss_power_reg_coef_sensitivity(
  true_var_Y = NULL,
  true_cov_YX = NULL,
  true_cov_XX = NULL,
  supposed_var_Y = NULL,
  supposed_cov_YX = NULL,
  supposed_cov_XX = NULL,
  specified_N = NULL,
  which_predictor = 1,
  desired_power = 0.85,
  alpha_level = 0.05,
  directional = FALSE,
  G = 1000,
  print_iter = TRUE,
  filename = NULL,
  seed = NULL
)
```

## Arguments

- true_var_Y:

  Population variance of the dependent variable (*Y*)

- true_cov_YX:

  Population covariance vector between the `p` predictor variables and
  the dependent variable (*Y*)

- true_cov_XX:

  Population covariance matrix of the `p` predictor variables

- supposed_var_Y:

  Supposed variance of the dependent variable (*Y*), the planning value
  used in sample size planning. Defaults to `true_var_Y`.

- supposed_cov_YX:

  Supposed covariance vector between the predictor variables and the
  dependent variable, the planning values used in sample size planning.
  Defaults to `true_cov_YX`.

- supposed_cov_XX:

  Supposed covariance matrix of the predictor variables, the planning
  values used in sample size planning. Defaults to `true_cov_XX`.

- specified_N:

  Directly specified sample size; if supplied, sample size planning is
  skipped.

- which_predictor:

  Index identifying which of the `p` predictors is the targeted
  predictor for the power test.

- desired_power:

  Desired degree of statistical power used for planning

- alpha_level:

  Type I error rate

- directional:

  Whether a one-sided or two-sided test is used

- G:

  Number of Monte Carlo replications

- print_iter:

  Whether to print the iteration number during the simulation

- filename:

  Optional path of a CSV file to receive the per-replication results
  (the coefficient estimate `reg_coef_j`, its standard error
  `se_reg_coef_j`, the *t* statistic `t_value`, and the observed \\R^2\\
  `R2`), overwriting any file already at that path; the default `NULL`
  writes nothing, and a throwaway run that wants the file should point
  it at `tempfile(fileext = ".csv")`.

- seed:

  `NULL` (the default) to draw from the current random number stream, or
  a single number that makes the call reproducible; a supplied seed is
  set for this call only, and the caller's random number state is
  restored when the function returns.

## Value

A `data.frame` with columns `term` and `value` summarizing the Monte
Carlo sensitivity analysis. The `term` entries are: `total_N` (the
sample size evaluated), `empirical_power`, `analytic_power` (computed
from
[`ss_power_reg_coef`](https://yelleknek.github.io/DMAR/reference/ss_power_reg_coef.md)),
the mean / median / SD of the realized \\\hat b_j\\ (`mean_reg_coef_j`,
`median_reg_coef_j`, `sd_reg_coef_j`), of its standard error
(`mean_se_reg_coef_j`, `median_se_reg_coef_j`, `sd_se_reg_coef_j`), of
the test statistic (`mean_t_value`, `median_t_value`, `sd_t_value`), and
of the squared multiple correlation coefficient (`mean_R2`, `median_R2`,
`sd_R2`), `t_crit` (the critical value), and the input echoes `p`,
`which_predictor`, `true_reg_coef_j` and `supposed_reg_coef_j` (the
population and planning values of the targeted coefficient implied by
the supplied covariance structures), `desired_power` (NA when
`specified_N` was supplied instead), and `alpha_level`. The result
carries the `dmar_ss_power_sensitivity` class, so
[`tidy`](https://generics.r-lib.org/reference/tidy.html) reports the
planned sample size beside the empirical and analytic power, and
[`glance`](https://generics.r-lib.org/reference/glance.html) adds the
simulated estimator distribution beside the echoed inputs.

## Details

When the supposed and true covariance structures are identical, the
function performs a Monte Carlo evaluation of the planning procedure (no
misspecification); when they differ, it performs a sensitivity analysis
on the consequences of misspecifying the population covariance structure
for the targeted coefficient's power. The planning step calls
[`ss_power_reg_coef`](https://yelleknek.github.io/DMAR/reference/ss_power_reg_coef.md)
with the supposed covariance structure; the simulation step generates
data from the true covariance structure.

The coefficient is the unstandardized one, \\b_j\\, the *j*th element of
\\\Sigma\_{XX}^{-1}\sigma\_{XY}\\, from `true_cov_XX` and `true_cov_YX`
for the true value and from their supposed counterparts for the supposed
value, and every replication is fit in the units of the data. The *t*
statistic of a coefficient is the same whether or not the data are
standardized, so the power is that of the test of the standardized
coefficient as well; the standardized coefficient itself, \\b_j
\sigma\_{X_j}/\sigma_Y\\ (Kelley & Maxwell, 2003, p. 316), and its
precision are the subject of
[`ss_aipe_src_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_src_sensitivity.md).

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

Kelley, K., & Maxwell, S. E. (2003). Sample size for multiple
regression: Obtaining regression coefficients that are accurate, not
simply significant. *Psychological Methods, 8*(3), 305–321.
[doi:10.1037/1082-989X.8.3.305](https://doi.org/10.1037/1082-989X.8.3.305)

Kelley, K., & Maxwell, S. E. (2008). Sample size planning with
applications to multiple regression: Power and accuracy for omnibus and
targeted effects. In P. Alasuutari, L. Bickman, & J. Brannen (Eds.),
*The Sage handbook of social research methods* (pp. 166–192). Sage.

Maxwell, S. E. (2000). Sample size and multiple regression analysis.
*Psychological Methods, 5*(4), 434–458.
[doi:10.1037/1082-989X.5.4.434](https://doi.org/10.1037/1082-989X.5.4.434)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 4 on individual comparisons of means and
Chapter 6 on trend analysis.)

## See also

[`ss_power_reg_coef`](https://yelleknek.github.io/DMAR/reference/ss_power_reg_coef.md),
[`ss_aipe_reg_coef_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reg_coef_sensitivity.md)

[`design_consequences`](https://yelleknek.github.io/DMAR/reference/design_consequences.md)
for what a chosen design delivers: power, the Type S (sign) and Type M
(exaggeration) errors of the significance filter, and the expected
confidence interval width.

Other sample size for power:
[`power_fisher_exact()`](https://yelleknek.github.io/DMAR/reference/power_fisher_exact.md),
[`ss_aipe_mixed_effects()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mixed_effects.md),
[`ss_power_R2()`](https://yelleknek.github.io/DMAR/reference/ss_power_R2.md),
[`ss_power_R2_change()`](https://yelleknek.github.io/DMAR/reference/ss_power_R2_change.md),
[`ss_power_R2_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_power_R2_sensitivity.md),
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
[`ss_power_rm_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_rm_anova.md),
[`ss_power_sc()`](https://yelleknek.github.io/DMAR/reference/ss_power_sc.md),
[`ss_power_sem()`](https://yelleknek.github.io/DMAR/reference/ss_power_sem.md),
[`ss_power_smd()`](https://yelleknek.github.io/DMAR/reference/ss_power_smd.md),
[`ss_power_split_plot_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_split_plot_anova.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Targeted coefficient power sensitivity with two predictors. The
# default G = 1000 replications is used in practice; G is reduced here
# so the example runs quickly.
set.seed(113)
cov_XX  <- matrix(c(1, 0.3, 0.3, 1), nrow = 2)
cov_YX  <- c(0.4, 0.3)
ss_power_reg_coef_sensitivity(
  true_var_Y = 1, true_cov_YX = cov_YX, true_cov_XX = cov_XX,
  which_predictor = 1, desired_power = 0.80,
  G = 100, print_iter = FALSE
)
#>  term                 value
#>  total_N              62   
#>  empirical_power      0.800
#>  analytic_power       0.801
#>  mean_reg_coef_j      0.351
#>  median_reg_coef_j    0.347
#>  sd_reg_coef_j        0.120
#>  mean_se_reg_coef_j   0.120
#>  median_se_reg_coef_j 0.119
#>  sd_se_reg_coef_j     0.016
#>  mean_t_value         2.939
#>  median_t_value       3.001
#>  sd_t_value           1.019
#>  mean_R2              0.222
#>  median_R2            0.210
#>  sd_R2                0.085
#>  t_crit               2.001
#>  p                    2    
#>  which_predictor      1    
#>  true_reg_coef_j      0.341
#>  supposed_reg_coef_j  0.341
#>  desired_power        0.800
#>  alpha_level          0.050
#> 
#> Replications: 100
```
