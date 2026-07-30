# Sample Size or Power for a Standardized Mean Difference (Two Independent Groups)

Determine the necessary per-group sample size to achieve a desired level
of statistical power for the two-sample (independent groups) *t*-test on
a standardized mean difference (Cohen's *d*; equivalently Hedges' *g*
and Glass's *g* for sample size purposes). Alternatively, given a
per-group sample size, return the realized statistical power.

## Usage

``` r
ss_power_smd(
  smd,
  desired_power = 0.85,
  alpha_level = 0.05,
  n_1 = NULL,
  n_2 = NULL,
  directional = FALSE
)
```

## Arguments

- smd:

  Supposed standardized mean difference (Cohen's *d*) the design is
  planned against: a value the researcher posits for the population,
  either a minimally important effect or a value believed to be true in
  the population, never a sample estimate. Echoed in the returned table
  as the `supposed_smd` row.

- desired_power:

  Desired statistical power (default 0.85)

- alpha_level:

  Type I error rate (default 0.05)

- n_1:

  Sample size for group 1 (if specified, the function returns the
  realized power; assumes `n_2 = n_1` unless `n_2` is also given)

- n_2:

  Sample size for group 2 (defaults to `n_1` when `n_1` is supplied)

- directional:

  Logical: `TRUE` for a one-sided test (in the same sign as `smd`),
  `FALSE` (default) for a two-sided test

## Value

A `data.frame` with `term` and `value` columns. The design result comes
first, followed by rows that echo the user-supplied planning inputs, so
the assumptions the power was evaluated under travel with the result.
The `supposed_smd` row is the supposed effect the plan is built on: a
value the researcher posits, either a minimally important effect or a
value believed to be true in the population, never a sample estimate.
The `tails` row is 2 for a nondirectional test and 1 for a directional
test.

- When `n_1` is `NULL`:

  Result rows `necessary_n_per_group`, `actual_power`, and
  `noncentral_t_parm`, then the planning inputs `supposed_smd`,
  `desired_power`, `alpha_level`, and `tails`.

- When `n_1` is specified:

  Result rows `specified_n_1`, `specified_n_2`, `actual_power`, and
  `noncentral_t_parm`, then the planning inputs `supposed_smd`,
  `alpha_level`, and `tails` (the supplied group sizes are the
  `specified_n_1` / `specified_n_2` rows).

## Details

The two-sample *t*-statistic with pooled standard deviation follows a
noncentral *t*-distribution with \\n_1 + n_2 - 2\\ degrees of freedom
and noncentrality parameter \\\lambda = \delta \sqrt{n_1 n_2 / (n_1 +
n_2)}\\, where \\\delta\\ is the population standardized mean
difference. For balanced designs (\\n_1 = n_2 = n\\) this simplifies to
\\\lambda = \delta \sqrt{n / 2}\\.

Power is computed as the probability that the absolute value of the test
statistic exceeds the critical value(s) under the alternative; the
function returns the per-group sample size for which power first reaches
`desired_power`.

Kelley and Rausch (2006) develop the accuracy in parameter estimation
approach to planning the sample size for the standardized mean
difference, implemented in
[`ss_aipe_smd`](https://yelleknek.github.io/DMAR/reference/ss_aipe_smd.md).

## References

Cohen, J. (1988). *Statistical power analysis for the behavioral
sciences* (2nd ed.). Hillsdale, NJ: Lawrence Erlbaum.

Kelley, K., Maxwell, S. E., & Rausch, J. R. (2003). Obtaining power or
obtaining precision: Delineating methods of sample size planning.
*Evaluation and the Health Professions, 26*(3), 258–287.
[doi:10.1177/0163278703255242](https://doi.org/10.1177/0163278703255242)

Kelley, K., & Rausch, J. R. (2006). Sample size planning for the
standardized mean difference: Accuracy in parameter estimation via
narrow confidence intervals. *Psychological Methods, 11*(4), 363–385.
[doi:10.1037/1082-989X.11.4.363](https://doi.org/10.1037/1082-989X.11.4.363)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.

Maxwell, S. E., Kelley, K., & Rausch, J. R. (2008). Sample size planning
for statistical power and accuracy in parameter estimation. *Annual
Review of Psychology, 59*, 537–563.
[doi:10.1146/annurev.psych.59.103006.093735](https://doi.org/10.1146/annurev.psych.59.103006.093735)

## See also

[`ss_aipe_smd`](https://yelleknek.github.io/DMAR/reference/ss_aipe_smd.md),
[`ci_smd`](https://yelleknek.github.io/DMAR/reference/ci_smd.md),
[`smd`](https://yelleknek.github.io/DMAR/reference/smd.md),
[`conf_limits_nct`](https://yelleknek.github.io/DMAR/reference/conf_limits_nct.md)

[`design_consequences`](https://yelleknek.github.io/DMAR/reference/design_consequences.md)
for what a chosen design delivers: power, the Type S (sign) and Type M
(exaggeration) errors of the significance filter, and the expected
confidence interval width.

Other sample size for power:
[`power_fisher_exact()`](https://yelleknek.github.io/DMAR/reference/power_fisher_exact.md),
[`ss_aipe_mixed_effects()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mixed_effects.md),
[`ss_aipe_tost_smd()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_tost_smd.md),
[`ss_power_R2()`](https://yelleknek.github.io/DMAR/reference/ss_power_R2.md),
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
[`ss_power_reg_coef_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_power_reg_coef_sensitivity.md),
[`ss_power_rm_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_rm_anova.md),
[`ss_power_sc()`](https://yelleknek.github.io/DMAR/reference/ss_power_sc.md),
[`ss_power_sem()`](https://yelleknek.github.io/DMAR/reference/ss_power_sem.md),
[`ss_power_split_plot_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_split_plot_anova.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Per-group sample size for d = 0.5, alpha = .05, power = .80, two-sided
ss_power_smd(smd = 0.5, desired_power = 0.80)
#>  term                  value
#>  necessary_n_per_group 64   
#>  actual_power          0.801
#>  noncentral_t_parm     2.83 
#>  supposed_smd          0.5  
#>  desired_power         0.8  
#>  alpha_level           0.05 
#>  tails                 2    

# Same with a directional (one-sided) test
ss_power_smd(smd = 0.5, desired_power = 0.80, directional = TRUE)
#>  term                  value
#>  necessary_n_per_group 51   
#>  actual_power          0.806
#>  noncentral_t_parm     2.52 
#>  supposed_smd          0.5  
#>  desired_power         0.8  
#>  alpha_level           0.05 
#>  tails                 1    

# Realized power given balanced n = 30 per group
ss_power_smd(smd = 0.5, n_1 = 30)
#>  term              value
#>  specified_n_1     30   
#>  specified_n_2     30   
#>  actual_power      0.478
#>  noncentral_t_parm 1.94 
#>  supposed_smd      0.5  
#>  alpha_level       0.05 
#>  tails             2    

# Realized power for unbalanced (n_1 = 30, n_2 = 50)
ss_power_smd(smd = 0.5, n_1 = 30, n_2 = 50)
#>  term              value
#>  specified_n_1     30   
#>  specified_n_2     50   
#>  actual_power      0.571
#>  noncentral_t_parm 2.17 
#>  supposed_smd      0.5  
#>  alpha_level       0.05 
#>  tails             2    
```
