# Sample Size or Power for a Factorial Between-Subjects ANOVA Effect

Determine the necessary per-cell sample size to achieve a desired level
of statistical power for a single *F* test (main effect or interaction)
in a between-subjects factorial ANOVA, or, given a per-cell sample size,
return the realized statistical power. The function handles two-way and
higher-order factorial designs.

## Usage

``` r
ss_power_factorial_anova(
  factor_levels,
  effect_indices,
  f = NULL,
  partial_eta_squared = NULL,
  desired_power = 0.85,
  alpha_level = 0.05,
  n_per_cell = NULL
)
```

## Arguments

- factor_levels:

  Integer vector giving the number of levels of each factor (e.g.,
  `c(2, 3)` for a 2 x 3 design, `c(2, 2, 4)` for a 2 x 2 x 4 design)

- effect_indices:

  Integer vector identifying the factors that define the effect of
  interest. For example, `1` requests the main effect of the first
  factor; `c(1, 2)` requests the AxB two-way interaction; `c(1, 2, 3)`
  requests the three-way interaction

- f:

  Cohen's *f* effect size for the chosen effect; supply this or
  `partial_eta_squared`, but not both

- partial_eta_squared:

  Partial eta squared for the chosen effect; supply this or `f`

- desired_power:

  Desired statistical power (default 0.85)

- alpha_level:

  Type I error rate (default 0.05)

- n_per_cell:

  Per-cell sample size; if specified, returns the realized power

## Value

A `data.frame` with rows for `necessary_n_per_cell` (or
`specified_n_per_cell`), `total_N`, `effect_df`, `error_df`,
`noncentrality`, and `actual_power`. The result carries the
`dmar_ss_power` class, so
[`tidy`](https://generics.r-lib.org/reference/tidy.html) and
[`glance`](https://generics.r-lib.org/reference/glance.html) summarize
it in broom convention (the reported size is the per-cell count).

## Details

For a between-subjects factorial design, the *F* statistic for the
chosen effect follows a noncentral *F* distribution under the
alternative with numerator degrees of freedom \\\prod\_{i \in S} (k_i -
1)\\ (where \\S\\ is `effect_indices` and \\k_i\\ is
`factor_levels[i]`), denominator degrees of freedom \\N - K\\ (where \\N
= n\_{cell} \prod k_i\\ and \\K = \prod k_i\\ is the number of cells),
and noncentrality parameter \\\lambda = N f^2\\. Cohen's *f* relates to
partial eta squared via \\f = \sqrt{\eta_p^2 / (1 - \eta_p^2)}\\.

The function searches over per-cell sample sizes until power reaches
`desired_power`; when `n_per_cell` is supplied it returns the realized
power.

For covariates, see
[`ss_power_factorial_ancova`](https://yelleknek.github.io/DMAR/reference/ss_power_factorial_ancova.md).
A complete worked three-factor example (a 2 x 4 x 3 design planned
effect by effect, simulated, analyzed with Type III sums of squares,
plotted, and followed up with focused contrasts) is the vignette
[`vignette("ancova_2x4x3_power", package = "DMAR")`](https://yelleknek.github.io/DMAR/articles/ancova_2x4x3_power.md).

## References

Cohen, J. (1988). *Statistical power analysis for the behavioral
sciences* (2nd ed.). Hillsdale, NJ: Lawrence Erlbaum.

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.

## See also

[`ss_power_one_way_anova`](https://yelleknek.github.io/DMAR/reference/ss_power_one_way_anova.md),
[`ss_power_c`](https://yelleknek.github.io/DMAR/reference/ss_power_c.md),
[`conf_limits_ncf`](https://yelleknek.github.io/DMAR/reference/conf_limits_ncf.md)

[`design_consequences`](https://yelleknek.github.io/DMAR/reference/design_consequences.md)
for what a chosen design delivers: power, the Type S (sign) and Type M
(exaggeration) errors of the significance filter, and the expected
confidence interval width.

Other sample size for power:
[`power_fisher_exact()`](https://yelleknek.github.io/DMAR/reference/power_fisher_exact.md),
[`ss_aipe_mixed_effects()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mixed_effects.md),
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
# 2 x 3 design, main effect of factor B (the 3-level factor), f = 0.25, power = .80
ss_power_factorial_anova(factor_levels = c(2, 3), effect_indices = 2,
                         f = 0.25, desired_power = 0.80)
#>  term                 value
#>  necessary_n_per_cell 27   
#>  total_N              162  
#>  effect_df            2    
#>  error_df             156  
#>  noncentrality        10.1 
#>  actual_power         0.813

# 2 x 2 design, AxB interaction, partial eta squared = 0.06, power = .80
ss_power_factorial_anova(factor_levels = c(2, 2), effect_indices = c(1, 2),
                         partial_eta_squared = 0.06, desired_power = 0.80)
#>  term                 value
#>  necessary_n_per_cell 32   
#>  total_N              128  
#>  effect_df            1    
#>  error_df             124  
#>  noncentrality        8.17 
#>  actual_power         0.81 

# 2 x 2 x 3 design, three-way interaction, f = 0.20
ss_power_factorial_anova(factor_levels = c(2, 2, 3), effect_indices = c(1, 2, 3),
                         f = 0.20, desired_power = 0.80)
#>  term                 value
#>  necessary_n_per_cell 21   
#>  total_N              252  
#>  effect_df            2    
#>  error_df             240  
#>  noncentrality        10.1 
#>  actual_power         0.814

# Realized power for n_per_cell = 20 in a 2x3 design, AxB interaction, f = 0.25
ss_power_factorial_anova(factor_levels = c(2, 3), effect_indices = c(1, 2),
                         f = 0.25, n_per_cell = 20)
#>  term                 value
#>  specified_n_per_cell 20   
#>  total_N              120  
#>  effect_df            2    
#>  error_df             114  
#>  noncentrality        7.5  
#>  actual_power         0.675
```
