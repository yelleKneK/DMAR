# Sample Size or Power for a One-Way Between-Subjects ANOVA Omnibus *F* Test

Determine the necessary total sample size to achieve a desired level of
statistical power for the omnibus *F* test in a one-way between-subjects
analysis of variance, or, given a total sample size, return the realized
statistical power.

## Usage

``` r
ss_power_one_way_anova(
  a,
  f = NULL,
  eta_squared = NULL,
  desired_power = 0.85,
  alpha_level = 0.05,
  N = NULL
)
```

## Arguments

- a:

  Number of groups (levels of the between-subjects factor)

- f:

  Cohen's *f* effect size (the population value); supply this or
  `eta_squared`, but not both

- eta_squared:

  Population eta squared (proportion of total variance accounted for by
  group membership); supply this or `f`

- desired_power:

  Desired statistical power (default 0.85)

- alpha_level:

  Type I error rate (default 0.05)

- N:

  Total sample size; if specified, the function returns the realized
  power (the ss_power\_\* family is not uniform here:
  [`ss_power_contrast`](https://yelleknek.github.io/DMAR/reference/ss_power_contrast.md)
  takes a *per-group* size)

## Value

A `data.frame`. When a sample size is being planned (`N` not supplied)
the rows are `necessary_N`, `n_per_group`, `a`, `noncentrality`, and
`actual_power`; the search constructs the total as a balanced design, so
`n_per_group` is a whole-number per-group count. When `N` is supplied,
power is evaluated at that total `N` directly and the rows are
`specified_N`, `a`, `noncentrality`, and `actual_power` (no
`n_per_group` row, since balance is not assumed). The result carries the
`dmar_ss_power` class, so
[`tidy`](https://generics.r-lib.org/reference/tidy.html) and
[`glance`](https://generics.r-lib.org/reference/glance.html) summarize
it in broom convention; the summarized sample size is the total `N`.

## Details

Under the alternative hypothesis, the omnibus *F* statistic follows a
noncentral *F* distribution with numerator df \\a - 1\\, denominator df
\\N - a\\, and noncentrality parameter \\\lambda = N f^2\\. Cohen's *f*
relates to eta squared via \\f = \sqrt{\eta^2 / (1 - \eta^2)}\\.

The function searches over total sample sizes \\N\\ (treating per-group
\\N/a\\ as balanced) until power first reaches `desired_power`. When `N`
is supplied it instead reports the realized power at that `N`.

## References

Cohen, J. (1988). *Statistical power analysis for the behavioral
sciences* (2nd ed.). Hillsdale, NJ: Lawrence Erlbaum.

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.

## See also

[`ss_power_factorial_anova`](https://yelleknek.github.io/DMAR/reference/ss_power_factorial_anova.md),
[`ss_power_c`](https://yelleknek.github.io/DMAR/reference/ss_power_c.md),
[`ss_power_sc`](https://yelleknek.github.io/DMAR/reference/ss_power_sc.md),
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
[`ss_power_factorial_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_factorial_anova.md),
[`ss_power_indirect_effect()`](https://yelleknek.github.io/DMAR/reference/ss_power_indirect_effect.md),
[`ss_power_mixed_effects()`](https://yelleknek.github.io/DMAR/reference/ss_power_mixed_effects.md),
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
# Three groups, f = 0.25, power = .80
ss_power_one_way_anova(a = 3, f = 0.25, desired_power = 0.80)
#>  term          value
#>  necessary_N   159  
#>  n_per_group   53   
#>  a             3    
#>  noncentrality 9.94 
#>  actual_power  0.805

# Same effect specified via eta squared
ss_power_one_way_anova(a = 3, eta_squared = 0.0588, desired_power = 0.80)
#>  term          value
#>  necessary_N   159  
#>  n_per_group   53   
#>  a             3    
#>  noncentrality 9.93 
#>  actual_power  0.805

# Realized power at N = 60 across 3 groups
ss_power_one_way_anova(a = 3, f = 0.25, N = 60)
#>  term          value
#>  specified_N   60   
#>  a             3    
#>  noncentrality 3.75 
#>  actual_power  0.374
```
