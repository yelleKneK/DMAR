# Sample Size or Power for a One-Way Repeated Measures ANOVA Omnibus *F* Test

Determine the necessary number of subjects (the total sample size \\N\\)
to achieve a desired level of statistical power for the omnibus *F* test
of the within-subjects factor in a one-way repeated measures ANOVA, or,
given a number of subjects, return the realized statistical power.

## Usage

``` r
ss_power_rm_anova(
  a,
  true_f = NULL,
  true_eta_squared = NULL,
  true_cor = 0,
  true_epsilon = 1,
  desired_power = 0.85,
  alpha_level = 0.05,
  specified_N = NULL
)
```

## Arguments

- a:

  Number of measurement occasions (levels of the within-subjects factor)

- true_f:

  Cohen's *f* effect size for the within-subjects factor (the population
  value); supply this or `true_eta_squared`, but not both

- true_eta_squared:

  Population eta squared (proportion of variance, on the relevant scale,
  accounted for by the within-subjects factor); supply this or `true_f`

- true_cor:

  Average correlation among the repeated measures (default 0). With
  `true_cor > 0`, the within-subjects test gains efficiency relative to
  a between-subjects analog

- true_epsilon:

  Greenhouse-Geisser / Huynh-Feldt sphericity adjustment in (0, 1\]
  (default 1, sphericity assumed). When `true_epsilon < 1`, both
  numerator and denominator degrees of freedom are multiplied by
  `true_epsilon`

- desired_power:

  Desired statistical power (default 0.85)

- alpha_level:

  Type I error rate (default 0.05)

- specified_N:

  Total sample size: the number of subjects, each measured at all `a`
  occasions; if specified, the function returns the realized power

## Value

A `data.frame` with rows for `necessary_N` (or `specified_N`), the
number of subjects, `a`, `df_effect`, `df_error`, `noncentral_f_parm`,
and `actual_power`. The result carries the `dmar_ss_power` class, so
[`tidy`](https://generics.r-lib.org/reference/tidy.html) and
[`glance`](https://generics.r-lib.org/reference/glance.html) summarize
it in broom convention (the reported size is the number of subjects).

## Details

Under the alternative hypothesis with sphericity (`true_epsilon = 1`),
the within-subjects *F* statistic follows a noncentral *F* distribution
with numerator df \\a - 1\\, denominator df \\(N - 1)(a - 1)\\, and
noncentrality parameter \\\lambda = N a f^2 / (1 - \rho)\\, where \\N\\
is the number of subjects, \\f\\ is Cohen's *f* (`true_f`) for the
within-subjects effect and \\\rho\\ is the average correlation across
the repeated measures (`true_cor`; Maxwell, Delaney, & Kelley, 2027).
Setting `true_cor = 0` reduces to the between-subjects expression.

When sphericity is violated, supplying `true_epsilon` (e.g., a
Greenhouse-Geisser estimate) rescales the test using the Muller-Barton
convention: both numerator and denominator degrees of freedom are
multiplied by `true_epsilon`, and the noncentrality parameter is
likewise multiplied by `true_epsilon`. Smaller `true_epsilon` therefore
reduces power and increases the necessary sample size.

## References

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.

Cohen, J. (1988). *Statistical power analysis for the behavioral
sciences* (2nd ed.). Hillsdale, NJ: Lawrence Erlbaum.

## See also

[`ss_power_one_way_anova`](https://yelleknek.github.io/DMAR/reference/ss_power_one_way_anova.md),
[`ss_power_pcm`](https://yelleknek.github.io/DMAR/reference/ss_power_pcm.md)

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
[`ss_power_reg_coef_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_power_reg_coef_sensitivity.md),
[`ss_power_sc()`](https://yelleknek.github.io/DMAR/reference/ss_power_sc.md),
[`ss_power_sem()`](https://yelleknek.github.io/DMAR/reference/ss_power_sem.md),
[`ss_power_smd()`](https://yelleknek.github.io/DMAR/reference/ss_power_smd.md),
[`ss_power_split_plot_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_split_plot_anova.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# 4 measurement occasions, f = 0.25, average within-subject correlation 0.5, power = .80
ss_power_rm_anova(a = 4, true_f = 0.25, true_cor = 0.5, desired_power = 0.80)
#>  term              value 
#>  necessary_N       24    
#>  a                 4     
#>  df_effect         3     
#>  df_error          69    
#>  noncentral_f_parm 12.000
#>  actual_power      0.817 

# Same but with Greenhouse-Geisser epsilon = 0.75
ss_power_rm_anova(a = 4, true_f = 0.25, true_cor = 0.5, true_epsilon = 0.75, desired_power = 0.80)
#>  term              value 
#>  necessary_N       29    
#>  a                 4     
#>  df_effect         2.250 
#>  df_error          63    
#>  noncentral_f_parm 10.875
#>  actual_power      0.814 

# Realized power at N = 20 subjects, a = 4
ss_power_rm_anova(a = 4, true_f = 0.25, true_cor = 0.5, specified_N = 20)
#>  term              value 
#>  specified_N       20    
#>  a                 4     
#>  df_effect         3     
#>  df_error          57    
#>  noncentral_f_parm 10.000
#>  actual_power      0.729 
```
