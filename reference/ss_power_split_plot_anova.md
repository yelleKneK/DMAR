# Sample Size or Power for a Mixed-Effects ANOVA (Between X Within Design)

Determine the necessary per-group sample size to achieve a desired level
of statistical power for one of the three *F* tests in a mixed-effects
ANOVA with one between-subjects factor and one within-subjects factor
(the between-subjects main effect, the within-subjects main effect, or
the between x within interaction), or, given a per-group sample size,
return the realized statistical power. (This design is also commonly
called a split-plot factorial.)

## Usage

``` r
ss_power_split_plot_anova(
  a,
  b,
  effect,
  true_f = NULL,
  true_eta_squared_partial = NULL,
  true_cor,
  true_epsilon = 1,
  desired_power = 0.85,
  alpha_level = 0.05,
  specified_n_per_group = NULL
)
```

## Arguments

- a:

  Number of levels of the between-subjects factor (i.e., number of
  groups)

- b:

  Number of levels of the within-subjects factor (i.e., number of
  measurement occasions)

- effect:

  Which *F* test to compute power for: `"between"`, `"within"`, or
  `"interaction"`

- true_f:

  Cohen's *f* effect size for the chosen effect (the population value);
  supply this or `true_eta_squared_partial`, but not both

- true_eta_squared_partial:

  Population partial eta squared for the chosen effect; supply this or
  `true_f`

- true_cor:

  Average correlation among the repeated measures within a subject (must
  lie in (-1, 1)). Higher `true_cor` reduces power for the
  between-subjects test (because subject means contain more redundant
  information) and increases power for the within-subjects and
  interaction tests

- true_epsilon:

  Greenhouse-Geisser / Huynh-Feldt sphericity adjustment in (0, 1\]
  (default 1, sphericity assumed). Applied to the within-subjects and
  interaction tests but not the between-subjects test. Both numerator
  and denominator df, and the noncentrality, are multiplied by
  `true_epsilon` (Muller-Barton convention)

- desired_power:

  Desired statistical power (default 0.85)

- alpha_level:

  Type I error rate (default 0.05)

- specified_n_per_group:

  Per-group (between-subjects) sample size; if specified, returns the
  realized power

## Value

A `data.frame` with rows for `necessary_n_per_group` (or
`specified_n_per_group`), `total_N`, `df_effect`, `df_error`,
`noncentral_f_parm`, and `actual_power`. The result carries the
`dmar_ss_power` class, so
[`tidy`](https://generics.r-lib.org/reference/tidy.html) and
[`glance`](https://generics.r-lib.org/reference/glance.html) summarize
it in broom convention (the reported size is the per-group count).

## Details

This is a two-factor mixed-effects design: one between-subjects factor
with \\a\\ levels and one within-subjects factor with \\b\\ levels;
\\n\\ subjects are randomly assigned to each between-subjects level and
each subject is measured at all \\b\\ within-subjects levels, for \\N =
na\\ subjects total and \\Nb\\ observations. The covariance among the
\\b\\ within-subject observations is summarized by `true_cor`, the
average pairwise correlation.

The three *F* tests have noncentrality parameters \$\$\lambda\_{B} = N b
f^2 / (1 + (b - 1) \rho)\$\$ for the between-subjects test (numerator df
\\a - 1\\, denominator df \\N - a\\), \$\$\lambda\_{W} = N b f^2 \\
\epsilon / (1 - \rho)\$\$ for the within-subjects test (numerator df
\\(b - 1)\epsilon\\, denominator df \\(N - a)(b - 1)\epsilon\\), and the
same form as \\\lambda_W\\ for the interaction (numerator df \\(a -
1)(b - 1)\epsilon\\, same denominator df). Cohen's *f* relates to
partial eta squared via \\f = \sqrt{\eta_p^2 / (1 - \eta_p^2)}\\. In
these formulas \\f\\ is `true_f` (or the value implied by
`true_eta_squared_partial`, which is \\\eta_p^2\\), \\\rho\\ is
`true_cor`, \\\epsilon\\ is `true_epsilon`, and \\n\\ is the per-group
size (`specified_n_per_group`, or the size the function plans).

This design is the compound-symmetry (random intercept) special case of
the two-level linear mixed-effects model: `true_cor` is the intraclass
correlation and the \\b\\ occasions are the observations within a
subject (i.e., the Level 1 units). For two between-subjects groups the
between-subjects *F*(1, .) test is therefore the two-level treatment *t*
test of
[`ss_power_mixed_effects`](https://yelleknek.github.io/DMAR/reference/ss_power_mixed_effects.md)
squared, so the two planners agree on that shared case.

## References

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.

Muller, K. E., & Barton, C. N. (1989). Approximate power for repeated
measures ANOVA lacking sphericity. *Journal of the American Statistical
Association, 84*, 549–555.

## See also

[`ss_power_one_way_anova`](https://yelleknek.github.io/DMAR/reference/ss_power_one_way_anova.md),
[`ss_power_factorial_anova`](https://yelleknek.github.io/DMAR/reference/ss_power_factorial_anova.md),
[`ss_power_rm_anova`](https://yelleknek.github.io/DMAR/reference/ss_power_rm_anova.md),
[`ss_power_mixed_effects`](https://yelleknek.github.io/DMAR/reference/ss_power_mixed_effects.md)

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
[`ss_power_rm_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_rm_anova.md),
[`ss_power_sc()`](https://yelleknek.github.io/DMAR/reference/ss_power_sc.md),
[`ss_power_sem()`](https://yelleknek.github.io/DMAR/reference/ss_power_sem.md),
[`ss_power_smd()`](https://yelleknek.github.io/DMAR/reference/ss_power_smd.md)

Other mixed models:
[`R2_mixed_effects()`](https://yelleknek.github.io/DMAR/reference/R2_mixed_effects.md),
[`R2_mixed_effects_decomposition()`](https://yelleknek.github.io/DMAR/reference/R2_mixed_effects_decomposition.md),
[`icc_lmer()`](https://yelleknek.github.io/DMAR/reference/icc_lmer.md),
[`manova_split_plot()`](https://yelleknek.github.io/DMAR/reference/manova_split_plot.md),
[`mixed_anova()`](https://yelleknek.github.io/DMAR/reference/mixed_anova.md),
[`ss_aipe_mixed_effects()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mixed_effects.md),
[`ss_aipe_mixed_effects_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mixed_effects_sensitivity.md),
[`ss_power_mixed_effects()`](https://yelleknek.github.io/DMAR/reference/ss_power_mixed_effects.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# 2 groups, 4 occasions, between-subjects effect, f = 0.25,
# average within-subject correlation 0.5, power = .80
ss_power_split_plot_anova(a = 2, b = 4, effect = "between", true_f = 0.25,
                    true_cor = 0.5, desired_power = 0.80)
#>  term                  value
#>  necessary_n_per_group 41   
#>  total_N               82   
#>  df_effect             1    
#>  df_error              80   
#>  noncentral_f_parm     8.200
#>  actual_power          0.808

# Same design, within-subjects (occasion) main effect, f = 0.25
ss_power_split_plot_anova(a = 2, b = 4, effect = "within", true_f = 0.25,
                    true_cor = 0.5, desired_power = 0.80)
#>  term                  value 
#>  necessary_n_per_group 12    
#>  total_N               24    
#>  df_effect             3     
#>  df_error              66    
#>  noncentral_f_parm     12.000
#>  actual_power          0.816 

# Same design, between x within interaction, f = 0.25
ss_power_split_plot_anova(a = 2, b = 4, effect = "interaction", true_f = 0.25,
                    true_cor = 0.5, desired_power = 0.80)
#>  term                  value 
#>  necessary_n_per_group 12    
#>  total_N               24    
#>  df_effect             3     
#>  df_error              66    
#>  noncentral_f_parm     12.000
#>  actual_power          0.816 

# Realized power for n = 25 per group on the interaction test, partial eta^2 = 0.06
ss_power_split_plot_anova(a = 2, b = 4, effect = "interaction",
                    true_eta_squared_partial = 0.06, true_cor = 0.5, specified_n_per_group = 25)
#>  term                  value 
#>  specified_n_per_group 25    
#>  total_N               50    
#>  df_effect             3     
#>  df_error              144   
#>  noncentral_f_parm     25.532
#>  actual_power          0.993 

# Greenhouse-Geisser correction with epsilon = 0.7 on the within-subjects test
ss_power_split_plot_anova(a = 2, b = 4, effect = "within", true_f = 0.25,
                    true_cor = 0.5, true_epsilon = 0.7, desired_power = 0.80)
#>  term                  value 
#>  necessary_n_per_group 15    
#>  total_N               30    
#>  df_effect             2.100 
#>  df_error              58.800
#>  noncentral_f_parm     10.500
#>  actual_power          0.808 
```
