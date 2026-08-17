# Sample Size Planning for Structural Equation Modeling From the Power Analysis Perspective

Calculate the necessary sample size for an SEM study, so as to have
enough power to reject the null hypothesis that (a) the model has
perfect fit, or (b) the difference in fit between two nested models
equal some specified amount.

## Usage

``` r
ss_power_sem(
  F_ML = NULL,
  df = NULL,
  RMSEA_null = NULL,
  RMSEA_true = NULL,
  F_full = NULL,
  F_res = NULL,
  RMSEA_full = NULL,
  RMSEA_res = NULL,
  df_full = NULL,
  df_res = NULL,
  alpha_level = 0.05,
  desired_power = 0.85
)
```

## Arguments

- F_ML:

  The true maximum likelihood fit function value in the population for
  the model of interest. Leave this argument NULL if you are doing
  nested model significance tests

- df:

  The degrees of freedom of the model of interest. Leave this argument
  NULL if you are doing nested model significance tests

- RMSEA_null:

  The model's population RMSEA under the null hypothesis. Leave this
  argument NULL if you are doing nested model significance tests

- RMSEA_true:

  The model's population RMSEA under the alternative hypothesis. This
  should be the model's true population RMSEA value. Leave this argument
  NULL if you are doing nested model significance tests

- F_full:

  The maximum likelihood fit function value for the full model

- F_res:

  The maximum likelihood fit function value for the restricted model

- RMSEA_full:

  The population RMSEA value for the full model

- RMSEA_res:

  The population RMSEA value for the restricted model

- df_full:

  The degrees of freedom for the full model

- df_res:

  The degrees of freedom for the restricted model

- alpha_level:

  The Type I error rate. Defaults to 0.05.

- desired_power:

  The desired power. Defaults to 0.85, matching the rest of the
  `ss_power_*` family.

## Value

A `data.frame` with a `necessary_N` row, the smallest integer *N* whose
power reaches `desired_power` under the supplied fit-function or RMSEA
alternative, and an `actual_power` row giving the realized power at that
*N*.

## References

MacCallum, R. C., Browne, M. W., & Sugawara, H. M. (1996). Power
analysis and determination of sample size for covariance structure
modeling. *Psychological Methods, 1*(2), 130–149.
[doi:10.1037/1082-989X.1.2.130](https://doi.org/10.1037/1082-989X.1.2.130)

Lai, K., & Kelley, K. (2011). Accuracy in parameter estimation for
targeted effects in structural equation modeling: Sample size planning
for narrow confidence intervals. *Psychological Methods, 16*(2),
127–148. [doi:10.1037/a0021764](https://doi.org/10.1037/a0021764)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.

## See also

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
[`ss_power_one_way_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_one_way_anova.md),
[`ss_power_pcm()`](https://yelleknek.github.io/DMAR/reference/ss_power_pcm.md),
[`ss_power_r()`](https://yelleknek.github.io/DMAR/reference/ss_power_r.md),
[`ss_power_rc()`](https://yelleknek.github.io/DMAR/reference/ss_power_rc.md),
[`ss_power_reg_coef()`](https://yelleknek.github.io/DMAR/reference/ss_power_reg_coef.md),
[`ss_power_reg_coef_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_power_reg_coef_sensitivity.md),
[`ss_power_rm_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_rm_anova.md),
[`ss_power_sc()`](https://yelleknek.github.io/DMAR/reference/ss_power_sc.md),
[`ss_power_smd()`](https://yelleknek.github.io/DMAR/reference/ss_power_smd.md),
[`ss_power_split_plot_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_split_plot_anova.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# One-model test: necessary N to reject H0: RMSEA = 0 in favor of
# a model whose population RMSEA is 0.05 at 80% power, alpha = .05,
# with df = 20.
ss_power_sem(RMSEA_null = 0, RMSEA_true = 0.05, df = 20,
             alpha_level = 0.05, desired_power = 0.80)
#>  term         value
#>  necessary_N  421  
#>  actual_power 0.801

# Equivalent input via the population fit function: F_ML = df * RMSEA^2.
ss_power_sem(F_ML = 20 * 0.05^2, df = 20, alpha_level = 0.05, desired_power = 0.80)
#>  term         value
#>  necessary_N  421  
#>  actual_power 0.801

# Two-model nested test: necessary N to detect the difference
# between a full model (RMSEA = 0.04, df = 18) and a restricted
# model (RMSEA = 0.06, df = 22) at 80% power.
ss_power_sem(RMSEA_full = 0.04, df_full = 18,
             RMSEA_res = 0.06, df_res = 22,
             alpha_level = 0.05, desired_power = 0.80)
#>  term         value
#>  necessary_N  238  
#>  actual_power 0.8  
```
