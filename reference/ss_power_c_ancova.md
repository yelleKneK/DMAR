# Sample Size or Power for an Unstandardized Contrast in a One-Way ANCOVA

Determine the necessary per-group sample size to achieve a desired level
of statistical power for the test of a single planned (unstandardized)
contrast on the adjusted means in a one-way analysis of covariance, or,
given a per-group sample size, return the realized statistical power.

## Usage

``` r
ss_power_c_ancova(
  psi,
  c_weights,
  sigma,
  rho,
  desired_power = 0.85,
  alpha_level = 0.05,
  n = NULL,
  directional = FALSE
)
```

## Arguments

- psi:

  The population unstandardized contrast effect on the adjusted means,
  \\\psi = \sum c_j \mu^{(adj)}\_j\\

- c_weights:

  Vector of contrast weights (must sum to zero); use fractional weights
  so the positive weights sum to 1

- sigma:

  Within-group population standard deviation of the response (the same
  \\\sigma\\ as in a one-way ANOVA on the response)

- rho:

  Within-group population correlation between the response and the
  covariate; must lie in (-1, 1)

- desired_power:

  Desired statistical power (default 0.85)

- alpha_level:

  Type I error rate (default 0.05)

- n:

  Per-group sample size (assumed balanced); if specified, returns the
  realized power

- directional:

  Logical: `TRUE` for a one-sided test (in the same sign as `psi`),
  `FALSE` (default) for a two-sided test

## Value

A `data.frame` with rows for `necessary_n_per_group` (or
`specified_n_per_group`), `actual_power`, and `noncentral_t_parm`. The
result carries the `dmar_ss_power` class, so
[`tidy`](https://generics.r-lib.org/reference/tidy.html) and
[`glance`](https://generics.r-lib.org/reference/glance.html) summarize
it in broom convention.

## Details

This function uses the standard large-sample formulation in which the
ANCOVA error variance is \\\sigma^2\_{adj} = \sigma^2 (1 - \rho^2)\\,
the contrast *t*-statistic has degrees of freedom \\N - J - 1\\ (one
less than the corresponding ANOVA contrast because of the covariate),
and the noncentrality parameter is \\\lambda = \psi / (\sigma \sqrt{1 -
\rho^2} \sqrt{\sum c_j^2 / n})\\. This assumes the covariate means are
equal across groups (the typical assumption under random assignment);
for designs with substantial group differences in the covariate, the
small-sample correction \\1 + (\bar X\_{j} - \bar X\_{\cdot})^2 /
SS^{(within)}\_X\\ would slightly inflate the standard error and reduce
power, an effect that is negligible for moderate or large \\n\\.

## References

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.

Lai, K., & Kelley, K. (2012). Accuracy in parameter estimation for
ANCOVA and ANOVA contrasts: Sample size planning via narrow confidence
intervals. *British Journal of Mathematical and Statistical Psychology,
65*, 350–370.
[doi:10.1111/j.2044-8317.2011.02029.x](https://doi.org/10.1111/j.2044-8317.2011.02029.x)

## See also

[`ss_power_c`](https://yelleknek.github.io/DMAR/reference/ss_power_c.md),
[`ci_c_ancova`](https://yelleknek.github.io/DMAR/reference/ci_c_ancova.md),
[`ss_aipe_c_ancova`](https://yelleknek.github.io/DMAR/reference/ss_aipe_c_ancova.md)

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
# Same population contrast as in the ANOVA example, with rho = 0.5 between
# outcome and covariate; ANCOVA is more efficient than ANOVA here.
ss_power_c_ancova(psi = 0.5, c_weights = c(0.5, 0.5, -0.5, -0.5),
                  sigma = 1, rho = 0.5, desired_power = 0.80)
#>  term                  value
#>  necessary_n_per_group 25   
#>  actual_power          0.815
#>  noncentral_t_parm     2.89 

# Realized power for n = 30 per group
ss_power_c_ancova(psi = 0.5, c_weights = c(0.5, 0.5, -0.5, -0.5),
                  sigma = 1, rho = 0.5, n = 30)
#>  term                  value
#>  specified_n_per_group 30   
#>  actual_power          0.88 
#>  noncentral_t_parm     3.16 
```
