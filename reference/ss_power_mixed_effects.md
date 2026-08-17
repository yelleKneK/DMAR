# Sample Size or Power for a Treatment Effect in a Two-Level Mixed-Effects Model

Determine the necessary number of level-2 units per arm to achieve a
desired level of statistical power for a treatment-versus-control
comparison in a two-level mixed-effects model with a random intercept
(e.g., individuals nested within clusters in a cluster-randomized trial,
or repeated measurements nested within subjects in a person-randomized
longitudinal study). Alternatively, given a number of level-2 units per
arm, return the realized statistical power.

## Usage

``` r
ss_power_mixed_effects(
  d,
  n,
  rho,
  J = NULL,
  desired_power = 0.85,
  alpha_level = 0.05,
  directional = FALSE
)
```

## Arguments

- d:

  Standardized treatment effect, defined as the population mean
  difference divided by the population standard deviation of the level-1
  outcome

- n:

  Number of level-1 units per level-2 unit (e.g., individuals per
  cluster, or measurements per subject); assumed equal across level-2
  units

- rho:

  Intra-class correlation (the proportion of total outcome variance
  attributable to differences between level-2 units); must be in \[0, 1)

- J:

  Number of level-2 units per arm (i.e., `J` treatment clusters and `J`
  control clusters); if specified, returns the realized power

- desired_power:

  Desired statistical power (default 0.85)

- alpha_level:

  Type I error rate (default 0.05)

- directional:

  Logical: `TRUE` for a one-sided test (in the same sign as `d`),
  `FALSE` (default) for a two-sided test

## Value

A `data.frame` with rows for `necessary_J_per_arm` (or
`specified_J_per_arm`), `total_N`, `noncentrality`, and `actual_power`.
The result carries the `dmar_ss_power` class, so
[`tidy`](https://generics.r-lib.org/reference/tidy.html) and
[`glance`](https://generics.r-lib.org/reference/glance.html) summarize
it in broom convention (the reported size is the number of clusters per
arm).

## Details

This function computes power for the fixed treatment effect at the
higher level of a two-level mixed-effects model with random intercept,
\$\$y\_{ij} = \beta_0 + \beta_1 T_j + u_j + \epsilon\_{ij},\$\$ where
\\u_j \sim N(0, \sigma_u^2)\\ is the level-2 random intercept and
\\\epsilon\_{ij} \sim N(0, \sigma_e^2)\\ is the level-1 residual. The
treatment indicator \\T_j\\ varies between level-2 units (i.e., entire
clusters or entire subjects are assigned to treatment or control). The
intra-class correlation is \\\rho = \sigma_u^2 / (\sigma_u^2 +
\sigma_e^2)\\ and the total outcome variance is \\\sigma_y^2 =
\sigma_u^2 + \sigma_e^2\\.

The standard error of the estimated treatment effect is
\\SE(\hat\beta_1) = \sigma_y \sqrt{2 (1 + (n - 1)\rho) / (J n)}\\,
giving a noncentrality parameter of \$\$\lambda = d \sqrt{J n / (2 (1 +
(n - 1)\rho))}\$\$ under a two-sample *t*-test with \\2J - 2\\ degrees
of freedom.

The factor \\1 + (n - 1)\rho\\ is the design effect: as the
within-cluster correlation grows, the effective information per level-2
unit shrinks, so more level-2 units are needed for a given level of
power.

## References

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.

Raudenbush, S. W. (1997). Statistical analysis and optimal design for
cluster randomized trials. *Psychological Methods, 2*, 173–185.
[doi:10.1037/1082-989X.2.2.173](https://doi.org/10.1037/1082-989X.2.2.173)

## See also

[`ss_power_split_plot_anova`](https://yelleknek.github.io/DMAR/reference/ss_power_split_plot_anova.md),
[`ss_power_smd`](https://yelleknek.github.io/DMAR/reference/ss_power_smd.md),
[`ss_power_pcm`](https://yelleknek.github.io/DMAR/reference/ss_power_pcm.md)

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

Other mixed models:
[`R2_mixed_effects()`](https://yelleknek.github.io/DMAR/reference/R2_mixed_effects.md),
[`R2_mixed_effects_decomposition()`](https://yelleknek.github.io/DMAR/reference/R2_mixed_effects_decomposition.md),
[`icc_lmer()`](https://yelleknek.github.io/DMAR/reference/icc_lmer.md),
[`manova_split_plot()`](https://yelleknek.github.io/DMAR/reference/manova_split_plot.md),
[`mixed_anova()`](https://yelleknek.github.io/DMAR/reference/mixed_anova.md),
[`ss_aipe_mixed_effects()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mixed_effects.md),
[`ss_aipe_mixed_effects_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mixed_effects_sensitivity.md),
[`ss_power_split_plot_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_split_plot_anova.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# 30 individuals per cluster, ICC = 0.05, standardized effect d = 0.30, power = .80
ss_power_mixed_effects(d = 0.30, n = 30, rho = 0.05, desired_power = 0.80)
#>  term                value
#>  necessary_J_per_arm 16   
#>  total_N             960  
#>  noncentrality       2.97 
#>  actual_power        0.819

# Same effect but with much higher ICC (e.g., schools or therapists)
ss_power_mixed_effects(d = 0.30, n = 30, rho = 0.20, desired_power = 0.80)
#>  term                value
#>  necessary_J_per_arm 41   
#>  total_N             2460 
#>  noncentrality       2.85 
#>  actual_power        0.805

# Realized power with 25 level-2 units per arm
ss_power_mixed_effects(d = 0.30, n = 30, rho = 0.05, J = 25)
#>  term                value
#>  specified_J_per_arm 25   
#>  total_N             1500 
#>  noncentrality       3.71 
#>  actual_power        0.953

# directional test
ss_power_mixed_effects(d = 0.30, n = 30, rho = 0.05, desired_power = 0.80,
                         directional = TRUE)
#>  term                value
#>  necessary_J_per_arm 12   
#>  total_N             720  
#>  noncentrality       2.57 
#>  actual_power        0.801
```
