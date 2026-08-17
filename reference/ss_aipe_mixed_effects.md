# AIPE Sample Size Planning for a Fixed Effect in a Two-Level Mixed-Effects Model

Computes the minimum number of clusters (level-2 units) needed so that
the confidence interval on a level-1 fixed-effect slope has expected
full width no larger than \\\omega\\ (Kelley, 2007; Raudenbush & Liu,
2001; Snijders & Bosker, 2012). The function inverts the closed-form
approximation for the variance of a fixed-effect slope in a balanced
two-level random-intercept model.

## Usage

``` r
ss_aipe_mixed_effects(
  sigma2_y,
  sigma2_x,
  icc,
  width,
  cluster_size = 20L,
  conf_level = 0.95
)
```

## Arguments

- sigma2_y:

  Total variance of the outcome variable.

- sigma2_x:

  Variance of the level-1 predictor (covariate).

- icc:

  Intraclass correlation of the outcome.

- width:

  Target full CI width on the slope.

- cluster_size:

  Per-cluster sample size (number of level-1 units per level-2 unit).
  Default `20L`.

- conf_level:

  Confidence level. Default `0.95`.

## Value

A `data.frame` with rows for the recommended number of clusters
`necessary_n_clusters`, the implied total sample size `total_N`
(`necessary_n_clusters * cluster_size`), the target `width`, the
intraclass correlation `icc`, the `cluster_size`, and the resulting
`ci_width_expected` (the expected full CI width at the recommended
size).

## Details

**Variance of the slope.** For a level-1 predictor centered within
cluster, the asymptotic variance of \\\hat\beta\\ is approximately
\$\$\mathrm{Var}(\hat\beta) \\\approx\\ \frac{\sigma^2_y (1 - \rho_I)}{N
\sigma^2_x},\$\$ where \\N = n\_{\mathrm{clusters}} \cdot m\\ is the
total number of level-1 units, \\m\\ is the cluster size, and \\\rho_I\\
the intraclass correlation. Because within-cluster centering removes the
cluster-level variation from the predictor, the design effect \\1 +
(m - 1) \rho_I\\ that inflates the variance of a cluster-level estimand
does not appear here; clustering enters only through the residual
variance \\\sigma^2_y (1 - \rho_I)\\. The function inverts this
expression for \\N\\. No anticipated slope value is needed: \\\beta\\
does not appear in the variance, so the recommended number of clusters
is the same whatever the slope.

**Scope.** Planning is for the most common single-level covariate case
(random intercept, fixed slope, level-1 predictor centered within
cluster). For cross-level interactions or random slopes, the variance
formula changes and a Monte Carlo planner should be used instead
(Schoemann, Boulton, & Short, 2017).

## References

Kelley, K. (2007). Confidence intervals for standardized effect sizes:
Theory, application, and implementation. *Journal of Statistical
Software, 20*(8), 1–24.
[doi:10.18637/jss.v020.i08](https://doi.org/10.18637/jss.v020.i08)

Raudenbush, S. W., & Liu, X.-F. (2001). Effects of study duration,
frequency of observation, and sample size on power in studies of group
differences in polynomial change. *Psychological Methods, 6*(4),
387–401.
[doi:10.1037/1082-989X.6.4.387](https://doi.org/10.1037/1082-989X.6.4.387)

Schoemann, A. M., Boulton, A. J., & Short, S. D. (2017). Determining
power and sample size for simple and complex mediation models. *Social
Psychological and Personality Science, 8*(4), 379–386.
[doi:10.1177/1948550617715068](https://doi.org/10.1177/1948550617715068)

Snijders, T. A. B., & Bosker, R. J. (2012). *Multilevel analysis: An
introduction to basic and advanced multilevel modeling* (2nd ed.). Sage.

## See also

[`ss_power_mixed_effects`](https://yelleknek.github.io/DMAR/reference/ss_power_mixed_effects.md),
[`var_icc`](https://yelleknek.github.io/DMAR/reference/var_icc.md),
[`ss_aipe_icc`](https://yelleknek.github.io/DMAR/reference/ss_aipe_icc.md)

[`design_consequences`](https://yelleknek.github.io/DMAR/reference/design_consequences.md)
for what a chosen design delivers: power, the Type S (sign) and Type M
(exaggeration) errors of the significance filter, and the expected
confidence interval width.

Other sample size for power:
[`power_fisher_exact()`](https://yelleknek.github.io/DMAR/reference/power_fisher_exact.md),
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
[`ss_power_smd()`](https://yelleknek.github.io/DMAR/reference/ss_power_smd.md),
[`ss_power_split_plot_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_split_plot_anova.md)

Other mixed models:
[`R2_mixed_effects()`](https://yelleknek.github.io/DMAR/reference/R2_mixed_effects.md),
[`R2_mixed_effects_decomposition()`](https://yelleknek.github.io/DMAR/reference/R2_mixed_effects_decomposition.md),
[`icc_lmer()`](https://yelleknek.github.io/DMAR/reference/icc_lmer.md),
[`manova_split_plot()`](https://yelleknek.github.io/DMAR/reference/manova_split_plot.md),
[`mixed_anova()`](https://yelleknek.github.io/DMAR/reference/mixed_anova.md),
[`ss_aipe_mixed_effects_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mixed_effects_sensitivity.md),
[`ss_power_mixed_effects()`](https://yelleknek.github.io/DMAR/reference/ss_power_mixed_effects.md),
[`ss_power_split_plot_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_split_plot_anova.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# 1. Plan a two-level study with cluster size 20, ICC = 0.10,
#        sigma_y = 1, sigma_x = 1, target CI width = 0.20:
ss_aipe_mixed_effects(sigma2_y = 1, sigma2_x = 1, icc = 0.10,
                           width = 0.20, cluster_size = 20)
#>  term                 value
#>  necessary_n_clusters 18   
#>  total_N              360  
#>  width                0.2  
#>  icc                  0.1  
#>  cluster_size         20   
#>  ci_width_expected    0.196
#> 
#> Confidence level: 95%

# 2. The same study with stronger clustering (ICC = 0.20):
ss_aipe_mixed_effects(sigma2_y = 1, sigma2_x = 1, icc = 0.20,
                           width = 0.20, cluster_size = 20)
#>  term                 value
#>  necessary_n_clusters 16   
#>  total_N              320  
#>  width                0.2  
#>  icc                  0.2  
#>  cluster_size         20   
#>  ci_width_expected    0.196
#> 
#> Confidence level: 95%
```
