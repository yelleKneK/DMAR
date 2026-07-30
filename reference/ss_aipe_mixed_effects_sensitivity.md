# Sensitivity Analysis for Sample Size Planning From the AIPE Perspective for a Mixed-Effects Fixed Effect

Quantifies how much misspecification of the variance components
(\\\sigma^2_Y\\, \\\sigma^2_X\\, and the intraclass correlation
\\\mathrm{icc}\\) distorts an AIPE-based sample size plan for a
cluster-level fixed effect under a two-level random-intercept model. On
each replication the function simulates *K* clusters of `cluster_size`
observations each from \$\$Y\_{ki} = \beta\\X_k + u_k +
\epsilon\_{ki},\$\$ with \\X_k \sim N(0, \sigma^2_X)\\, \\u_k \sim N(0,
\mathrm{icc}\cdot\sigma^2_Y)\\, and \\\epsilon\_{ki} \sim N(0, (1 -
\mathrm{icc})\sigma^2_Y)\\. The model is then refit by either
[`lme4::lmer`](https://rdrr.io/pkg/lme4/man/lmer.html) (if available) or
by GLS-by-cluster aggregation, and a Wald CI on the fixed effect is
recorded.

## Usage

``` r
ss_aipe_mixed_effects_sensitivity(
  true_sigma2_y = NULL,
  true_sigma2_x = NULL,
  true_icc = NULL,
  true_beta = 0,
  estimated_sigma2_y = NULL,
  estimated_sigma2_x = NULL,
  estimated_icc = NULL,
  width,
  cluster_size = 20L,
  specified_K = NULL,
  conf_level = 0.95,
  G = 1000,
  print_iter = FALSE,
  save = FALSE,
  filename = "ss_aipe_mixed_effects_sensitivity_result.csv"
)
```

## Arguments

- true_sigma2_y:

  Population total variance of *Y*.

- true_sigma2_x:

  Population variance of the cluster-level predictor.

- true_icc:

  Population intraclass correlation (between-cluster share of total
  variance).

- true_beta:

  Population fixed-effect slope (default `0`).

- estimated_sigma2_y, estimated_sigma2_x, estimated_icc:

  Planning values passed to
  [`ss_aipe_mixed_effects`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mixed_effects.md);
  supply all three or `specified_K` but not both.

- width:

  Desired full width of the CI on the fixed effect.

- cluster_size:

  Number of observations per cluster (assumed balanced).

- specified_K:

  Number of clusters to evaluate.

- conf_level:

  Confidence level (default `0.95`).

- G:

  Number of Monte Carlo replications.

- print_iter:

  Logical.

- save:

  Logical. Save per-replication CSV.

- filename:

  Path used when `save = TRUE`.

## Value

A `data.frame` with rows for mean / median / SD of the realized
fixed-effect estimate and CI width, the proportion of intervals at or
below `width`, tail-specific and overall non-coverage of `true_beta`,
and the input echoes.

## References

McNeish, D., & Kelley, K. (2019). Fixed effects versus mixed effects
models for clustered data: Reviewing the approaches, disentangling the
differences, and making recommendations. *Psychological Methods, 24*,
20–35. [doi:10.1037/met0000182](https://doi.org/10.1037/met0000182)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapters 15 and 16 on mixed-effects models.)

## See also

[`ss_aipe_mixed_effects`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mixed_effects.md),
[`icc_lmer`](https://yelleknek.github.io/DMAR/reference/icc_lmer.md)

[`design_consequences`](https://yelleknek.github.io/DMAR/reference/design_consequences.md)
for what a chosen design delivers: power, the Type S (sign) and Type M
(exaggeration) errors of the significance filter, and the expected
confidence interval width.

Other AIPE sample size planning:
[`ss_aipe_c_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_c_sensitivity.md),
[`ss_aipe_cliff_delta()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cliff_delta.md),
[`ss_aipe_cliff_delta_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cliff_delta_sensitivity.md),
[`ss_aipe_composite_sem()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_composite_sem.md),
[`ss_aipe_icc()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_icc.md),
[`ss_aipe_icc_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_icc_sensitivity.md),
[`ss_aipe_indirect_effect()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_indirect_effect.md),
[`ss_aipe_indirect_effect_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_indirect_effect_sensitivity.md),
[`ss_aipe_omega_squared()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_omega_squared.md),
[`ss_aipe_omega_squared_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_omega_squared_sensitivity.md),
[`ss_aipe_partial_r()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_partial_r.md),
[`ss_aipe_partial_r_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_partial_r_sensitivity.md),
[`ss_aipe_pcm_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_pcm_sensitivity.md),
[`ss_aipe_reliability_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reliability_sensitivity.md),
[`ss_aipe_semipartial_r()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_semipartial_r.md),
[`ss_aipe_semipartial_r_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_semipartial_r_sensitivity.md),
[`ss_aipe_tost_smd_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_tost_smd_sensitivity.md)

Other mixed models:
[`R2_mixed_effects()`](https://yelleknek.github.io/DMAR/reference/R2_mixed_effects.md),
[`R2_mixed_effects_decomposition()`](https://yelleknek.github.io/DMAR/reference/R2_mixed_effects_decomposition.md),
[`icc_lmer()`](https://yelleknek.github.io/DMAR/reference/icc_lmer.md),
[`manova_split_plot()`](https://yelleknek.github.io/DMAR/reference/manova_split_plot.md),
[`mixed_anova()`](https://yelleknek.github.io/DMAR/reference/mixed_anova.md),
[`ss_aipe_mixed_effects()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mixed_effects.md),
[`ss_power_mixed_effects()`](https://yelleknek.github.io/DMAR/reference/ss_power_mixed_effects.md),
[`ss_power_split_plot_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_split_plot_anova.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# \donttest{
# Monte Carlo sensitivity check, reduced sizes for a fast example.
set.seed(113)
ss_aipe_mixed_effects_sensitivity(
  true_sigma2_y = 1, true_sigma2_x = 1, true_icc = 0.10,
  true_beta = 0.30,
  specified_K = 25, cluster_size = 10,
  width = 0.40,
  G = 25, print_iter = FALSE
)
#>  term               value 
#>  mean_beta          0.325 
#>  median_beta        0.322 
#>  sd_beta            0.108 
#>  mean_ci_width      0.384 
#>  median_ci_width    0.383 
#>  sd_ci_width        0.0694
#>  pct_ci_less_w      0.64  
#>  pct_ci_miss_low    0     
#>  pct_ci_miss_high   0.08  
#>  total_type_I_error 0.08  
#>  n_clusters         25    
#>  cluster_size       10    
#>  total_N            250   
#>  true_sigma2_y      1     
#>  true_sigma2_x      1     
#>  true_icc           0.1   
#>  true_beta          0.3   
#>  estimated_sigma2_y <NA>  
#>  estimated_sigma2_x <NA>  
#>  estimated_icc      <NA>  
#>  width              0.4   
#>  conf_level         0.95  
#> 
#> Confidence level: 95%
# }
```
