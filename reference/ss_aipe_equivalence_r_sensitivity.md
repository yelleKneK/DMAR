# Sensitivity Analysis for Sample Size Planning From the AIPE Perspective for an Equivalence-Test Correlation

Quantifies how much misspecification of the population correlation
distorts an AIPE-based sample size plan for the two-one-sided-tests
(TOST) confidence interval on the Pearson correlation. On each
replication the function simulates *N* bivariate normal pairs with
population correlation `true_r`, computes the sample correlation and its
Fisher's \\Z\\ confidence interval via
[`ci_r`](https://yelleknek.github.io/DMAR/reference/ci_correlation.md),
and summarizes the realized widths and the proportion of replications in
which the computed interval falls entirely inside (`equivalent`) the
specified equivalence bounds.

## Usage

``` r
ss_aipe_equivalence_r_sensitivity(
  true_r = 0,
  estimated_r = NULL,
  width,
  rho_lower = NULL,
  rho_upper = NULL,
  specified_N = NULL,
  conf_level = 0.95,
  assurance = NULL,
  G = 1000,
  print_iter = FALSE,
  save = FALSE,
  filename = "ss_aipe_equivalence_r_sensitivity_result.csv"
)
```

## Arguments

- true_r:

  Population correlation (the data generating value). Defaults to `0`
  (no association, the exact equivalence case).

- estimated_r:

  Planning value of the population correlation passed to
  [`ss_aipe_equivalence_r`](https://yelleknek.github.io/DMAR/reference/ss_aipe_equivalence_r.md);
  supply this or `specified_N` but not both.

- width:

  Desired full width of the two-sided CI on the correlation.

- rho_lower, rho_upper:

  Equivalence bounds on the correlation, as positive magnitudes with the
  same meaning as in
  [`equivalence_r`](https://yelleknek.github.io/DMAR/reference/equivalence_r.md):
  the region is \\(-\rho_L, +\rho_U)\\. `rho_upper` is required;
  `rho_lower` defaults to `rho_upper` (a symmetric region). The
  simulator records whether the realized CI falls entirely inside the
  region.

- specified_N:

  Sample size to evaluate.

- conf_level:

  Confidence level (default `0.95`).

- assurance:

  Optional assurance probability.

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
correlation and CI width, the proportion of intervals at or below
`width`, tail-specific and overall non-coverage of `true_r`, the
proportion of intervals classified as `equivalent` (CI fully inside the
bounds), and the input echoes, including `assurance` (present only when
an assurance was supplied).

## References

Counsell, A., & Cribbie, R. A. (2015). Equivalence tests for comparing
correlation and regression coefficients. *British Journal of
Mathematical and Statistical Psychology, 68*(2), 292–309.
[doi:10.1111/bmsp.12045](https://doi.org/10.1111/bmsp.12045)

Goertzen, J. R., & Cribbie, R. A. (2010). Detecting a lack of
association: An equivalence testing approach. *British Journal of
Mathematical and Statistical Psychology, 63*(3), 527–537.
[doi:10.1348/000711009X475853](https://doi.org/10.1348/000711009X475853)

## See also

[`ss_aipe_equivalence_r`](https://yelleknek.github.io/DMAR/reference/ss_aipe_equivalence_r.md),
[`equivalence_r`](https://yelleknek.github.io/DMAR/reference/equivalence_r.md),
[`ss_aipe_r_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_r_sensitivity.md),
[`ss_aipe_equivalence_smd_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_equivalence_smd_sensitivity.md)

Other AIPE sample size planning:
[`ss_aipe_c_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_c_sensitivity.md),
[`ss_aipe_cliff_delta()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cliff_delta.md),
[`ss_aipe_cliff_delta_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cliff_delta_sensitivity.md),
[`ss_aipe_composite_sem()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_composite_sem.md),
[`ss_aipe_equivalence_r()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_equivalence_r.md),
[`ss_aipe_equivalence_smd()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_equivalence_smd.md),
[`ss_aipe_equivalence_smd_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_equivalence_smd_sensitivity.md),
[`ss_aipe_icc()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_icc.md),
[`ss_aipe_icc_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_icc_sensitivity.md),
[`ss_aipe_indirect_effect()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_indirect_effect.md),
[`ss_aipe_indirect_effect_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_indirect_effect_sensitivity.md),
[`ss_aipe_mixed_effects_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mixed_effects_sensitivity.md),
[`ss_aipe_omega_squared()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_omega_squared.md),
[`ss_aipe_omega_squared_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_omega_squared_sensitivity.md),
[`ss_aipe_partial_r()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_partial_r.md),
[`ss_aipe_partial_r_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_partial_r_sensitivity.md),
[`ss_aipe_pcm_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_pcm_sensitivity.md),
[`ss_aipe_r()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_r.md),
[`ss_aipe_r_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_r_sensitivity.md),
[`ss_aipe_reliability_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reliability_sensitivity.md),
[`ss_aipe_semipartial_r()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_semipartial_r.md),
[`ss_aipe_semipartial_r_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_semipartial_r_sensitivity.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Reduced Monte Carlo sweep (small G) for a fast, illustrative run.
set.seed(113)
ss_aipe_equivalence_r_sensitivity(
  true_r      = 0.0,
  estimated_r = 0.0,
  width       = 0.30,
  rho_upper   = 0.20,
  G = 50, print_iter = FALSE
)
#>  term               value  
#>  mean_r             0.00852
#>  median_r           0.00378
#>  sd_r               0.0758 
#>  mean_ci_width      0.298  
#>  median_ci_width    0.298  
#>  sd_ci_width        0.00206
#>  pct_ci_less_w      1      
#>  pct_equivalent     0.5    
#>  pct_ci_miss_low    0.02   
#>  pct_ci_miss_high   0      
#>  total_type_I_error 0.02   
#>  total_N            172    
#>  true_r             0      
#>  estimated_r        0      
#>  width              0.3    
#>  conf_level         0.95   
#>  rho_lower          -0.2   
#>  rho_upper          0.2    
#> 
#> Confidence level: 95%
```
