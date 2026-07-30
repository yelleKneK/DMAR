# Sensitivity Analysis for Sample Size Planning From the AIPE Perspective for an Equivalence-Test SMD

Quantifies how much misspecification of the population standardized mean
difference distorts an AIPE-based sample size plan for the
two-one-sided-tests (TOST) confidence interval on the SMD. On each
replication the function simulates two normal groups of size *n* per
group with population standardized mean difference `true_smd`, computes
the SMD and its noncentral *t* confidence interval via
[`ci_smd`](https://yelleknek.github.io/DMAR/reference/ci_smd.md), and
summarizes the realized widths and the proportion of replications in
which the computed interval falls entirely inside (`equivalent`) the
specified equivalence bounds.

## Usage

``` r
ss_aipe_tost_smd_sensitivity(
  true_smd = 0,
  estimated_smd = NULL,
  width,
  equivalence_lower = NULL,
  equivalence_upper = NULL,
  n_per_group = NULL,
  conf_level = 0.95,
  assurance = NULL,
  G = 1000,
  print_iter = FALSE,
  save = FALSE,
  filename = "ss_aipe_tost_smd_sensitivity_result.csv"
)
```

## Arguments

- true_smd:

  Population standardized mean difference (the data generating value).
  Defaults to `0` (perfect equivalence).

- estimated_smd:

  Planning value of the population SMD passed to
  [`ss_aipe_tost_smd`](https://yelleknek.github.io/DMAR/reference/ss_aipe_tost_smd.md);
  supply this or `n_per_group` but not both.

- width:

  Desired full width of the two-sided CI on the SMD.

- equivalence_lower, equivalence_upper:

  Equivalence bounds on the SMD. When supplied, the simulator records
  whether the realized CI falls entirely inside
  `[equivalence_lower, equivalence_upper]`; if both are `NULL` the
  bounds default to symmetric \\\pm 0.20\\.

- n_per_group:

  Per-group sample size to evaluate.

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

A `data.frame` with rows for mean / median / SD of the realized SMD and
CI width, the proportion of intervals at or below `width`, tail-specific
and overall non-coverage of `true_smd`, the proportion of intervals
classified as `equivalent` (CI fully inside the bounds), and the input
echoes.

## References

Kelley, K., & Rausch, J. R. (2006). Sample size planning for the
standardized mean difference: Accuracy in parameter estimation via
narrow confidence intervals. *Psychological Methods, 11*, 363–385.
[doi:10.1037/1082-989X.11.4.363](https://doi.org/10.1037/1082-989X.11.4.363)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.

## See also

[`ss_aipe_tost_smd`](https://yelleknek.github.io/DMAR/reference/ss_aipe_tost_smd.md),
[`tost_smd`](https://yelleknek.github.io/DMAR/reference/tost_smd.md),
[`ss_aipe_smd_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_smd_sensitivity.md)

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
[`ss_aipe_mixed_effects_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mixed_effects_sensitivity.md),
[`ss_aipe_omega_squared()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_omega_squared.md),
[`ss_aipe_omega_squared_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_omega_squared_sensitivity.md),
[`ss_aipe_partial_r()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_partial_r.md),
[`ss_aipe_partial_r_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_partial_r_sensitivity.md),
[`ss_aipe_pcm_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_pcm_sensitivity.md),
[`ss_aipe_reliability_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reliability_sensitivity.md),
[`ss_aipe_semipartial_r()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_semipartial_r.md),
[`ss_aipe_semipartial_r_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_semipartial_r_sensitivity.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# \donttest{
# Reduced Monte Carlo sweep (small G) for a fast, illustrative run.
set.seed(113)
ss_aipe_tost_smd_sensitivity(
  true_smd      = 0.0,
  estimated_smd = 0.0,
  width         = 0.30,
  equivalence_lower = -0.20, equivalence_upper = 0.20,
  G = 50, print_iter = FALSE
)
#>  term               value   
#>  mean_smd           -0.0106 
#>  median_smd         -0.0104 
#>  sd_smd             0.0795  
#>  mean_ci_width      0.3     
#>  median_ci_width    0.3     
#>  sd_ci_width        0.000153
#>  pct_ci_less_w      0.86    
#>  pct_equivalent     0.5     
#>  pct_ci_miss_low    0.02    
#>  pct_ci_miss_high   0.04    
#>  total_type_I_error 0.06    
#>  n_per_group        342     
#>  total_N            684     
#>  true_smd           0       
#>  estimated_smd      0       
#>  width              0.3     
#>  conf_level         0.95    
#>  equivalence_lower  -0.2    
#>  equivalence_upper  0.2     
#> 
#> Confidence level: 95%
# }
```
