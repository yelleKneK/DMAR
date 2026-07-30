# Sensitivity Analysis for Sample Size Planning From the AIPE Perspective for an Unstandardized Contrast

Quantifies how much misspecification of the population error variance
distorts an AIPE-based sample size plan for an unstandardized contrast
of means. Because the half-width of the confidence interval on \\\psi =
\sum_j c_j \mu_j\\ depends on the error variance, the contrast weights,
and the per-group sample size, but not on the value of \\\psi\\ itself,
this sensitivity analysis varies the planning value of the error
variance. On each replication the function simulates \\n\\ observations
per group from a normal population with variance `true_error_variance`,
builds the confidence interval via
[`ci_c`](https://yelleknek.github.io/DMAR/reference/ci_c.md), and
summarizes the realized widths and coverage of `true_psi`.

## Usage

``` r
ss_aipe_c_sensitivity(
  true_error_variance = NULL,
  estimated_error_variance = NULL,
  c_weights,
  width,
  true_psi = 0,
  n_per_group = NULL,
  conf_level = 0.95,
  assurance = NULL,
  G = 1000,
  print_iter = FALSE,
  save = FALSE,
  filename = "ss_aipe_c_sensitivity_result.csv"
)
```

## Arguments

- true_error_variance:

  Population error variance (the data generating value). Must be
  positive.

- estimated_error_variance:

  Error variance used to plan the study (the value passed to
  [`ss_aipe_c`](https://yelleknek.github.io/DMAR/reference/ss_aipe_c.md)).
  Supply this or `n_per_group` but not both.

- c_weights:

  Contrast weight vector. Must sum to zero.

- width:

  Desired full width of the confidence interval on the unstandardized
  contrast.

- true_psi:

  Population value of the contrast; the simulator places group means
  such that \\\sum_j c_j \mu_j = \\`true_psi`. The width of the interval
  does not depend on this value but the realized coverage of `true_psi`
  does. Default `0`.

- n_per_group:

  Per-group sample size to evaluate (incompatible with
  `estimated_error_variance`); when used, the planner is bypassed.

- conf_level:

  Confidence level (default `0.95`).

- assurance:

  Optional probability that the realized interval is no wider than
  `width`; passed to
  [`ss_aipe_c`](https://yelleknek.github.io/DMAR/reference/ss_aipe_c.md)
  when resolving the planned sample size.

- G:

  Number of Monte Carlo replications (default 1000).

- print_iter:

  Logical. Print the iteration index after each replication (helpful for
  long runs); default `FALSE`.

- save:

  Logical. If `TRUE` the per-replication results are appended to
  `filename`; default `FALSE`.

- filename:

  Path used when `save = TRUE`.

## Value

A `data.frame` with rows for mean / median / SD of the realized
estimator and interval width, the proportion of intervals at or below
`width`, the tail-specific and overall empirical non-coverage of
`true_psi`, and the input echoes (per-group sample size, total sample
size, true and estimated error variances, width, confidence level).

## References

Lai, K., & Kelley, K. (2012). Accuracy in parameter estimation for
ANCOVA and ANOVA contrasts: Sample size planning via narrow confidence
intervals. *British Journal of Mathematical and Statistical Psychology,
65*, 350–370.
[doi:10.1111/j.2044-8317.2011.02029.x](https://doi.org/10.1111/j.2044-8317.2011.02029.x)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 4 on individual comparisons.)

## See also

[`ss_aipe_c`](https://yelleknek.github.io/DMAR/reference/ss_aipe_c.md),
[`ci_c`](https://yelleknek.github.io/DMAR/reference/ci_c.md),
[`ss_aipe_sc_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_sc_sensitivity.md)

[`design_consequences`](https://yelleknek.github.io/DMAR/reference/design_consequences.md)
for what a chosen design delivers: power, the Type S (sign) and Type M
(exaggeration) errors of the significance filter, and the expected
confidence interval width.

Other AIPE sample size planning:
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
[`ss_aipe_semipartial_r_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_semipartial_r_sensitivity.md),
[`ss_aipe_tost_smd_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_tost_smd_sensitivity.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# \donttest{
# Monte Carlo sweep; G is small here so the example runs quickly.
# Well-specified: planner used error_variance = 4, truth is 4.
set.seed(113)
ss_aipe_c_sensitivity(
  true_error_variance      = 4,
  estimated_error_variance = 4,
  c_weights = c(-1, 0, 1),
  width = 1, G = 50, print_iter = FALSE
)
#>  term                     value   
#>  mean_psi                 -0.01   
#>  median_psi               -0.00595
#>  sd_psi                   0.213   
#>  mean_ci_width            0.995   
#>  median_ci_width          1       
#>  sd_ci_width              0.0321  
#>  pct_ci_less_w            0.52    
#>  pct_ci_miss_low          0       
#>  pct_ci_miss_high         0.02    
#>  total_type_I_error       0.02    
#>  n_per_group              124     
#>  total_N                  372     
#>  true_error_variance      4       
#>  estimated_error_variance 4       
#>  true_psi                 0       
#>  width                    1       
#>  conf_level               0.95    
#> 
#> Confidence level: 95%

# Misspecified: planner used 4, truth is 9. Realized widths inflate.
set.seed(113)
ss_aipe_c_sensitivity(
  true_error_variance      = 9,
  estimated_error_variance = 4,
  c_weights = c(-1, 0, 1),
  width = 1, G = 50, print_iter = FALSE
)
#>  term                     value   
#>  mean_psi                 -0.0151 
#>  median_psi               -0.00892
#>  sd_psi                   0.32    
#>  mean_ci_width            1.49    
#>  median_ci_width          1.5     
#>  sd_ci_width              0.0481  
#>  pct_ci_less_w            0       
#>  pct_ci_miss_low          0       
#>  pct_ci_miss_high         0.02    
#>  total_type_I_error       0.02    
#>  n_per_group              124     
#>  total_N                  372     
#>  true_error_variance      9       
#>  estimated_error_variance 4       
#>  true_psi                 0       
#>  width                    1       
#>  conf_level               0.95    
#> 
#> Confidence level: 95%
# }
```
