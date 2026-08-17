# Sensitivity Analysis for Sample Size Planning From the AIPE Perspective for a Partial Correlation

Quantifies how much misspecification of the population partial
correlation distorts an AIPE-based sample size plan. The function
constructs an \\(J + 1) \times (J + 1)\\ population covariance matrix
whose implied partial correlation between \\Y\\ and \\X_1\\ (controlling
for \\X_2, \ldots, X_J\\) equals `true_rho`, then on each replication
draws an *n*-row sample from the corresponding multivariate normal
distribution and computes the sample partial correlation and its
Fisher's \\Z\\ CI.

## Usage

``` r
ss_aipe_partial_r_sensitivity(
  true_rho = NULL,
  estimated_rho = NULL,
  J,
  width,
  specified_N = NULL,
  conf_level = 0.95,
  assurance = NULL,
  G = 1000,
  print_iter = FALSE,
  save = FALSE,
  filename = "ss_aipe_partial_r_sensitivity_result.csv"
)
```

## Arguments

- true_rho:

  Population partial correlation between *Y* and \\X_1\\ controlling for
  \\X_2, \ldots, X_J\\; must lie in \\(-1, 1)\\.

- estimated_rho:

  Planning value of the partial correlation passed to
  [`ss_aipe_partial_r`](https://yelleknek.github.io/DMAR/reference/ss_aipe_partial_r.md);
  supply this or `specified_N` but not both.

- J:

  Total number of predictors (so the partial correlation is between *Y*
  and one of the *J* predictors, partialing out the other \\J - 1\\).
  Must be at least 1.

- width:

  Desired full width of the CI on the partial correlation.

- specified_N:

  Sample size to evaluate (incompatible with `estimated_rho`).

- conf_level:

  Confidence level (default `0.95`).

- assurance:

  Optional assurance probability passed to
  [`ss_aipe_partial_r`](https://yelleknek.github.io/DMAR/reference/ss_aipe_partial_r.md).

- G:

  Number of Monte Carlo replications (default 1000).

- print_iter:

  Logical. Print iteration index per replication.

- save:

  Logical. If `TRUE` write per-replication results to `filename`.

- filename:

  Path used when `save = TRUE`.

## Value

A `data.frame` with rows for the realized partial correlation, the
interval width, the proportion of intervals at or below `width`,
tail-specific and overall non-coverage of `true_rho`, and the input
echoes, including `assurance` (present only when an assurance was
supplied).

## References

Kelley, K., & Maxwell, S. E. (2003). Sample size for multiple
regression: Obtaining regression coefficients that are accurate, not
simply significant. *Psychological Methods, 8*(3), 305–321.
[doi:10.1037/1082-989X.8.3.305](https://doi.org/10.1037/1082-989X.8.3.305)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.

## See also

[`ss_aipe_partial_r`](https://yelleknek.github.io/DMAR/reference/ss_aipe_partial_r.md),
[`ss_aipe_semipartial_r_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_semipartial_r_sensitivity.md)

[`design_consequences`](https://yelleknek.github.io/DMAR/reference/design_consequences.md)
for what a chosen design delivers: power, the Type S (sign) and Type M
(exaggeration) errors of the significance filter, and the expected
confidence interval width.

Other AIPE sample size planning:
[`ss_aipe_c_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_c_sensitivity.md),
[`ss_aipe_cliff_delta()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cliff_delta.md),
[`ss_aipe_cliff_delta_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cliff_delta_sensitivity.md),
[`ss_aipe_composite_sem()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_composite_sem.md),
[`ss_aipe_equivalence_r()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_equivalence_r.md),
[`ss_aipe_equivalence_r_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_equivalence_r_sensitivity.md),
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
# Reduced replications and a wide target interval keep this fast.
set.seed(113)
ss_aipe_partial_r_sensitivity(
  true_rho = 0.40, estimated_rho = 0.40, J = 3, width = 0.40,
  G = 50, print_iter = FALSE
)
#>  term               value 
#>  mean_partial_r     0.413 
#>  median_partial_r   0.429 
#>  sd_partial_r       0.108 
#>  mean_ci_width      0.391 
#>  median_ci_width    0.39  
#>  sd_ci_width        0.0389
#>  pct_ci_less_w      0.58  
#>  pct_ci_miss_low    0.04  
#>  pct_ci_miss_high   0.02  
#>  total_type_I_error 0.06  
#>  total_N            72    
#>  J                  3     
#>  true_rho           0.4   
#>  estimated_rho      0.4   
#>  width              0.4   
#>  conf_level         0.95  
#> 
#> Confidence level: 95%
```
