# Sensitivity Analysis for Sample Size Planning From the AIPE Perspective for a Reliability Coefficient

Quantifies how much misspecification of the population reliability
coefficient distorts an AIPE-based sample size plan for the
composite-score reliability. On each replication the function simulates
an *n* \\\times\\ *i* item-by-subject data matrix from a single-factor
parallel-tests model whose population reliability of the sum score
equals `true_reliability`, fits the requested estimator (alpha or omega)
via the corresponding `reliability_*` function with the supplied
`ci_method`, and records the realized reliability estimate and its
confidence interval.

**Population model.** Each item has a single common-factor loading and
uncorrelated unique error. With per-item variance normalized to 1, the
loading and unique variance are chosen so the Cronbach-style sum-score
reliability equals `true_reliability`: \$\$\lambda^2 \\=\\
\frac{\rho}{i(1 - \rho) + \rho}, \qquad \psi^2 \\=\\ 1 - \lambda^2,\$\$
where \\\rho = \\`true_reliability` and \\i\\ is the item count. Item
scores are \\y\_{ij} = \lambda T_i + e\_{ij}\\, with \\T_i \sim N(0,
1)\\ and \\e\_{ij} \sim N(0, \psi^2)\\.

## Usage

``` r
ss_aipe_reliability_sensitivity(
  true_reliability = NULL,
  estimated_reliability = NULL,
  i,
  width,
  specified_N = NULL,
  estimator = c("alpha", "omega"),
  ci_method = NULL,
  conf_level = 0.95,
  assurance = NULL,
  G = 1000,
  print_iter = FALSE,
  save = FALSE,
  filename = "ss_aipe_reliability_sensitivity_result.csv"
)
```

## Arguments

- true_reliability:

  Population reliability coefficient (in \\\[0, 1)\\).

- estimated_reliability:

  Reliability used to plan the study; the function passes the implied
  lambda / psi^2 to
  [`ss_aipe_reliability`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reliability.md).

- i:

  Number of items in the composite.

- width:

  Desired full width of the CI on reliability.

- specified_N:

  Sample size to evaluate (incompatible with `estimated_reliability`).

- estimator:

  One of `"alpha"` (default; coefficient alpha via
  [`reliability_alpha`](https://yelleknek.github.io/DMAR/reference/reliability_alpha.md))
  or `"omega"` (composite reliability via
  [`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md)).
  For a parallel-tests population the two coincide; differences in
  sample estimates reflect estimator-specific finite-sample bias and CI
  behavior.

- ci_method:

  CI method passed to the estimator. Default `"bonett"` for alpha and
  `"mlr"` for omega.

- conf_level:

  Confidence level (default `0.95`).

- assurance:

  Optional assurance probability passed to the planner.

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
reliability and CI width, the proportion of intervals at or below
`width`, tail-specific and overall non-coverage of `true_reliability`,
and the input echoes.

## References

Kelley, K., & Pornprasertmanit, S. (2016). Confidence intervals for
population reliability coefficients: Evaluation of methods,
recommendations, and software for composite measures. *Psychological
Methods, 21*, 69–92.
[doi:10.1037/a0040086](https://doi.org/10.1037/a0040086)

Terry, L. J., & Kelley, K. (2012). Sample size planning for composite
reliability coefficients: Accuracy in parameter estimation via narrow
confidence intervals. *British Journal of Mathematical and Statistical
Psychology, 65*, 371–401.
[doi:10.1111/j.2044-8317.2011.02030.x](https://doi.org/10.1111/j.2044-8317.2011.02030.x)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.

## See also

[`ss_aipe_reliability`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reliability.md),
[`reliability_alpha`](https://yelleknek.github.io/DMAR/reference/reliability_alpha.md),
[`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md)

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
[`ss_aipe_semipartial_r()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_semipartial_r.md),
[`ss_aipe_semipartial_r_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_semipartial_r_sensitivity.md),
[`ss_aipe_tost_smd_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_tost_smd_sensitivity.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# \donttest{
# Reduced Monte Carlo sweep (small G) so the example runs quickly;
# raise G for a production sensitivity analysis.
set.seed(113)
ss_aipe_reliability_sensitivity(
  true_reliability      = 0.80,
  estimated_reliability = 0.80,
  i = 4, width = 0.15,
  estimator = "alpha",
  G = 20, print_iter = FALSE
)
#>  term                  value 
#>  mean_reliability      0.782 
#>  median_reliability    0.792 
#>  sd_reliability        0.0412
#>  mean_ci_width         0.169 
#>  median_ci_width       0.161 
#>  sd_ci_width           0.0319
#>  pct_ci_less_w         0.4   
#>  pct_ci_miss_low       0     
#>  pct_ci_miss_high      0.05  
#>  total_type_I_error    0.05  
#>  total_N               74    
#>  items                 4     
#>  true_reliability      0.8   
#>  estimated_reliability 0.8   
#>  width                 0.15  
#>  conf_level            0.95  
#> 
#> Confidence level: 95%
# }
```
