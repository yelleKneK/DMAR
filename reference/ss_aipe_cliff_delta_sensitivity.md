# Sensitivity Analysis for Sample Size Planning From the AIPE Perspective for Cliff's Delta

Quantifies how much misspecification of the population Cliff's delta
(\\\delta = \Pr(X \> Y) - \Pr(X \< Y)\\) distorts an AIPE-based sample
size plan. On each replication the function simulates two independent
samples whose population Cliff's delta equals `true_delta`, computes the
sample
[`cliff_delta`](https://yelleknek.github.io/DMAR/reference/cliff_delta.md)
and its CI, and summarizes the realized widths and coverage.

**Data generating mechanism.** The simulator draws each sample from a
normal distribution and chooses the mean shift so that the implied
Cliff's delta equals `true_delta`. For normal samples \\\delta = 2
\Phi(\Delta/\sqrt{2}) - 1\\ where \\\Delta\\ is the standardized mean
difference, so the simulator sets \\\Delta = \sqrt{2} \cdot
\Phi^{-1}((1 + \delta)/2)\\.

## Usage

``` r
ss_aipe_cliff_delta_sensitivity(
  true_delta = NULL,
  estimated_delta = NULL,
  ratio = 1,
  width,
  specified_N = NULL,
  conf_level = 0.95,
  assurance = NULL,
  G = 1000,
  print_iter = FALSE,
  save = FALSE,
  filename = "ss_aipe_cliff_delta_sensitivity_result.csv"
)
```

## Arguments

- true_delta:

  Population Cliff's delta (the data generating value); in \\(-1, 1)\\.

- estimated_delta:

  Planning value passed to
  [`ss_aipe_cliff_delta`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cliff_delta.md);
  supply this or `specified_N` but not both.

- ratio:

  Allocation ratio \\n_1 / n_2\\ (default 1).

- width:

  Desired full width of the CI on Cliff's delta.

- specified_N:

  Total sample size to evaluate (split per `ratio`).

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

A `data.frame` with rows for mean / median / SD of the realized Cliff's
delta and CI width, the proportion of intervals at or below `width`,
tail-specific and overall non-coverage of `true_delta`, and the input
echoes, including `assurance` (present only when an assurance was
supplied).

## References

Cliff, N. (1993). Dominance statistics: Ordinal analyses to answer
ordinal questions. *Psychological Bulletin, 114*(3), 494–509.
[doi:10.1037/0033-2909.114.3.494](https://doi.org/10.1037/0033-2909.114.3.494)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.

## See also

[`ss_aipe_cliff_delta`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cliff_delta.md),
[`cliff_delta`](https://yelleknek.github.io/DMAR/reference/cliff_delta.md)

[`design_consequences`](https://yelleknek.github.io/DMAR/reference/design_consequences.md)
for what a chosen design delivers: power, the Type S (sign) and Type M
(exaggeration) errors of the significance filter, and the expected
confidence interval width.

Other AIPE sample size planning:
[`ss_aipe_c_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_c_sensitivity.md),
[`ss_aipe_cliff_delta()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cliff_delta.md),
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
set.seed(113)
# Small G keeps the Monte Carlo sweep fast; raise G for a real plan.
ss_aipe_cliff_delta_sensitivity(
  true_delta = 0.30, estimated_delta = 0.30,
  width = 0.30, G = 25, print_iter = FALSE
)
#>  term               value  
#>  mean_cliff_delta   0.31   
#>  median_cliff_delta 0.303  
#>  sd_cliff_delta     0.068  
#>  mean_ci_width      0.24   
#>  median_ci_width    0.241  
#>  sd_ci_width        0.00786
#>  pct_ci_less_w      1      
#>  pct_ci_miss_low    0.08   
#>  pct_ci_miss_high   0      
#>  total_type_I_error 0.08   
#>  n_1                156    
#>  n_2                156    
#>  total_N            312    
#>  true_delta         0.3    
#>  estimated_delta    0.3    
#>  ratio              1      
#>  width              0.3    
#>  conf_level         0.95   
#> 
#> Confidence level: 95%
```
