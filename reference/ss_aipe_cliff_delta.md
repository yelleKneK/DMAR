# Sample Size for AIPE on Cliff's \\\delta\\

Determines the sample size needed for the confidence interval on Cliff's
(1993) \\\delta\\ (and equivalently Vargha-Delaney's \\A = (\delta + 1)
/ 2\\) to have a desired width, using the large-sample U-statistic
variance of DeLong, DeLong, & Clarke-Pearson (1988) under
exchangeability.

## Usage

``` r
ss_aipe_cliff_delta(
  delta,
  width,
  which_width = c("Full", "Lower", "Upper"),
  conf_level = 0.95,
  ratio = 1,
  assurance = NULL
)
```

## Arguments

- delta:

  Anticipated population Cliff's \\\delta\\; numeric scalar in \\(-1,
  1)\\.

- width:

  Desired full width of the CI on \\\delta\\.

- which_width:

  `"Full"` (default), `"Lower"`, or `"Upper"`.

- conf_level:

  Desired confidence level. Default `0.95`.

- ratio:

  Ratio \\n_1 / n_2\\ of the two group sample sizes. Default `1`
  (balanced).

- assurance:

  Optional. Probability that the realized CI is no wider than `width`.

## Value

A `data.frame` with rows for the recommended group sample sizes \\n_1,
n_2\\, the expected CI width, and the inputs echoed back.

## Details

**DeLong-DeLong-Clarke-Pearson variance.** Under exchangeability of
paired observations within and across groups, the asymptotic variance of
\\\hat\delta\\ for groups of sizes \\n_1, n_2\\ is approximated by the
upper bound (DeLong et al., 1988): \$\$\mathrm{Var}(\hat\delta) \\\le\\
\frac{(1 - \delta^2)}{\min(n_1, n_2)}.\$\$ Setting the half-width of a
Wald-style CI \\z\_{1-\alpha/2} \sqrt{\mathrm{Var}(\hat\delta)}\\ equal
to the target half-width and solving gives the recommended per-group
sample size. The bound is conservative; the realized CI is generally
narrower than the target.

**Allocation.** The bound is dominated by \\\min(n_1, n_2)\\, so
balanced allocation (`ratio = 1`) is approximately optimal under
standard conditions; unbalanced allocations require the larger total *N*
to achieve the same precision.

## References

Cliff, N. (1993). Dominance statistics: Ordinal analyses to answer
ordinal questions. *Psychological Bulletin, 114*(3), 494–509.
[doi:10.1037/0033-2909.114.3.494](https://doi.org/10.1037/0033-2909.114.3.494)

DeLong, E. R., DeLong, D. M., & Clarke-Pearson, D. L. (1988). Comparing
the areas under two or more correlated receiver operating characteristic
curves: A nonparametric approach. *Biometrics, 44*(3), 837–845.

Vargha, A., & Delaney, H. D. (2000). A critique and improvement of the
CL common language effect size statistics of McGraw and Wong. *Journal
of Educational and Behavioral Statistics, 25*(2), 101–132.
[doi:10.3102/10769986025002101](https://doi.org/10.3102/10769986025002101)

## See also

[`cliff_delta`](https://yelleknek.github.io/DMAR/reference/cliff_delta.md),
[`ss_aipe_partial_r`](https://yelleknek.github.io/DMAR/reference/ss_aipe_partial_r.md),
[`ss_aipe_semipartial_r`](https://yelleknek.github.io/DMAR/reference/ss_aipe_semipartial_r.md)

[`design_consequences`](https://yelleknek.github.io/DMAR/reference/design_consequences.md)
for what a chosen design delivers: power, the Type S (sign) and Type M
(exaggeration) errors of the significance filter, and the expected
confidence interval width.

Other AIPE sample size planning:
[`ss_aipe_c_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_c_sensitivity.md),
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
# 1. Plan a balanced design so the 95% CI on delta has full width
#        <= 0.20 when anticipating delta = 0.30.
ss_aipe_cliff_delta(delta = 0.30, width = 0.20)
#>  term           value
#>  n_1            350  
#>  n_2            350  
#>  necessary_N    700  
#>  expected_width 0.2  
#>  delta          0.3  
#>  ratio          1    
#>  width_target   0.2  
#>  conf_level     0.95 
#> 
#> Confidence level: 95%

# 2. Unbalanced: twice as many in group 1.
ss_aipe_cliff_delta(delta = 0.30, width = 0.20, ratio = 2)
#>  term           value
#>  n_1            700  
#>  n_2            350  
#>  necessary_N    1050 
#>  expected_width 0.2  
#>  delta          0.3  
#>  ratio          2    
#>  width_target   0.2  
#>  conf_level     0.95 
#> 
#> Confidence level: 95%
```
