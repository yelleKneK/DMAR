# Sensitivity Analysis for Sample Size Planning From the AIPE Perspective for an Indirect Effect

Quantifies how much misspecification of the population mediation path
coefficients \\a\\ and \\b\\ distorts an AIPE-based sample size plan for
the indirect effect \\ab\\. On each replication the function simulates a
three-variable mediation system \\X \to M \to Y\\ of size *n* with
population path coefficients `true_a` and `true_b`, fits the two
regressions of *M* on *X* and *Y* on *M* and *X*, computes the sample
indirect effect \\\hat a \hat b\\, and forms the interval the plan
targeted: the symmetric Wald interval from the delta method standard
error (`method = "closed_form"`) or the Monte Carlo interval
(`method = "monte_carlo"`), matching
[`ss_aipe_indirect_effect`](https://yelleknek.github.io/DMAR/reference/ss_aipe_indirect_effect.md).

## Usage

``` r
ss_aipe_indirect_effect_sensitivity(
  true_a = NULL,
  true_b = NULL,
  estimated_a = NULL,
  estimated_b = NULL,
  width,
  specified_N = NULL,
  method = c("closed_form", "monte_carlo"),
  conf_level = 0.95,
  B = 5000L,
  G = 1000,
  print_iter = FALSE,
  save = FALSE,
  filename = "ss_aipe_indirect_effect_sensitivity_result.csv"
)
```

## Arguments

- true_a:

  Population path coefficient *a* (from *X* to *M*); the data generating
  value.

- true_b:

  Population path coefficient *b* (from *M* to *Y* after controlling for
  *X*); the data generating value.

- estimated_a, estimated_b:

  Path coefficients used to plan the study (passed to
  [`ss_aipe_indirect_effect`](https://yelleknek.github.io/DMAR/reference/ss_aipe_indirect_effect.md)).
  Supply both or neither (if neither, supply `specified_N`).

- width:

  Desired full width of the CI on \\ab\\.

- specified_N:

  Sample size to evaluate (incompatible with `estimated_a` /
  `estimated_b`).

- method:

  One of `"closed_form"` (default) or `"monte_carlo"`; the interval
  computed on each replication, also forwarded to the planner when the
  sample size is planned from `estimated_a` and `estimated_b`. A
  planning call with `method = "monte_carlo"` runs the planner's a
  priori Monte Carlo search at its default `G`, so it takes a few
  seconds.

- conf_level:

  Confidence level (default `0.95`).

- B:

  Number of Monte Carlo draws used for the indirect-effect CI when
  `method = "monte_carlo"` (default 5000).

- G:

  Number of outer simulation replications (default 1000).

- print_iter:

  Logical.

- save:

  Logical. Save per-replication CSV.

- filename:

  Path used when `save = TRUE`.

## Value

A `data.frame` with rows for mean / median / SD of \\\hat a \hat b\\ and
the CI width, the proportion of intervals at or below `width`,
tail-specific and overall non-coverage of the population value
`true_a * true_b`, and the input echoes.

## References

Preacher, K. J., & Kelley, K. (2011). Effect size measures for mediation
models: Quantitative strategies for communicating indirect effects.
*Psychological Methods, 16*(2), 93–115.
[doi:10.1037/a0022658](https://doi.org/10.1037/a0022658)

Tofighi, D., & Kelley, K. (2020). Improved inference in mediation
analysis: Introducing the model-based constrained optimization
procedure. *Psychological Methods, 25*, 496–515.
[doi:10.1037/met0000259](https://doi.org/10.1037/met0000259)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.

## See also

[`ss_aipe_indirect_effect`](https://yelleknek.github.io/DMAR/reference/ss_aipe_indirect_effect.md),
[`var_indirect_effect`](https://yelleknek.github.io/DMAR/reference/var_indirect_effect.md)

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
# Reduced replications and a wide target width keep this fast.
set.seed(113)
ss_aipe_indirect_effect_sensitivity(
  true_a = 0.4, true_b = 0.3,
  estimated_a = 0.4, estimated_b = 0.3,
  width = 0.40, method = "closed_form",
  G = 50, print_iter = FALSE
)
#>  term               value
#>  mean_ab            0.128
#>  median_ab          0.134
#>  sd_ab              0.107
#>  mean_ci_width      0.435
#>  median_ci_width    0.434
#>  sd_ci_width        0.149
#>  pct_ci_less_w      0.42 
#>  pct_ci_miss_low    0    
#>  pct_ci_miss_high   0.1  
#>  total_type_I_error 0.1  
#>  total_N            27   
#>  true_a             0.4  
#>  true_b             0.3  
#>  true_ab            0.12 
#>  estimated_a        0.4  
#>  estimated_b        0.3  
#>  width              0.4  
#>  conf_level         0.95 
#> 
#> Confidence level: 95%
```
