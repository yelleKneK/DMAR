# Sample Size for AIPE on Omega Squared (ANOVA Effect Size)

Determines the sample size needed for the noncentral *F* confidence
interval on the population omega squared (\\\omega^2\\) to have a
desired width (Accuracy in Parameter Estimation; Kelley, 2008; Steiger,
2004). The function uses the same noncentral *F* machinery as
[`ci_omega_squared`](https://yelleknek.github.io/DMAR/reference/ci_omega_squared.md):
at each candidate \\N\\ it computes the expected CI width by inverting
the noncentral *F* distribution and stops at the smallest \\N\\ that
achieves the target width.

## Usage

``` r
ss_aipe_omega_squared(
  population_omega_squared,
  df_effect,
  width,
  which_width = c("Full", "Lower", "Upper"),
  conf_level = 0.95,
  assurance = NULL
)
```

## Arguments

- population_omega_squared:

  Anticipated population \\\omega^2\\. Must lie in \\\[0, 1)\\.

- df_effect:

  Numerator degrees of freedom for the effect (e.g., \\a - 1\\ for an
  \\a\\-group one-way ANOVA, or the appropriate per-effect numerator df
  in a factorial design).

- width:

  Desired full width of the CI on \\\omega^2\\.

- which_width:

  `"Full"` (default) or `"Lower"` / `"Upper"` half-width.

- conf_level:

  Desired confidence level (default `0.95`).

- assurance:

  Optional. Probability that the realized CI is no wider than `width`;
  when supplied, the sample size is inflated by the standard chi squared
  correction (Kelley, 2008).

## Value

A `data.frame` with rows for the recommended *total* sample size
`necessary_N`, the expected CI width at that sample size, and the inputs
echoed back.

## Details

**Connection to noncentral *F* machinery.** The CI on \\\omega^2\\ is
built by inverting the noncentral *F* sampling distribution of the
observed *F* statistic, following Steiger (2004) and Kelley (2007); see
[`ci_omega_squared`](https://yelleknek.github.io/DMAR/reference/ci_omega_squared.md).
To plan a sample size, we iterate: for each candidate \\N\\, compute the
*F* the analyst would observe *at the population effect size*, build its
CI on \\\omega^2\\, and stop at the smallest \\N\\ whose CI width is
below the target.

**Population-effect-to-*F* mapping.** Given a target \\\omega^2\\, the
expected sample *F* that yields exactly that \\\omega^2\\ as the point
estimate from \\\hat\omega^2 = df\_{\text{eff}}(F - 1) /
\[df\_{\text{eff}}(F - 1) + N\]\\ is \\F = 1 + \omega^2 N /
\[df\_{\text{eff}} (1 - \omega^2)\]\\. This is the *F* value used at
each iteration of the search.

**Tolerance behavior at small *N*.** For small candidate *N* the
noncentral *F* lower limit is often clamped to zero (see
[`?conf_limits_ncf`](https://yelleknek.github.io/DMAR/reference/conf_limits_ncf.md)).
The search ignores these clamps in the iteration and reports the final
clamp count, if any, as an informational message; this matches the
convention in
[`ss_aipe_R2`](https://yelleknek.github.io/DMAR/reference/ss_aipe_R2.md).

## References

Algina, J., Moulder, B. C., & Moser, B. K. (2002). Sample size
requirements for accurate estimation of squared semi-partial correlation
coefficients. *Multivariate Behavioral Research, 37*(1), 37–57.
[doi:10.1207/s15327906mbr3701_02](https://doi.org/10.1207/s15327906mbr3701_02)

Kelley, K. (2007). Confidence intervals for standardized effect sizes:
Theory, application, and implementation. *Journal of Statistical
Software, 20*(8), 1–24.
[doi:10.18637/jss.v020.i08](https://doi.org/10.18637/jss.v020.i08)

Kelley, K. (2008). Sample size planning for the squared multiple
correlation coefficient: Accuracy in parameter estimation via narrow
confidence intervals. *Multivariate Behavioral Research, 43*(4),
524–555.
[doi:10.1080/00273170802490632](https://doi.org/10.1080/00273170802490632)

Kelley, K., & Preacher, K. J. (2012). On effect size. *Psychological
Methods, 17*, 137–152.
[doi:10.1037/a0028086](https://doi.org/10.1037/a0028086)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 3 on \\\eta^2\\, Chapter 7 on factorial
designs, and Chapter 11 on generalized \\\eta^2\\ for within-subjects
designs.)

Steiger, J. H. (2004). Beyond the *F* test: Effect size confidence
intervals and tests of close fit in the analysis of variance and
contrast analysis. *Psychological Methods, 9*(2), 164–182.
[doi:10.1037/1082-989X.9.2.164](https://doi.org/10.1037/1082-989X.9.2.164)

## See also

[`ci_omega_squared`](https://yelleknek.github.io/DMAR/reference/ci_omega_squared.md),
[`omega_squared`](https://yelleknek.github.io/DMAR/reference/omega_squared.md),
[`omega_squared_partial`](https://yelleknek.github.io/DMAR/reference/omega_squared_partial.md),
[`ss_aipe_R2`](https://yelleknek.github.io/DMAR/reference/ss_aipe_R2.md)

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
# 1. Plan total N so the 95% CI on omega^2 has full width <= 0.10
#        in a 3-group one-way ANOVA (df_effect = 2), anticipated
#        omega^2 = 0.10.
ss_aipe_omega_squared(population_omega_squared = 0.10,
                      df_effect = 2,
                      width = 0.10)
#> During the iterative sample size search, the noncentral F lower-limit clamp in conf_limits_ncf() fired in 10 intermediate evaluations.
#>  term                     value 
#>  necessary_N              473   
#>  expected_width           0.0999
#>  population_omega_squared 0.1   
#>  df_effect                2     
#>  width_target             0.1   
#>  conf_level               0.95  
#> 
#> Confidence level: 95%

# 2. Same problem with 80% assurance:
ss_aipe_omega_squared(population_omega_squared = 0.10,
                      df_effect = 2,
                      width = 0.10,
                      assurance = 0.80)
#> During the iterative sample size search, the noncentral F lower-limit clamp in conf_limits_ncf() fired in 10 intermediate evaluations.
#>  term                     value 
#>  necessary_N              499   
#>  expected_width           0.0973
#>  population_omega_squared 0.1   
#>  df_effect                2     
#>  width_target             0.1   
#>  conf_level               0.95  
#> 
#> Confidence level: 95%
```
