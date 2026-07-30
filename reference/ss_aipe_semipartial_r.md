# Sample Size for AIPE on a Semipartial (Part) Correlation

Determines the sample size needed for a confidence interval on a
population semipartial correlation \\r\_{Y(X \cdot Z_1 \cdots Z_J)}\\
(the unique contribution of \\X\\ to \\Y\\ after controlling for \\Z_1,
\ldots, Z_J\\, with \\Y\\ *not* residualized) to have a desired width,
using the Olkin-Finn (1995) / Algina-Olejnik (2003) asymptotic variance
and the AIPE framework of Kelley & Maxwell (2003).

## Usage

``` r
ss_aipe_semipartial_r(
  r_sp,
  J,
  width,
  which_width = c("Full", "Lower", "Upper"),
  conf_level = 0.95,
  assurance = NULL
)
```

## Arguments

- r_sp:

  Anticipated population semipartial correlation, in \\(-1, 1)\\.

- J:

  Number of variables partialled out of \\X\\ (count of \\Z_1, \ldots,
  Z_J\\); must be at least 1.

- width:

  Desired full width of the confidence interval on the semipartial
  correlation.

- which_width:

  Whether `width` refers to the `"Full"` width (default) or the
  `"Lower"`/`"Upper"` half-width.

- conf_level:

  Desired confidence level (default `0.95`).

- assurance:

  Optional. Probability that the realized CI is no wider than `width`.
  When supplied, the sample size is inflated using the standard chi
  squared correction (Kelley & Maxwell, 2003).

## Value

A `data.frame` with rows for the recommended sample size, the expected
CI width at that sample size, and the inputs echoed back.

## Details

**Asymptotic variance of the semipartial.** Under multivariate
normality, the sample semipartial correlation \\r\_{Y(X \cdot Z)}\\ has
asymptotic variance \$\$\mathrm{Var}(\hat r\_{Y(X \cdot Z)}) \\\approx\\
\frac{(1 - r\_{Y(X \cdot Z)}^2)^2}{n - J - 1}\$\$ (Olkin & Finn, 1995,
with the partial-correlation degrees-of-freedom correction). Inverting
for the sample size needed to achieve a target half-width \\w\_{1/2}\\
at confidence level \\1 - \alpha\\: \$\$n \\=\\ J + 1 + \Big\lceil
z\_{1 - \alpha/2}^{2} \cdot (1 - r\_{Y(X \cdot Z)}^{2})^2 / w\_{1/2}^{2}
\Big\rceil.\$\$

**Comparison with partial-r planning.** The partial correlation \\r\_{XY
\cdot Z}\\ divides the covariance after residualizing both \\X\\ and
\\Y\\ on \\Z\\; the semipartial divides after residualizing only \\X\\.
The semipartial is the natural effect size companion to a standardized
regression coefficient: its square equals the \\\Delta R^2\\ contributed
by \\X\\ above and beyond the controls. See
[`var_semipartial_r`](https://yelleknek.github.io/DMAR/reference/var_semipartial_r.md)
for the asymptotic variance, and
[`ss_aipe_partial_r`](https://yelleknek.github.io/DMAR/reference/ss_aipe_partial_r.md)
for the partial-correlation analog of this function.

**Note on conservatism of the assurance plan.** The empirical simulation
reported in
[`vignette("aipe_simulation_study", package = "DMAR")`](https://yelleknek.github.io/DMAR/articles/aipe_simulation_study.md)
finds that `ss_aipe_semipartial_r()` is on the boundary of its valid
range at 80% assurance and modestly conservative at 99% assurance. At
\\\gamma = 0.80\\, the realized assurance at the recommended sample size
is within Monte-Carlo error of the target, that is, the bound is
operating at the edge of its validity. At \\\gamma = 0.99\\, the ideal
sample size is about 15 to 20 subjects smaller than the recommended
sample size, reflecting the looser upper-tail bound at the 99% level.
The recommended sample size is therefore a sufficient sample size rather
than the smallest possible sample size. A small safety margin (5 to 10
subjects) is advisable when planning at \\\gamma = 0.80\\. See the
simulation vignette for the per-condition overshoot.

## References

Algina, J., & Olejnik, S. (2003). Sample size tables for correlation
analysis with applications in partial correlation and multiple
regression analysis. *Multivariate Behavioral Research, 38*(3), 309–323.
[doi:10.1207/s15327906mbr3803_02](https://doi.org/10.1207/s15327906mbr3803_02)

Cohen, J., Cohen, P., West, S. G., & Aiken, L. S. (2003). *Applied
multiple regression/correlation analysis for the behavioral sciences*
(3rd ed.). Lawrence Erlbaum.

Kelley, K., & Maxwell, S. E. (2003). Sample size for multiple
regression: Obtaining regression coefficients that are accurate, not
simply significant. *Psychological Methods, 8*(3), 305–321.
[doi:10.1037/1082-989X.8.3.305](https://doi.org/10.1037/1082-989X.8.3.305)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 3 on the one-way ANOVA and Chapter 4 on
contrasts.)

Olkin, I., & Finn, J. D. (1995). Correlations redux. *Psychological
Bulletin, 118*(1), 155–164.
[doi:10.1037/0033-2909.118.1.155](https://doi.org/10.1037/0033-2909.118.1.155)

## See also

[`var_semipartial_r`](https://yelleknek.github.io/DMAR/reference/var_semipartial_r.md),
[`ss_aipe_partial_r`](https://yelleknek.github.io/DMAR/reference/ss_aipe_partial_r.md),
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
[`ss_aipe_omega_squared()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_omega_squared.md),
[`ss_aipe_omega_squared_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_omega_squared_sensitivity.md),
[`ss_aipe_partial_r()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_partial_r.md),
[`ss_aipe_partial_r_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_partial_r_sensitivity.md),
[`ss_aipe_pcm_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_pcm_sensitivity.md),
[`ss_aipe_reliability_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reliability_sensitivity.md),
[`ss_aipe_semipartial_r_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_semipartial_r_sensitivity.md),
[`ss_aipe_tost_smd_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_tost_smd_sensitivity.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# 1. Plan n so the 95% CI on r_sp (J = 3) has full width <= 0.15
#        when the anticipated semipartial is 0.25.
ss_aipe_semipartial_r(r_sp = 0.25, J = 3, width = 0.15)
#>  term           value
#>  necessary_N    605  
#>  expected_width 0.15 
#>  r_sp           0.25 
#>  J              3    
#>  width_target   0.15 
#>  conf_level     0.95 
#> 
#> Confidence level: 95%

# 2. With 80% assurance:
ss_aipe_semipartial_r(r_sp = 0.25, J = 3, width = 0.15,
                      assurance = 0.80)
#>  term           value
#>  necessary_N    635  
#>  expected_width 0.146
#>  r_sp           0.25 
#>  J              3    
#>  width_target   0.15 
#>  conf_level     0.95 
#> 
#> Confidence level: 95%
```
