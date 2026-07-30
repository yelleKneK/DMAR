# Sample Size for AIPE on a Partial Correlation

Determines the sample size needed for a confidence interval on a
population partial correlation \\\rho\_{XY \cdot Z_1 \cdots Z_J}\\ to
have a desired width (Accuracy in Parameter Estimation; Kelley, 2008).
The function inverts the asymptotic variance of the partial Pearson
correlation, either on the raw scale (Olkin & Finn, 1995) or on the
Fisher *z*-transformed scale (Hotelling, 1953; Bonett, 2008), and solves
for the smallest \\n\\ that achieves the target half-width (or full
width).

## Usage

``` r
ss_aipe_partial_r(
  rho,
  J,
  width,
  which_width = c("Full", "Lower", "Upper"),
  conf_level = 0.95,
  fisher_z = FALSE,
  assurance = NULL
)
```

## Arguments

- rho:

  Anticipated population partial correlation, in \\(-1, 1)\\.

- J:

  Number of variables partialled out (count of \\Z_1, \ldots, Z_J\\);
  must be at least 1.

- width:

  Desired full width of the confidence interval on the partial
  correlation.

- which_width:

  Whether `width` refers to the `"Full"` width (default) or the
  `"Lower"` or `"Upper"` half-width.

- conf_level:

  Desired confidence level (default `0.95`).

- fisher_z:

  Logical. If `TRUE`, the half-width target is applied on the Fisher-*z*
  scale (variance \\1/(n - J - 3)\\; Bonett, 2008) and the resulting CI
  is back-transformed via \\\tanh\\. If `FALSE` (default), uses the
  raw-scale Olkin-Finn (1995) asymptotic variance.

- assurance:

  Optional. Probability that the realized CI is no wider than `width`
  (\\1 - \gamma\\). When supplied, the sample size is inflated using the
  standard chi squared correction (Kelley, 2008); when `NULL`, the
  assurance is fixed at 0.5.

## Value

A `data.frame` with the rows `necessary_N` (the recommended total sample
size, rounded up), `expected_width` at that sample size, and the inputs
echoed back.

## Details

**Raw-scale Olkin-Finn asymptotic variance.** The half-width of a
\\100(1 - \alpha)\\\\ CI on the partial Pearson correlation is
approximately \$\$w\_{1/2} \\\approx\\ z\_{1 - \alpha/2} \cdot
\sqrt{\\\frac{(1 - \rho\_{XY \cdot Z}^{\\2})^2}{n - J - 1}\\}.\$\$
Solving for \\n\\: \$\$n \\=\\ J + 1 + \Big\lceil (z\_{1 - \alpha/2})^2
\cdot (1 - \rho\_{XY \cdot Z}^{\\2})^2 / w\_{1/2}^{2} \Big\rceil.\$\$
This is the planning analog of the half-width of
[`ci_r`](https://yelleknek.github.io/DMAR/reference/ci_r.md) applied to
a partial correlation.

**Fisher-*z* scale (recommended for small \\\rho\\, near boundary, or
small \\n - J\\).** Bonett (2008) advocates planning on the
variance-stabilized Fisher-*z* scale and back-transforming the bounds.
On the \\z\\ scale, the asymptotic half-width is \$\$w^{(z)}\_{1/2}
\\\approx\\ z\_{1 - \alpha/2} / \sqrt{n - J - 3}.\$\$ Solving for the
\\n\\ that achieves a given back-transformed \\w\_{1/2}\\ is done by a
1-D search; this is generally the more accurate route when \\n\\ is
small or \\\|\rho\|\\ is large.

**When to use partial vs. simple correlation planning.** Use this
function when the inferential target is the population correlation
between \\X\\ and \\Y\\ *after* statistically controlling for \\Z_1,
\ldots, Z_J\\. For the simple Pearson correlation, see
[`ss_aipe_rc`](https://yelleknek.github.io/DMAR/reference/ss_aipe_rc.md)
or equivalent.

**Note on conservatism of the assurance plan.** The empirical simulation
reported in
[`vignette("aipe_simulation_study", package = "DMAR")`](https://yelleknek.github.io/DMAR/articles/aipe_simulation_study.md)
finds that `ss_aipe_partial_r()` is tight (zero overshoot) at 80%
assurance but modestly conservative at 99% assurance, with an empirical
ideal sample size of about 5 to 10 subjects smaller than the recommended
sample size. The mechanism is the usual one for AIPE assurance plans:
the Olkin-Finn (1995) Wald- style upper bound on \\\Pr(\widehat W \>
\omega)\\ that the planner inverts is not tight at the recommended
sample size, especially at the 99% level where the inversion has to push
further into the upper tail of \\\widehat W\\. The recommended sample
size is a sufficient sample size rather than the smallest possible
sample size. See the simulation vignette for the per-condition
overshoot.

## References

Algina, J., & Olejnik, S. (2003). Sample size tables for correlation
analysis with applications in partial correlation and multiple
regression analysis. *Multivariate Behavioral Research, 38*(3), 309–323.
[doi:10.1207/s15327906mbr3803_02](https://doi.org/10.1207/s15327906mbr3803_02)

Bonett, D. G. (2008). Confidence intervals for standardized linear
contrasts of means. *Psychological Methods, 13*(2), 99–109.
[doi:10.1037/1082-989X.13.2.99](https://doi.org/10.1037/1082-989X.13.2.99)

Hotelling, H. (1953). New light on the correlation coefficient and its
transforms. *Journal of the Royal Statistical Society, Series B, 15*(2),
193–232.

Kelley, K. (2008). Sample size planning for the squared multiple
correlation coefficient: Accuracy in parameter estimation via narrow
confidence intervals. *Multivariate Behavioral Research, 43*(4),
524–555.
[doi:10.1080/00273170802490632](https://doi.org/10.1080/00273170802490632)

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

[`var_partial_r`](https://yelleknek.github.io/DMAR/reference/var_partial_r.md),
[`expected_partial_r`](https://yelleknek.github.io/DMAR/reference/expected_partial_r.md),
[`ss_aipe_semipartial_r`](https://yelleknek.github.io/DMAR/reference/ss_aipe_semipartial_r.md),
[`ss_aipe_rc`](https://yelleknek.github.io/DMAR/reference/ss_aipe_rc.md)

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
# 1. Plan n so the 95% CI on rho_XY.Z (J = 2 controls) has
#        full width <= 0.20, when the anticipated partial r is 0.30.
ss_aipe_partial_r(rho = 0.30, J = 2, width = 0.20)
#>  term           value
#>  necessary_N    322  
#>  expected_width 0.2  
#>  rho            0.3  
#>  J              2    
#>  width_target   0.2  
#>  conf_level     0.95 
#> 
#> Confidence level: 95%

# 2. Same problem on the Fisher-z scale (Bonett 2008):
ss_aipe_partial_r(rho = 0.30, J = 2, width = 0.20, fisher_z = TRUE)
#>  term           value
#>  necessary_N    322  
#>  expected_width 0.2  
#>  rho            0.3  
#>  J              2    
#>  width_target   0.2  
#>  conf_level     0.95 
#> 
#> Confidence level: 95%

# 3. With 80% assurance (Kelley 2008):
ss_aipe_partial_r(rho = 0.30, J = 2, width = 0.20, assurance = 0.80)
#>  term           value
#>  necessary_N    344  
#>  expected_width 0.193
#>  rho            0.3  
#>  J              2    
#>  width_target   0.2  
#>  conf_level     0.95 
#> 
#> Confidence level: 95%
```
