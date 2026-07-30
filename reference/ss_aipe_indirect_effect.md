# Sample Size for AIPE on a Mediated (Indirect) Effect \\ab\\

Determines the sample size needed for the confidence interval on a
mediated effect \\ab\\ (the product of the \\X \to M\\ and \\M \to Y\\
coefficients in a simple three-variable mediation model) to have a
desired full width. Two CI methods are supported: `"sobel"` (asymptotic
/ delta method standard error; MacKinnon et al., 2002) and
`"monte_carlo"` (Monte Carlo confidence intervals on the distribution of
the product; Tofighi & MacKinnon, 2011; Schoemann, Boulton, & Short,
2017). The asymptotic method is essentially closed-form and very fast;
Monte Carlo is more accurate near the null but slower.

## Usage

``` r
ss_aipe_indirect_effect(
  a,
  b,
  width,
  method = c("sobel", "monte_carlo"),
  conf_level = 0.95,
  n_max = 10000L,
  B = 5000L
)
```

## Arguments

- a:

  Anticipated population coefficient for \\X \to M\\, typically on the
  standardized scale. Numeric scalar.

- b:

  Anticipated population coefficient for \\M \to Y\\ controlling for
  \\X\\, typically standardized. Numeric scalar.

- width:

  Desired full width of the confidence interval on \\ab\\.

- method:

  One of `"sobel"` (default) or `"monte_carlo"`.

- conf_level:

  Desired confidence level (default `0.95`).

- n_max:

  Upper bound on the search; default `10000`.

- B:

  Monte Carlo replications per candidate \\n\\ when
  `method = "monte_carlo"`; default `5000`. Larger values give a more
  stable estimate of CI width at the cost of speed.

## Value

A `data.frame` with rows for the recommended sample size, the expected
CI width at that size, and the inputs echoed back. The CI method is
carried on the returned object as the `ci_method` attribute.

## Details

**The mediation model.** The simple mediator model is \$\$M = \alpha_1 +
a X + \varepsilon_M,\$\$ \$\$Y = \alpha_2 + c' X + b M +
\varepsilon_Y,\$\$ with the indirect (mediated) effect of \\X\\ on \\Y\\
through \\M\\ equal to \\ab\\ (MacKinnon, Lockwood, Hoffman, West, &
Sheets, 2002).

**Sobel asymptotic standard error.** Under joint normality of the
regression estimators and *independent* \\\hat a\\ and \\\hat b\\, the
delta method standard error of \\\hat a \hat b\\ is \$\$\mathrm{SE}(\hat
a \hat b) \\=\\ \sqrt{\\a^2 \sigma_b^2 + b^2 \sigma_a^2\\},\$\$ where
for standardized predictors with no other covariates, \\\sigma_a^2 =
(1 - a^2)/(n - 2)\\ and \\\sigma_b^2 = (1 - r\_{YM \cdot X}^{2})/(n -
3)\\ with \\r\_{YM \cdot X}\\ the partial correlation of \\Y\\ and \\M\\
controlling for \\X\\. This implementation uses the simplifying
approximation \\\sigma_b^2 \approx (1 - b^2)/(n - 3)\\ when the \\c'\\
(direct effect) is unspecified, which is conservative for small \\b\\.
The CI is \\\hat a \hat b \pm z\_{1 - \alpha/2} \mathrm{SE}(\hat a \hat
b)\\.

**Monte Carlo CI.** Tofighi & MacKinnon (2011) and Schoemann, Boulton, &
Short (2017) recommend Monte Carlo CIs, which sample \\(\tilde a, \tilde
b)\\ pairs from independent normals centered at the point estimates with
the Sobel standard errors, multiply, and read off the empirical
\\(\alpha/2, 1 - \alpha/2)\\ quantiles. The resulting CI accommodates
the skewness of the product distribution and is more accurate near \\a =
0\\ or \\b = 0\\ than the symmetric Sobel interval.

**Caveats.** Both methods assume standardized variables and independent
\\\hat a\\, \\\hat b\\. For models with covariates, moderators, or
non-standard scaling, more detailed planning (e.g., via the lavaan or
PowerAnalysis frameworks) is recommended.

## References

Fritz, M. S., & MacKinnon, D. P. (2007). Required sample size to detect
the mediated effect. *Psychological Science, 18*(3), 233–239.
[doi:10.1111/j.1467-9280.2007.01882.x](https://doi.org/10.1111/j.1467-9280.2007.01882.x)

Lachowicz, M. J., Preacher, K. J., & Kelley, K. (2018). A novel measure
of effect size for mediation analysis. *Psychological Methods, 23*,
244–261. [doi:10.1037/met0000165](https://doi.org/10.1037/met0000165)

MacKinnon, D. P., Lockwood, C. M., Hoffman, J. M., West, S. G., &
Sheets, V. (2002). A comparison of methods to test mediation and other
intervening variable effects. *Psychological Methods, 7*(1), 83–104.
[doi:10.1037/1082-989X.7.1.83](https://doi.org/10.1037/1082-989X.7.1.83)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.

Preacher, K. J., & Kelley, K. (2011). Effect size measures for mediation
models: Quantitative strategies for communicating indirect effects.
*Psychological Methods, 16*(2), 93–115.
[doi:10.1037/a0022658](https://doi.org/10.1037/a0022658)

Schoemann, A. M., Boulton, A. J., & Short, S. D. (2017). Determining
power and sample size for simple and complex mediation models. *Social
Psychological and Personality Science, 8*(4), 379–386.
[doi:10.1177/1948550617715068](https://doi.org/10.1177/1948550617715068)

Tofighi, D., & Kelley, K. (2020). Improved inference in mediation
analysis: Introducing the model-based constrained optimization
procedure. *Psychological Methods, 25*, 496–515.
[doi:10.1037/met0000259](https://doi.org/10.1037/met0000259)

Tofighi, D., & MacKinnon, D. P. (2011). RMediation: An R package for
mediation analysis confidence intervals. *Behavior Research Methods,
43*(3), 692–700.
[doi:10.3758/s13428-011-0076-x](https://doi.org/10.3758/s13428-011-0076-x)

## See also

[`ss_aipe_partial_r`](https://yelleknek.github.io/DMAR/reference/ss_aipe_partial_r.md),
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
# 1. Plan n so the 95% CI on ab has full width <= 0.20, with
#        anticipated standardized a = 0.40 and b = 0.40. Sobel.
ss_aipe_indirect_effect(a = 0.40, b = 0.40, width = 0.20,
                         method = "sobel")
#>  term           value
#>  necessary_N    106  
#>  expected_width 0.2  
#>  a              0.4  
#>  b              0.4  
#>  ab             0.16 
#>  width_target   0.2  
#>  conf_level     0.95 
#> 
#> Confidence level: 95%

# 2. Monte Carlo CI (more accurate for small ab); B reduced here
#        to keep the example fast.
# \donttest{
set.seed(113)
ss_aipe_indirect_effect(a = 0.40, b = 0.40, width = 0.20,
                         method = "monte_carlo", B = 1000)
#>  term           value
#>  necessary_N    118  
#>  expected_width 0.194
#>  a              0.4  
#>  b              0.4  
#>  ab             0.16 
#>  width_target   0.2  
#>  conf_level     0.95 
#> 
#> Confidence level: 95%
# }
```
