# Confidence Interval for the Proportion of Variance Accounted for (in the Dependent Variable by Knowing the Levels of the Factor)

Computes the exact confidence limits for the proportion of variance in
the dependent variable accounted for by knowing the levels of the factor
(group status in a single factor design) in a fixed effects analysis of
variance, so an omnibus *F*-test is accompanied by an effect size with a
statement of its precision.

## Usage

``` r
ci_pvaf(
  F_value = NULL,
  df_1 = NULL,
  df_2 = NULL,
  N = NULL,
  conf_level = 0.95,
  alpha_lower = NULL,
  alpha_upper = NULL,
  ...
)
```

## Arguments

- F_value:

  Observed *F*-value from fixed effects analysis of variance

- df_1:

  Numerator degrees of freedom

- df_2:

  Denominator degrees of freedom

- N:

  Sample size

- conf_level:

  Confidence interval coverage (i.e., 1-Type I error rate); default is
  .95

- alpha_lower:

  Type I error for the lower confidence limit

- alpha_upper:

  Type I error for the upper confidence limit

- ...:

  Allows one to potentially include parameter values for inner functions

## Value

A 4-row `data.frame` with columns `term`, `value`, `prob_less`, and
`prob_greater`. The `term` values are `"lower_limit"` (the lower
confidence limit on the proportion of variance accounted for, on the
\[0, 1\] scale), `"pvaf"` (the sample proportion of variance accounted
for, `df_1 * F_value / (df_1 * F_value + df_2)`, the same value that eta
squared reports, so the point estimate sits between its confidence
limits), `"upper_limit"` (the upper confidence limit), and
`"actual_coverage"` (the achieved coverage probability, which equals
`conf_level` when both tail targets are met). The `prob_less` and
`prob_greater` columns report the achieved tail-error probabilities at
the two limits; `NA` on the `"pvaf"` and `"actual_coverage"` rows.

## Details

The confidence level must be specified in one of following two ways:
using confidence interval coverage (`conf_level`), or lower and upper
confidence limits (`alpha_lower` and `alpha_upper`).

This function uses the confidence interval transformation principle
(Steiger, 2004) to transform the confidence limits for the noncentrality
parameter to the confidence limits for the population proportion of
variance accounted for by knowing the group status. The confidence
interval for the noncentral *F* parameter can be obtained from the
function
[`conf_limits_ncf`](https://yelleknek.github.io/DMAR/reference/conf_limits_ncf.md),
which is used within this function.

## Note

This function can be used for single or factorial ANOVA designs.

## References

Fleishman, A. I. (1980). Confidence intervals for correlation ratios.
*Educational and Psychological Measurement, 40*(3), 659–670.

Kelley, K. (2007). Confidence intervals for standardized effect sizes:
Theory, application, and implementation. *Journal of Statistical
Software, 20*(8), 1–24.
[doi:10.18637/jss.v020.i08](https://doi.org/10.18637/jss.v020.i08)

Kelley, K. (2008). Sample size planning for the squared multiple
correlation coefficient: Accuracy in parameter estimation via narrow
confidence intervals. *Multivariate Behavioral Research, 43*, 524–555.
[doi:10.1080/00273170802490632](https://doi.org/10.1080/00273170802490632)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 3 on \\R^2\\ as a model comparison effect
size.)

Steiger, J. H. (2004). Beyond the *F* Test: Effect size confidence
intervals and tests of close fit in the Analysis of Variance and
Contrast Analysis. *Psychological Methods, 9*(2), 164–182.
[doi:10.1037/1082-989X.9.2.164](https://doi.org/10.1037/1082-989X.9.2.164)

## See also

[`conf_limits_ncf`](https://yelleknek.github.io/DMAR/reference/conf_limits_ncf.md)

Other confidence intervals for effect sizes:
[`ci_R2()`](https://yelleknek.github.io/DMAR/reference/ci_R2.md),
[`ci_c()`](https://yelleknek.github.io/DMAR/reference/ci_c.md),
[`ci_c_ancova()`](https://yelleknek.github.io/DMAR/reference/ci_c_ancova.md),
[`ci_c_ancova_bp()`](https://yelleknek.github.io/DMAR/reference/ci_c_ancova_bp.md),
[`ci_correlation`](https://yelleknek.github.io/DMAR/reference/ci_correlation.md),
[`ci_cv()`](https://yelleknek.github.io/DMAR/reference/ci_cv.md),
[`ci_eta_squared()`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared.md),
[`ci_eta_squared_generalized()`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared_generalized.md),
[`ci_eta_squared_partial()`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared_partial.md),
[`ci_mahalanobis()`](https://yelleknek.github.io/DMAR/reference/ci_mahalanobis.md),
[`ci_omega_squared()`](https://yelleknek.github.io/DMAR/reference/ci_omega_squared.md),
[`ci_rc()`](https://yelleknek.github.io/DMAR/reference/ci_rc.md),
[`ci_reg_coef()`](https://yelleknek.github.io/DMAR/reference/ci_reg_coef.md),
[`ci_rmsea()`](https://yelleknek.github.io/DMAR/reference/ci_rmsea.md),
[`ci_sc()`](https://yelleknek.github.io/DMAR/reference/ci_sc.md),
[`ci_sc_ancova()`](https://yelleknek.github.io/DMAR/reference/ci_sc_ancova.md),
[`ci_sm()`](https://yelleknek.github.io/DMAR/reference/ci_sm.md),
[`ci_smd()`](https://yelleknek.github.io/DMAR/reference/ci_smd.md),
[`ci_smd_c()`](https://yelleknek.github.io/DMAR/reference/ci_smd_c.md),
[`ci_snr()`](https://yelleknek.github.io/DMAR/reference/ci_snr.md),
[`ci_src()`](https://yelleknek.github.io/DMAR/reference/ci_src.md),
[`ci_srsnr()`](https://yelleknek.github.io/DMAR/reference/ci_srsnr.md),
[`contrast_adjusted()`](https://yelleknek.github.io/DMAR/reference/contrast_adjusted.md),
[`plot_smd()`](https://yelleknek.github.io/DMAR/reference/plot_smd.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
## Bargman (1970) gave an example in which a 5-group ANOVA with 11 subjects in each
## group is conducted and the observed F value is 11.221. This example was used
## in Venables (1975),  Fleishman (1980), and Steiger (2004). If one wants to calculate the
## exact confidence interval for the proportion of variance accounted for in that example,
## this function can be used.
ci_pvaf(F_value = 11.221, df_1 = 4, df_2 = 50, N = 55)
#>  term            value prob_less prob_greater
#>  lower_limit     0.226 0.025     0.975       
#>  pvaf            0.473 <NA>      <NA>        
#>  upper_limit     0.587 0.975     0.025       
#>  actual_coverage 0.95  <NA>      <NA>        
#> 
#> Confidence level: 95%

ci_pvaf(F_value = 11.221, df_1 = 4, df_2 = 50, N = 55, conf_level = .90)
#>  term            value prob_less prob_greater
#>  lower_limit     0.261 0.05      0.95        
#>  pvaf            0.473 <NA>      <NA>        
#>  upper_limit     0.565 0.95      0.05        
#>  actual_coverage 0.9   <NA>      <NA>        
#> 
#> Confidence level: 90%

ci_pvaf(F_value = 11.221, df_1 = 4, df_2 = 50, N = 55, alpha_lower = 0, alpha_upper = .05)
#>  term            value prob_less prob_greater
#>  lower_limit     0     0         1           
#>  pvaf            0.473 <NA>      <NA>        
#>  upper_limit     0.565 0.95      0.05        
#>  actual_coverage 0.95  <NA>      <NA>        
#> 
#> Confidence level: 95%
```
