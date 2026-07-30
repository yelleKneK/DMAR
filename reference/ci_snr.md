# Confidence Interval for the Signal-to-Noise Ratio

Computes the exact confidence interval for the signal-to-noise ratio in
a fixed effects analysis of variance, the variance due to the factor of
interest divided by the error variance, expressing the magnitude of an
effect relative to the unexplained variability.

## Usage

``` r
ci_snr(
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

  Observed *F*-value from the analysis of variance

- df_1:

  Numerator degrees of freedom

- df_2:

  Denominator degrees of freedom

- N:

  Sample size

- conf_level:

  Confidence interval coverage (i.e., 1 - Type I error rate), default is
  .95

- alpha_lower:

  Type I error for the lower confidence limit

- alpha_upper:

  Type I error for the upper confidence limit

- ...:

  Allows one to potentially include parameter values for inner functions

## Value

A 2-row `data.frame` with columns `term` and `value`. The `term` values
are `"lower_limit"` and `"upper_limit"`, giving the lower and upper
confidence limits on the signal-to-noise ratio.

## Details

The confidence level must be specified in one of following two ways:
using confidence interval coverage (`conf_level`), or lower and upper
confidence limits (`alpha_lower` and `alpha_upper`).

This function uses the confidence interval transformation principle
(Steiger, 2004) to transform the confidence limits for the noncentrality
parameter to the confidence limits for the population's signal-to-noise
ratio. The confidence interval for noncentral *F* parameter can be
obtained from the `conf_limits_ncf` function in DMAR, which is used
internally within this function.

## Note

The signal to noise ratio is defined as the variance due to the
particular factor over the error variance (i.e., the mean square error).

## References

Kelley, K. (2007). Confidence intervals for standardized effect sizes:
Theory, application, and implementation. *Journal of Statistical
Software, 20*(8), 1–24.
[doi:10.18637/jss.v020.i08](https://doi.org/10.18637/jss.v020.i08)

Fleishman, A. I. (1980). Confidence intervals for correlation ratios.
*Educational and Psychological Measurement, 40*(3), 659–670.

Steiger, J. H. (2004). Beyond the *F* Test: Effect size confidence
intervals and tests of close fit in the Analysis of Variance and
Contrast Analysis. *Psychological Methods, 9*(2), 164–182.
[doi:10.1037/1082-989X.9.2.164](https://doi.org/10.1037/1082-989X.9.2.164)

## See also

[`ci_srsnr`](https://yelleknek.github.io/DMAR/reference/ci_srsnr.md),
[`conf_limits_ncf`](https://yelleknek.github.io/DMAR/reference/conf_limits_ncf.md)

Other confidence intervals for effect sizes:
[`ci_R2()`](https://yelleknek.github.io/DMAR/reference/ci_R2.md),
[`ci_c()`](https://yelleknek.github.io/DMAR/reference/ci_c.md),
[`ci_c_ancova()`](https://yelleknek.github.io/DMAR/reference/ci_c_ancova.md),
[`ci_c_ancova_bp()`](https://yelleknek.github.io/DMAR/reference/ci_c_ancova_bp.md),
[`ci_cc()`](https://yelleknek.github.io/DMAR/reference/ci_cc.md),
[`ci_cv()`](https://yelleknek.github.io/DMAR/reference/ci_cv.md),
[`ci_eta_squared()`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared.md),
[`ci_eta_squared_generalized()`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared_generalized.md),
[`ci_eta_squared_partial()`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared_partial.md),
[`ci_mahalanobis()`](https://yelleknek.github.io/DMAR/reference/ci_mahalanobis.md),
[`ci_omega_squared()`](https://yelleknek.github.io/DMAR/reference/ci_omega_squared.md),
[`ci_pvaf()`](https://yelleknek.github.io/DMAR/reference/ci_pvaf.md),
[`ci_r()`](https://yelleknek.github.io/DMAR/reference/ci_r.md),
[`ci_rc()`](https://yelleknek.github.io/DMAR/reference/ci_rc.md),
[`ci_reg_coef()`](https://yelleknek.github.io/DMAR/reference/ci_reg_coef.md),
[`ci_rmsea()`](https://yelleknek.github.io/DMAR/reference/ci_rmsea.md),
[`ci_sc()`](https://yelleknek.github.io/DMAR/reference/ci_sc.md),
[`ci_sc_ancova()`](https://yelleknek.github.io/DMAR/reference/ci_sc_ancova.md),
[`ci_sm()`](https://yelleknek.github.io/DMAR/reference/ci_sm.md),
[`ci_smd()`](https://yelleknek.github.io/DMAR/reference/ci_smd.md),
[`ci_smd_c()`](https://yelleknek.github.io/DMAR/reference/ci_smd_c.md),
[`ci_src()`](https://yelleknek.github.io/DMAR/reference/ci_src.md),
[`ci_srsnr()`](https://yelleknek.github.io/DMAR/reference/ci_srsnr.md),
[`contrast_adjusted()`](https://yelleknek.github.io/DMAR/reference/contrast_adjusted.md),
[`plot_smd()`](https://yelleknek.github.io/DMAR/reference/plot_smd.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
## Bargman (1970) gave an example in which a 5-group ANOVA with 11 subjects in each
## group is conducted and the observed F value is 11.2213. This example was
## used in Venables (1975),  Fleishman (1980), and Steiger (2004). If one wants to calculate
## the exact confidence interval for the signal-to-noise ratio of that example, this
## function can be used.

ci_snr(F_value = 11.221, df_1 = 4, df_2 = 50, N = 55)
#>  term        value
#>  lower_limit 0.293
#>  upper_limit 1.42 
#> 
#> Confidence level: 95%

ci_snr(F_value = 11.221, df_1 = 4, df_2 = 50, N = 55, conf_level = .90)
#>  term        value
#>  lower_limit 0.352
#>  upper_limit 1.3  
#> 
#> Confidence level: 90%

ci_snr(F_value = 11.221, df_1 = 4, df_2 = 50, N = 55, alpha_lower = .02, alpha_upper = .03)
#>  term        value
#>  lower_limit 0.276
#>  upper_limit 1.39 
#> 
#> Confidence level: 95%
```
