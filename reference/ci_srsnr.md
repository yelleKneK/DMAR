# Confidence Interval for the Square Root of the Signal-to-Noise Ratio

Computes the exact confidence interval for the square root of the
signal-to-noise ratio, the standard deviation of the group means
relative to the error standard deviation. On this root scale the
quantity is an effect size in standard deviation units, the multi-group
analog of the standardized mean difference.

## Usage

``` r
ci_srsnr(
  F_value = NULL,
  df_1 = NULL,
  df_2 = NULL,
  N = NULL,
  means = NULL,
  sigma_squared = NULL,
  n_per_group = NULL,
  conf_level = 0.95,
  alpha_lower = NULL,
  alpha_upper = NULL,
  ...
)
```

## Arguments

- F_value:

  Observed *F*-value from the analysis of variance. Use this argument
  when re-analyzing existing data.

- df_1:

  Numerator degrees of freedom

- df_2:

  Denominator degrees of freedom

- N:

  Sample size

- means:

  Numeric vector of population or hypothesized group means. Supply
  together with `sigma_squared` and `n_per_group` as a design-stage
  alternative to `F_value`: the function then computes the *F*-value
  implied by these design parameters and proceeds with the same
  noncentral F machinery.

- sigma_squared:

  The within-group variance. Used with `means`.

- n_per_group:

  A single per-group sample size, or a vector of per-group sample sizes
  the same length as `means`. Used with `means`.

- conf_level:

  Confidence interval coverage (i.e., 1 - Type I error rate); default is
  .95

- alpha_lower:

  Type I error for the lower confidence limit

- alpha_upper:

  Type I error for the upper confidence limit

- ...:

  Allows one to potentially include parameter values for inner functions

## Value

A 2-row `data.frame` with columns `term` and `value`. The `term` values
are `"lower_limit"` and `"upper_limit"`, giving the square roots of the
corresponding signal-to-noise-ratio confidence limits.

## Details

The confidence level must be specified in one of following two ways:
using confidence interval coverage (`conf_level`), or lower and upper
confidence limits (`alpha_lower` and `alpha_upper`).

The square root of the signal-to-noise ratio is defined as the standard
deviation due to the particular factor over the standard deviation of
the error (i.e., the square root of the mean square error). This
function uses the confidence interval transformation principle (Steiger,
2004) to transform the confidence limits for the noncentrality parameter
to the confidence limits for square root of signal-to-noise ratio. The
confidence interval for noncentral *F* parameter can be obtained from
function `conf_limits_ncf` in DMAR.

## References

Fleishman, A. I. (1980). Confidence intervals for correlation ratios.
*Educational and Psychological Measurement, 40*(3), 659–670.

Kelley, K. (2007). Confidence intervals for standardized effect sizes:
Theory, application, and implementation. *Journal of Statistical
Software, 20*(8), 1–24.
[doi:10.18637/jss.v020.i08](https://doi.org/10.18637/jss.v020.i08)

Steiger, J. H. (2004). Beyond the *F* Test: Effect size confidence
intervals and tests of close fit in the Analysis of Variance and
Contrast Analysis. *Psychological Methods, 9*(2), 164–182.
[doi:10.1037/1082-989X.9.2.164](https://doi.org/10.1037/1082-989X.9.2.164)

## See also

[`ci_snr`](https://yelleknek.github.io/DMAR/reference/ci_snr.md),
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
[`ci_snr()`](https://yelleknek.github.io/DMAR/reference/ci_snr.md),
[`ci_src()`](https://yelleknek.github.io/DMAR/reference/ci_src.md),
[`contrast_adjusted()`](https://yelleknek.github.io/DMAR/reference/contrast_adjusted.md),
[`plot_smd()`](https://yelleknek.github.io/DMAR/reference/plot_smd.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
## To illustrate the calculation of the confidence interval for noncentral
## F parameter,Bargman (1970) gave an example in which a 5-group ANOVA with
## 11 subjects in each group is conducted and the observed F value is 11.2213.
## This example continued to be used in Venables (1975),  Fleishman (1980),
## and Steiger (2004). If one wants to calculate the exact confidence interval
## for square root of the signal-to-noise ratio of that example, this
## function can be used.

ci_srsnr(F_value = 11.221, df_1 = 4, df_2 = 50, N = 55)
#>  term        value
#>  lower_limit 0.541
#>  upper_limit 1.19 
#> 
#> Confidence level: 95%

ci_srsnr(F_value = 11.221, df_1 = 4, df_2 = 50, N = 55, conf_level = .90)
#>  term        value
#>  lower_limit 0.594
#>  upper_limit 1.14 
#> 
#> Confidence level: 90%

ci_srsnr(F_value = 11.221, df_1 = 4, df_2 = 50, N = 55, alpha_lower = .02, alpha_upper = .03)
#>  term        value
#>  lower_limit 0.525
#>  upper_limit 1.18 
#> 
#> Confidence level: 95%

# Design-stage call with population means + within-group variance + n.
# Useful when planning a study, before data are observed: derives the
# implied F-value internally and returns the resulting CI on the square
# root of the signal-to-noise ratio.
ci_srsnr(means = c(94, 91, 92, 83), sigma_squared = 67.375, n_per_group = 6)
#> Warning: The observed F_value is below the alpha_lower critical value of the central F-distribution; the lower noncentrality limit has been clamped to 0 and the reported 'prob_greater' on the lower_limit row reflects the actual upper-tail probability at lambda = 0.
#>  term        value
#>  lower_limit 0    
#>  upper_limit 0.874
#> 
#> Confidence level: 95%
```
