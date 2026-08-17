# Confidence Interval for the Standardized Mean

Computes the exact confidence interval for the standardized mean, the
mean divided by the standard deviation, by inverting the noncentral *t*
distribution. The standardized mean is the one-sample analog of the
standardized mean difference and shares its noncentral interval theory.

## Usage

``` r
ci_sm(
  sm = NULL,
  mean = NULL,
  sd = NULL,
  ncp = NULL,
  N = NULL,
  conf_level = 0.95,
  alpha_lower = NULL,
  alpha_upper = NULL,
  ...
)
```

## Arguments

- sm:

  Standardized mean

- mean:

  Mean

- sd:

  Standard deviation

- ncp:

  Noncentral parameter

- N:

  Sample size

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

A 3-row `data.frame` with columns `term` and `value`. The `term` values
are `"lower_limit"` (the lower confidence limit on the standardized
mean), `"std_mean"` (the standardized mean), and `"upper_limit"` (the
upper confidence limit).

## Details

The user must specify the standardized mean in one and only one of the
three ways: a) mean and standard deviation (`mean` and `sd`), b)
standardized mean (`sm`), and c) noncentral parameter (`ncp`). The
confidence level must be specified in one of following two ways: using
confidence interval coverage (`conf_level`), or lower and upper
confidence limits (`alpha_lower` and `alpha_upper`). This function uses
the exact confidence interval method based on noncentral
*t*-distributions. The confidence interval for noncentral *t*-parameter
can be obtained from the `conf_limits_nct` function in DMAR.

## Note

The standardized mean is the mean divided by the standard deviation.

## References

Kelley, K. (2007). Confidence intervals for standardized effect sizes:
Theory, application, and implementation. *Journal of Statistical
Software, 20*(8), 1–24.
[doi:10.18637/jss.v020.i08](https://doi.org/10.18637/jss.v020.i08)

Steiger, J. H., & Fouladi, R. T. (1997). Noncentrality interval
estimation and the evaluation of statistical methods. In L. L. Harlow,
S. A. Mulaik, & J. H. Steiger (Eds.), *What if there were no
significance tests?* (pp. 221–257). Mahwah, NJ: Lawrence Erlbaum.

## See also

[`conf_limits_nct`](https://yelleknek.github.io/DMAR/reference/conf_limits_nct.md)

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
[`ci_pvaf()`](https://yelleknek.github.io/DMAR/reference/ci_pvaf.md),
[`ci_rc()`](https://yelleknek.github.io/DMAR/reference/ci_rc.md),
[`ci_reg_coef()`](https://yelleknek.github.io/DMAR/reference/ci_reg_coef.md),
[`ci_rmsea()`](https://yelleknek.github.io/DMAR/reference/ci_rmsea.md),
[`ci_sc()`](https://yelleknek.github.io/DMAR/reference/ci_sc.md),
[`ci_sc_ancova()`](https://yelleknek.github.io/DMAR/reference/ci_sc_ancova.md),
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
ci_sm(sm = 2.037905, N = 13, conf_level = .95)
#>  term        value
#>  lower_limit 1.05 
#>  std_mean    2.04 
#>  upper_limit 3    
#> 
#> Confidence level: 95%
ci_sm(mean = 30, sd = 14.721, N = 13, conf_level = .95)
#>  term        value
#>  lower_limit 1.05 
#>  std_mean    2.04 
#>  upper_limit 3    
#> 
#> Confidence level: 95%
ci_sm(ncp = 7.347771, N = 13, conf_level = .95)
#>  term        value
#>  lower_limit 1.05 
#>  std_mean    2.04 
#>  upper_limit 3    
#> 
#> Confidence level: 95%
ci_sm(sm = 2.037905, N = 13, alpha_lower = .05, alpha_upper = 0)
#>  term        value
#>  lower_limit 1.2  
#>  std_mean    2.04 
#>  upper_limit Inf  
ci_sm(mean = 50, sd = 10, N = 25, conf_level = .95)
#>  term        value
#>  lower_limit 3.54 
#>  std_mean    5    
#>  upper_limit 6.45 
#> 
#> Confidence level: 95%
```
