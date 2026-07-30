# Confidence Interval for the Population Correlation Coefficient

This function forms a confidence interval for the population correlation
coefficient \\\rho\\. The confidence interval is for the population
value \\\rho\\; the required input is the corresponding sample value,
the observed sample correlation coefficient *r*. This approach assumes
that the two variables on which the correlation is based are bivariate
normally distributed (e.g., Hays, 1994, Chapter 14).

## Usage

``` r
ci_cc(r, n, conf_level = 0.95, alpha_lower = NULL, alpha_upper = NULL)
```

## Arguments

- r:

  Observed value of the sample correlation coefficient (specifically the
  zero-order Pearson's product-moment correlation coefficient)

- n:

  sample size

- conf_level:

  Desired confidence level, where the error rate is the same on each
  side

- alpha_lower:

  The Type I error rate for the lower confidence interval limit

- alpha_upper:

  The Type I error rate for the upper confidence interval limit

## Value

A 3-row `data.frame` with columns `term` and `value`. The `term` values
are `"lower_limit"` (the lower confidence limit on the population
correlation \\\rho\\), `"est_cor"` (the observed sample correlation
coefficient), and `"upper_limit"` (the upper limit on \\\rho\\).

## Details

Note that this approach to confidence intervals will not generally lead
to a symmetric confidence interval. The function first transforms \\r\\
into *Z*', forms a confidence interval for the population value (i.e.,
\\\zeta\\), and then transforms the confidence limits for \\\zeta\\ into
the scale of the correlation coefficient.

## Note

This confidence interval assumes that the two variables the correlation
is based are bivariate normal. See Hays (1994, Chapter 14) for details.

## References

Kelley, K. (2007). Confidence intervals for standardized effect sizes:
Theory, application, and implementation. *Journal of Statistical
Software, 20*(8), 1–24.
[doi:10.18637/jss.v020.i08](https://doi.org/10.18637/jss.v020.i08)

Hays, W. L. (1994). *Statistics* (5th ed.). Fort Worth, TX: Harcourt
Brace College Publishers.

## See also

[`convert_z_r`](https://yelleknek.github.io/DMAR/reference/convert_Z_r.md),
[`convert_r_z`](https://yelleknek.github.io/DMAR/reference/convert_r_Z.md)

Other confidence intervals for effect sizes:
[`ci_R2()`](https://yelleknek.github.io/DMAR/reference/ci_R2.md),
[`ci_c()`](https://yelleknek.github.io/DMAR/reference/ci_c.md),
[`ci_c_ancova()`](https://yelleknek.github.io/DMAR/reference/ci_c_ancova.md),
[`ci_c_ancova_bp()`](https://yelleknek.github.io/DMAR/reference/ci_c_ancova_bp.md),
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
[`ci_srsnr()`](https://yelleknek.github.io/DMAR/reference/ci_srsnr.md),
[`contrast_adjusted()`](https://yelleknek.github.io/DMAR/reference/contrast_adjusted.md),
[`plot_smd()`](https://yelleknek.github.io/DMAR/reference/plot_smd.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Example, from Hays. Suppose n=100 and r=.35.
ci_cc(r = .35, n = 100, conf_level = .95)
#>  term        value
#>  lower_limit 0.165
#>  est_cor     0.35 
#>  upper_limit 0.511
#> 
#> Confidence level: 95%

# Here is another way to enter the above example.
ci_cc(r = .35, n = 100, conf_level = NULL, alpha_lower = .025, alpha_upper = .025)
#>  term        value
#>  lower_limit 0.165
#>  est_cor     0.35 
#>  upper_limit 0.511

# Here are examples of one-sided confidence intervals.
ci_cc(r = .35, n = 100, conf_level = NULL, alpha_lower = 0, alpha_upper = .05)
#>  term        value
#>  lower_limit -1   
#>  est_cor     0.35 
#>  upper_limit 0.487
ci_cc(r = .35, n = 100, conf_level = NULL, alpha_lower = .05, alpha_upper = 0)
#>  term        value
#>  lower_limit 0.196
#>  est_cor     0.35 
#>  upper_limit 1    
```
