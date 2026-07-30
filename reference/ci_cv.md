# Confidence Interval for the Coefficient of Variation

Computes the noncentral *t*-based confidence interval for the population
coefficient of variation, the standard deviation relative to the mean,
so variability can be reported on a scale that is free of the
measurement units.

## Usage

``` r
ci_cv(
  cv = NULL,
  mean = NULL,
  sd = NULL,
  n = NULL,
  data = NULL,
  conf_level = 0.95,
  alpha_lower = NULL,
  alpha_upper = NULL,
  ...
)
```

## Arguments

- cv:

  Coefficient of variation

- mean:

  Sample mean

- sd:

  Sample standard deviation (square root of the unbiased estimate of the
  variance

- n:

  Sample size

- data:

  Vector of data for which the confidence interval for the coefficient
  of variation is to be calculated

- conf_level:

  Desired confidence level (1-Type I error rate)

- alpha_lower:

  The proportion of values beyond the lower limit of the confidence
  interval (cannot be used with `conf_level`)

- alpha_upper:

  The proportion of values beyond the upper limit of the confidence
  interval (cannot be used with `conf_level`)

- ...:

  Allows one to potentially include parameter values for inner functions

## Value

A 4-row `data.frame` with columns `term`, `value`, `prob_less`, and
`prob_greater`. The rows are ordered so the two point estimates sit
between the confidence limits: `"lower_limit"` (lower confidence limit
on the coefficient of variation), `"c_of_v"` (the sample coefficient of
variation), `"c_of_v_unbiased"` (the unbiased estimator), and
`"upper_limit"` (upper confidence limit). The `prob_less` and
`prob_greater` columns report the achieved tail probabilities of the
noncentral t search at the limit values; they are `NA` for the
point-estimate rows.

## Details

Uses the noncentral *t*-distribution to calculate the confidence
interval for the population coefficient of variation.

## References

Johnson, N. L., & Welch, B. L. (1940). Applications of the non-central
*t*-distribution. *Biometrika, 31*(3–4), 362–389.
[doi:10.1093/biomet/31.3-4.362](https://doi.org/10.1093/biomet/31.3-4.362)

Kelley, K. (2007). Sample size planning for the coefficient of variation
from the accuracy in parameter estimation approach. *Behavior Research
Methods, 39*(4), 755–766.
[doi:10.3758/BF03192966](https://doi.org/10.3758/BF03192966)

Kelley, K. (2007). Confidence intervals for standardized effect sizes:
Theory, application, and implementation. *Journal of Statistical
Software, 20*(8), 1–24.
[doi:10.18637/jss.v020.i08](https://doi.org/10.18637/jss.v020.i08)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 3.)

McKay, A. T. (1932). Distribution of the coefficient of variation and
the extended *t* distribution. *Journal of the Royal Statistical
Society, 95*(4), 695–698.

## See also

[`cv`](https://yelleknek.github.io/DMAR/reference/cv.md)

Other confidence intervals for effect sizes:
[`ci_R2()`](https://yelleknek.github.io/DMAR/reference/ci_R2.md),
[`ci_c()`](https://yelleknek.github.io/DMAR/reference/ci_c.md),
[`ci_c_ancova()`](https://yelleknek.github.io/DMAR/reference/ci_c_ancova.md),
[`ci_c_ancova_bp()`](https://yelleknek.github.io/DMAR/reference/ci_c_ancova_bp.md),
[`ci_cc()`](https://yelleknek.github.io/DMAR/reference/ci_cc.md),
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
set.seed(113)
N <- 15
X <- rnorm(N, 5, 1)
mean.X <- mean(X)
sd.X <- var(X)^.5

ci_cv(mean = mean.X, sd = sd.X, n = N, alpha_lower = .025,
alpha_upper = .025, conf_level = NULL)
#>  term            value prob_less prob_greater
#>  lower_limit     0.15  0.025     0.975       
#>  c_of_v          0.207 <NA>      <NA>        
#>  c_of_v_unbiased 0.21  <NA>      <NA>        
#>  upper_limit     0.335 0.975     0.025       
ci_cv(data = X, conf_level = .95)
#>  term            value prob_less prob_greater
#>  lower_limit     0.15  0.025     0.975       
#>  c_of_v          0.207 <NA>      <NA>        
#>  c_of_v_unbiased 0.21  <NA>      <NA>        
#>  upper_limit     0.335 0.975     0.025       
#> 
#> Confidence level: 95%
ci_cv(cv = sd.X / mean.X, n = N, conf_level = .95)
#>  term            value prob_less prob_greater
#>  lower_limit     0.15  0.025     0.975       
#>  c_of_v          0.207 <NA>      <NA>        
#>  c_of_v_unbiased 0.21  <NA>      <NA>        
#>  upper_limit     0.335 0.975     0.025       
#> 
#> Confidence level: 95%
```
