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

  Type I error for the upper confidence limit. The table records the
  coverage, `1 - alpha_lower - alpha_upper`, as its level and, when the
  two tail areas differ, the areas themselves, which the line beneath
  the printed table,
  [`results_sentence`](https://yelleknek.github.io/DMAR/reference/dmar_output_helpers.md),
  and the `alpha_lower` and `alpha_upper` columns of `tidy()` state;
  with a tail area of 0 the interval is a one-sided confidence bound.

- ...:

  Allows one to potentially include parameter values for inner functions

## Value

A 2-row `data.frame` with columns `term` and `value`. The `term` values
are `"ci_lower"` and `"ci_upper"`, giving the lower and upper confidence
limits on the signal-to-noise ratio.

## Details

The confidence level must be specified in one of following two ways:
using confidence interval coverage (`conf_level`), or lower and upper
confidence limits (`alpha_lower` and `alpha_upper`).

This function uses the confidence interval transformation principle
(Steiger, 2004) to transform the confidence limits for the noncentrality
parameter to the confidence limits for the population's signal-to-noise
ratio. The confidence interval for noncentral *F* parameter can be
obtained from the `ci_nc_F` function in DMAR, which is used internally
within this function.

## Note

The signal to noise ratio is defined as the variance due to the
particular factor over the error variance (i.e., the mean square error).

## References

Bargmann, R. E. (1970). Interpretation and use of a generalized
discriminant function. In R. C. Bose, I. M. Chakravarti, P. C.
Mahalanobis, C. R. Rao, & K. J. C. Smith (Eds.), *Essays in probability
and statistics* (pp. 35–60). University of North Carolina Press.

Fleishman, A. I. (1980). Confidence intervals for correlation ratios.
*Educational and Psychological Measurement, 40*(3), 659–670.
[doi:10.1177/001316448004000309](https://doi.org/10.1177/001316448004000309)

Kelley, K. (2007). Confidence intervals for standardized effect sizes:
Theory, application, and implementation. *Journal of Statistical
Software, 20*(8), 1–24.
[doi:10.18637/jss.v020.i08](https://doi.org/10.18637/jss.v020.i08)

Steiger, J. H. (2004). Beyond the *F* Test: Effect size confidence
intervals and tests of close fit in the Analysis of Variance and
Contrast Analysis. *Psychological Methods, 9*(2), 164–182.
[doi:10.1037/1082-989X.9.2.164](https://doi.org/10.1037/1082-989X.9.2.164)

Venables, W. (1975). Calculation of confidence intervals for
noncentrality parameters. *Journal of the Royal Statistical Society,
Series B, 37*(3), 406–412.
[doi:10.1111/j.2517-6161.1975.tb01554.x](https://doi.org/10.1111/j.2517-6161.1975.tb01554.x)

## See also

[`ci_srsnr`](https://yelleknek.github.io/DMAR/reference/ci_srsnr.md),
[`ci_nc_F`](https://yelleknek.github.io/DMAR/reference/ci_nc_F.md)

Other confidence intervals for effect sizes:
[`ci_R2()`](https://yelleknek.github.io/DMAR/reference/ci_R2.md),
[`ci_c()`](https://yelleknek.github.io/DMAR/reference/ci_c.md),
[`ci_c_ancova()`](https://yelleknek.github.io/DMAR/reference/ci_c_ancova.md),
[`ci_c_ancova_bp()`](https://yelleknek.github.io/DMAR/reference/ci_c_ancova_bp.md),
[`ci_correlation`](https://yelleknek.github.io/DMAR/reference/ci_correlation.md),
[`ci_cv()`](https://yelleknek.github.io/DMAR/reference/ci_cv.md),
[`ci_eta_squared()`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared.md),
[`ci_eta_squared_generalized()`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared_generalized.md),
[`ci_mahalanobis()`](https://yelleknek.github.io/DMAR/reference/ci_mahalanobis.md),
[`ci_omega_squared()`](https://yelleknek.github.io/DMAR/reference/ci_omega_squared.md),
[`ci_pvaf()`](https://yelleknek.github.io/DMAR/reference/ci_pvaf.md),
[`ci_ratio_of_means()`](https://yelleknek.github.io/DMAR/reference/ci_ratio_of_means.md),
[`ci_ratio_of_means_paired()`](https://yelleknek.github.io/DMAR/reference/ci_ratio_of_means_paired.md),
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
## Bargmann (1970) gave an example of a one-way ANOVA with an observed
## F value of 11.2213 on 4 and 50 degrees of freedom, which Venables (1975),
## Fleishman (1980), and Steiger (2004) used in turn. Steiger describes it
## as five groups of 11 and reports F as 11.221, the value used here.
## The exact confidence interval for the signal-to-noise ratio of that example:

ci_snr(F_value = 11.221, df_1 = 4, df_2 = 50, N = 55)
#>  term     value
#>  ci_lower 0.293
#>  ci_upper 1.423
#> 
#> Confidence level: 95%

ci_snr(F_value = 11.221, df_1 = 4, df_2 = 50, N = 55, conf_level = .90)
#>  term     value
#>  ci_lower 0.352
#>  ci_upper 1.301
#> 
#> Confidence level: 90%

ci_snr(F_value = 11.221, df_1 = 4, df_2 = 50, N = 55, alpha_lower = .02, alpha_upper = .03)
#>  term     value
#>  ci_lower 0.276
#>  ci_upper 1.392
#> 
#> Confidence level: 95% (unequal tail areas: 0.02 below the lower limit, 0.03 above the upper limit)
```
