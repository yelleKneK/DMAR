# Confidence Interval for an Unstandardized Regression Coefficient

Computes a confidence interval for a population regression coefficient
in its raw (unstandardized) metric, by the standard *t*-based approach
or the noncentral *t* approach. A thin convenience wrapper around
[`ci_reg_coef`](https://yelleknek.github.io/DMAR/reference/ci_reg_coef.md),
which is the general engine; for the standardized coefficient use
[`ci_src`](https://yelleknek.github.io/DMAR/reference/ci_src.md).

## Usage

``` r
ci_rc(
  reg_coef_j = NULL,
  se_reg_coef_j = NULL,
  sd_Y = NULL,
  sd_X = NULL,
  N,
  p,
  R2_Y_X = NULL,
  R2_j_X_without_j = NULL,
  conf_level = 0.95,
  R2_Y_X_without_j = NULL,
  t_value = NULL,
  alpha_lower = NULL,
  alpha_upper = NULL,
  noncentral = FALSE
)
```

## Arguments

- reg_coef_j:

  Value of the regression coefficient for the *j*th predictor variable,
  or `NULL` (the default) when the other inputs determine it: `t_value`
  with `R2_Y_X` and `R2_j_X_without_j`, or the three squared multiple
  correlations, which fix its magnitude but not its sign (see
  [`ci_reg_coef`](https://yelleknek.github.io/DMAR/reference/ci_reg_coef.md))

- se_reg_coef_j:

  Standard error for the *j*th predictor variable

- sd_Y:

  Standard deviation of *Y*, the dependent variable

- sd_X:

  Standard deviation of *X*, the predictor variable of interest

- N:

  Sample size

- p:

  The number of predictors

- R2_Y_X:

  The squared multiple correlation coefficient predicting *Y* from the
  *p* predictor variables

- R2_j_X_without_j:

  The squared multiple correlation coefficient predicting the *j*th
  predictor variable (i.e., the predictor of interest) from the
  remaining *p*-1 predictor variables

- conf_level:

  Desired level of confidence for the computed interval (i.e., 1 - the
  Type I error rate). A level above .999999998, whose tail areas are
  below 1e-9, is refused

- R2_Y_X_without_j:

  The squared multiple correlation coefficient predicting *Y* from the
  *p*-1 predictor variable with the *j*th predictor of interest excluded

- t_value:

  The *t*-value evaluating the null hypothesis that the population
  regression coefficient for the *j*th predictor equals zero

- alpha_lower:

  The Type I error rate for the lower confidence interval limit

- alpha_upper:

  The Type I error rate for the upper confidence interval limit. The
  table records the coverage, `1 - alpha_lower - alpha_upper`, as its
  level and, when the two tail areas differ, the areas themselves, which
  the line beneath the printed table,
  [`results_sentence`](https://yelleknek.github.io/DMAR/reference/dmar_output_helpers.md),
  and the `alpha_lower` and `alpha_upper` columns of `tidy()` state;
  with a tail area of 0 the interval is a one-sided confidence bound. A
  nonzero tail area below 1e-9 is refused (see "Small tail areas" in
  [`ci_nc_t`](https://yelleknek.github.io/DMAR/reference/ci_nc_t.md)).

- noncentral:

  `TRUE` or `FALSE` statement specifying whether or not the noncentral
  approach to confidence intervals should be used

## Value

A 3-row `data.frame` with columns `term`, `value`, `prob_less`, and
`prob_greater`. The `term` values are `"ci_lower"`, `"reg_coef"` (the
coefficient as given, or as the other inputs determine it), and
`"ci_upper"`, and `value` holds the coefficient and its confidence
limits in the raw metric. The `prob_less` and `prob_greater` columns
report the tail probabilities below and above each limit; when the
noncentral *t* approach is used they are the achieved tail
probabilities; they are `NA` on the row of the estimate. The estimate
sits between its limits, as in
[`ci_src`](https://yelleknek.github.io/DMAR/reference/ci_src.md) and
[`ci_reg_coef`](https://yelleknek.github.io/DMAR/reference/ci_reg_coef.md).

## Details

Returns the confidence limits for the regression coefficient of interest
from the standard approach to confidence interval formation or from the
noncentral approach to confidence interval formation using the
noncentral *t*-distribution.

## Note

Not all of the values need to be specified, only those that contain all
of the necessary information in order to compute the confidence interval
(options are thus given for the values that need to be specified).

## References

Kelley, K. (2007). Confidence intervals for standardized effect sizes:
Theory, application, and implementation. *Journal of Statistical
Software, 20*(8), 1–24.
[doi:10.18637/jss.v020.i08](https://doi.org/10.18637/jss.v020.i08)

Kelley, K., & Maxwell, S. E. (2003). Sample size for multiple
regression: Obtaining regression coefficients that are accurate, not
simply significant. *Psychological Methods, 8*(3), 305–321.
[doi:10.1037/1082-989X.8.3.305](https://doi.org/10.1037/1082-989X.8.3.305)

Kelley, K., & Maxwell, S. E. (2008). Sample size planning with
applications to multiple regression: Power and accuracy for omnibus and
targeted effects. In P. Alasuutari, L. Bickman, & J. Brannen (Eds.),
*The Sage handbook of social research methods* (pp. 166–192). Sage.

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 4 on individual comparisons of means and
Chapter 6 on trend analysis.)

Smithson, M. (2003). *Confidence intervals*. Thousand Oaks, CA: Sage
Publications.

Steiger, J. H. (2004). Beyond the *F* Test: Effect size confidence
intervals and tests of close fit in the Analysis of Variance and
Contrast Analysis. *Psychological Methods, 9*(2), 164–182.
[doi:10.1037/1082-989X.9.2.164](https://doi.org/10.1037/1082-989X.9.2.164)

## See also

[`ss_aipe_reg_coef`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reg_coef.md),
[`ci_nc_t`](https://yelleknek.github.io/DMAR/reference/ci_nc_t.md),
[`ci_reg_coef`](https://yelleknek.github.io/DMAR/reference/ci_reg_coef.md),
[`ci_src`](https://yelleknek.github.io/DMAR/reference/ci_src.md)

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
ci_rc(reg_coef_j = 0.61319, se_reg_coef_j = 0.16098, N = 30, p = 6, conf_level = 0.95)
#>  term     value prob_less prob_greater
#>  ci_lower 0.280 0.025     0.975       
#>  reg_coef 0.613 <NA>      <NA>        
#>  ci_upper 0.946 0.975     0.025       
#> 
#> Confidence level: 95%
```
